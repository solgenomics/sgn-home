#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <getopt.h>

#include "kp_types.h"
#include "log_message.h"

static uchar *output_basename = NULL;
int verbosity_level = 0;

#define ALLOC_STEP (128)

inline int white_space(uchar c) {

  if (c == '\n' || c == '\t' || c == '\r' || c == ' ') return 1;
  return 0;
}

inline int nucleotide(uchar c) {
  
  /* Lower case letter to simply if statement */
  if (c<97) c+=32;

  if (c == 'a' || c == 'c' || c == 't' || c == 'g' || c == 'n' || c == 'x')
    return 1;
  return 0;
}

/* openfile()
   Input:  filename to open and open (read/write) mode.
   Output: Returns a FILE * pointer to the open file.

   Purpose: Eliminate the if (f == NULL) statement from the main code. All
   file open errors are fatal to this program, so this function will exit
   the program via logmsg(MSG_FATAL,...) if fopen() fails. <filetype> 
   parameter allows more verbose error message to the user indicating
   (hopefully) what we thought we were going to do with the file. */
static FILE *openfile(uchar *filename, uchar *mode, uchar *filetype) {
  FILE *f;

  f = fopen(filename, mode);
  if (f == NULL) {
    logmsg(MSG_FATAL,"Can't open %s \"%s\" (%s)\n",
	   filetype, filename, strerror(errno));
  }

  return f;
}


static inline int base_toint(uchar base) {

  switch(base) {
  case 'A':
    return 0;
  case 'C':
    return 1;
  case 'G':
    return 2;
  case 'T':
    return 3;
  default:
    return 0;
  }
}

static uint read_fullline(uchar **rline, uint *space, FILE *f) {
  uchar inputline[1024];
  int length, i;

  length = 0;
  while(fgets(inputline, 1024, f) != NULL) {
    i = 0;
    while(inputline[i]!=0 && inputline[i]!='\n') {
      if (*space <= length) {
	*rline = realloc(*rline, *space + ALLOC_STEP);
	if (*rline == NULL) {
	  logmsg(MSG_FATAL,"Failed allocating memory to read input file "
		 "(%d bytes)\n",*space + ALLOC_STEP);
	}
	*space += ALLOC_STEP;
      }
      (*rline)[length++] = inputline[i++];
    }
    if (inputline[i] == '\n') break;
  }

  if (*space <= length) {
    *rline = realloc(*rline, *space + ALLOC_STEP);
    if (*rline == NULL) {
      logmsg(MSG_FATAL,"Failed allocating memory to read input file "
	     "(%d bytes)\n",*space + ALLOC_STEP);
    }
    *space += ALLOC_STEP;
  }
  (*rline)[length] = 0;

  return length;
}

static uint format_input(uchar *input_seqfile, FILE *indfile, FILE *strfile,
			 FILE *binfile) {
  FILE *sf;
  uint line_no, input_length, name_start, name_end, i, j;
  uint input_allocsize;
  uchar *inputline;
  seqmeta_t *seqmeta;
  seqmeta_t *seq;
  uchar *sequence, *binseq, *comp;
  uchar **seq_names;
  int binsize, seqsize;
  int seq_length;
  int seq_id;

  int name_ptr, binfile_ptr, strfile_ptr;

  inputline = NULL;
  input_allocsize = 0;
  binfile_ptr = strfile_ptr = 4;
  name_ptr = 0;

  /* This holds metadata for output when we are finished. */
  MA(seq, sizeof(seqmeta_t *));
  seqmeta = NULL; 
  seq_names = NULL;

  /* This holds sequence as its read in */
  sequence = NULL;
  binseq = NULL;
  binsize = 0;
  seqsize = 0;

  sf = openfile(input_seqfile, "r", "FASTA sequence file");

  input_length = read_fullline(&inputline, &input_allocsize, sf);
  line_no = 1;
  seq_id = 0;
  while(input_length>0) {
    PUSH(seqmeta, seq_id, sizeof(seqmeta_t));
    seq = seqmeta + seq_id;

    /* Assumption: inputline at this place always contains a FASTA 
       header line */
    if (inputline[0]!='>') {
      logmsg(MSG_FATAL,"FASTA parse error at line %u in file %s: header"
	     " line expected, beginning with '>'",line_no, input_seqfile);
    }
    
    /* Recover the sequence name */
    name_start = 1;
    while(name_start < input_length && white_space(inputline[name_start])) 
      name_start++;
    name_end = name_start + 1;
    while(name_end < input_length && !white_space(inputline[name_end]))
      name_end++;

    if (name_start == input_length) {
      logmsg(MSG_FATAL,"FASTA parse error at line %u in file %s:\nsequence "
	     "name not found in header: %s\n",line_no, input_seqfile, 
	     inputline);
    }

    inputline[name_end] = 0;

    PUSH(seq_names, seq_id, sizeof(uchar *));
    seq_names[seq_id] = strdup(inputline + name_start);

    seq->name_length = name_end - name_start;
    seq->name_pos = name_ptr;
    name_ptr += seq->name_length + 1;

    /* Read in sequence data */
    seq_length = 0;
    input_length = read_fullline(&inputline, &input_allocsize, sf);
    line_no++;
    i = 0;
    while(input_length>0 && inputline[0]!='>') {
      j=0;
      if (seq_length + input_length >= seqsize) {
	seqsize += input_length;
	RA(sequence, sizeof(uchar), seqsize);
      }	
      while(j<input_length) {
	if (nucleotide(inputline[j])) {
	  sequence[seq_length++] = inputline[j];
	}
	j++;
      }
      input_length = read_fullline(&inputline, &input_allocsize, sf);
      line_no++;
    }

    seq->seq_length = seq_length;
    seq->seqstr_pos = strfile_ptr;
    seq->seqbin_pos = binfile_ptr;

    /* Output the sequence in text (string) form to the string file */
    fwrite(sequence, sizeof(uchar), seq_length, strfile);

    if (seq_length > binsize) {
      binsize = seq_length;
      RA(binseq, sizeof(uchar), binsize);
    }

    /* Convert the sequence to binary form (not packed) */
    for(j=0;j<seq_length;j++) {
      binseq[j] = base_toint(sequence[j]);
    }
    fwrite(binseq, sizeof(uchar), seq_length, binfile);

    binfile_ptr += seq_length;
    strfile_ptr += seq_length;
    seq_id++;

#if 0
    /* Generate a reverse complement meta record */
    PUSH(seq_names, seq_id, sizeof(uchar *));
    MA(seq_names[seq_id], seq->name_length+2);
    strcpy(seq_names[seq_id], seq_names[seq_id-1]);
    strcat(seq_names[seq_id], "-");

    PUSH(seqmeta, seq_id, sizeof(seqmeta_t));
    seq = seqmeta + seq_id;

    seq->name_length = strlen(seq_names[seq_id]);
    seq->name_pos = name_ptr;
    seq->seq_length = seq_length;
    seq->seqstr_pos = strfile_ptr;
    seq->seqbin_pos = binfile_ptr;

    /* Generate the actual reverse complement text sequence from the 
       "template strand" */
    MA(comp, sizeof(uchar)*seq_length);
    for(j=0;j<seq_length;j++) {
      switch(sequence[seq_length-j-1]) {
      case 'A':
	comp[j] = 'T';
	break;
      case 'C':
	comp[j] = 'G';
	break;
      case 'G':
	comp[j] = 'C';
	break;
      case 'T':
	comp[j] = 'A';
	break;
      case 'N':
      case 'X':
	comp[j] = sequence[seq_length-j-1];
	break;
      default:
	logmsg(MSG_ERROR,"Unknown nucleotide '%c' in sequence %s at position "
	       "%d\n",sequence[seq_length-j-1],seq_names[seq_id-1],
	       seq_length-j-1);
	break;
      }
    }

    /* Write out complement sequence to string file */
    fwrite(comp, sizeof(uchar), seq_length, strfile);

    /* Build complement sequence in binary form and write out to binary file */
    for(j=0;j<seq_length;j++) {
      binseq[j] = base_toint(comp[j]);
    }
    fwrite(binseq, sizeof(uchar), seq_length, binfile);

    name_ptr += seq->name_length;
    strfile_ptr += seq->seq_length;
    binfile_ptr += seq->seq_length;
    seq_id++;

    free(comp);
#endif
  }

  fclose(sf);

  free(inputline);
  free(binseq);
  free(sequence);

  fwrite(&seq_id, sizeof(uint), 1, indfile);
  fwrite(seqmeta, sizeof(seqmeta_t), seq_id, indfile);
  for(i=0;i<seq_id;i++) {
    fwrite(seq_names[i], sizeof(uchar), seqmeta[i].name_length + 1, indfile);
  }

  return seq_id;
}

#if 0
void generate_reverse_complement(void) {
  int i;

  logmsg(MSG_INFO,"Generating corresponding reverse complement sequences..\n");

  sequences = realloc(sequences, sizeof(seq_t *)*n_seq*2);
  if (sequences == NULL) {
    logmsg(MSG_FATAL,"Failed allocating memory for reverse complement "
	   "sequences (%d bytes)\n",sizeof(seq_t *)*n_seq*2);
  }

  for(i=0;i<n_seq;i++) {
    seq_t *comp, *source;
    int j,k;

    source = sequences[i];
    MA(comp, sizeof(seq_t));
    sequences[i + n_seq] = comp;
    comp->seq = NULL;

    MA(comp->label, sizeof(uchar*)*(strlen(source->label)+2));
    strcpy(comp->label, source->label);
    comp->label[strlen(source->label)] = '-';
    comp->label[strlen(source->label)+1] = 0;
    comp->length = source->length;
    MA(comp->seqstr, (source->length+1)*sizeof(uchar)); 
    MA(comp->qual, source->length*sizeof(int));

    for(j=0;j<source->length;j++) {
      k = source->length - j - 1;
      comp->qual[j] = source->qual[k];
      switch(source->seqstr[k]) {
      case 'A':
	comp->seqstr[j] = 'T';
	break;
      case 'C':
	comp->seqstr[j] = 'G';
	break;
      case 'G':
	comp->seqstr[j] = 'C';
	break;
      case 'T':
	comp->seqstr[j] = 'A';
	break;
      case 'N':
	comp->seqstr[j] = 'N';
	break;
      default:
	logmsg(MSG_ERROR,"Unknown nucleotide '%c' in sequence %s at position "
	       "%d\n",source->seqstr[k], source->label, k);
	break;
      }
    }
    comp->seqstr[j] = 0;    
  }
  
  n_seq = n_seq*2;
  compute_intsequence();
  logmsg(MSG_INFO,"Finished generating complement sequences\n");
}
#endif

static void usage(char *program_name) {

  fprintf(stderr,"\n\n%s:\n\n"
"Program to read FASTA format sequence data and format it for rapid approx.       \n"
"comparison for overlap detection and preclustering. Output of this program	  \n"
"is binary files meant to be read by other programs. This program does no real	  \n"
"work other than translation.							  \n"
"										  \n"
"Options:									  \n"
"--seqfile=<filename> (-s) (required)						  \n"
"    FASTA format input file to be translated/formatted				  \n"
"--basename=<string> (-o)							  \n"
"    String to use as prefix for output files created. Uses input filename by 	  \n"
"    default.									  \n"
"--verbose=<integer> (-v)                                                         \n"
"    Verbosity level. 0 (normal) by default. Negative enables debugging messages  \n"
"    Positive makes program quieter.						  \n"
"--help (-h)									  \n"
"    Prints this message.							  \n"
"										  \n"
,program_name);								  
										  
}										  
static void parse_arguments(uchar **seqfilename, int argc, char *argv[]) {	  
  int option_index, commandline_error, rval;					  
  struct option longopts[] = {
    { "seqfile", 1, NULL, 's'},
    { "basename", 1, NULL, 'o'},  
    { "verbose", 1, NULL, 'v'},
    { "help", 1, NULL, 'h'},
    { NULL, 0, NULL, 0}
  };
  char *optstring = "s:v:o:";

  *seqfilename = output_basename = NULL;
  commandline_error = 0;
  while((rval = getopt_long(argc, argv, optstring, longopts, &option_index))
	!= -1) {
    switch(rval) {
    case ':':
      logmsg(MSG_ERROR,"! Option \"%s\" requires an argument.\n",
	     longopts[option_index].name);
      commandline_error = 1;
      break;
    case 's':
      *seqfilename = strdup(optarg);
      break;
    case 'v':
      verbosity_level = atoi(optarg);
      break;
    case 'o':
      output_basename = strdup(optarg);
      break;
    case 'h':
      usage(argv[0]);
      exit(0);
      break;
    case '?':
    default:
      logmsg(MSG_ERROR,"! Option \"%c\" unknown.\n",optopt);
      commandline_error = 1;
      break;
    }
  }

  if (*seqfilename == NULL) {
    logmsg(MSG_ERROR,"! Input sequence file in FASTA format must be "
	   "specified with -s <filename> or --seqfile=<filename> option\n");
    commandline_error = 1;
  }
  
  if (commandline_error) {
    logmsg(MSG_ERROR,"! Program halted due to command line option errors\n");
    usage(argv[0]);
    exit(-1);
  }

  if (output_basename == NULL) {
    output_basename = strdup(*seqfilename);
  }
}

int main(int argc, char *argv[]) {
  uchar *input_seqfile;
  uchar *temp;
  uint n_seq, x;
  int l;
  FILE *indfile, *strfile, *binfile;

  configure_logmsg(MSG_DEBUG1);
  parse_arguments(&input_seqfile, argc, argv);
  configure_logmsg(verbosity_level);

  logmsg(MSG_INFO,"Output basename set to %s\n",output_basename);

  l = strlen(output_basename);
  MA(temp, (l+6)*sizeof(char));

  strcpy(temp, output_basename);
  strcat(temp, ".ind");
  indfile = openfile(temp, "w", "sequence index file");
  x = INDFILE_MAGIC;
  fwrite(&x, sizeof(uint), 1, indfile);
  
  strcpy(temp, output_basename);
  strcat(temp, ".seq");
  strfile = openfile(temp, "w", "sequence string file");
  x = STRFILE_MAGIC;
  fwrite(&x, sizeof(uint), 1, strfile);

  strcpy(temp, output_basename);
  strcat(temp, ".sbin");
  binfile = openfile(temp, "w", "sequence binary file");
  x = BINFILE_MAGIC;
  fwrite(&x, sizeof(uint), 1, binfile);
  
  n_seq = format_input(input_seqfile, indfile, strfile, binfile);
  fclose(indfile);
  fclose(strfile);
  fclose(binfile);
  logmsg(MSG_INFO,"%d sequences formatted\n",n_seq);

  return 0;
}
