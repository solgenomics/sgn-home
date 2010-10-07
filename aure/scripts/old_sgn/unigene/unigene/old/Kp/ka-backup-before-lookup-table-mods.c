#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>
#include <getopt.h>
#include <limits.h>

#include "ka.h"
#include "log_message.h"

void test_shit();

typedef struct {
  int s1;
  int s2;
  int s1_ltrim;
  int s2_ltrim;
  int s1_rtrim;
  int s2_rtrim;
  int length;
  float mismatch_score;
  float trim_score;
  float align_score;
} align_t;


sequence_t *sequences;
uchar **seq_names;
int n_seq;

static align_t **overlap;
static int *n_overlap;

static uint **nolist;
static int *n_nolist;

static uchar *output_basename = NULL;

static int verbosity_level = 0;

static void usage() {
  fprintf(stderr,"\nUsage:\n
--seqfile=<fasta input sequence file> (-s)
--qualfile=<fasta input phred quality file> (-q)
--help (-h) prints this message

\n");
}

static void parse_arguments(uchar **seqfilename, uchar **qualfilename,
			    int argc, char *argv[]) {
  int option_index, commandline_error, rval;
  struct option longopts[] = {
    { "seqfile", 1, NULL, 's' },
    { "qualfile", 1, NULL, 'q' },
    { "verbosity", 1, NULL, 'v' },
    { "output-basename", 1, NULL, 'o' },
    { "help", 1, NULL, 'h'},
    { NULL, 0, NULL, 0 }
  };

  char *optstring = "s:q:hv:o:";

  *seqfilename = *qualfilename = NULL;
  commandline_error = 0;
  while((rval = getopt_long(argc, argv, optstring, longopts, &option_index)) != -1) {
    switch(rval) {
    case ':':
      logmsg(MSG_ERROR,"! Option \"%s\" requires an argument.\n",longopts[option_index].name);
      commandline_error = 1;
      break;
    case 's':
      *seqfilename = strdup(optarg);
      break;
    case 'q':
      *qualfilename = strdup(optarg);
      break;
    case 'h':
      usage();
      exit(0);
      break;
    case 'v':
      verbosity_level = atoi(optarg);
      break;
    case 'o':
      output_basename = strdup(optarg);
      break;
    case '?':
    default:
      logmsg(MSG_ERROR,"! Option \"%c\" unknown.\n",optopt);
      commandline_error = 1;
      break;
    }
  }

  if (*seqfilename == NULL || *qualfilename == NULL) {
    logmsg(MSG_ERROR,"! You must specify FASTA files containing sequence and corresponding quality\n  information.\n");
    commandline_error = 1;
  }

  if (commandline_error) {
    logmsg(MSG_ERROR,"! Program halted due to command line option errors\n");
    usage();
    exit(-1);
  }

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

static void *local_calloc(size_t n, size_t size, uchar *description) {
  void *q;

  q = calloc(n, size);
  if (q == NULL) {
    logmsg(MSG_FATAL,"Failed allocating memory for %s (%d bytes)\n",
	   description, size);
  }
  return q;
}


static void *local_realloc(void *p, size_t size, uchar *description) {
  void *q;
  
  q = realloc(p, size);
  if (q == NULL) {
    logmsg(MSG_FATAL,"Failed allocating memory for %s (%d bytes)\n",
	   description, size);
  }
  return q;
}

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

#define ALLOC_STEP (128)
static uint read_fullline(uchar **rline, uint *space, FILE *f) {
  uchar inputline[1024];
  int length, i;

  length = 0;
  while(fgets(inputline, 1024, f) != NULL) {
    i = 0;
    while(inputline[i]!=0 && inputline[i]!='\n') {
      if (*space <= length) {
	*rline = local_realloc(*rline, *space + ALLOC_STEP, "loading input");
	*space += ALLOC_STEP;
      }
      (*rline)[length++] = inputline[i++];
    }
    if (inputline[i] == '\n') break;
  }

  if (*space <= length) {
    *rline = local_realloc(*rline, *space + ALLOC_STEP, "loading input");
    *space += ALLOC_STEP;
  }
  (*rline)[length] = 0;

  return length;
}

static void load_inputsequence(uchar *input_seqfile, uchar *input_qualfile) {
  FILE *sf, *qf;
  uint line_no, input_length, name_start, name_end, i, j;
  uint seqname_size, input_allocsize;
  uchar *inputline;
  uchar *seq;

  n_seq = seqname_size = input_allocsize = 0;
  seq_names = NULL;
  sequences = NULL;
  inputline = NULL;


  sf = openfile(input_seqfile, "r", "FASTA sequence file");
  qf = openfile(input_qualfile, "r", "FASTA quality file");

  input_length = read_fullline(&inputline, &input_allocsize, sf);
  line_no = 1;
  while(input_length>0) {
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
	     "name not found in header: %s\n",inputline);
    }

    inputline[name_end] = 0;

    if (seqname_size % ALLOC_STEP == 0) {
      seq_names = local_realloc(seq_names, 
				sizeof(uchar *)*(seqname_size + ALLOC_STEP),
				"index of sequence names");
    }
    seq_names[seqname_size] = strdup(inputline + name_start);

    input_length = read_fullline(&inputline, &input_allocsize, sf);
    line_no++;
    i = 0; seq = NULL;
    while(input_length>0 && inputline[0]!='>') {
      j=0;
      while(j<input_length) {
	if (nucleotide(inputline[j])) {
	  if (i % ALLOC_STEP == 0) {
	    seq = local_realloc(seq, sizeof(uchar)*(i+ALLOC_STEP+1),
				"loading input sequence");
	  }
	  seq[i++] = inputline[j];
	}
	j++;
      }
      input_length = read_fullline(&inputline, &input_allocsize, sf);
      line_no++;
    }
    seq[i] = 0;

    if (n_seq % ALLOC_STEP == 0) {
      sequences = local_realloc(sequences, 
				sizeof(sequence_t)*(n_seq+ALLOC_STEP),
				"storing table of input sequences");
    }
    sequences[n_seq].readname_index = seqname_size++;
    sequences[n_seq].length = i;
    seq = local_realloc(seq, i+1, "");
    sequences[n_seq].sequence = seq;
    sequences[n_seq].quality = NULL;
    n_seq++;
  }

  logmsg(MSG_INFO,"Loaded %d sequences from %s\n",n_seq, input_seqfile);
  fclose(sf);

  input_length = read_fullline(&inputline, &input_allocsize, qf);
  if (input_length == 0 || inputline[0]!='>') {
    logmsg(MSG_FATAL,"FASTA parse error at line %u in file %s: header"
	   " line expected, beginning with '>'",line_no, input_seqfile);
  }
  line_no = 1;
  while(input_length>0) {
    uint seqindex;

    /* Recover the sequence name */
    name_start = 1;
    while(name_start < input_length && white_space(inputline[name_start])) 
      name_start++;
    name_end = name_start + 1;
    while(name_end < input_length && !white_space(inputline[name_end]))
      name_end++;

    if (name_start == input_length) {
      logmsg(MSG_FATAL,"FASTA parse error at line %u in file %s:\nsequence "
	     "name not found in header: %s\n",inputline);
    }

    inputline[name_end] = 0;
    /* Yes, linear search is stupid... */
    i = 0;
    while(i<seqname_size) {
      if (strcmp(inputline + name_start, seq_names[i]) == 0) break;
      i++;
    }
    if (i == seqname_size) {
      logmsg(MSG_WARNING,"Sequence %s in quality file was not found in "
	     "FASTA sequence input file.\n",inputline + name_start);
      do {
	input_length = read_fullline(&inputline, &input_allocsize, qf);
	line_no++;
      } while(input_length>0 && inputline[0]!='>');
      continue;
    }
    seqindex = i;
    if (sequences[seqindex].quality != NULL) {
      logmsg(MSG_FATAL,"Sequence %s has more than one entry in quality file\n",
	     seq_names[seqindex]);
    }
    sequences[seqindex].quality = local_calloc(sequences[seqindex].length,
					       sizeof(uint),
					       "sequence quality values");

    input_length = read_fullline(&inputline, &input_allocsize, qf);
    line_no++;
    i = 0;
    while(input_length>0 && inputline[0]!='>') {
      j = 0;
      while(j<input_length) {
	uint score, k;

	while(j<input_length && white_space(inputline[j])) j++;
	k = j;
	while(k<input_length && inputline[k]>='0' && inputline[k]<='9') k++;
	if (k<input_length && !white_space(inputline[k])) {
	  logmsg(MSG_FATAL,"FASTA quality parse error at line %u: non-numeric"
		 " characters found where phred quality values expected\n",
		 line_no);
	}
	score = 0;
	while(j<k) {
	  score = score*10 + (inputline[j] - '0');
	  j++;
	}
	if (i >= sequences[seqindex].length) {
	  logmsg(MSG_FATAL,"FASTA quality parse error for sequence %s:\n"
		 "more quality values found than sequence letters\n",
		 seq_names[seqindex]);
	}
	sequences[seqindex].quality[i++] = score;
      }
      input_length = read_fullline(&inputline, &input_allocsize, qf);
      line_no++;
    }

    if (i != sequences[seqindex].length) {
      logmsg(MSG_FATAL,"FASTA quality parse error for sequence %s:\n",
	     "less quality values found than sequence letters\n",
	     seq_names[seqindex]);
    }
  }
  for(i=0;i<n_seq;i++) {
    if (sequences[i].quality == NULL) {
      logmsg(MSG_ERROR,"Sequence %s has no entry in quality file\n",
	     seq_names[i]);
    }
  }
  logmsg(MSG_INFO,"Loaded %d corresponding quality scores from %s\n",
	 n_seq, input_qualfile);
  fclose(qf);

  free(inputline);

}

static void debug_sequence(uchar *name, sequence_t *seq) {
  int i;

  logmsg(MSG_DEBUG4,">%s\n%s\n",name,seq->sequence);
  for(i=0;i<seq->length;i++) {
    logmsg(MSG_DEBUG4,"%d ",seq->quality[i]);
  }
  logmsg(MSG_DEBUG4,"\n");

}

static void polya_truncate() {
  int i,j,k;
  int na, qual, pasta_qual;
  sequence_t *seqobj;

  typedef struct {
    uint truncate_pos;
    uint a_length;
    uint a_qual;
    uint pasta_qual;
  } candidate_t;

  candidate_t *candidates;
  int cand, cand_allocate;

  candidates = NULL;
  cand_allocate = 0;
  for(i=0;i<n_seq;i++) {
    j = 0;
    seqobj = sequences + i;

    cand = 0;
    while(j<seqobj->length) {
      if (seqobj->sequence[j++] == 'A') {
	na = 1; qual = seqobj->quality[j-1];
	while(j<seqobj->length && seqobj->sequence[j] == 'A') {
	  na++;
	  qual += seqobj->quality[j++];
	}
	if (na > 11) {
	  k = j;
	  pasta_qual = 0;
	  while(k<seqobj->length && (k-j)<na) {
	    pasta_qual += seqobj->quality[k++];
	  }
	  if (k>j && qual/na > (pasta_qual/(k-j))*1.5 && 
	      seqobj->length-j<seqobj->length/3) {
	    if (cand >= cand_allocate) {
	      cand_allocate += 10;
	      candidates = local_realloc(candidates, 
					 sizeof(candidate_t)*cand_allocate,
					 "candidates for poly-A truncation");
	    }
	    candidates[cand].truncate_pos = j;
	    candidates[cand].a_length = na;
	    candidates[cand].a_qual = qual/na;
	    candidates[cand].pasta_qual = pasta_qual/(k-j);
	    cand++;
	  }
	}
      }
    }

    if (cand > 0) {
      int select;
      /* In preliminary testing of this simple algorithm, I have never seen
	 more than one candidate presented. However, if more than one is
	 found, we will truncate at the one which has the longer A
	 run. Other alternatives might be always truncate at the position
	 where we truncate the least amount of sequence, or at the position
	 where the quality difference is sharpest. */
      select = 0;
      if (cand > 1) {
	for(j=0;j<cand;j++) {
	  if (candidates[j].a_length > candidates[select].a_length)
	    select = j;
	}
      }
      logmsg(MSG_DEBUG1,"Truncating post poly-A noise for sequence %s at "
	     "position %d\n", seq_names[i], candidates[select].truncate_pos);
      seqobj->sequence = local_realloc(seqobj->sequence,
					    candidates[select].truncate_pos + 1,
					    "PolyA Truncation");
      seqobj->sequence[candidates[select].truncate_pos] = 0;
      seqobj->quality = local_realloc(seqobj->quality,
				      candidates[select].truncate_pos + 1,
				      "PolyA Quality Truncation");
      seqobj->length = candidates[select].truncate_pos;
    }
  }

  for(i=0;i<n_seq;i++) {
    seqobj = sequences + i;
    j = seqobj->length;

    cand = 0;
    while(j>=0) {
      if (seqobj->sequence[j--] == 'T') {
	na = 1; qual = seqobj->quality[j+1];
	while(j>=0 && seqobj->sequence[j] == 'T') {
	  na++;
	  qual += seqobj->quality[j--];
	}
	if (na > 11) {
	  k = j;
	  pasta_qual = 0;
	  while(k>=0 && (j-k)<na) {
	    pasta_qual += seqobj->quality[k--];
	  }
	  if (k<j && qual/na > (pasta_qual/(j-k))*1.5 && j<seqobj->length/3) {
	    if (cand >= cand_allocate) {
	      cand_allocate += 10;
	      candidates = local_realloc(candidates, 
					 sizeof(candidate_t)*cand_allocate,
					 "candidates for poly-T truncation");
	    }
	    candidates[cand].truncate_pos = j + 1;
	    candidates[cand].a_length = na;
	    candidates[cand].a_qual = qual/na;
	    candidates[cand].pasta_qual = pasta_qual/(j-k);
	    cand++;
	  }
	}
      }
    }

    if (cand > 0) {
      int select;
      select = 0;
      if (cand > 1) {
	for(j=0;j<cand;j++) {
	  if (candidates[j].a_length > candidates[select].a_length)
	    select = j;
	}
      }
      logmsg(MSG_DEBUG1,"Truncating leading pre-poly-T noise for sequence %s at "
	     "position %d\n", seq_names[i], candidates[select].truncate_pos);
      seqobj->length -= candidates[select].truncate_pos + 1;
      memmove(seqobj->sequence, 
	      seqobj->sequence + candidates[select].truncate_pos + 1,
	      seqobj->length);
      seqobj->sequence = local_realloc(seqobj->sequence, seqobj->length + 1,
				       "Poly-T truncation");
      memmove(seqobj->quality, 
	      seqobj->quality + candidates[select].truncate_pos,
	      seqobj->length);
      seqobj->quality = local_realloc(seqobj->quality, 
				      seqobj->length*(sizeof(uint)), 
				      "Poly-T Truncation");
    }
  }
  free(candidates);
}

static void generate_reverse_complement() {
  int i;

  logmsg(MSG_INFO,"Generating corresponding reverse complement sequences...\n");
  seq_names = local_realloc(seq_names, sizeof(uchar *)*n_seq*2,
			     "generate_reverse_complement()");
  sequences = local_realloc(sequences, sizeof(sequence_t)*n_seq*2,
			    "generate_reverse_complement()");
  for(i=0;i<n_seq;i++) {
    sequence_t *comp, *source;
    int j,k;

    source = sequences + i;
    comp = sequences + i + n_seq;

    comp->readname_index = source->readname_index + n_seq;
    comp->length = source->length;
    comp->sequence = local_calloc(source->length+1, sizeof(uchar), 
				  "generate_reverse_complement()");
    comp->quality = local_calloc(source->length, sizeof(int),
				 "generate_reverse_complement()");
    for(j=0;j<source->length;j++) {
      k = source->length - j - 1;
      comp->quality[j] = source->quality[k];
      switch(source->sequence[k]) {
      case 'A':
	comp->sequence[j] = 'T';
	break;
      case 'C':
	comp->sequence[j] = 'G';
	break;
      case 'G':
	comp->sequence[j] = 'C';
	break;
      case 'T':
	comp->sequence[j] = 'A';
	break;
      case 'N':
	comp->sequence[j] = 'N';
	break;
      default:
	logmsg(MSG_ERROR,"Unknown nucleotide '%c' in sequence %s at position "
	       "%d\n",source->sequence[k], seq_names[i], k);
	break;
      }
    }

    
    seq_names[i + n_seq] = malloc(strlen(seq_names[i])+2);
    strcpy(seq_names[i + n_seq], seq_names[i]);
    strcat(seq_names[i + n_seq], "-");
  }
  
  n_seq = n_seq*2;
  logmsg(MSG_INFO,"Finished generating complement sequences\n");
}


typedef struct {
  int s;
  int start;
} lookup_t;

#undef ALLOC_STEP
#define ALLOC_STEP (128)
#define WORDSIZE (9)
inline int dna_to_int(uchar *dna_string) {
  int i;
  int key;

  key = 0;
  for(i=0;i<WORDSIZE;i++) {
    switch(dna_string[i]) {
    case 'A':
      key = (key << 2);
      break;
    case 'C':
      key = (key << 2) + 1;
      break;
    case 'G':
      key = (key << 2) + 2;
      break;
    case 'T':
      key = (key << 2) + 3;
      break;
    case 'N':
    default:
      return -1;
    }
  }    

  return key;
}

static void build_subsequence_lookup_table(lookup_t ***rval_lookup_table,
					   int **rval_lookup_table_lengths) {
  int i,j;
  int key;
  lookup_t **lookup_table;
  int *lookup_table_lengths;

  lookup_table = local_calloc(1<<(WORDSIZE*2), sizeof(lookup_t *),
			      "subsequence lookup table");
  lookup_table_lengths = local_calloc(1<<(WORDSIZE*2), sizeof(int),
				      "subsequence lookup table");

  for(i=0;i<n_seq;i++) {
    for(j=0;j<sequences[i].length-WORDSIZE;j++) {
      key = dna_to_int(sequences[i].sequence + j);
      if (key != -1) {
	if (lookup_table_lengths[key] % ALLOC_STEP == 0) {
	  lookup_table[key] = 
	    local_realloc(lookup_table[key], sizeof(lookup_t)*
			  (lookup_table_lengths[key] + ALLOC_STEP),
			  "subsequence lookup table");
	}
	lookup_table[key][lookup_table_lengths[key]].s = i;
	lookup_table[key][lookup_table_lengths[key]].start = j;
	lookup_table_lengths[key]++;
      }
    }
  }
  for(i=0;i<(0x1 << WORDSIZE*2);i++) {
    if (lookup_table_lengths[i]) 
      lookup_table[i] = 
	realloc(lookup_table[i], sizeof(lookup_t)*(lookup_table_lengths[i]));
  }
  
  *rval_lookup_table = lookup_table;
  *rval_lookup_table_lengths = lookup_table_lengths;
}

typedef struct {
  int diagonal;
  int s1_start;
  int s2_start;
  int length;
  int score;
} match_t;

static void find_word_matches(match_t ***rval_matches, int **rval_n_matches,
			     lookup_t **lookup_table,
			     int *lookup_lengths, int seq_index) {
  int i,j;
  match_t **matches;
  int *n_matches;
  int key;
  sequence_t *seqobj;

  seqobj = sequences + seq_index;

  matches = local_calloc(n_seq, sizeof(match_t *), "find_word_matches");
  n_matches = local_calloc(n_seq, sizeof(int), "find_word_matches");

  for(i=0;i<seqobj->length-WORDSIZE;i++) {
    key = dna_to_int(seqobj->sequence + i);
    
    if (key == -1 || lookup_lengths[key] == 0) continue;

    for(j=0;j<lookup_lengths[key];j++) {
      match_t *m;
      int match_seqid;
      match_seqid = lookup_table[key][j].s;

      if (match_seqid == seq_index) continue;

      if (n_matches[match_seqid] % ALLOC_STEP == 0) {
	matches[match_seqid] = 
	  local_realloc(matches[match_seqid],
			sizeof(match_t)*(n_matches[match_seqid] + ALLOC_STEP),
			"find_word_matches");
      }
      
      m = matches[match_seqid] + n_matches[match_seqid]++;
      m->s1_start = i;
      m->s2_start = lookup_table[key][j].start;
      m->diagonal = m->s1_start - m->s2_start;
    }
  }

  *rval_matches = matches;
  *rval_n_matches = n_matches;
}

#define MISMATCH_WEIGHT (1.0) 
#define QC(x) (pow(10,(x)/-10.0))

static void banded_smith_waterman(align_t *al, sequence_t *a, sequence_t *b, 
				  int diagonal, int bandwidth) {
  float **mat;
  uchar **bt;
  int i,j,k, maxi, maxj, move, mini, minj, alignment_length;
  float s[3];
  float max;
  int p;
  int *query_pos, *sbjct_pos;
  uchar *query, *sbjct, *align;
  int match, mismatch, gaps, l_trim_length, r_trim_length;
  float mismatch_score, l_trim, r_trim;

  mat = local_calloc(a->length+1, sizeof(float *), "banded_smith_waterman()");
  bt = local_calloc(a->length+1, sizeof(uchar *), "banded_smith_waterman()");
  for(i=0;i<=a->length;i++) {
    mat[i] = local_calloc(b->length+1, sizeof(float), 
			  "banded_smith_waterman()");
    bt[i] = local_calloc(b->length+1, sizeof(uchar), 
			 "banded_smith_waterman()");
  }

  maxi = maxj = 0;
  for(i=1;i<a->length+1;i++) {
    for(j=i + diagonal - bandwidth; j< i + diagonal + bandwidth; j++) {
      if (j < 1) continue;
      if (j > b->length) break;
      if (a->sequence[i-1] == b->sequence[j-1])
	s[0] = mat[i-1][j-1] + 2.0;
      else
	s[0] = mat[i-1][j-1] - 5.00;
      s[1] = 0.0;
      for(k=1;k<=i;k++) {
	float val;
	val = mat[i-k][j] - (6.0 + 2*k);
	if (val > s[1]) s[1] = val;
      }
      s[2] = 0.0;
      for(k=1;k<=j;k++) {
	float val;
	val = mat[i][j-k] - (6.0 + 2*k);
	if (val > s[2]) s[2] = val;
      }
      max = 0.0;
      move = -1;
      for(k=0;k<=2;k++) {
	if (s[k] > max) {
	  max = s[k];
	  move = k;
	}
      }
      mat[i][j] = max;
      if (max > mat[maxi][maxj]) {
	maxi = i; maxj = j;
      }
      bt[i][j] = move;
    }
  }


  logmsg(MSG_DEBUG1,"maxi = %d maxj = %d length_a = %d length_b = %d\n",
	  maxi, maxj, a->length, b->length);
  i = maxi;
  j = maxj;
  mismatch = match = gaps = 0;
  mismatch_score = 0.0;
  
  alignment_length = 0;
  while(mat[i][j]>0.0) {
    int divisor;
    float s_qual, q_qual;
    alignment_length++;
    switch(bt[i][j]) {
    case 0:
      if (a->sequence[--i] == b->sequence[--j]) {
	match++;
      } else {
	mismatch++;
	mismatch_score += MISMATCH_WEIGHT*pow((1.0-QC(a->quality[i]))*(1.0-QC(b->quality[j])),0.5);
      }
      
      break;
    case 1:
    case 2:
      q_qual = 0;
      s_qual = 0;
      divisor = 0;
      for(k=-2;k<1;k++) {
	if (i+k<0 || i+k>=a->length) continue;
	q_qual += QC(a->quality[i+k]);
	divisor++;
      }
      q_qual /= divisor;

      divisor = 0;
      for(k=-2;k<1;k++) {
	if (j+k<0 || i+k>=b->length) continue;
	s_qual += QC(a->quality[j+k]);
	divisor++;
      }
      s_qual /= divisor;
      mismatch_score += MISMATCH_WEIGHT*sqrt((1.0-q_qual)*(1.0-s_qual));      
      gaps += 1;
      if (bt[i][j] == 1) i--;
      else j--;
      break;
    default:
      abort();
    }
  }
  mini = i; minj = j;

  l_trim = 0.0;
  l_trim_length = 0;
  while(i>0 && j>0) {
    logmsg(MSG_DEBUG4,"LTRIM: %c:%c %d:%d\n",a->sequence[i-1],b->sequence[j-1],a->quality[i-1],b->quality[j-1]);
    l_trim += pow(((1.0-QC(a->quality[--i]))*(1.0-QC(b->quality[--j]))),0.2);
    l_trim_length++;
  }
  i = maxi; j = maxj;
  r_trim = 0.0;
  r_trim_length = 0;
  while(i<a->length && j<b->length) {
    logmsg(MSG_DEBUG4,"RTRIM: %c:%c %d:%d\n",a->sequence[i],b->sequence[j],a->quality[i],b->quality[j]);
    r_trim += pow(((1.0-QC(a->quality[i++]))*(1.0-QC(b->quality[j++]))),0.2);
    r_trim_length++;
  }

  logmsg(MSG_DEBUG1,"Max score for %s vs %s: %f\n",seq_names[a->readname_index],
	  seq_names[b->readname_index], mat[maxi][maxj]);
  logmsg(MSG_DEBUG1,"Mismatches %d Matches %d Gaps %d Mismatch score %f\n",
	  mismatch, match, gaps, mismatch_score/mismatch);
  logmsg(MSG_DEBUG1,"l_trim %f l_trim length %d\n",l_trim, l_trim_length);
  logmsg(MSG_DEBUG1,"r_trim %f r_trim length %d\n",r_trim, r_trim_length);

  /* Skip this time consuming debuging output generation and display if
     the user has opt'd not to see it. This code generates a BLAST-like
     display of the alignment we just found */
  if (verbosity_level <= -2) {
    query = local_calloc(maxi+maxj+2, sizeof(uchar), "");
    sbjct = local_calloc(maxi+maxj+2, sizeof(uchar), "");
    align = local_calloc(maxi+maxj+2, sizeof(uchar), "");
    query_pos = local_calloc(maxi+maxj+2, sizeof(int), "");
    sbjct_pos = local_calloc(maxi+maxj+2, sizeof(int), "");
    p = maxi+maxj;
    i = maxi;
    j = maxj;
    while(mat[i][j]>0.0) {
      switch(bt[i][j]) {
      case 0:
	if (a->sequence[i-1] == b->sequence[j-1]) {
	  align[p] = '|';
	} else {
	  align[p] = ' ';
	}
	query[p] = a->sequence[--i];
	sbjct[p] = b->sequence[--j];
	break;
      case 1:
	query[p] = a->sequence[--i];
	sbjct[p] = '-';
	align[p] = ' ';
	break;
      case 2:
	query[p] = '-';
	sbjct[p] = b->sequence[--j];
	align[p] = ' ';
	break;
      default:
	abort();
      }
      p--;
      query_pos[p] = i;
      sbjct_pos[p] = j;
    }
    
    p += 1;
    while(p<(maxi+maxj)) {
      logmsg(MSG_DEBUG3,"Query %3.3d %-65.65s\n", query_pos[p], query + p);
      logmsg(MSG_DEBUG3,"          %-65.65s\n",align + p);
      logmsg(MSG_DEBUG3,"Sbjct %3.3d %-65.65s\n\n", sbjct_pos[p], sbjct + p);
      for(k=0;k<65;k++) {
	if (k + p > (maxi+maxj)) break;
	if (align[k + p] !='|') {
	  logmsg(MSG_DEBUG2,"Mismatch at position %d-%d %c:%c\t ",query_pos[p+k],
		 sbjct_pos[p+k], query[p+k], sbjct[p+k]);
	  logmsg(MSG_DEBUG2,"Quality values %d:%d\n",a->quality[query_pos[p+k]],
		 b->quality[sbjct_pos[p+k]]);
	  logmsg(MSG_DEBUG2,"Error probability: %4.3f:%4.3f\n",
		pow(10,(a->quality[query_pos[p+k]]/-10.0)),
		pow(10,(b->quality[sbjct_pos[p+k]]/-10.0)));
	}
      }
      p += 65;
    }
    free(query);
    free(sbjct);
    free(align);
    free(query_pos);
    free(sbjct_pos);
  }

  al->s1 = a->readname_index;
  al->s2 = b->readname_index;
  al->s1_ltrim = mini + 1;
  al->s2_ltrim = minj + 1;
  al->s1_rtrim = a->length - maxi;
  al->s2_rtrim = b->length - maxj;
  al->length = alignment_length;
  al->mismatch_score = mismatch_score;
  al->trim_score = l_trim + r_trim;
  al->align_score = mat[maxi][maxj];
  for(i=0;i<=a->length;i++) {
    free(mat[i]);
    free(bt[i]);
  }
  free(mat);
  free(bt);

}


static int match_compare(const void *c, const void *d) {
  const match_t *a, *b;
  int x;
  a = c;
  b = d;
  
  x = a->diagonal - b->diagonal;
  if (x != 0) return x;
  return (a->s1_start - b->s1_start);
}

static void add_overlap(int seq_index, align_t *al) {

  if (n_overlap[seq_index] % ALLOC_STEP == 0) {
    overlap[seq_index] = local_realloc(overlap[seq_index],sizeof(align_t)*
				 (n_overlap[seq_index] + ALLOC_STEP),
				 "add_overlap()");
  }
  memcpy(overlap[seq_index] + n_overlap[seq_index]++, al, sizeof(align_t));
}

static void add_nolist(uint s1, uint s2) {

  if (n_nolist[s1] % ALLOC_STEP == 0) {
    nolist[s1] = local_realloc(nolist[s1], 
			       sizeof(uint)*(n_nolist[s1] + ALLOC_STEP), 
			       "add_nolist()");
  }
  nolist[s1][n_nolist[s1]++] = s2;
}

static int match_compare_top(const void *a, const void *b) {
  return ((match_t *) b)->score - ((match_t *) a)->score;
}

static void combine_consecutive_matches(int *n_matches, match_t **match,
					int seq_index) {
  int i, j, k, p;

  for(i=0;i<n_seq;i++) {
    if (n_matches[i] == 0) continue;

    p = 0; j = 0;
    while(j<n_matches[i]) {
      k = j+1;
      while(k<n_matches[i] && (match[i][j].diagonal == match[i][k].diagonal) &&
	    ((match[i][k].s1_start - match[i][j].s1_start) <= (k-j)+WORDSIZE*2)) k++;
      if (k - j == 1 && p == j) {
	j = k;
	match[i][p].length = WORDSIZE;
	match[i][p].score = WORDSIZE*2;
	p++;
	continue;
      }

      match[i][p].diagonal = match[i][j].diagonal;
      match[i][p].s1_start = match[i][j].s1_start;
      match[i][p].s2_start = match[i][j].s2_start;
      match[i][p].length = match[i][k-1].s1_start - match[i][j].s1_start 
	+ WORDSIZE;
      match[i][p].score = 0;
      for(;j<k-1;j++) {
	if (match[i][j].s1_start == match[i][j+1].s1_start - 1) {
	  match[i][p].score += 2;
	} else {
	  match[i][p].score += WORDSIZE*2 - 5;
	}
      }
      match[i][p].score += WORDSIZE*2;
      p++;
      j = k;
    }
    /* Take the top ten matches
    if (p>20) {
      qsort(match[i], p, sizeof(match_t), match_compare_top);
      p = 20;
    } */
    n_matches[i] = p;
    match[i] = local_realloc(match[i], p*sizeof(match_t), 
			     "combine_consecutive_matches()");
  }
}

/*
  NOTES:

  After above function, should have consecutive matches combined, but any
  mismatch will break the consecutive hits so each diagonal may have several
  runs of consecutive hits. Need to below collect the best of these
  consecutive runs (or diagonals) and store for fast finishing of the
  alignment with gaps/mismatches by dynamic programming. */

typedef struct {
  int weight;
  int length;
  int *out_edges;
  int n_outedge;
  int s1_start;
  int s2_start;
  int diagonal;
} node_t;

typedef struct {
  int s_node;
  int e_node;
  int weight;
} edge_t;


fasta_t **fasta_scores;
int *n_fasta;

static int match_compare2(const void *a, const void *b) {
  
  return ((match_t *) a)->s1_start - ((match_t *) b)->s1_start;
}

#undef ALLOC_STEP
#define ALLOC_STEP (10)
static void build_graph(int *n_matches, match_t **match,
			int seq_index) {
  int i,j,k;
  int n_edges, n_nodes;
  edge_t *edges;
  node_t *nodes;
  int *pred;
  int *score;
  int max;
  int start, end;
  
  /* Resort by starting point rather than diagonal */
  for(i=0;i<n_seq;i++)
    if (n_matches[i])
      qsort(match[i], n_matches[i], sizeof(match_t), match_compare2);

  for(i=0;i<n_seq;i++) {

    if (n_matches[i] == 0) continue;
    /* Nodes graph should be initialized here !! */
    n_nodes = 0;
    nodes = local_calloc(ALLOC_STEP, sizeof(node_t), "build_graph()");
    nodes[0].weight = 0;
    nodes[0].out_edges = NULL;
    nodes[0].n_outedge = 0;
    nodes[0].s1_start = -1;
    nodes[0].s2_start = -1;
    nodes[0].diagonal = 0;
    n_nodes = 1;

    n_edges = 0;
    edges = NULL;

    for(j=0;j<n_matches[i];j++) {
      if (n_nodes % ALLOC_STEP == 0)
	nodes = local_realloc(nodes, sizeof(node_t)*(n_nodes + ALLOC_STEP),
			      "build_graph()");
      nodes[n_nodes].weight = match[i][j].score;
      nodes[n_nodes].length = match[i][j].length;
      nodes[n_nodes].out_edges = NULL;
      nodes[n_nodes].n_outedge = 0;
      nodes[n_nodes].s1_start = match[i][j].s1_start;
      nodes[n_nodes].s2_start = match[i][j].s2_start;
      nodes[n_nodes].diagonal = match[i][j].diagonal;
      for(k=0;k<n_nodes;k++) {
	if (nodes[k].s1_start < match[i][j].s1_start) {
	  if (n_edges % ALLOC_STEP == 0) 
	    edges = local_realloc(edges, sizeof(edge_t)*(ALLOC_STEP + n_edges),
				  "build_graph()");

	  edges[n_edges].s_node = k;
	  edges[n_edges].e_node = n_nodes;
	  /* Weight of edge is cost of connecting these nodes. It should be
	     the gap cost (difference in diagonals) plus an adjustment if the 
	     match runs overlap. If the edge is to the start node, set the
	     weight to the trim cost.  */
	  if (k == 0) {
	    if (nodes[n_nodes].s1_start < nodes[n_nodes].s2_start) {
	      edges[n_edges].weight = nodes[n_nodes].s1_start * 5;
	    } else {
	      edges[n_edges].weight = nodes[n_nodes].s2_start * 5;
	    }
	  } else {
	    edges[n_edges].weight = 
	      abs(nodes[k].diagonal-match[i][j].diagonal)*6;
	    if (match[i][j].s1_start < (nodes[k].s1_start + nodes[k].length)) {
	      edges[n_edges].weight += (nodes[k].s1_start + nodes[k].length 
		- match[i][k].s1_start)*2;
	    }
	  }

	  if (nodes[k].n_outedge % ALLOC_STEP == 0)
	    nodes[k].out_edges = 
	      local_realloc(nodes[k].out_edges, 
			    sizeof(int)*(nodes[k].n_outedge + ALLOC_STEP),
			    "build_graph()");
	  nodes[k].out_edges[nodes[k].n_outedge++] = n_edges;
	  n_edges ++;
	}
      }
      n_nodes++;
    }
    if (n_nodes % ALLOC_STEP == 0)
      nodes = local_realloc(nodes, sizeof(node_t)*(n_nodes + ALLOC_STEP),
			    "build_graph()");
    nodes[n_nodes].weight = 0;
    nodes[n_nodes].length = 0;
    nodes[n_nodes].out_edges = NULL;
    nodes[n_nodes].n_outedge = 0;
    nodes[n_nodes].s1_start = -1;
    nodes[n_nodes].s2_start = -1;
    nodes[n_nodes].diagonal = 0;
    for(k=1;k<n_nodes;k++) {
      int x, y;
      if (n_edges % ALLOC_STEP == 0) {
	edges = local_realloc(edges, sizeof(edge_t)*(ALLOC_STEP + n_edges),
			      "build_graph()");
      }
      edges[n_edges].s_node = k;
      edges[n_edges].e_node = n_nodes;
      x = nodes[k].s1_start + nodes[k].length;
      y = nodes[k].s2_start + nodes[k].length;
      if ( sequences[seq_index].length - x < sequences[i].length - y) {
	edges[n_edges].weight = (sequences[seq_index].length - x) * 5;
      } else {
	edges[n_edges].weight = (sequences[i].length - y) * 5;
      }
      
      if (nodes[k].n_outedge % ALLOC_STEP == 0)
	nodes[k].out_edges = 
	  local_realloc(nodes[k].out_edges, 
			sizeof(int)*(nodes[k].n_outedge + ALLOC_STEP),
			"build_graph()");
      nodes[k].out_edges[nodes[k].n_outedge++] = n_edges;
      n_edges ++;
    }
    n_nodes ++;

    /* Do something with graph here */
    logmsg(MSG_DEBUG4,"Graph has %d nodes and %d edges\n",n_nodes, n_edges);

    /* Single-source shortest path (sort of) */
    pred = local_calloc(n_nodes, sizeof(int), "build_graph()");
    score = local_calloc(n_nodes, sizeof(int), "build_graph()");
    for(k=0;k<n_nodes;k++) {
      pred[k] = -1;
      score[k] = INT_MIN;
    }
    score[0] = 0;
    for(k=0;k<n_nodes;k++) {
      int l;
      for(l=0;l<nodes[k].n_outedge;l++) {
	int edge, s;
	edge = nodes[k].out_edges[l];
	s = score[k] - edges[edge].weight + nodes[edges[edge].e_node].weight;
	if (s > score[edges[edge].e_node]) {
	  pred[edges[edge].e_node] = k;
	  score[edges[edge].e_node] = s + nodes[edges[edge].e_node].weight;
	}
      }
    }
    max = 0;
    for(k=1;k<n_nodes;k++) {
      if (score[k] > score[max]) max = k;
    }
    logmsg(MSG_DEBUG4,"Best path score = %d (%d)\n",score[max],score[n_nodes-1]);
    k = pred[n_nodes-1];
    end = nodes[k].s1_start + nodes[k].length;
    while(k!=0) {
      start = nodes[k].s1_start;
      k = pred[k];
    }

    if (score[n_nodes-1] > 200) {
      fasta_t *f;
      if (n_fasta[seq_index] % 10 == 0)
	fasta_scores[seq_index] = 
	  local_realloc(fasta_scores[seq_index], 
			sizeof(fasta_t)*(10+n_fasta[seq_index]), 
			"build_graph()");
      
      /* Record this score for later use */
      f = fasta_scores[seq_index] + n_fasta[seq_index];
      f->s1 = seq_index;
      f->s2 = i;
      f->score = score[n_nodes-1];
      f->start = start;
      f->end = end;
      n_fasta[seq_index]++;
    }

    free(pred);
    free(score);
    for(k=0;k<n_nodes;k++) {
      if (nodes[k].out_edges) free(nodes[k].out_edges);
    }
    free(nodes);
    if (edges) free(edges);
    nodes = NULL;
    edges = NULL;
    n_nodes = 0;
    n_edges = 0;
  }
}


#define SWAP(x,y) ((x)^=(y),(y)^=(x),(x)^=(y))
static void find_hits(lookup_t **lookup_table, int *lookup_lengths, 
		      int seq_index) {
  int *n_matches;
  int max_matches, favorite_diagonal;
  match_t **matches;
  int i,j,k;
  sequence_t *seqobj;
  align_t al;

  seqobj = sequences + seq_index;
  logmsg(MSG_DEBUG0,"Searching for word matches for sequence %s\n",
	  seq_names[seqobj->readname_index]);
  find_word_matches(&matches, &n_matches, lookup_table, lookup_lengths, 
		    seq_index);
  for(i=0;i<n_seq;i++) 
    if (n_matches[i])
      qsort(matches[i], n_matches[i], sizeof(match_t), match_compare);

  combine_consecutive_matches(n_matches, matches, seq_index);

  build_graph(n_matches, matches, seq_index);

#if 0
  for(i=0;i<n_matches;i++) {
    j = i;
    while(j<n_matches && matches[i].s2_sequence == matches[j].s2_sequence) j++;
    /* The comparision is symmetric, so sequences with matches to indices
       less than their own have already been run previously. */
    if (matches[i].s2_sequence > seq_index) {
      k = i;
      max_matches = 0;
      favorite_diagonal = -999;
      while(i<j) {
	while(k<j && (matches[i].s2_start - matches[i].s1_start) == 
	      (matches[k].s2_start - matches[k].s1_start)) k++;
	if ( k-i > max_matches) {
	  favorite_diagonal = matches[i].s2_start - matches[i].s1_start;
	  max_matches = k-i;
	}
	i = k;
      }
      logmsg(MSG_DEBUG1,"%s best diagonal with %s is %d (%d word matches)\n",
	     seq_names[seqobj->readname_index], 
	     seq_names[matches[i-1].s2_sequence], favorite_diagonal, 
	     max_matches);
      if (max_matches > 30) {
	banded_smith_waterman(&al, seqobj, sequences + matches[i-1].s2_sequence, favorite_diagonal, 20);
	if ((al.mismatch_score + al.trim_score)/ (float) al.length < 0.10) {
	  add_overlap(al.s1, &al);
	  SWAP(al.s1,al.s2);
	  SWAP(al.s1_ltrim, al.s2_ltrim);
	  SWAP(al.s1_rtrim, al.s2_rtrim);
	  add_overlap(al.s1, &al);
	} else {
	  logmsg(MSG_DEBUG1,"Discarding overlap between %s and %s because mismatches and\ntrimming differences are larger than 5%% of the alignment found.\n",
		 seq_names[al.s1],seq_names[al.s2]);
	  logmsg(MSG_DEBUG1,"Difference score: %f\n",
		 (al.mismatch_score + al.trim_score)/al.length);
	  add_nolist(al.s1,al.s2);
	  add_nolist(al.s2,al.s1);
	}
      }
    } else {
      i = j;
    }
  }
#endif
  for(i=0;i<n_seq;i++) {
    free(matches[i]);
  }
  free(matches);
  free(n_matches);
}

static void find_overlaps() {
  int i;
  lookup_t **lookup_table;
  int *lookup_table_lengths;

  fasta_scores = calloc(n_seq, sizeof(fasta_t *));
  n_fasta = calloc(n_seq, sizeof(int));
  build_subsequence_lookup_table(&lookup_table, &lookup_table_lengths);

  for(i=0;i<n_seq;i++) {      
    find_hits(lookup_table, lookup_table_lengths, i);
  }

}

static int alignment_compare(const void *c, const void *d) {
  const align_t *a;
  const align_t *b;

  a = c;
  b = d;

  if (a->align_score > b->align_score) {
    return -1;
  } else if (a->align_score < b->align_score) {
    return 1;
  } else {
    if (a->mismatch_score < b->mismatch_score) {
      return -1;
    } else if (a->mismatch_score > b->mismatch_score) {
      return 1;
    } else {
      return 0;
    }
  }
}

static void make_contigs() {
  int i,j,k,n_used,n_contigs;
  int *used;
  uint **contig;
  int *contig_size;
  int *contig_nolist;
  float max;

  used = local_calloc(n_seq, sizeof(int), "make_contigs()");
  n_used = 0;

  contig_nolist = local_calloc(n_seq, sizeof(int), "make_contigs()");

  for(i=0;i<n_seq;i++) {
    if (n_overlap[i]) 
      qsort(overlap[i], n_overlap[i], sizeof(align_t), alignment_compare);
  }

  contig = NULL;
  contig_size = NULL;
  n_contigs = 0;
  while(n_used < n_seq) {
    max = 0.0;
    j = -1;
    for(i=0;i<n_seq;i++) {
      if (! used[i] && n_overlap[i]>0 && overlap[i][0].align_score>max) {
	j = i;
	max = overlap[i][0].align_score;
      }
    }
    if (max == 0.0) break;
    if (n_contigs % ALLOC_STEP == 0) {
      contig = local_realloc(contig, sizeof(uint *)*(n_contigs + ALLOC_STEP),
			     "make_contigs()");
      contig_size = local_realloc(contig_size, 
				  sizeof(int)*(n_contigs + ALLOC_STEP),
				  "make_contigs()");
      for(i=n_contigs;i<n_contigs+ALLOC_STEP;i++) {
	contig[i] = NULL;
	contig_size[i] = 0;
      }
    }
    contig[n_contigs] = local_calloc(ALLOC_STEP, sizeof(uint), 
				     "make_contigs()");
    contig_size[n_contigs] = 1;
    contig[n_contigs][0] = j;
    used[j] = 1;
    memset(contig_nolist, 0x0, sizeof(uint)*n_seq);
    for(i=0;i<n_nolist[j];i++)
      contig_nolist[i] = 1;
    if (j > n_seq/2) {
      used[j - n_seq/2] = 1;
    } else {
      used[j + n_seq/2] = 1;
    }
    n_used += 2;
    j = 0;
    /* Take transitive closure */
    while(j<contig_size[n_contigs]) {
      int cc;

      cc = contig[n_contigs][j];
      for(i=0;i<n_overlap[cc];i++) {
	if (! used[overlap[cc][i].s2] && ! contig_nolist[overlap[cc][i].s2]) {
	  int dont_add;
	  
	  dont_add = 0;
	  fprintf(stderr,"Adding sequence %s because of overlap with sequence %s\n",seq_names[overlap[cc][i].s2],seq_names[cc]);
	  for(k=0;k<n_nolist[overlap[cc][i].s2];k++) {
	    int l;
	    for(l=0;l<contig_size[n_contigs];l++) {
	      if (nolist[overlap[cc][i].s2][k] == contig[n_contigs][l]) {
		fprintf(stderr,"Yomama: sequence %s is in nolist of sequence %s which is being added to contig!\n",seq_names[contig[n_contigs][l]], seq_names[overlap[cc][i].s2]);
		dont_add = 1;
	      }
	    }
	  }
	  
	  if (! dont_add) {
	    if (contig_size[n_contigs] % ALLOC_STEP == 0) {
	      contig[n_contigs] = local_realloc(contig[n_contigs], 
						(sizeof(uint)*(contig_size[n_contigs] + ALLOC_STEP)),
						"make_contigs()");
	    }
	    contig[n_contigs][contig_size[n_contigs]++] = overlap[cc][i].s2;
	    used[overlap[cc][i].s2] = 1;
	    for(k=0;k<n_nolist[overlap[cc][i].s2];k++)
	      contig_nolist[nolist[overlap[cc][i].s2][k]] = 1;
	    if (overlap[cc][i].s2 > n_seq/2) {
	      used[overlap[cc][i].s2 - n_seq/2] = 1;
	    } else {
	      used[overlap[cc][i].s2 + n_seq/2] = 1;
	    }
	    n_used += 2;
	  }
	}
      }
      j++;
    }
    n_contigs++;
  }
  logmsg(MSG_DEBUG1,"Found %d contigs\n",n_contigs);
  for(i=0;i<n_contigs;i++) {
    logmsg(MSG_DEBUG1,"Contig %d\n",i);
    for(j=0;j<contig_size[i];j++) {
      logmsg(MSG_DEBUG1,"\t%s\n",seq_names[contig[i][j]]);
    }
  }
  logmsg(MSG_DEBUG1,"Assembled %d out of %d\n",n_used/2,n_seq/2);
  
  /* Temporary code! -- output FASTA files of the members of each contig
     if output_basename is set */
  if (output_basename) {
    uchar *fname;
    FILE *f;

    fname = malloc(strlen(output_basename) + 256);
    for(i=0;i<n_contigs;i++) {
      sprintf(fname,"%s-contig-%d.seq",output_basename,i);
      f = openfile(fname, "w", "output contig components FASTA file");
      for(j=0;j<contig_size[i];j++) {
	fprintf(f,">%s\n",seq_names[contig[i][j]]);
	fprintf(f,"%s\n",sequences[contig[i][j]].sequence);
      }
    }
    free(fname);
  }
}

static void assess_fiveprime_threeprime() {
  int i, real_nseq, j, comp_matches;

  real_nseq = n_seq/2;
  for(i=0;i<real_nseq;i++) {
    comp_matches = 0;
    for(j=0;j<n_overlap[i];i++) {
      if (overlap[i][j].s2 >= real_nseq) {
	comp_matches++;
      }
    }
    if (comp_matches > n_overlap[i]/2) {
      fprintf(stderr,"%s matches more complement sequences (%d)\n",
	      seq_names[i],comp_matches);
    }
  } 
}

int main(int argc, char *argv[]) {
  uchar *input_seqfile, *input_qualfile;

  configure_logmsg(0);

  parse_arguments(&input_seqfile, &input_qualfile, argc, argv);
  configure_logmsg(verbosity_level);
  logmsg(MSG_DEBUG0,"Inputfile = %s\tQualfile = %s\n",input_seqfile, 
	  input_qualfile);

  load_inputsequence(input_seqfile, input_qualfile);

  polya_truncate(sequences, seq_names, n_seq);

  generate_reverse_complement();

  n_overlap = local_calloc(n_seq, sizeof(int), "main()");
  overlap = local_calloc(n_seq, sizeof(align_t *), "main()");

  n_nolist = local_calloc(n_seq, sizeof(int), "main()");
  nolist = local_calloc(n_seq, sizeof(uint *), "main()");

  find_overlaps();

  test_shit();

#if 0  
  make_contigs();

  assess_fiveprime_threeprime();
#endif

  return 0;
}
