#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <getopt.h>
#include <math.h>

#include "kp_types.h"
#include "log_message.h"

static uchar *database_basename = NULL;
static uchar *output_basename = NULL;
static int verbosity_level = 0; 
static int mem_coresize = 192;
static uint wordsize = 9;
static uint n_seq = -1;
static int forward_only = 0;

static seqmeta_t *seqmeta = NULL;

static void usage(char *program_name) {

  fprintf(stderr,"\n\n%s:\n\n"
"Program to create lookup tables from formatted sequence database.             \n"
"This program does no real work, except for creating lookup table files for    \n"
"partitions of the sequence database. Assumed available memory should to       \n"
"hold individual lookup table resident should be specified on command line.    \n"
"                                                                              \n"
"Options:								       \n"
"--database=<database basename> (-d) (required)				       \n"
"    Preformatted binary database to be word-indexed into lookup tables	       \n"
"--basename=<string> (-o)                                                      \n"
"    Filename prefix for lookup tables. Lookup tables will be created as       \n"
"    <string>.N where N is an integer					       \n"
"--memsize=<integer> (-m)						       \n"
"    Assumed available core RAM size. Lookup tables will be made not much      \n"
"    larger than this size. Value is in megabytes (MB)			       \n"
"--verbose=<integer> (-v)						       \n"
"    Verbosity level. 0 (normal) by default. Negative enables debugging messages\n"
"    Positive makes program quieter.                                           \n"
"--forward-only (-f)							       \n"
"    Skip every other sequence from the input. Useful to prevent reverse       \n"
"    complement data in the lookup table. Note that usually this utility is    \n"
"    used with a preformatted sequence database file which automatically       \n"
"    puts the reverse complement of the sequence after each input sequence.    \n"
"    It is not necessary to compare a reverse complement with another 	       \n"
"    reverse complement, that is the same as forward vs. forward.	       \n"
"--help (-h)                                                                   \n"
"    Prints this message.						       \n"
"									       \n"
,program_name);							       

}


static void parse_arguments(int argc, char *argv[]) {
  int option_index, commandline_error, rval;
  struct option longopts[] = {
    { "database", 1, NULL, 'd'},
    { "basename", 1, NULL, 'o'},
    { "memsize", 1, NULL, 'm'},
    { "verbose", 1, NULL, 'v'},
    { "forward-only", 0, NULL, 'f'},
    { "help", 0, NULL, 'h'},
    { NULL, 0, NULL, 0}
  };
  char *optstring = "d:v:o:m:hf";

  database_basename = output_basename = NULL;
  commandline_error = 0;
  while((rval = getopt_long(argc, argv, optstring, longopts, &option_index))
	!= -1) {
    switch(rval) {
    case ':':
      logmsg(MSG_ERROR,"\n! Option \"%s\" requires an argument.\n",
	     longopts[option_index].name);
      commandline_error = 1;
      break;
    case 'd':
      database_basename = strdup(optarg);
      break;
    case 'v':
      verbosity_level = atoi(optarg);
      break;
    case 'o':
      output_basename = strdup(optarg);
      break;
    case 'm':
      mem_coresize = atoi(optarg);
      break;
    case 'f':
      forward_only = 1;
      break;
    case 'h':
      usage(argv[0]);
      exit(0);
      break;
    case '?':
    default:
      logmsg(MSG_ERROR,"\n! Option \"%c\" unknown.\n",optopt);
      commandline_error = 1;
      break;
    }
  }

  if (database_basename == NULL) {
    logmsg(MSG_ERROR,"! Formatted sequence database basename must be "
	   "specified with -d <basename> or --database=<basename> option\n");
    commandline_error = 1;
  }

  if (mem_coresize <= 0) {
    logmsg(MSG_ERROR,"! Specified RAM size assumption must be larger than "
	   "0\n");
    commandline_error = 1;
  }
  
  if (commandline_error) {
    logmsg(MSG_ERROR,"! Program halted due to command line option errors\n");
    usage(argv[0]);
    exit(-1);
  }

  if (output_basename == NULL) {
    output_basename = strdup(database_basename);
  }
}

static void open_databasefiles(FILE **indfile, FILE **binfile) {
  int l;
  uchar *temp;
  uint x;
  FILE *f;

  l = strlen(database_basename) + 6;
  MA(temp, l);
  strcpy(temp, database_basename);
  strcat(temp, ".ind");
  f = fopen(temp, "r");
  if (f == NULL) {
    logmsg(MSG_FATAL,"! Failed opening database index file %s (%s)\n",
	   temp, strerror(errno));
  }
  fread(&x, sizeof(uint), 1, f);
  if (x != INDFILE_MAGIC) {
    logmsg(MSG_FATAL,"! Database index file does not appear to be properly formatted\n");
  }
  fread(&n_seq, sizeof(uint), 1, f);
  MA(seqmeta, sizeof(seqmeta_t)*n_seq);
  fread(seqmeta, sizeof(seqmeta_t), n_seq, f);
  *indfile = f;
  
  strcpy(temp, database_basename);
  strcat(temp, ".sbin");
  f = fopen(temp, "r");
  if (f == NULL) {
    logmsg(MSG_FATAL,"! Failed opening database binary file %s (%s)\n",
	   temp, strerror(errno));
  }
  fread(&x, sizeof(uint), 1, f);
  if (x != BINFILE_MAGIC) {
    logmsg(MSG_FATAL,"! Database binary file does not appear to be properly formatted\n");
  }
  *binfile = f;

  free(temp);
}

static void write_header(FILE *f, uint start, uint stop, uint table_size, int tn) {
  uint x;
  
  x = LOOKUP_MAGIC;
  fwrite(&x, sizeof(uint), 1, f);

  fwrite(&wordsize, sizeof(uint), 1, f);
  fwrite(&start, sizeof(uint), 1, f);
  fwrite(&stop, sizeof(uint), 1, f);
  fwrite(&tn, sizeof(int), 1, f);
  fwrite(&table_size, sizeof(uint), 1, f);

}

static uint build_lookuptable(lookupmeta_t *lookup_meta, word_t **ld,
			      uint *total_words, uint start_seq, 
			      FILE *binfile) {
  uint word, mask;
  uint limit;
  int length, j, total, seqsize;
  uchar *seq;
  uint seq_id, end_seq;
  int *fill;
  word_t *lookup_data;
  double p, var, expect;
  uint censored;
  
  limit = (mem_coresize*1024*1024)/sizeof(word_t);

  mask = (0x1 << wordsize*2) - 1;

  seq = NULL;
  seqsize = 0;
  total = 0;
  seq_id = start_seq;
  fseek(binfile, seqmeta[start_seq].seqbin_pos, SEEK_SET);
  while(seq_id<n_seq && total < limit) {
    length = seqmeta[seq_id].seq_length;
    if (seqsize < length) {
      seqsize = length;
      RA(seq, seqsize, sizeof(uchar));
    }

    fread(seq, sizeof(uchar), length, binfile);
    /* Cheap hack to prevent cataloging of reverse complement sequences which
       ought to be odd numbered sequence ids */
    if (forward_only && (seq_id & 0x1)) {
      seq_id++;
      continue;
    }
    word = 0;
    for(j=0;j<wordsize;j++) {
      word = (word << 2) | seq[j];
    }
    lookup_meta[word].n_words++;
    total++;
    for(;j<length;j++) {
      word = ((word << 2) & mask) | seq[j];
      lookup_meta[word].n_words++;
      total++;
    }

    seq_id++;
  }

  p = 1.0/(double) mask;
  expect = p*total;
  var = p*(1.0-p);
  censored = 0;
  for(word=0;word<=mask;word++) {
    if (lookup_meta[word].n_words > expect*50) {
      fprintf(stderr,"Censoring word: %0X (%d obs out of %d total, expect = %5.2f)\n",word,
	      lookup_meta[word].n_words, total, expect);
      censored += lookup_meta[word].n_words;
      lookup_meta[word].n_words = 0;
    }
  }
  //  total -= censored;

  MA(fill, sizeof(int)*(mask+1));
  MA(lookup_data, total*sizeof(word_t));
  fill[0] = 0;
  lookup_meta[0].start_pos = 5*sizeof(uint) + sizeof(lookupmeta_t)*(mask+1);
  for(word=1;word<=mask;word++) {
    fill[word] = fill[word-1] + lookup_meta[word-1].n_words;
    lookup_meta[word].start_pos = lookup_meta[word-1].start_pos + 
      lookup_meta[word-1].n_words*sizeof(word_t);
  }

  end_seq = seq_id;
  seq_id = start_seq;
  fseek(binfile, seqmeta[start_seq].seqbin_pos, SEEK_SET);
  while(seq_id < end_seq) {
    length = seqmeta[seq_id].seq_length;
    fread(seq, sizeof(uchar), length, binfile);
    /* Cheap hack to prevent cataloging of reverse complement sequences which
       ought to be odd numbered sequence ids */
    if (forward_only && (seq_id & 0x1)) {
      seq_id++;
      continue;
    }
    word = 0;
    for(j=0;j<wordsize;j++) {
      word = (word << 2) | seq[j];
    }
    lookup_data[fill[word]].seq_id = seq_id;
    lookup_data[fill[word]].seq_pos = 0;
    fill[word]++;
    for(;j<length;j++) {
      word = ((word << 2) & mask) | seq[j];
      lookup_data[fill[word]].seq_id = seq_id;
      lookup_data[fill[word]].seq_pos = j - wordsize;
      fill[word]++;
    }

    seq_id++;
  }

  *ld = lookup_data;
  *total_words = total;
  return (end_seq - start_seq - 1);
}

static void create_lookup_tables(void) {
  uint n_words;
  int i,j, table_number, l;
  uchar *lookup_filename;
  uint n, total;
  lookupmeta_t *lookup_meta;
  word_t *lookup_data;
  FILE *lf;
  FILE *indfile, *binfile;

  /* n_seq is read out of index file header */
  open_databasefiles(&indfile, &binfile);

  n_words = 0x1 << (wordsize*2);

  MA(lookup_meta, sizeof(lookupmeta_t)*n_words);
  MA(lookup_data, sizeof(word_t *)*n_words);
  l = strlen(output_basename) + 32;
  MA(lookup_filename, l);

  i = 0;
  table_number = 0;
  while(i < n_seq) {
    for(j=0;j<n_words;j++) {
      lookup_meta[j].n_words = 0;
      lookup_meta[j].start_pos = 0;
    }
    sprintf(lookup_filename,"%s.lt.%d",output_basename,table_number);
    n = build_lookuptable(lookup_meta, &lookup_data, &total, i, binfile);
    lf = fopen(lookup_filename, "w");
    if (lf == NULL) {
      logmsg(MSG_FATAL,"! Failed opening output file %s (%s)\n",
	     lookup_filename, strerror(errno));
    }
    logmsg(MSG_INFO,"Writing lookup table %d spanning sequences %u - %u\n",
	   table_number, i, i + n);
    write_header(lf, i, i + n, total, table_number);
    fwrite(lookup_meta, sizeof(lookupmeta_t), n_words, lf);
    fwrite(lookup_data, sizeof(word_t), total, lf);
    free(lookup_data);
    fclose(lf);
    i += n + 1;
    table_number++;
  }

  free(lookup_meta);
  free(lookup_filename);
}

int main(int argc, char *argv[]) {


  configure_logmsg(MSG_DEBUG1);
  parse_arguments(argc, argv);
  configure_logmsg(verbosity_level);

  logmsg(MSG_INFO,"Output basename set to %s\n",output_basename);

  create_lookup_tables();

  return 0;
}
