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

void pairwise_prescan(void);
void polya_truncate(void);
void generate_reverse_complement(void);
void load_inputsequence(uchar *, uchar *);
void connected_components(void);
void spanning_tree(void);

seq_t *sequences;
int n_seq;

static uchar *output_basename = NULL;

static int verbosity_level = 0;
static int memsize = 256;

static void usage() {
  fprintf(stderr,"\nUsage:\n
--seqfile=<fasta input sequence file> (-s)
--qualfile=<fasta input phred quality file> (-q)
--memsize=<megabytes of RAM to use> (-m)
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
    { "memsize", 1, NULL, 'm' },
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
    case 'm':
      memsize = atoi(optarg);
      if (memsize <= 0) {
	logmsg(MSG_ERROR,"! Option memsize should be larger than 1 MB\n");
	commandline_error = 1;
      }
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


int main(int argc, char *argv[]) {
  uchar *input_seqfile, *input_qualfile;

  configure_logmsg(MSG_DEBUG1);

  parse_arguments(&input_seqfile, &input_qualfile, argc, argv);
  configure_logmsg(verbosity_level);
  logmsg(MSG_DEBUG0,"Inputfile = %s\tQualfile = %s\n",input_seqfile, 
	  input_qualfile);

  load_inputsequence(input_seqfile, input_qualfile);

  polya_truncate();

  generate_reverse_complement();

  pairwise_prescan();

  connected_components();

  spanning_tree();

  //test_shit();

  return 0;
}
