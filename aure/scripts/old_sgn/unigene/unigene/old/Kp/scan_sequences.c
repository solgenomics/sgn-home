#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <getopt.h>
#include <limits.h>

#include "kp_types.h"
#include "log_message.h"

#define SCORE_THRESHOLD (75)

static uchar *lookup_filename = NULL;
static uchar *seq_filename = NULL;
static uint verbosity_level = 0;
static uint wordsize;
static uint mask;

static lookupmeta_t *lookup_meta;
static word_t **lookup_data;
static word_t *lookup;

static uint n_seq = -1;
static seqmeta_t *seqmeta = NULL;
static uint ltable_start, ltable_end;


typedef struct {
  int db_seq;
  int di;
  int pos;
  int length;
} wordhit_t;

inline int mergesort_compare(wordhit_t *a, wordhit_t *b) {
  if (a->db_seq < b->db_seq) return -1;
  else if (a->db_seq > b->db_seq) return 1;
  else {
    if (a->di < b->di) return -1;
    else if (a->di > b->di) return 1;
    else {
      if (a->pos < b->pos) return -1;
      else if (a->pos > b->pos) return 1;
      else return 0;
    }
  }
}


#define SWAP(x,y) ((x)^=(y),(y)^=(x),(x)^=(y))
static void wordhit_mergesort_actual(wordhit_t *hits, int s, int e, 
				     wordhit_t *tspace) {
  int i,j,half,f,x;

  if (e - s > 7) {
    half = s + (e-s)/2;
    wordhit_mergesort_actual(hits, s, half, tspace);
    wordhit_mergesort_actual(hits, half, e, tspace);
    i = s;
    j = half;
    f = 0;
    while(i < half && j < e) {
      x = mergesort_compare(hits + i, hits + j);
      if (x == -1) {
	tspace[f++] = hits[i++];
      } else {
	tspace[f++] = hits[j++];
      }
    }
    assert((i < half || j < e) && !(i<half && j<e));
    if (j < e) {
      memcpy(hits + s, tspace, (j-s)*sizeof(wordhit_t));
    } else {
      memcpy(tspace + f, hits + i, (half-i)*sizeof(wordhit_t));
      memcpy(hits + s, tspace, (e-s)*sizeof(wordhit_t));
    }
  } else {
    for(i=s;i<e;i++) {
      int min;
      min = i;
      for(j=i+1;j<e;j++)
	if (mergesort_compare(hits + min, hits + j) == 1) min = j;
      if (min != i) {
	wordhit_t t;
	t = hits[min];
	hits[min] = hits[i];
	hits[i] = t;
      }
    }
  }
}

static void mergesort_bypos_actual(wordhit_t *hits, int s, int e, 
				   wordhit_t *tspace) {
  int i,j,half,f;

  if (e - s > 7) {
    half = s + (e-s)/2;
    wordhit_mergesort_actual(hits, s, half, tspace);
    wordhit_mergesort_actual(hits, half, e, tspace);
    i = s;
    j = half;
    f = 0;
    while(i < half && j < e) {
      if (hits[i].pos < hits[j].pos) {
	tspace[f++] = hits[i++];
      } else {
	tspace[f++] = hits[j++];
      }
    }
    assert((i < half || j < e) && !(i<half && j<e));
    if (j < e) {
      memcpy(hits + s, tspace, (j-s)*sizeof(wordhit_t));
    } else {
      memcpy(tspace + f, hits + i, (half-i)*sizeof(wordhit_t));
      memcpy(hits + s, tspace, (e-s)*sizeof(wordhit_t));
    }
  } else {
    for(i=s;i<e;i++) {
      int min;
      min = i;
      for(j=i+1;j<e;j++)
	if (hits[j].pos < hits[min].pos) min = j;
      if (min != i) {
	wordhit_t t;
	t = hits[min];
	hits[min] = hits[i];
	hits[i] = t;
      }
    }
  }
}

static void wordhit_mergesort(wordhit_t *hits, int s, int e) {
  wordhit_t *tspace;
  
  MA(tspace, (e-s)*sizeof(wordhit_t));
  wordhit_mergesort_actual(hits, s, e, tspace);
  free(tspace);
}

static void mergesort_bypos(wordhit_t *hits, int s, int e) {
  wordhit_t *tspace;

  MA(tspace, (e-s)*sizeof(wordhit_t));
  mergesort_bypos_actual(hits, s, e, tspace);
  free(tspace);

}


static int wordhit_compare(const void *c, const void *d) {
  wordhit_t *a, *b;

  a = (wordhit_t *) c;
  b = (wordhit_t *) d;
  
  if (a->db_seq < b->db_seq) return -1;
  else if (a->db_seq > b->db_seq) return 1;
  else {
    if (a->di < b->di) return -1;
    else if (a->di > b->di) return 1;
    else {
      if (a->pos < b->pos) return -1;
      else if (a->pos > b->pos) return 1;
      else return 0;
    }
  }
}

typedef struct {
  int db_seq;
  int min_di;
  int max_di;
  int score;
  int start;
  int end;
  int s_start;
  int s_end;
  int length;
} hit_report_t;

/* This is allocated dynamically by main() after the lookup table is loaded,
   for efficiency, since it is always needed by the below function, which is
   called for every sequence comparison, but the size of the array is not
   known until the lookup table is loaded */
static int *hits_byseq = NULL;

static wordhit_t *find_wordmatches(uchar *seq, uint seq_id, int length, int *return_nhits) {
  int n_hits;
  int i,j,t;
  uint word;
  wordhit_t *hits;
  
  /* This implements a censoring technique to speed the execution of the 
     program, by excluding spurious word matches to sequences (which would
     trigger a single-source shortest path analysis) that can not in total
     reach the output threshold. Also, by excluding these words from the
     list of word hits, the sorting time for combining the word hits is also
     reduced. */
  for(j=0;j<(ltable_end - ltable_start);j++)
    hits_byseq[j] = 0;
  
  /* Count the word hits */
  word = 0;
  for(i=0;i<wordsize;i++) 
    word = (word << 2) + seq[i];
  for(j=0;j<lookup_meta[word].n_words;j++) {
    hits_byseq[lookup_data[word][j].seq_id - ltable_start]++;
  }
  for(;i<length;i++) {
    word = ((word << 2) & mask) + seq[i];
    for(j=0;j<lookup_meta[word].n_words;j++) {
      hits_byseq[lookup_data[word][j].seq_id - ltable_start]++;
    }
  }

  n_hits = 0;
  for(j=0;j<(ltable_end - ltable_start);j++) {
    if ((j+ltable_start)>=seq_id && hits_byseq[j]*2 >= SCORE_THRESHOLD) {
      n_hits += hits_byseq[j];
    } else {
      hits_byseq[j] = 0;
    }
  }
  if (n_hits == 0) {
    *return_nhits = 0;
    return NULL;
  }
  MA(hits, n_hits*sizeof(wordhit_t));
  
  /* Allocate the memory needed and then compile the records for each word */
  t = 0;
  word = 0;
  for(i=0;i<wordsize;i++) 
    word = (word << 2) + seq[i];
  for(j=0;j<lookup_meta[word].n_words;j++) {
    if (hits_byseq[lookup_data[word][j].seq_id - ltable_start] > 0) {
      hits[t].db_seq = lookup_data[word][j].seq_id;
      hits[t].di = lookup_data[word][j].seq_pos - (i - wordsize);
      hits[t].pos = (i - wordsize);
      t++;
    }
  }
  for(;i<length;i++) {
    word = ((word << 2) & mask) + seq[i];
    for(j=0;j<lookup_meta[word].n_words;j++) {
      if (hits_byseq[lookup_data[word][j].seq_id - ltable_start] > 0) {
	hits[t].db_seq = lookup_data[word][j].seq_id;
	hits[t].di = lookup_data[word][j].seq_pos - (i - wordsize);
	hits[t].pos = (i - wordsize);
	t++;
      }
    }
  }
  assert(t == n_hits);

  *return_nhits = n_hits;
  return hits;
}

static int combine_hits(wordhit_t *hits, int n_hits) {
  int i,j,f;

  f = i = j = 0;
  while(i < n_hits) {
    while(j < n_hits                       && 
	  hits[i].db_seq == hits[j].db_seq &&
	  hits[i].di == hits[j].di         && 
	  hits[j].pos - hits[i].pos == j - i) j++;
    hits[f].di = hits[i].di;
    hits[f].db_seq = hits[i].db_seq;
    hits[f].pos = hits[i].pos;
    hits[f].length = j - i + wordsize - 1;
    f++;
    i = j;
  }


  return f;
}

static int single_source_shortest_path(int *adjmatrix, int *pred, int *score,
					wordhit_t *hits, int i, int j) {
  int n_nodes;
  int k, l, s;

  /* Single source shortest path algorithm */
  /* Selection sort on position */
  if (j-i > 10) {
    mergesort_bypos(hits, i, j);
  } else {
    for(k=i;k<j;k++) {
      int min;
      min = k;
      for(l=k+1;l<j;l++)
	if (hits[l].pos < hits[min].pos) min = l;
      if (min != k) {
	SWAP(hits[min].db_seq,hits[k].db_seq);
	SWAP(hits[min].di,hits[k].di);
	SWAP(hits[min].pos,hits[k].pos);
	SWAP(hits[min].length,hits[k].length);
      }
    }
  }

  /* Adj matrix initialization */
  n_nodes = j - i + 2;
  memset(adjmatrix, 0x0, sizeof(int)*n_nodes*n_nodes);
  for(k=i;k<j;k++) {
    for(l=i;l<j;l++) {
      if (hits[k].pos < hits[l].pos) {
	adjmatrix[(k-i+1)*n_nodes + (l-i+1)] = abs(hits[k].di - hits[l].di)
	  + abs(hits[k].pos + hits[k].length - hits[l].pos) + 1;
      }	    
    }
  }
  /* Boundary conditions: source has edge to every node (except sink node)
     Every node (except source) has edge to sink node */
  for(l=1;l<n_nodes-1;l++) {
    adjmatrix[l] = 1;
    adjmatrix[l*n_nodes + (n_nodes-1)] = 1;
  }
  
  /* State initialization for single-source shortest path */
  for(k=0;k<n_nodes;k++) {
    pred[k] = -1;
    score[k] = INT_MIN;
  }
  pred[0] = -1;
  score[0] = 0;
  
  /* SSSP proper */
  for(k=0;k<n_nodes;k++) {
    for(l=0;l<n_nodes;l++) {
      if (adjmatrix[k*n_nodes + l]) {
	if (l>0 && l<n_nodes-1) {
	  s = score[k] - adjmatrix[k*n_nodes + l] + hits[l + i - 1].length;
	} else {
	  s = score[k];
	}
	if (s > score[l]) {
	  pred[l] = k;
	  score[l] = s;
	}
      }
    }
  }

  return n_nodes;
}

static int fasta_scan(uchar *seq, uint seq_id, int length, 
		      hit_report_t *report_hits) {
  int i, j, k, f;
  int n_hits, n_nodes, max, max_span;
  int min_di, max_di, total_length;
  wordhit_t *hits;
  int *adjmatrix;
  int *pred, *score;
  int end, start, s_start, s_end;

  
  hits = find_wordmatches(seq, seq_id, length, &n_hits);
  if (hits == NULL) return 0;

  wordhit_mergesort(hits, 0, n_hits);

#if 0
  for(i=0;i<n_hits-1;i++) {
    if (mergesort_compare(hits + i, hits + i + 1) == 1) {
      fprintf(stderr,"Out of order!\n");
      fprintf(stderr,"%d %d %d %d\n",hits[i].db_seq,hits[i].di, hits[i].pos,
	      hits[i].length);
      fprintf(stderr,"%d %d %d %d\n",hits[i+1].db_seq,hits[i+1].di, 
	      hits[i+1].pos, hits[i+1].length);

    }
  }
#endif
  //qsort(hits, n_hits, sizeof(wordhit_t), wordhit_compare);

  f = combine_hits(hits, n_hits);

  max_span = i = j = 0;
  while(i < f) {
    while(j < f && hits[i].db_seq == hits[j].db_seq) j++;
    if ((j - i) > max_span)
      max_span = j - i;
    i = j;
  }

  MA(adjmatrix, (max_span + 3)*(max_span + 3)*sizeof(int));
  MA(pred, sizeof(int)*(max_span + 3));
  MA(score, sizeof(int)*(max_span + 3));

  n_hits = 0;
  i = j = 0;
  while(i < f) {
    while(j < f && hits[j].db_seq == hits[i].db_seq) j++;

    n_nodes = single_source_shortest_path(adjmatrix, pred, score, hits, i, j);

    /* Recovery of results from SSSP */
    max = 0;
    for(k=1;k<n_nodes;k++) {
      if (score[k] > score[max]) max = k;
    }
    k = max;
    min_di = hits[i+k-1].di;
    max_di = hits[i+k-1].di;
    start = end = hits[i+k-1].pos + hits[i+k-1].length;
    s_start = s_end = hits[i+k-1].pos + hits[i+k-1].di + hits[i+k-1].length;
    while(k!=0) {
      total_length += hits[i+k-1].length;
      if (hits[i+k-1].di < min_di) {
	min_di = hits[i+k-1].di;
      } else if (hits[i+k-1].di > max_di) {
	max_di = hits[i+k-1].di;
      }
      if (pred[k] == 0) {
	start = hits[i+k-1].pos;
	s_start = hits[i+k-1].pos + hits[i+k-1].di;
      }
      k = pred[k];
    }

    /* Arbitrary selection of results to report */
    if (score[max] >= SCORE_THRESHOLD) {
      report_hits[n_hits].db_seq = hits[i].db_seq;
      report_hits[n_hits].min_di = min_di - 5;
      report_hits[n_hits].max_di = max_di + 5;
      report_hits[n_hits].score = score[max];
      report_hits[n_hits].start = start;
      report_hits[n_hits].end = end;
      report_hits[n_hits].s_start = s_start;
      report_hits[n_hits].s_end = s_end;
      n_hits++;
    }

    i = j;
  }
  free(adjmatrix);
  free(pred);
  free(score);
  free(hits);
  return n_hits;
}

static void usage(char *program_name) {

  fprintf(stderr,"\n\n%s:\n\n"
"Quick program to scan formatted sequence file against a pre-formatted \n"
"database of words (sub-sequence), to approximate alignment by linking \n"
"together consecutive sequences of matching words.\n"
"\n"
"Options:\n"
"--seqfile=<basename> (-s) (required)\n"
"    Basename of preformatted sequence 'database'\n"
"--lookupfile=<lookup file> (-l) (required)\n"
"    Preformatted lookup table\n"
"--verbose=<integer> (-v)\n"
"    Verbosity level. 0 (normal) by default. Negative enables debugging messages\n"
"    Positive makes program quieter.\n"
"--help (-h)\n"
"    Prints this message.\n"
,program_name);

}


static void parse_arguments(int argc, char *argv[]) {
  int option_index, commandline_error, rval;
  struct option longopts[] = {
    { "seqfile", 1, NULL, 's'},
    { "lookupfile", 1, NULL, 'l'},
    { "verbose", 1, NULL, 'v'},
    { "help", 1, NULL, 'h'},
    { NULL, 0, NULL, 0}
  };
  char *optstring = "s:l:v:h";

  commandline_error = 0;
  while((rval = getopt_long(argc, argv, optstring, longopts, &option_index))
	!= -1) {
    switch(rval) {
    case ':':
      logmsg(MSG_ERROR,"\n! Option \"%s\" requires an argument.\n",
	     longopts[option_index].name);
      commandline_error = 1;
      break;
    case 'l':
      lookup_filename = strdup(optarg);
      break;
    case 'v':
      verbosity_level = atoi(optarg);
      break;
    case 's':
      seq_filename = strdup(optarg);
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

  if (lookup_filename == NULL) {
    logmsg(MSG_ERROR,"! Formatted lookup file must be "
	   "specified with -l <lookup file> or --seqfile=<lookup file> "
	   "option\n");
    commandline_error = 1;
  }

  if (seq_filename == NULL) {
    logmsg(MSG_ERROR,"! Formatted sequence database basename must be "
	   "specified with -s <basename> or --seqfile=<basename> "
	   "option\n");
    commandline_error = 1;
  }
  
  if (commandline_error) {
    logmsg(MSG_ERROR,"! Program halted due to command line option errors\n");
    usage(argv[0]);
    exit(-1);
  }

}

static void open_databasefiles(FILE **indfile, FILE **binfile) {
  int l;
  uchar *temp;
  uint x;
  FILE *f;

  l = strlen(seq_filename) + 6;
  MA(temp, l);
  strcpy(temp, seq_filename);
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
  
  strcpy(temp, seq_filename);
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

static void open_lookupfile(FILE **lookupfile) {
  FILE *lf;
  uint x;
  int table_index;
  uint table_size;
  uint i;

  lf = fopen(lookup_filename, "r");
  if (lf == NULL) {
    logmsg(MSG_FATAL,"! Failed opening lookup file %s (%s)\n",
	   lookup_filename, strerror(errno));
  }
  
  fread(&x, sizeof(uint), 1, lf);
  if (x != LOOKUP_MAGIC) {
    logmsg(MSG_FATAL,"! Lookup file does not appear to be properly formatted\n");
  }

  fread(&wordsize, sizeof(uint), 1, lf);
  if (wordsize < 2 || wordsize > 24) {
    logmsg(MSG_FATAL,"! Lookup file does not appear to be properly formatted\n");
  }
  mask = (0x1 << (wordsize*2)) - 1;

  fread(&ltable_start, sizeof(uint), 1, lf);
  fread(&ltable_end, sizeof(uint), 1, lf);
  fread(&table_index, sizeof(int), 1, lf);

  logmsg(MSG_INFO,"Loading lookup table file %d: covering sequences "
	 "%lu - %lu\n",table_index, ltable_start, ltable_end);
  ltable_end += 1;

  fread(&table_size, sizeof(uint), 1, lf);
  MA(lookup_meta, sizeof(lookupmeta_t)*(mask+1));
  MA(lookup_data, sizeof(word_t *)*mask+1);
  MA(lookup, sizeof(word_t)*table_size);
  
  fread(lookup_meta, sizeof(lookupmeta_t), mask+1, lf);
  fread(lookup, sizeof(word_t), table_size, lf);
  
  lookup_data[0] = lookup;
  for(i=1;i<=mask;i++) {
    lookup_data[i] = lookup_data[i-1] + lookup_meta[i-1].n_words;
  }

  *lookupfile = lf;
}

static void reverse_complement(uchar *seq, int length) {
  int i,j;
  uchar *comp;

  comp = alloca(sizeof(uchar)*length);
  for(i=0;i<length;i++) {
    j = length - i - 1;
    comp[i] = seq[j] ^ 0x3;
  }
  memcpy(seq, comp, length);
}

#define MIN(x,y) ((x)<(y)?(x):(y))
int main(int argc, char *argv[]) {
  FILE *indfile, *binfile;
  FILE *lookupfile;
  hit_report_t *report_hits;
  uint i, j, n_hits;
  int seqsize, length;
  uchar *seq;

  configure_logmsg(MSG_DEBUG1);
  parse_arguments(argc, argv);
  configure_logmsg(verbosity_level);

  logmsg(MSG_INFO,"Input database basename set to %s\n",seq_filename);
  open_databasefiles(&indfile, &binfile);
  open_lookupfile(&lookupfile);
  MA(hits_byseq, sizeof(int)*ltable_end);

  MA(report_hits, sizeof(hit_report_t)*ltable_end);
  seq = NULL;
  seqsize = 0;
  for(i=0;i<n_seq;i++) {
    length = seqmeta[i].seq_length;
    if (length > seqsize) {
      seqsize = length;
      RA(seq, seqsize, sizeof(uchar));
    }

    fread(seq, sizeof(uchar), length, binfile);
    n_hits = fasta_scan(seq, i, length, report_hits);
    for(j=0;j<n_hits;j++) {
      int db_seq, start, end, s_start, s_end, s_length, discount, score;

      db_seq = report_hits[j].db_seq;
      start = report_hits[j].start;
      end = report_hits[j].end;
      s_start = report_hits[j].s_start;
      s_end = report_hits[j].s_end;
      s_length = seqmeta[db_seq].seq_length;
      score = report_hits[j].score;

      discount = MIN(start, s_start) + MIN(length - end - 1, s_length - s_end - 1);
      fprintf(stdout,"%u %u %d %d %d %d %d %d %d %d %d\n",i,db_seq,score,
	      discount,score-discount,length,s_length, start, end, 
	      s_start, s_end);
    }

    reverse_complement(seq, length);
    n_hits = fasta_scan(seq, i, length, report_hits);
    for(j=0;j<n_hits;j++) {
      int db_seq, start, end, s_start, s_end, s_length, discount, score;

      db_seq = report_hits[j].db_seq;
      start = report_hits[j].start;
      end = report_hits[j].end;
      s_start = report_hits[j].s_start;
      s_end = report_hits[j].s_end;
      s_length = seqmeta[db_seq].seq_length;
      score = report_hits[j].score;

      discount = MIN(start, s_start) + MIN(length - end - 1, s_length - s_end - 1);
      fprintf(stdout,"%u %u %d %d %d %d %d %d %d %d %d RC\n",i,db_seq,score,
	      discount,score-discount,length,s_length, start, end, 
	      s_start, s_end);
    }


  }

  return 0;
}
