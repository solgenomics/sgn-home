#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <math.h>

#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>

#include "ka.h"
#include "log_message.h"

/* Quick C program to find a local alignment and compute edit distance.

   Purpose: For processing a bunch of ESTs with homology to the RBCS tomato 
   genes as detected by BLAST. Here, we wish to compute the exact edit
   distance and determine which RBCS gene this EST is most likely to be
   a match with. BLAST scores/alignments/evalues are too rough for this
   fine level of discrimination */


typedef unsigned int uint32;

/* To speed the development of this "quick" program, I'm not coding for
   FASTA input. I am assuming the following input is given in a single file:
  
   type, length, value encoding (binary file format)

   TYPES:

   SEQ_HEADER (flags: Quality present or not?)
   LABEL
   SEQUENCE
   QUALITY

*/


int match_reward = 2;
int mismatch_penalty = -5;
int gap_penalty = -6;


typedef struct {
  int db_seq;
  int pos;
} word_t;

extern uint n_seq;
extern seq_t **sequences;

static double low_threshold, high_threshold;

static word_t *lookup_table;
static int *lookup_size;
static word_t **lookup;

#define WORDSIZE (9)
#define MASK (0x1FFFF)
static void build_wordlookup(seq_t **seq, int n_seq) {
  uint32 n_words;
  uint32 word;
  uint32 i, j;
  int total_size, length;
  int *fill_ptr;

  logmsg(MSG_INFO,"Building subsequence lookup table...\n");
  n_words = 0x1 << (WORDSIZE*2);
  CA(lookup_size, n_words, sizeof(int));
  CA(lookup, n_words, sizeof(word_t *));

  for(i=0;i<n_seq;i++) {
    if (sequences[i]->expected_error < low_threshold || 
	sequences[i]->expected_error > high_threshold) {
      continue;
    }
    length = seq[i]->length - WORDSIZE;
    word = 0;
    for(j=0;j<WORDSIZE;j++)
      word = (word << 2) | seq[i]->seq[j];
    word &= MASK;
    lookup_size[word]++;
    for(;j<length;j++) {
      word = ((word << 2) & MASK) | seq[i]->seq[j];
      lookup_size[word]++;
    }  
  }
  total_size = 0;
  for(i=0;i<n_words;i++)
    total_size += lookup_size[i];
  logmsg(MSG_INFO,"Finished counting words... sequence lookup table "
	 "requires %d bytes", total_size*sizeof(word_t));
  MA(lookup_table, sizeof(word_t)*total_size);
  total_size = 0;
  for(i=0;i<n_words;i++) {
    if (lookup_size[i] == 0) continue;
    lookup[i] = lookup_table + total_size;
    total_size += lookup_size[i];
  }

  CA(fill_ptr, n_words, sizeof(int));
  for(i=0;i<n_seq;i++) {
    if (sequences[i]->expected_error < low_threshold || 
	sequences[i]->expected_error > high_threshold) {
      continue;
    }
    length = seq[i]->length - WORDSIZE;
    word = 0;
    for(j=0;j<WORDSIZE;j++)
      word = (word << 2) | seq[i]->seq[j];
    word &= MASK;
    lookup[word][fill_ptr[word]].db_seq = i;
    lookup[word][fill_ptr[word]].pos = 0;
    fill_ptr[word]++;
    for(;j<length;j++) {
      word = ((word << 2) & MASK) | seq[i]->seq[j];
      lookup[word][fill_ptr[word]].db_seq = i;
      lookup[word][fill_ptr[word]].pos = j - WORDSIZE;
      fill_ptr[word]++;
    }
  }
  free(fill_ptr);

  logmsg(MSG_INFO,"Finished building lookup table.\n");
}

typedef struct {
  int db_seq;
  int di;
  int pos;
  int length;
} wordhit_t;

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
} hit_report_t;

#define SWAP(x,y) ((x)^=(y),(y)^=(x),(x)^=(y))
static int fasta_scan(seq_t *seq, seq_t **db, int db_size, 
		      hit_report_t *report_hits) {
  int i, j, k, l, f, t;
  uint32 word;
  int n_hits, length, n_nodes, s, max, max_span;
  int min_di, max_di, total_length;
  wordhit_t *hits;
  int *adjmatrix;
  int *pred, *score;

  n_hits = 0;
  length = seq->length - WORDSIZE;
  word = 0;
  for(i=0;i<WORDSIZE;i++) 
    word = (word << 2) + seq->seq[i];
  n_hits += lookup_size[word];
  for(;i<length;i++) {
    word = ((word << 2) & MASK) + seq->seq[i];
    n_hits += lookup_size[word];
  }
  if (n_hits == 0) return 0;
  
  t = 0;
  MA(hits, n_hits*sizeof(wordhit_t));
  word = 0;
  for(i=0;i<WORDSIZE;i++) 
    word = (word << 2) + seq->seq[i];
  for(j=0;j<lookup_size[word];j++) {
    hits[t].db_seq = lookup[word][j].db_seq;
    hits[t].di = lookup[word][j].pos - (i - WORDSIZE);
    hits[t].pos = (i - WORDSIZE);
    t++;
  }
  for(;i<length;i++) {
    word = ((word << 2) & MASK) + seq->seq[i];
    for(j=0;j<lookup_size[word];j++) {
      hits[t].db_seq = lookup[word][j].db_seq;
      hits[t].di = lookup[word][j].pos - (i - WORDSIZE);
      hits[t].pos = (i - WORDSIZE);
      t++;
    }
  }
  assert(t == n_hits);
  
  qsort(hits, n_hits, sizeof(wordhit_t), wordhit_compare);
  f = i = j = 0;
  while(i < n_hits) {
    while(j < n_hits                       && 
	  hits[i].db_seq == hits[j].db_seq &&
	  hits[i].di == hits[j].di         && 
	  hits[j].pos - hits[i].pos == j - i) j++;
    hits[f].di = hits[i].di;
    hits[f].db_seq = hits[i].db_seq;
    hits[f].pos = hits[i].pos;
    hits[f].length = j-i + WORDSIZE;
    f++;
    i = j;
  }

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
    for(l=1;l<n_nodes-1;l++) {
      adjmatrix[l] = 1;
      adjmatrix[l*n_nodes + (n_nodes-1)] = 1;
    }
    for(k=0;k<n_nodes;k++) {
      pred[k] = -1;
      score[k] = INT_MIN;
    }
    pred[0] = -1;
    score[0] = 0;
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
    max = 0;
    for(k=1;k<n_nodes;k++) {
      if (score[k] > score[max]) max = k;
    }
    k = max;
    min_di = hits[i+k-1].di;
    max_di = hits[i+k-1].di;
    total_length = 0;
    while(k!=0) {
      total_length += hits[i+k-1].length;
      if (hits[i+k-1].di < min_di) {
	min_di = hits[i+k-1].di;
      } else if (hits[i+k-1].di > max_di) {
	max_di = hits[i+k-1].di;
      }
      k = pred[k];
    }
    if (total_length >= 75 && abs(max_di - min_di) < 50) {
      report_hits[n_hits].db_seq = hits[i].db_seq;
      report_hits[n_hits].min_di = min_di - 5;
      report_hits[n_hits].max_di = max_di + 5;
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


typedef struct {
  int matches;
  int mis_matches;
  int insertions;
  int deletions;
  seq_t *query;
  seq_t *db;
  int q_start;
  int q_end;
  int s_start;
  int s_end;
  double expected_error;
} align_t;


typedef struct {
  int score;
  int bt;
} align_matrix_t;

static align_t align_sequences(seq_t *q, seq_t *s, int min_di, int max_di) {
  int n_rows, n_cols, i, j;
  align_matrix_t *matrix;
  align_matrix_t **m;
  uchar *qstr, *sstr, *astr;
  int max_score, score, dir;
  int p, max_p;
  int matches, mis_matches;
  int insertions, deletions;
  int start_row, start_col, di_span;
  int best_i, best_j;
  align_t r;
  double expected_error;

  n_rows = q->length + 1;
  n_cols = s->length + 1;

  di_span = abs(min_di - max_di);

  matrix = malloc(n_cols*n_rows*sizeof(align_matrix_t));
  m = alloca(n_rows*sizeof(align_matrix_t *));
  //  gap_q = malloc(n_rows*sizeof(int *));

  for(i=0;i<n_rows;i++) {
    m[i] = matrix + i*n_cols;
    //gap_q[i] = calloc(n_cols, sizeof(int));
  }

  if (max_di >= 0) {
    start_row = 1;
  } else {
    start_row = abs(max_di) + 1;
  }

  start_col = min_di + start_row - 1;

  for(j=start_col-1;j<start_col+di_span+4;j++) {
    if (j>=0 && j<n_cols) {
      m[start_row-1][j].score = 0;
      m[start_row-1][j].bt = 0;
    }
  }

  for(i=start_row;i<n_rows;i++) {
    j = min_di + i - 1;
    if (j<0) j = 0;
    else if (j>=n_cols) j= n_cols-1;
    m[i][j].score = 0;
    m[i][j].bt = 0;
    j = max_di + i + 1;
    if (j<0) j = 0;
    else if (j>=n_cols) j= n_cols-1;
    m[i][j].score= 0;
    m[i][j].bt = 0;
  }

  best_i = best_j = 0;
  m[0][0].score = 0;
  m[0][0].bt = 0;
  for(i=start_row;i<n_rows;i++) {
    int stop;
    start_col = min_di + i;
    if (start_col <= 0) start_col = 1;
    stop = max_di + i + 1;
    if (stop >= n_cols) stop = n_cols;
    for(j=start_col;j<stop;j++) {
      dir = 0;
      max_score = 0;

      if (q->seq[i-1] == s->seq[j-1]) {
	score = m[i-1][j-1].score + match_reward;
      } else {
	score = m[i-1][j-1].score + mismatch_penalty;
      }
      if (score > max_score) {
	max_score = score;
	dir = 1;      
      }
      
      if (j < n_cols - 1)
	score = m[i-1][j].score + gap_penalty;
      else
	score = m[i-1][j].score;

      if (score > max_score) {
	max_score = score;
	dir = 2;
      }

      if (i < n_rows - 1)
	score = m[i][j-1].score + gap_penalty;
      else
	score = m[i][j-1].score;
      if (score > max_score) {
	max_score = score;
	dir = 3;
      }

      m[i][j].bt = dir;
      m[i][j].score = max_score;
      if (max_score > m[best_i][best_j].score) {
	best_i = i;
	best_j = j;
      }
    }
  }

  i = best_i;
  j = best_j;
  
#if 0
  qstr = calloc(1, n_rows + n_cols + 1);
  sstr = calloc(1, n_rows + n_cols + 1);
  astr = calloc(1, n_rows + n_cols + 1);

  p = best_i + best_j + 1;
  while(i<n_rows-1 || j<n_cols-1) {
    if (i<n_rows-1) {
      qstr[p] = q->seqstr[i++];
    } else {
      qstr[p] = '-';
    }
    if (j<n_cols-1) {
      sstr[p] = s->seqstr[j++];
    } else {
      sstr[p] = '-';
    }
    if (sstr[p] == qstr[p] && qstr[p]!='-') {
      astr[p] = '!';
    } else {
      astr[p] = ' ';
    }
    p++;
  }
  max_p = p;

  i = best_i;
  j = best_j;
  p = best_i + best_j;
  matches = mis_matches = insertions = deletions = 0;
  while(m[i][j].bt != 0) {
    switch(m[i][j].bt) {
    case 1:
      i--;
      j--;
      if (q->seq[i] == s->seq[j]) {
	astr[p] = '|';
	matches++;
      } else {
	astr[p] = ' ';
	mis_matches++;
      }
      qstr[p] = q->seqstr[i];
      sstr[p] = s->seqstr[j];
      break;
    case 2:
      if (matches > 0) insertions++;
      astr[p] = ' ';
      i--;
      qstr[p] = q->seqstr[i];
      sstr[p] = '-';
      break;
    case 3:
      if (matches > 0) deletions++;
      astr[p] = ' ';
      j--;
      qstr[p] = '-';
      sstr[p] = s->seqstr[j];
      break;
    }
    p--;
  }
  while(i>0 || j>0) {
    if (i>0) {
      i--;
      qstr[p] = q->seqstr[i];
    } else {
      qstr[p] = '-';
    }
    if (j>0) {
      j--;
      sstr[p] = s->seqstr[j];
    } else {
      sstr[p] = '-';
    }
    if (sstr[p] == qstr[p]) {
      astr[p] = '!';
    } else {
      astr[p] = ' ';
    }
    p--;
  }

  fprintf(stdout,"Diagonal span (%d,%d)\n",min_di, max_di);
  p+=1;
  while(p<max_p) {
    fprintf(stdout,"%-20.20s   %-120.120s\n",q->label, qstr + p);
    fprintf(stdout,"                       %-120.120s\n",astr + p);
    fprintf(stdout,"%-20.20s   %-120.120s\n",s->label, sstr + p);
    fprintf(stdout,"\n");
    p += 120;
  }
  fprintf(stdout,"\n\n");
  fprintf(stdout,"Matches = %d\tMismatches = %d\tInsertions = %d\tDeletions = %d\n",matches,mis_matches,insertions,deletions);
  free(astr);
  free(qstr);
  free(sstr);
#else
  i = best_i;
  j = best_j;
  r.q_end = best_i;
  r.s_end = best_j;
  matches = mis_matches = insertions = deletions = 0;
  expected_error = 0.0;
  while(m[i][j].bt != 0) {
    double a,b;
    switch(m[i][j].bt) {
    case 1:
      i--;
      j--;
      if (q->seq[i] == s->seq[j]) {
	matches++;
      } else {
	mis_matches++;
	a = pow(10.0, q->qual[i]/-10.0);
	b = pow(10.0, s->qual[j]/-10.0);
	expected_error += a + b - (a*b);
	  
      }
      break;
    case 2:
      i--;
      insertions++;
      a = pow(10.0, q->qual[i]/-10.0);
      b = pow(10.0, s->qual[j]/-10.0);
      expected_error += a + b - (a*b);
      a = pow(10.0, q->qual[i-1]/-10.0);
      b = pow(10.0, s->qual[j-1]/-10.0);
      expected_error += a + b - (a*b);
      a = pow(10.0, q->qual[i+1]/-10.0);
      b = pow(10.0, s->qual[j+1]/-10.0);
      expected_error += a + b - (a*b);
      break;
    case 3:
      j--;
      deletions++;
      a = pow(10.0, q->qual[i]/-10.0);
      b = pow(10.0, s->qual[j]/-10.0);
      expected_error += a + b - (a*b);
      a = pow(10.0, q->qual[i-1]/-10.0);
      b = pow(10.0, s->qual[j-1]/-10.0);
      expected_error += a + b - (a*b);
      a = pow(10.0, q->qual[i+1]/-10.0);
      b = pow(10.0, s->qual[j+1]/-10.0);
      expected_error += a + b - (a*b);
      break;
    }
  }

#endif

  free(matrix);

  r.q_start = i;
  r.s_start = j;
  r.expected_error = expected_error;
  r.matches = matches;
  r.mis_matches = mis_matches;
  r.insertions = insertions;
  r.deletions = deletions;
  
  return r;
}

#if 0
int main(int argc, char *argv[]) {
  seq_t *seq;
  seq_t **db;
  int db_size, i, j, n_hits, x, length;
  int database_size;
  FILE *f;
  struct stat stat_buf;
  align_t al;
  hit_report_t *hit_report;
  double expected_error;

  if (argc!=2) {
    fprintf(stderr,"Usage: %s <database filename>\n",argv[0]);
    fprintf(stderr,"Send query data as standard input.\n");
    exit(-1);
  }

  if (stat(argv[1],&stat_buf)) {
    fprintf(stderr,"Can't stat() file %s (%s)\n",argv[1],strerror(errno));
    fprintf(stderr,"Usage: %s <database filename>\n",argv[0]);
    fprintf(stderr,"Send query data as standard input.\n");
    exit(-1);
  }
  database_size = stat_buf.st_size;
  f = fopen(argv[1],"r");
  if (f == NULL) {
    fprintf(stderr,"Failed opening file %s (%s)\n",argv[1],strerror(errno));
    fprintf(stderr,"Usage: %s <database filename>\n",argv[0]);
    fprintf(stderr,"Send query data as standard input.\n");
    exit(-1);
  }

  db_size = 0;
  db = NULL;
  while(ftell(f) < database_size) {
    seq = read_sequence(f);
    if (seq == NULL) {
      fprintf(stderr,"Failed reading sequence at position %ld\n",
	      ftell(f));
      continue;
    }
    db = realloc(db, (sizeof(seq_t *)*(db_size+2)));
    db[db_size++] = seq;
    db[db_size++] = reverse_complement(seq);
  }
  fclose(f);
  fprintf(stderr,"Building lookup table of database sequence words..");
  build_wordlookup(db, db_size);
  MA(hit_report, sizeof(hit_report_t)*db_size);
  fprintf(stderr,"Done.\n");

  for(i=0;i<db_size;i++) {
    fprintf(stderr,"Scanning sequence %s\n",db[i]->label);
    n_hits = fasta_scan(db[i], db, db_size, hit_report);
    for(j=0;j<n_hits;j++) {
      //      if (hit_report[j].db_seq == i) continue;
      al = align_sequences(db[i], db[hit_report[j].db_seq], 
			   hit_report[j].min_di, hit_report[j].max_di);
      x = al.mis_matches + al.insertions + al.deletions;      
      if (al.matches > 100) {
	length = x + al.matches;
	expected_error = (db[i]->expected_error + 
			  db[hit_report[j].db_seq]->expected_error)*0.75;
	fprintf(stdout,"%s\t%s\t%30.30f\t%30.30f\t%30.30f\t (%d)\n",
		db[i]->label, db[hit_report[j].db_seq]->label,
		x/(double) length, expected_error, 
		((x)/(double) length) - expected_error, x);
      } 
    }

  } 
  

  return 0;
}
#endif


fasta_t **fasta_scores;
int *n_fasta;

static int fasta_obsdiff_compare(const void *a, const void *b) {

  if (((fasta_t *) a)->obs_diff < ((fasta_t *) b)->obs_diff)
    return -1;
  else if (((fasta_t *) a)->obs_diff > ((fasta_t *) b)->obs_diff)
    return 1;
  else
    return 0;
}

static double select_parameter(void) {
  int i, j, n_bins, t, bin;
  int total_obs;
  double bin_size, top, bot;
  double *bins;
  

  bot = 1000.0;
  top = -1000.0;
  total_obs = 0;
  for(i=0;i<n_seq;i++) {
    for(j=0;j<n_fasta[i];j++) {
      if (fasta_scores[i][j].obs_diff < bot) 
	bot = fasta_scores[i][j].obs_diff;
      if (fasta_scores[i][j].obs_diff > top)
	top = fasta_scores[i][j].obs_diff;
    }
    total_obs += n_fasta[i];
  }
  fprintf(stderr,"Total observations %d\n", total_obs);

  n_bins = n_seq;
  bins = alloca(sizeof(double)*n_bins);
  bin_size = (top - bot + 0.001)/(double) n_bins;
  
  for(i=0;i<n_bins;i++) {
    bins[i] = 0.0;
  }
  for(i=0;i<n_seq;i++) {
    for(j=0;j<n_fasta[i];j++) {
      bin = (int) ((fasta_scores[i][j].obs_diff - bot)/(double) bin_size);
      bins[bin] += 1.0;
    }
  }

  for(i=0;i<n_bins;i++) {
    bins[i] /= (double) total_obs;
    fprintf(stderr,"bin %d (%6.6f , %6.6f)\t (%6.6f%%)\n",i,
	    bot + i*bin_size, bot + (i+1)*bin_size, bins[i]);
  }
  for(i=0;i<n_bins;i++) {
    if (bins[i] == 0.0) break;
  }
  while(i<n_bins && bins[i] == 0.0) i++;
  if (i < n_bins) {
    fprintf(stderr,"Global parameter selected: %30.30f\n", bot + i*bin_size);
    return bot + i*bin_size;
  } else {
    fprintf(stderr,"Global paramter selection failed: %30.30f\n",
	    (top - bot) / 2.0);
    return (top - bot) / 2.0;
  }
}


static void analyze_fasta_scores(int seq_id, double global_parameter) {
  int i, f;

  f = 0;
  if (sequences[seq_id]->expected_error < low_threshold ||
      sequences[seq_id]->expected_error > high_threshold) {
    qsort(fasta_scores[seq_id], n_fasta[seq_id], sizeof(fasta_t),
	  fasta_obsdiff_compare);
    n_fasta[seq_id] = 3;
    fprintf(stderr,"Sequence %s censored due to significant differences in quality with sample (%8.8f %8.8f - %8.8f)\n",sequences[seq_id]->label,
	    sequences[seq_id]->expected_error, low_threshold, high_threshold);
  }
  else {
    int i, j, max, maxi;
    int select_bin;
    double top, bot, bin_size, parm, k;
    volatile int n_bins, pow2, bin; 
    double sum;
    double *bins;

    select_bin = 0;
    parm = -1.0;
    qsort(fasta_scores[seq_id], n_fasta[seq_id], sizeof(fasta_t),
	  fasta_obsdiff_compare);
    top = fasta_scores[seq_id][n_fasta[seq_id]-1].obs_diff + 0.01;
    bot = fasta_scores[seq_id][0].obs_diff;
    n_bins = n_fasta[seq_id];
    pow2 = 1;
    while(n_bins > 0) {
      pow2 <<= 1;
      n_bins >>= 1;
    }
    assert((pow2 & (pow2 - 1)) == 0);
    bins = alloca(sizeof(double)*(pow2));
    n_bins = 2;
    while(n_bins <= pow2) {
      bin_size = (top - bot) /(double) n_bins;
      for(i=0;i<n_bins;i++) {
	bins[i] = 0.0;
      }
      for(i=0;i<n_fasta[seq_id];i++) {
	bin = (int) ((fasta_scores[seq_id][i].obs_diff - bot)/bin_size);
	assert(bin>=0 && bin < n_bins);
	bins[bin]++;
      }
#if 1
      fprintf(stderr,"Bin size = %6.6f\t N = %d\n",bin_size,n_bins);
      for(i=0;i<n_bins;i++) {
	fprintf(stderr,"bin %d (%6.6f,%6.6f) %30.30f\n",i,
	 	i*bin_size + bot, (i+1)*bin_size + bot,
		bins[i]/n_fasta[seq_id]);
	
      }
#endif
#if 0
      i = 0;
      maxi = 0;
      max = 0;
      while(i<n_bins) {
	if (bins[i] == 0.0) {
	  j = i+1;
	  while(j<n_bins && bins[j] == 0.0) j++;
	  if (j-i > max) {
	    max = j-i;
	    maxi = i;
	  }
	}
	i++;
      }
      parm = bot + bin_size*maxi;
#else
      k = n_fasta[seq_id]/(double) n_bins;
      //      while(i<n_bins && bins[i] < (2.0/(double) n_fasta[seq_id])) i++;
      select_bin = 0;
      for(i=0;i<n_bins;i++) {
	bins[i] /= (double) n_fasta[seq_id];
      }
      for(i=1;i<n_bins;i++) {
	if (bins[i] == 0.0) {
	  sum = 0.0;
	  for(j=i+1;j<n_bins;j++)
	    sum += bins[j];
	  if (sum >= (1.0/(double) n_bins)) { 
	    select_bin = i;
	    break;
	  }
	}
      }
      if (select_bin > 0) {
	parm = bot + select_bin*bin_size;
	break;
      }  
#endif
      n_bins <<= 1;
    }
    i = 0;
    while(i<n_fasta[seq_id] && fasta_scores[seq_id][i].obs_diff < parm)
      i++;
    fprintf(stderr,"%s Selected parameter = %30.30f (%5.3f%%)\n",sequences[seq_id]->label, parm,
	    i/(double) n_fasta[seq_id] * 100.0);
    n_fasta[seq_id] = i;
  }
}


static void identify_badsequences(void) {
  int i;
  double sum_sq, sum, var, stdev, mean;

  sum_sq = sum = 0.0;
  for(i=0;i<n_seq;i++) {
    sum_sq += sequences[i]->expected_error*sequences[i]->expected_error;
    sum += sequences[i]->expected_error;
  }
  
  mean = sum/n_seq;
  var = (sum_sq/n_seq - mean*mean)*((n_seq + 1.0)/(double) n_seq);
  stdev = sqrt(var);

  low_threshold = mean - 3*stdev;
  high_threshold = mean + 3*stdev;
}

void pairwise_prescan(void) {
  int i,j;
  int n_hits;
  align_t al;
  int x, length;
  hit_report_t *hit_report;
  double expected_error, global_parameter;
  int stupid = 0;

  CA(fasta_scores, n_seq, sizeof(fasta_t *));
  CA(n_fasta, n_seq, sizeof(int));

  identify_badsequences();

  build_wordlookup(sequences, n_seq);
  MA(hit_report, sizeof(hit_report_t)*n_seq);
  for(i=0;i<n_seq/2;i++) {
    fprintf(stderr,"Scanning sequence %s....",sequences[i]->label);
    n_hits = fasta_scan(sequences[i], sequences, n_seq, hit_report);
    fprintf(stderr,"%d hits\n",n_hits);
    for(j=0;j<n_hits;j++) {
      int sbjct, complement;
      if (hit_report[j].db_seq <= i) continue;
      sbjct = hit_report[j].db_seq;
      al = align_sequences(sequences[i], sequences[sbjct], 
			   hit_report[j].min_di, hit_report[j].max_di);
      x = al.mis_matches + al.insertions + al.deletions;      
      length = x + al.matches;
      if (length >= 50) {
	double query_ee, sbjct_ee;
	int k;

#if 0
	query_ee = sequences[i]->expected_error;
	for(k=0;k<al.q_start;k++) {
	  query_ee -= pow(10.0, sequences[i]->qual[k]/-10.0);
	}
	for(k=al.q_end;k<sequences[i]->length;k++) {
	  query_ee -= pow(10.0, sequences[i]->qual[k]/-10.0);
	}
	query_ee /= (double) length;
	sbjct_ee = sequences[sbjct]->expected_error;
	for(k=0;k<al.s_start;k++) {
	  sbjct_ee -= pow(10.0, sequences[sbjct]->qual[k]/-10.0);
	}
	for(k=al.s_end;k<sequences[sbjct]->length;k++) {
	  sbjct_ee -= pow(10.0, sequences[sbjct]->qual[k]/-10.0);
	}
	sbjct_ee /= (double) length;
	expected_error = (sequences[i]->expected_error + 
			  sequences[sbjct]->expected_error)*0.75;
	/*	fprintf(stdout,"%s\t%s\t%30.30f\t%30.30f\t%30.30f\t (%d)\n",
		sequences[i]->label, sequences[sbjct]->label,
		x/(double) length, expected_error, 
		((x)/(double) length) - expected_error, x);*/
#endif
	PUSH(fasta_scores[i], n_fasta[i], sizeof(fasta_t));
	fasta_scores[i][n_fasta[i]].s2 = sbjct;
	fasta_scores[i][n_fasta[i]].obs_diff = (x - al.expected_error)/(double) length;
	fasta_scores[i][n_fasta[i]].score = al.matches - x;
	n_fasta[i]++;
	
	if (sequences[i]->expected_error < high_threshold) {
	  PUSH(fasta_scores[sbjct], n_fasta[sbjct], sizeof(fasta_t));
	  fasta_scores[sbjct][n_fasta[sbjct]].s2 = i;
	  fasta_scores[sbjct][n_fasta[sbjct]].obs_diff = (x - al.expected_error)/(double) length;
	  fasta_scores[sbjct][n_fasta[sbjct]].score = al.matches - x;
	  n_fasta[sbjct]++;
	}
	
	complement = i + n_seq/2;
	if (sbjct >= n_seq/2) {
	  sbjct -= n_seq/2;
	} else {
	  sbjct += n_seq/2;
	}
	PUSH(fasta_scores[complement], n_fasta[complement], sizeof(fasta_t));
	fasta_scores[complement][n_fasta[complement]].s2 = sbjct;
	fasta_scores[complement][n_fasta[complement]].obs_diff = (x - al.expected_error)/(double) length;
	fasta_scores[complement][n_fasta[complement]].score = al.matches - x;
	n_fasta[complement]++;
	
	if (sequences[complement]->expected_error < high_threshold) {
	  PUSH(fasta_scores[sbjct], n_fasta[sbjct], sizeof(fasta_t));
	  fasta_scores[sbjct][n_fasta[sbjct]].s2 = complement;
	  fasta_scores[sbjct][n_fasta[sbjct]].obs_diff = (x - al.expected_error)/(double) length;
	  fasta_scores[sbjct][n_fasta[sbjct]].score = al.matches - x;
	  n_fasta[sbjct]++;
	}
      }  else {
	fprintf(stdout,"Stupid comparision: %d %d %d %d\n",al.matches,al.mis_matches,al.insertions,al.deletions);
	stupid++;
      }
    }
  }
  free(hit_report);  

  global_parameter = select_parameter();
  for(i=0;i<n_seq;i++) {
    for(j=0;j<n_fasta[i];j++) {
      assert(sequences[fasta_scores[i][j].s2]->expected_error <= high_threshold);
    }
    if (n_fasta[i] > 1) {
      analyze_fasta_scores(i, global_parameter);
    }
    for(j=0;j<n_fasta[i];j++) {
      assert(sequences[fasta_scores[i][j].s2]->expected_error <= high_threshold);
    }
  }

#if 1
  for(i=0;i<n_seq;i++) {
    int j;
    for(j=0;j<n_fasta[i];j++) {
      int target, k;
      target = fasta_scores[i][j].s2;
      for(k=0;k<n_fasta[target];k++) {
	if (fasta_scores[target][k].s2 == i) break;
      }
      if (k == n_fasta[target]) {
	assert(sequences[target]->expected_error <= high_threshold);
	fasta_scores[target][n_fasta[target]].s2 = i;
	fasta_scores[target][n_fasta[target]].obs_diff = 
	  fasta_scores[i][j].obs_diff;
	fasta_scores[target][n_fasta[target]].score = 
	  fasta_scores[i][j].score;
	n_fasta[target]++;
      }
    }
  }
#endif 
  fprintf(stderr,"Total stupidity: %d\n",stupid);
}
