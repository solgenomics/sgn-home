#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <limits.h>

#include "ka.h"
#include "log_message.h"

extern sequence_t *sequences;
extern int n_seq;
extern uchar **seq_names;

static int wordsize;

#define GAP (6)
#define MISMATCH (3)
#define MATCH (2)

#define ALLOC_STEP (128)

#define PUSH(a,l,t) \
if ((l) % ALLOC_STEP == 0)   { \
   (a) = realloc((a), (t)*((l) + ALLOC_STEP)); \
   if ((a) == NULL) {                     \
      logmsg(MSG_FATAL,"Failed allocating memory at %s:%d (%d bytes)\n", \
             (t)*(l + ALLOC_STEP),__FILE__,__LINE__);  \
} \
}

#define CA(p,n,s)     \
(p) = calloc((n),(s)); \
if ((p) == NULL) {    \
   logmsg(MSG_FATAL,"Failed allocating memory at %s:%d (%d bytes)\n", \
          __FILE__,__LINE__,(n)*(s)); \
}

#define MA(p,s)     \
(p) = malloc((s));   \
if ((p) == NULL) {  \
  logmsg(MSG_FATAL,"Failed allocating memory at %s:%d (%d bytes)\n", \
	 __FILE__,__LINE__,(s));  \
}

typedef struct {
  int s;
  int start;
} lookup_t;

static lookup_t **lookup_table;
static int *lookup_lengths;

typedef struct {
  int diagonal;
  int s1_start;
  int s2_start;
  int length;
  int score;
} match_t;

static match_t **match;
static int *n_matches;

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

inline int dna_to_int(uchar *dna_string) {
  int i;
  int key;

  key = 0;
  for(i=0;i<wordsize;i++) {
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
      /* Cheat a little and treat N like an A */
      key = (key << 2);
      break;
    default:
      abort();
    }
  }    

  return key;
}

static void build_subsequence_lookup_table(void) {
  int i,j, key, total;
  lookup_t *table;
  int *temp_lengths;
  int table_size;

  logmsg(MSG_INFO,"Building subsequence lookup table...\n");
  table_size = 1 << (wordsize*2);

  CA(lookup_table, table_size, sizeof(lookup_t *));
  CA(lookup_lengths, table_size, sizeof(int));
  total = 0;
  for(i=0;i<n_seq;i++) {
    for(j=0;j<sequences[i].length-wordsize;j++) {
      key = dna_to_int(sequences[i].sequence + j);
      lookup_lengths[key]++;
      total++;
    }
  }
  
  logmsg(MSG_INFO,"Finished counting words... sequence lookup table "
	 "requires %d bytes", total*sizeof(lookup_t));
  
  MA(table, (total*sizeof(lookup_t)));
  CA(temp_lengths, table_size, sizeof(int));
  total = 0;
  for(i=0;i<(0x1 << wordsize*2);i++) {
    lookup_table[i] = table + total;
    total += lookup_lengths[i];
  }

  for(i=0;i<n_seq;i++) {
    for(j=0;j<sequences[i].length-wordsize;j++) {
      key = dna_to_int(sequences[i].sequence + j);
      lookup_table[key][temp_lengths[key]].s = i;
      lookup_table[key][temp_lengths[key]].start = j;
      temp_lengths[key]++;
    }
  }

  free(temp_lengths);
  logmsg(MSG_INFO,"Finished building lookup table.\n");
}

static void find_word_matches(int seq_index) {
  int i,j;
  int key;
  sequence_t *seqobj;

  seqobj = sequences + seq_index;

  CA(match, n_seq, sizeof(match_t *));
  CA(n_matches, n_seq, sizeof(int));

  for(i=0;i<seqobj->length-wordsize;i++) {
    key = dna_to_int(seqobj->sequence + i);
    
    if (lookup_lengths[key] == 0) continue;

    for(j=0;j<lookup_lengths[key];j++) {
      match_t *m;
      int match_seqid;
      match_seqid = lookup_table[key][j].s;

      if (match_seqid == seq_index) continue;

      PUSH(match[match_seqid], n_matches[match_seqid], sizeof(match_t));
      m = match[match_seqid] + n_matches[match_seqid]++;
      m->s1_start = i;
      m->s2_start = lookup_table[key][j].start;
      m->diagonal = m->s1_start - m->s2_start;
    }
  }
}

static int match_compare_length(const void *a, const void *b) {

  return ((match_t *) b)->length - ((match_t *) a)->length;
}

static void combine_consecutive_matches(int seq_index) {
  int i, j, k, p;

  for(i=0;i<n_seq;i++) {
    if (n_matches[i] == 0) continue;

    p = 0; j = 0;
    while(j<n_matches[i]) {
      k = j+1;
      while(k<n_matches[i] && (match[i][j].diagonal == match[i][k].diagonal) &&
	    ((match[i][k].s1_start-match[i][j].s1_start)<=(k-j)+wordsize*2)) 
	k++;

      if (k - j == 1 && p == j) {
	j = k;
	match[i][p].length = wordsize;
	match[i][p].score = wordsize*2;
	p++;
	continue;
      }

      match[i][p].diagonal = match[i][j].diagonal;
      match[i][p].s1_start = match[i][j].s1_start;
      match[i][p].s2_start = match[i][j].s2_start;
      match[i][p].length = match[i][k-1].s1_start - match[i][j].s1_start 
	+ wordsize;
      match[i][p].score = 0;
      for(;j<k-1;j++) {
	if (match[i][j].s1_start == match[i][j+1].s1_start - 1) {
	  match[i][p].score += 2;
	} else {
	  match[i][p].score += wordsize*2 - 5;
	}
      }
      match[i][p].score += wordsize*2;
      p++;
      j = k;
    }

    if (p > 50) {
      qsort(match[i], p, sizeof(match_t), match_compare_length);
      p = 50;
    }
    n_matches[i] = p;
    //    match[i] = realloc(match[i], p*sizeof(match_t)); 
  }
}

static int match_compare2(const void *a, const void *b) {
  
  return ((match_t *) a)->s1_start - ((match_t *) b)->s1_start;
}

/* NOTES:

   After above function, should have consecutive matches combined, but any
   mismatch will break the consecutive hits so each diagonal may have several
   runs of consecutive hits. Need to below collect the best of these
   consecutive runs (or diagonals) and store for fast finishing of the
   alignment with gaps/mismatches by dynamic programming. */
static void compute_fasta_init1_scores(int seq_index) {
  int i,j,k;
  int n_edges, n_nodes;
  edge_t *edges;
  node_t *nodes;
  int *pred;
  int *score;
  int max;
  int start, end;
  
  /* Resort by starting point rather than diagonal */
  max = 0;
  for(i=0;i<n_seq;i++) {
    if (n_matches[i])
      qsort(match[i], n_matches[i], sizeof(match_t), match_compare2);
    if (n_matches[i] > n_matches[max]) max = i;
  }

  /* Allocate the maximum needed upfront and recycle as we look through the
     sequences. This saves time managing allocation and deallocation. */
  n_nodes = n_matches[max]+2;
  MA(nodes, n_nodes*sizeof(node_t));
  MA(edges, n_nodes*(n_nodes)*sizeof(edge_t));
  MA(pred, n_nodes*sizeof(int));
  MA(score, n_nodes*sizeof(int));
  for(i=0;i<n_nodes;i++) {
    MA(nodes[i].out_edges, n_nodes*sizeof(int));
  }

  for(i=0;i<n_seq;i++) {
    if (n_matches[i] == 0) continue;
    n_nodes = 0;
    nodes[0].weight = 0;
    nodes[0].n_outedge = 0;
    nodes[0].s1_start = -1;
    nodes[0].s2_start = -1;
    nodes[0].diagonal = 0;
    n_nodes = 1;
    n_edges = 0;

    for(j=0;j<n_matches[i];j++) {
      nodes[n_nodes].weight = match[i][j].score;
      nodes[n_nodes].length = match[i][j].length;
      nodes[n_nodes].n_outedge = 0;
      nodes[n_nodes].s1_start = match[i][j].s1_start;
      nodes[n_nodes].s2_start = match[i][j].s2_start;
      nodes[n_nodes].diagonal = match[i][j].diagonal;

      /* Loop over existing nodes, adding edges if the new node starting
         position comes after a existing node's starting position */
      for(k=0;k<n_nodes;k++) {
	if (nodes[k].s1_start < match[i][j].s1_start) {
	  edges[n_edges].s_node = k;
	  edges[n_edges].e_node = n_nodes;

	  /* Weight of edge is cost of connecting these nodes. It should be
	     the gap cost (difference in diagonals) plus an adjustment if the 
	     match runs overlap. If the edge is to the start node, set the
	     weight to the trim cost.  */
	  if (k == 0) {
	    /* No gap cost at the begining -- only a charge for unaligned
	       overlapping sequence (a mismatch) */
	    if (nodes[n_nodes].s1_start < nodes[n_nodes].s2_start) {
	      edges[n_edges].weight = nodes[n_nodes].s1_start * MISMATCH;
	    } else {
	      edges[n_edges].weight = nodes[n_nodes].s2_start * MISMATCH;
	    }
	  } else {
	    /* GAP cost to get from prior nodes' diagonal to this one */
	    edges[n_edges].weight = 
	      abs(nodes[k].diagonal-match[i][j].diagonal) * GAP;

	    /* *Subtract* (positve edge weight) the portion which overlaps
	       in both matches since the node weights will include each of
	       these portions respectively, and thus result in double counting
	       of the overlap portion */
	    if (match[i][j].s1_start < (nodes[k].s1_start + nodes[k].length)) {
	      edges[n_edges].weight += (nodes[k].s1_start + nodes[k].length 
		- match[i][j].s1_start) * MATCH;
	    }
	  }

	  nodes[k].out_edges[nodes[k].n_outedge++] = n_edges;
	  n_edges ++;
	}
      }
      n_nodes++;
    }
    nodes[n_nodes].weight = 0;
    nodes[n_nodes].length = 0;
    nodes[n_nodes].n_outedge = 0;
    nodes[n_nodes].s1_start = -1;
    nodes[n_nodes].s2_start = -1;
    nodes[n_nodes].diagonal = 0;
    for(k=1;k<n_nodes;k++) {
      int x, y;
      edges[n_edges].s_node = k;
      edges[n_edges].e_node = n_nodes;
      x = nodes[k].s1_start + nodes[k].length;
      y = nodes[k].s2_start + nodes[k].length;
      if ( sequences[seq_index].length - x < sequences[i].length - y) {
	edges[n_edges].weight = (sequences[seq_index].length - x) * MISMATCH;
      } else {
	edges[n_edges].weight = (sequences[i].length - y) * MISMATCH;
      }
      
      nodes[k].out_edges[nodes[k].n_outedge++] = n_edges;
      n_edges ++;
    }
    n_nodes ++;

    logmsg(MSG_DEBUG4,"Graph has %d nodes and %d edges\n",n_nodes, n_edges);

    /* Single-source shortest path (sort of) */
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
    logmsg(MSG_DEBUG4,"Best path score = %d (%d)\n",score[max],
	   score[n_nodes-1]);
    k = pred[n_nodes-1];
    end = nodes[k].s1_start + nodes[k].length;
    start = nodes[k].s1_start;
    k = pred[k];
    while(k!=0) {
      start = nodes[k].s1_start;
      k = pred[k];
    }

    if (score[n_nodes-1] > 200) {
      fasta_t *f;

      /* Record this score for later use */
      PUSH(fasta_scores[seq_index], n_fasta[seq_index], sizeof(fasta_t));      
      f = fasta_scores[seq_index] + n_fasta[seq_index];
      f->s1 = seq_index;
      f->s2 = i;
      f->score = score[n_nodes-1];
      f->start = start;
      f->end = end;
      n_fasta[seq_index]++;
    }
  }

  free(nodes);
  free(edges);
  free(pred);
  free(score);
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

void pairwise_prescan(void) {
  int i,j;

  wordsize = 12;

  CA(fasta_scores, n_seq, sizeof(fasta_t *));
  CA(n_fasta, n_seq, sizeof(int));

  build_subsequence_lookup_table();

  for(i=0;i<n_seq;i++) {
    logmsg(MSG_DEBUG0,"Searching for pairwise alignments with sequence %s\n",
	   seq_names[i]);
    find_word_matches(i);
    
    for(j=0;j<n_seq;j++) {
      if (n_matches[i])
	qsort(match[i], n_matches[i], sizeof(match_t), match_compare);
    }

    combine_consecutive_matches(i);

    compute_fasta_init1_scores(i);

    for(j=0;j<n_seq;j++) {
      free(match[j]);
    }
    free(n_matches);
    free(match);
  }

  free(lookup_table[0]);
  free(lookup_lengths);
  free(lookup_table);
}


