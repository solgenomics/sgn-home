#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <limits.h>

#include "ka.h"
#include "log_message.h"

#define ALLOC_STEP (32)

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

extern seq_t **sequences;
extern int n_seq;

extern fasta_t **fasta_scores;
extern int *n_fasta;


static int **tree_edges;
static int *n_tedges;
static int **back_edges;
static int *n_bedges;

static void cc_dfs_visit(int *color, int **component, int *component_size, 
			 int n) {
  int i;

  color[n] = 1;
  if (n > n_seq/2) {
    color[n - n_seq/2] = 1;
  } else {
    color[n + n_seq/2] = 1;
  }
  PUSH(*component, *component_size, sizeof(int));
  (*component)[(*component_size)++] = n;
  for(i=0;i<n_fasta[n];i++) {
    if (color[fasta_scores[n][i].s2] == 0) {
      cc_dfs_visit(color, component, component_size, fasta_scores[n][i].s2);
      PUSH(tree_edges[n], n_tedges[n], sizeof(int));
      tree_edges[n][n_tedges[n]++] = fasta_scores[n][i].s2;
    } else {
      PUSH(back_edges[n], n_bedges[n], sizeof(int));
      back_edges[n][n_bedges[n]++] = fasta_scores[n][i].s2;
    }
  }
}

static int scan_forward(int candidate, int ancestor) {
  int i;

  if (ancestor == candidate) return 1;
  for(i=0;i<n_tedges[ancestor];i++) {
    if (scan_forward(candidate, tree_edges[ancestor][i])) return 1;
  }

  return 0;
}

static int scan_descendant(int candidate, int descendant) {
  int i;

  for(i=0;i<n_bedges[descendant];i++) {
    /* We want to find a *proper* ancestor of candidate */
    if (back_edges[descendant][i] == candidate) continue;
    if (scan_forward(candidate, back_edges[descendant][i])) return 1;
  }
  /* Check all descendants for a back edges */
  for(i=0;i<n_tedges[descendant];i++)
    if (scan_descendant(candidate, tree_edges[descendant][i])) return 1;

  /* If this node has no backedges, and no descendants have back edges to
     proper ancestors of <candidate>, then candidate is an articulation point */
  return 0;
}

static int scan_atri_point(int candidate) {
  int i;

  /* For <candidate> to be a articulation point, it must have a child s.t. 
     no descendant of the child (including the child itself) has a back-edge
     to a proper anscestor of candidate. So, finding one such child indicates
     an articulation point. */
  for(i=0;i<n_tedges[candidate];i++) {
    if (scan_descendant(candidate, tree_edges[candidate][i]) == 0) return 1;
  }

  return 0;
}

int **components;
int *component_size;
int n_components;

void connected_components() {
  int *color;
  int i;
  int rc_count, divide_point, j;

  CA(color, n_seq, sizeof(int));
  CA(back_edges, n_seq, sizeof(int *));
  CA(tree_edges, n_seq, sizeof(int *));
  CA(n_tedges, n_seq, sizeof(int));
  CA(n_bedges, n_seq, sizeof(int));


  n_components = 0;
  components = NULL;
  component_size = NULL;
  for(i=0;i<n_seq;i++) {
    if (color[i]) continue;
    PUSH(components, n_components, sizeof(int *));
    PUSH(component_size, n_components, sizeof(int));
    components[n_components] = NULL;
    component_size[n_components] = 0;
    cc_dfs_visit(color, components + n_components, 
		 component_size + n_components, i);
    n_components++;
  }

  divide_point = n_seq/2;
  for(i=0;i<n_components;i++) {
    rc_count = 0;
    for(j=0;j<component_size[i];j++) {
      if (components[i][j] > divide_point)
	rc_count++;
    }
    if (rc_count > component_size[i]/2) {
      for(j=0;j<component_size[i];j++)
	if (components[i][j] > divide_point) {
	  components[i][j] -= divide_point;
	} else {
	  components[i][j] += divide_point;
	}
    }
  }

  fprintf(stderr,"Found %d connected components\n", n_components);
  for(i=0;i<n_components;i++) {
    int j;
    fprintf(stderr,"Component %d\n",i);
    if (n_tedges[components[i][0]] > 1) {
      fprintf(stderr,"Articulation point %d %s\n",components[i][0],
	      sequences[components[i][0]]->label);
    }
    for(j=0;j<component_size[i];j++) {
      fprintf(stderr,"\t%s\n",sequences[components[i][j]]->label);
      if (j!=0 && scan_atri_point(components[i][j])) {
	fprintf(stderr,"%s (%d) is an atriculation point\n",
		sequences[components[i][j]]->label,components[i][j]);
	/* Scan hit list for chimera potential */
      }
    }
    fprintf(stderr,"\n\n");
  }

  free(color);
  free(back_edges);
  free(tree_edges);
  free(n_tedges);
  free(n_bedges);

}

typedef struct {
  int node;
  int key;
} pri_queue_t;

asm_order_t **assembly_order;

static int asm_order_compare(const void *a, const void *b) {
  
  return ((asm_order_t *) b)->score - ((asm_order_t *) a)->score;
}

void spanning_tree(void) {
  pri_queue_t *q;
  int *p, *used, *allowed;
  int root, max_score;
  int i, j, k, n_used;

  MA(allowed, (n_seq*sizeof(int)));
  MA(used, (n_seq*sizeof(int)));
  MA(q, (n_seq*sizeof(pri_queue_t)));
  MA(p, (n_seq*sizeof(int)));
  CA(assembly_order, n_components, sizeof(asm_order_t *));
  for(i=0;i<n_components;i++) {
    if (component_size[i] == 1) continue;
    /* Initialize state... empty queue, nothing allowed, nothing used, 
       all predecessors are invalid: */
    for(j=0;j<n_seq;j++) {
      allowed[j] = 0;
      used[j] = 0;
      q[j].node = j;
      q[j].key = INT_MIN;
      p[j] = -1;
    }

    root = -1;
    max_score = -1;
    for(j=0;j<component_size[i];j++) {
      allowed[components[i][j]] = 1;
      for(k=0;k<n_fasta[components[i][j]];k++) {
	if (fasta_scores[components[i][j]][k].score > max_score) {
	  root = components[i][j];
	  max_score = fasta_scores[components[i][j]][k].score;
	}
      }
    }

    assert(root >= 0);
    q[root].key = 0;

    n_used = component_size[i];;
    while(n_used > 0) {

      /* Find the best current edge in the graph, connecting an unused
         node... this is a lazy and unefficient way to extract the maximum
         from the "priority queue" -- it should be implemented as a heap. */
      max_score = -1;
      for(j=0;j<component_size[i];j++) {
	if (used[components[i][j]]) continue;
	if (max_score == -1 || q[components[i][j]].key > q[max_score].key) {
	  max_score = components[i][j];
	}
      }
      used[max_score] = 1;
      n_used--;

      /* Scan the edges from our new node and update the priority queue */
      for(j=0;j<n_fasta[max_score];j++) {
	if (! allowed[fasta_scores[max_score][j].s2]) continue;
	if (used[fasta_scores[max_score][j].s2]) continue;
	if (fasta_scores[max_score][j].score > 
	    q[fasta_scores[max_score][j].s2].key) {
	  p[fasta_scores[max_score][j].s2] = max_score;
	  q[fasta_scores[max_score][j].s2].key = 
	    fasta_scores[max_score][j].score;
	}
      }
    }

    fprintf(stderr,"Spanning tree for connected component %d (%d nodes)\n", 
	    i, component_size[i]);
    CA(assembly_order[i], component_size[i], sizeof(asm_order_t));
    for(j=0;j<component_size[i];j++) {
      assembly_order[i][j].s1 = p[components[i][j]];
      assembly_order[i][j].s2 = components[i][j];
      assembly_order[i][j].score = q[components[i][j]].key;
    }
    qsort(assembly_order[i], component_size[i], sizeof(asm_order_t),
	  asm_order_compare);
  }

  free(allowed);
  free(used);
  free(q);
  free(p);
}
