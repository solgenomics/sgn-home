#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <limits.h>
#include <time.h>
#include <math.h>

#include "ka.h"
#include "log_message.h"


extern seq_t **sequences;
extern int n_seq;
extern uchar **seq_names;

extern fasta_t **fasta_scores;
extern int *n_fasta;

extern int **components;
extern int n_components;
extern int *component_size;

extern asm_order_t **assembly_order;

#undef PUSH
#define PUSH(a,l,t) \
if ((l) % 10 == 0)   { \
   (a) = realloc((a), (t)*((l) + 10)); \
   if ((a) == NULL) {                     \
      logmsg(MSG_FATAL,"Failed allocating memory (%d bytes)\n", \
             (t)*(l + 10));  \
} \
memset((void *)(a) + (t)*(l), 0x0, 10*(t)); \
}


typedef struct {
  uchar base;
  int n_align;
  int *aligned_nodes;
  int n_seq;
  int *seq;
  int *seq_pos;

  /* These items used to optimize alignemnt algorithm. They do not provide
     information about the structure PO-MSA. */
  int gap_score;
  int gap_node;
} node_t; 

typedef struct edge {
  int node;
  struct edge *next;
} edge_t;

typedef struct {
  int n_nodes;
  int n_seq;
  int *seq;
  /* Array <n_nodes> of nodes */
  node_t *nodes;

  /* Array <n_nodes> of pointers to edge linked lists. NULL in this array
     means no edges. */
  edge_t **out_edges;
  edge_t **in_edges;
} po_t;


typedef struct list_element {
  int value;
  struct list_element *next;
} list_element_t;

typedef struct {
  int size;
  list_element_t *head;
  list_element_t *tail;
} list_t;


inline list_t *new_list(void) {
  list_t *l;

  l = malloc(sizeof(list_t));
  l->size = 0;
  l->head = l->tail = NULL;

  return l;
}

inline void push(list_t *list, int value) {
  list_element_t *l;

  if (list->size == 0) {
    l = list->head = list->tail = malloc(sizeof(list_element_t));
  } else {
    l = list->tail->next = malloc(sizeof(list_element_t));
  }
  l->value = value;
  l->next = NULL;
  list->tail = l;
  list->size++;
}

inline int shift(list_t *list) {
  int value;
  list_element_t *t;

  if (list->size == 0) {
    logmsg(MSG_ERROR,"Can't shift on empty list.\n");
    return 0;
  }
  value = list->head->value;
  t = list->head;
  list->head = list->head->next;
  list->size--;
  if (list->size == 0) {
    list->tail = NULL;
  }

  free(t);
  return value;
}

inline void delete_list(list_t *list) {
  list_element_t *last_p, *p;

  p = list->head;
  while(p!=NULL) {
    last_p = p;
    p = p->next;
    free(last_p);
  }

  free(list);
}

static void debug_alignlists(po_t *po) {
  int i,j,k;
  int *un;
  
  CA(un, po->n_nodes, sizeof(int));
  for(i=0;i<po->n_nodes;i++) {
    if (po->nodes[i].n_align > 0) {
      for(j=0;j<po->n_nodes;j++) {
	un[j] = 0;
      }
      for(j=0;j<po->nodes[i].n_align;j++) {
	un[po->nodes[i].aligned_nodes[j]] = 1;
      }
      un[i] = 1;
      for(j=0;j<po->nodes[i].n_align;j++) {
	int al_node;

	al_node = po->nodes[i].aligned_nodes[j];
	for(k=0;k<po->nodes[al_node].n_align;k++) {
	  if (un[po->nodes[al_node].aligned_nodes[k]] == 0) {
	    fprintf(stderr,"Aligned node %d does not have %d\n",
		    al_node, po->nodes[al_node].aligned_nodes[k]);
	  }
	}
      }
    }
  }
  free(un);
}

typedef struct {
  int bt_row;
  int bt_col;
  int dir;
} bt_t;

typedef struct {
  int gap_node;
  int gap_score;
} gap_t;

po_t *poa_align(po_t *a, po_t *b) {
  po_t *r;
  bt_t **bt;
  bt_t *backtrace;
  int **m;
  int n_rows, n_cols, i, j;
  int max_score;
  int *aligned_nodes, *matched_nodes, *translate;
  int n;
  int max_gapscore, gap_node;
  gap_t **gapa_matrix;
  gap_t *gapa_mat;
  clock_t start_time, tclock, gapa_clock, gapb_clock;
  int maxi, maxj;
  int *matrix;

  n_rows  = a->n_nodes + 1;
  n_cols   = b->n_nodes + 1;


  CA(matrix, n_rows*n_cols, sizeof(int));
  CA(backtrace, n_rows*n_cols, sizeof(bt_t));
  CA(gapa_mat, n_rows*n_cols, sizeof(gap_t));
  m = alloca(n_rows*sizeof(int *));
  bt = alloca(n_rows*sizeof(int *));
  gapa_matrix = alloca(n_rows*sizeof(int *));
  for(i=0;i<n_rows;i++) {
    //    MA(m[i], n_cols, sizeof(int));
    m[i] = matrix + i*n_cols;
    bt[i] = backtrace + i*n_cols;
    gapa_matrix[i] = gapa_mat + i*n_cols;
    //CA(bt[i], n_cols, sizeof(bt_t));
    //CA(gapa_matrix[i], n_cols, sizeof(gap_t));
  }
  for(j=1;j<n_cols;j++) {
    bt[0][j].dir = 2;
    /*bt[0][j].bt_row = 0;
    bt[0][j].bt_col = 0;
    m[0][j] = 0;
    gapa_matrix[0][j].gap_score = 0;
    gapa_matrix[0][j].gap_node = 0;*/
  }
  for(i=1;i<n_rows;i++) {
    bt[i][0].dir = 1;
    /*    bt[i][0].bt_row = 0;
    bt[i][0].bt_col = 0;
    m[i][0] = 0;*/
  }
  bt[0][0].dir = 0;
  m[0][0] = 0;

  start_time = clock();
  gapa_clock = gapb_clock = 0;
  fprintf(stderr,"Starting PO Alignment\n");
  for(i=1;i<n_rows;i++) {
    for(j=1;j<n_cols;j++) {
      int d_move;
      int bt_row, bt_col, bt_dir;
      edge_t *a_inedge, *b_inedge, *edge;
      list_t *gray_list, *next_list;
      int gap_length, node;
      int gap_score;
      
      /* For an on diagonal move */
      if (a->nodes[i-1].base == b->nodes[j-1].base) {
	d_move = 2;
      } else {
	d_move = -4;
      }

      /* Check all possible paths leading to an on diagonal match or mismatch
	 We keep the maximum for comparison with gap options */
      max_score = INT_MIN;
      bt_dir = 3;
      if (a->in_edges[i-1] == NULL || b->in_edges[j-1] == NULL) {
	max_score = d_move;
	bt_dir = 0;
	if (a->in_edges[i-1] == NULL) {
	  bt_row = 0;
	} else {
	  bt_row = a->in_edges[i-1]->node;
	}
	if (b->in_edges[j-1] == NULL) {
	  bt_col = 0;
	} else {
	  bt_col = b->in_edges[j-1]->node;
	}
      } else {
	a_inedge = a->in_edges[i-1];
	while(a_inedge != NULL) {
	  b_inedge = b->in_edges[j-1];
	  while(b_inedge != NULL) {
	    int node_row, node_col;
	    node_row = a_inedge->node + 1;
	    node_col = b_inedge->node + 1;

#ifdef DEBUG
#warning debug code is defined
	    if (node_row > i || node_col > j ) {
	      fprintf(stderr,"%d %d %d %d\n",i, node_row,j, node_col);
	    }
#endif
	    if (m[node_row][node_col] > max_score) {
	      bt_row = node_row;
	      bt_col = node_col;
	      max_score = m[node_row][node_col];
	      bt_dir = 0;
	    }
	    b_inedge = b_inedge->next;
	  }
	  a_inedge = a_inedge->next;
	}
	max_score += d_move;
      }

      tclock = clock();
      edge = a->in_edges[i-1];
      max_gapscore = 0;
      gap_node = -1;
      if (b->out_edges[j-1] != NULL) {
	while(edge!=NULL) {
	  if (gapa_matrix[edge->node][j].gap_score > max_gapscore) {
	    max_gapscore = gapa_matrix[edge->node][j].gap_score;
	    gap_node = edge->node;
	  }
	  edge = edge->next;
	}
	gapa_matrix[i-1][j].gap_node = gap_node;
	gapa_matrix[i-1][j].gap_score = max_gapscore - 4;
	max_gapscore -= 16;
	if (max_gapscore > max_score) {
	  max_score = max_gapscore;
	  bt_dir = 1;
	  bt_row = gap_node + 1;
	  bt_col = j;
	}
      } else {
	while(edge!=NULL) {
	  if (m[edge->node + 1][j] > max_gapscore) {
	    max_gapscore = m[edge->node + 1][j];
	    gap_node = edge->node;
	  }
	  edge = edge->next;
	}
	if (max_gapscore > max_score) {
	  bt_dir = 1;
	  bt_row = gap_node + 1;
	  bt_col = j;
	  max_score = max_gapscore;
	}
      }
      gapa_clock += (clock() - tclock);
      tclock = clock();
      edge = b->in_edges[j-1];
      gap_node = -1;
      max_gapscore = 0;
      /* No gap penalty for gaps introduced where sequences don't overlap */
      if (a->out_edges[i-1] != NULL) {
	while(edge!=NULL) {
	  if (b->nodes[edge->node].gap_score > max_gapscore) {
	    max_gapscore = b->nodes[edge->node].gap_score;
	    gap_node = edge->node;
	  }
	  edge = edge->next;
	}
	b->nodes[j-1].gap_node = gap_node;
	b->nodes[j-1].gap_score = max_gapscore - 4;
	max_gapscore -= 16;
	if (max_gapscore > max_score) {
	  max_score = max_gapscore;
	  bt_dir = 2;
	  bt_row = i;
	  bt_col = gap_node + 1;
	}
      } else {
	while(edge!=NULL) {
	  if (m[i][edge->node + 1] > max_gapscore) {
	    max_gapscore = m[i][edge->node + 1];
	    gap_node = edge->node;
	  }
	  edge = edge->next;
	}
	if (max_gapscore > max_score) {
	  bt_dir = 2;
	  bt_row = i;
	  bt_col = gap_node + 1;
	  max_score = max_gapscore;
	}
      }
      gapb_clock += (clock() - tclock);

#ifdef DEBUG
      if (bt_dir == 3) {
	abort();
	bt_row = i-1;
	bt_col = j-1;
	bt_dir = 0;
      }
#endif

      m[i][j] = max_score;
      bt[i][j].bt_row = bt_row;
      bt[i][j].bt_col = bt_col;
      bt[i][j].dir = bt_dir;

      if (bt_dir != 2) {
	b->nodes[j-1].gap_score = max_score;
	b->nodes[j-1].gap_node = j-1;
      }
      if (bt_dir != 1) {
	gapa_matrix[i-1][j].gap_score = max_score;
	gapa_matrix[i-1][j].gap_node = i-1;
      }
    }
  }
  fprintf(stderr,"Alignment finished (%1.2f seconds)\n",(clock()-start_time)/ (float) CLOCKS_PER_SEC);
  fprintf(stderr,"GapA (%1.2f seconds) GapB (%1.2f seconds)\n",
	  gapa_clock/(float) CLOCKS_PER_SEC, gapb_clock/(float) CLOCKS_PER_SEC);
  start_time = clock();
  /* Back track here */
  fprintf(stderr,"Alignment score: %d\n",m[n_rows-1][n_cols-1]);

  if (0) {
    int p,k;
    int maxi, maxj;
    uchar *query, *sbjct, *align;
    query = calloc(n_rows+n_cols+2, sizeof(uchar));
    sbjct = calloc(n_rows+n_cols+2, sizeof(uchar));
    align = calloc(n_rows+n_cols+2, sizeof(uchar));
    p = n_rows+n_cols;
    i = n_rows-1;
    j = n_cols-1;
#if 1
    maxi = n_rows - 1;
    maxj = n_cols - 1;
    for(i=0;i<n_rows;i++) {
      if (m[i][n_cols-1] > m[maxi][maxj]) {
	maxi = i;
	maxj = n_cols-1;
      }
    }
    for(j=0;j<n_cols;j++) {
      if (m[n_rows-1][j] >= m[maxi][maxj]) {
	maxi = n_rows-1;
	maxj = j;
      }
    }
    i = maxi;
    j = maxj;
#endif
    fprintf(stderr,"Starting back-track: score = %d %d %d\n",m[i][j],i,j);
    while(i>0 || j>0) {
      int t, nexti, nextj;
      switch(bt[i][j].dir) {
      case 0:
	assert(p>=0);
	if(a->nodes[i-1].base == b->nodes[j-1].base) {
	  align[p] = '|';
	} else {
	  align[p] = ' ';
	}
	query[p] = a->nodes[i-1].base;
	sbjct[p] = b->nodes[j-1].base;
	nexti = bt[i][j].bt_row;
	nextj = bt[i][j].bt_col;
	i = nexti;
	j = nextj;
	p--;
	break;
      case 1:
	t = bt[i][j].bt_row;
	while(i>t) {
	  assert(p>=0);
	  query[p] = a->nodes[--i].base;
	  sbjct[p] = '-';
	  align[p] = ' ';
	  p--;
	}
	break;
      case 2:
	t = bt[i][j].bt_col;
	while(j>t) {
	  assert(p>=0);
	  query[p] = '-';
	  sbjct[p] = b->nodes[--j].base;
	  align[p] = ' ';
	  p--;
	}
	break;
      case 3:
	fprintf(stderr,"Alignment end at %d %d\n",i,j);
	break;
      default:
	fprintf(stderr,"Unknown back-track code %d\n",bt[i][j].dir);
	break;
      }
    }
    p += 1;
    fprintf(stderr,"Alignment stoped at %d %d\n",i,j);
    while(p<(n_rows + n_cols)) {
      fprintf(stderr,"Query %-65.65s\n",query + p);
      fprintf(stderr,"      %-65.65s\n",align + p);
      fprintf(stderr,"Sbjct %-65.65s\n\n",sbjct + p);
      p += 65;
    }
    free(query);
    free(align);
    free(sbjct);
  }
 
  MA(r, sizeof(po_t));
  r->n_nodes = a->n_nodes + b->n_nodes;
  r->nodes = NULL;
  r->n_seq = 0;
  r->seq = NULL;
  MA(r->seq, (a->n_seq+b->n_seq)*sizeof(int));
  for(i=0;i<a->n_seq;i++) {
    r->seq[r->n_seq++] = a->seq[i];
  }
  for(i=0;i<b->n_seq;i++) {
    r->seq[r->n_seq++] = b->seq[i];
  }


  /* Discover the number of nodes we'll have in the resulting MSA:

     This should be the sum of the nodes in the two input MSAs, minus the
     nodes which are "aligned" with identical bases. */    
  aligned_nodes = malloc(sizeof(int)*b->n_nodes);
  matched_nodes = malloc(sizeof(int)*b->n_nodes);
  translate = malloc(sizeof(int)*b->n_nodes);
  for(i=0;i<b->n_nodes;i++) {
    aligned_nodes[i] = matched_nodes[i] = translate[i] = -1;
  }

  maxi = n_rows - 1;
  maxj = n_cols - 1;
  for(i=0;i<n_rows;i++) {
    if (m[i][n_cols-1] > m[maxi][maxj]) {
      maxi = i;
      maxj = n_cols-1;
    }
  }
  for(j=0;j<n_cols;j++) {
    if (m[n_rows-1][j] >= m[maxi][maxj]) {
      maxi = n_rows-1;
      maxj = j;
    }
  }
  i = maxi;
  j = maxj;
  while(i > 0 || j > 0) {
    int nexti, nextj;
    nexti = bt[i][j].bt_row;
    nextj = bt[i][j].bt_col;
    if (bt[i][j].dir == 0) {
      aligned_nodes[j-1] = i-1;
      if (a->nodes[i-1].base == b->nodes[j-1].base) {
	matched_nodes[j-1] = i-1;
	r->n_nodes --;
      }
    } 
    i = nexti;
    j = nextj;
  }
  n = a->n_nodes;
  for(i=0;i<b->n_nodes;i++) {
    if (matched_nodes[i] == -1) {
      translate[i] = n++;
    } else {
      translate[i] = matched_nodes[i];
    }
  }


  fprintf(stderr,"%d nodes needed for merged MSA\n",r->n_nodes);
  r->nodes = calloc(r->n_nodes, sizeof(node_t));
  r->in_edges = calloc(r->n_nodes, sizeof(edge_t *));
  r->out_edges = calloc(r->n_nodes, sizeof(edge_t *));

  /* Create of deep copy of PO a */
  for(i=0;i<a->n_nodes;i++) {
    node_t *copy, *src;
    edge_t *edge;
    copy = r->nodes + i;
    src = a->nodes + i;

    copy->base = src->base;
    for(j=0;j<src->n_align;j++) {
      PUSH(copy->aligned_nodes, copy->n_align, sizeof(int));
      copy->aligned_nodes[copy->n_align++] = src->aligned_nodes[j];
    }
    for(j=0;j<src->n_seq;j++) {
      PUSH(copy->seq, copy->n_seq, sizeof(int));
      copy->seq[copy->n_seq] = src->seq[j];
      PUSH(copy->seq_pos, copy->n_seq, sizeof(int));
      copy->seq_pos[copy->n_seq] = src->seq_pos[j];
      copy->n_seq++;
    }

    edge = a->in_edges[i];
    while(edge!=NULL) {
      edge_t *t;
      t = malloc(sizeof(edge_t));
      t->node = edge->node;
      t->next = r->in_edges[i];
      r->in_edges[i] = t;
      edge = edge->next;
    }
    
    edge = a->out_edges[i];
    while(edge!=NULL) {
      edge_t *t;
      t = malloc(sizeof(edge_t));
      t->node = edge->node;
      t->next = r->out_edges[i];
      r->out_edges[i] = t;
      edge = edge->next;
    }
  }

  for(i=0;i<b->n_nodes;i++) {
    node_t *node;
    edge_t *edge;

    node = r->nodes + translate[i];
    if (matched_nodes[i] == -1) {
      node->base = b->nodes[i].base;
    }

    /* Append list of sequences containing this node in PO b  */
    for(j=0;j<b->nodes[i].n_seq;j++) {
      PUSH(node->seq, node->n_seq, sizeof(int));
      PUSH(node->seq_pos, node->n_seq, sizeof(int));
      node->seq[node->n_seq] = b->nodes[i].seq[j];
      node->seq_pos[node->n_seq] = b->nodes[i].seq_pos[j];
      node->n_seq++;
    }

    // Should be superceded by below code
    /* Append nodes which are aligned in PO b to this matched node */
    for(j=0;j<b->nodes[i].n_align;j++) {
      PUSH(node->aligned_nodes, node->n_align, sizeof(int));
      node->aligned_nodes[node->n_align++] = 
	translate[b->nodes[i].aligned_nodes[j]];
    }
    
    /* Append the in/out edges for this node in PO b, translated to the
       new node indicies for the resultant PO */
    edge = b->in_edges[i];
    while(edge!=NULL) {
      edge_t *t;

      /* Check if the edge already exists. Only add an edge if it doesn't
         exist. */
      t = r->in_edges[translate[i]];
      while(t!=NULL) {
	if (t->node == translate[edge->node]) break;
	t = t->next;
      }
      if (t == NULL) {
	t = malloc(sizeof(edge_t));
	t->node = translate[edge->node];
	t->next = r->in_edges[translate[i]];
	r->in_edges[translate[i]] = t;
      }
      edge = edge->next;
    }

    edge = b->out_edges[i];
    while(edge!=NULL) {
      edge_t *t;
      t = r->out_edges[translate[i]];
      while(t!=NULL) {
	if (t->node == translate[edge->node]) break;
	t = t->next;
      }
      if (t == NULL) {
	t = malloc(sizeof(edge_t));
	t->node = translate[edge->node];
	t->next = r->out_edges[translate[i]];
	r->out_edges[translate[i]] = t;
      }
      edge = edge->next;
    }

    /* In this case, we have aligned two nodes from PO <a> and PO <b> but
       they are not identical (mismatch) so the nodes were not fused. 
       
       Here we must find the union of aligned nodes and rebuild the 
       aligned_nodes list for each member */
    if (1 || aligned_nodes[i] != -1) {
      node_t *node_a, *node_b, *t_node;
      int *un;
      int *al_list;
      int k;
      int n_align;

 
      un = calloc(r->n_nodes, sizeof(int));
      if (aligned_nodes[i] != -1) {
	node_a = r->nodes + aligned_nodes[i];
	for(j=0;j<node_a->n_align;j++) {
	  un[node_a->aligned_nodes[j]] = 1;
	}
	un[aligned_nodes[i]] = 1;
      }
      node_b = r->nodes + translate[i];


      for(j=0;j<node_b->n_align;j++) {
	un[node_b->aligned_nodes[j]] = 1;
      }
      un[translate[i]] = 1;
      
      al_list = NULL;
      n_align = 0;
      for(j=0;j<r->n_nodes;j++) {
	if (un[j] == 0) continue;
	PUSH(al_list, n_align, sizeof(int));
	al_list[n_align++] = j;
      }

      for(j=0;j<n_align;j++) {
	t_node = r->nodes + al_list[j];
	for(k=0;k<t_node->n_align;k++) {
	  if (un[t_node->aligned_nodes[k]] == 0) {
	    un[t_node->aligned_nodes[k]] = 1;
	    PUSH(al_list, n_align, sizeof(int));
	    al_list[n_align++] = t_node->aligned_nodes[k];
	  }
	}
      }

      for(j=0;j<n_align;j++) {
#if 0
	if (n_align - r->nodes[al_list[j]].n_align == 1) {
	  int xor;
	  /* Tricks: al_list and r->nodes[al_list[j]].aligned_nodes should
	     be the same, except r->nodes[al_list[j]] should be missing
	     al_list[j] (itself). Thus, the cummulative XOR of both lists
	     should result in al_list[j] if this is true. It is possible
	     xor = al_list[j] by chance if the lists are different, but that
	     is minimized by requring the lists to be the same length. The
	     chance of of it happening is 1/r->n_nodes otherwise. */
	  xor = 0;
	  for(k=0;k<r->nodes[al_list[j]].n_align;k++) {
	    xor ^= r->nodes[al_list[j]].aligned_nodes[k];
	  }
	  for(k=0;k<n_align;k++) {
	    xor ^= al_list[k];
	  }
	  if (xor == al_list[j]) { 
	    //	    fprintf(stderr,"XOR = 0, skipping aligned list rebuild\n");
	    continue;
	  }
	}
#endif
	if (r->nodes[al_list[j]].aligned_nodes)
	  free(r->nodes[al_list[j]].aligned_nodes);
	r->nodes[al_list[j]].n_align = 0;
	r->nodes[al_list[j]].aligned_nodes = NULL;

	t_node = r->nodes+ al_list[j];
	for(k=0;k<n_align;k++) {
	  if (k == j) continue;
	  PUSH(t_node->aligned_nodes, t_node->n_align, sizeof(int));
	  t_node->aligned_nodes[t_node->n_align++] = al_list[k];
	}
      }
      free(al_list);
      free(un);
    }
  }
  debug_alignlists(r);
  free(aligned_nodes);
  free(matched_nodes);
  free(translate);

  /*
  for(i=0;i<n_rows;i++) {
    //free(m[i]);
    free(bt[i]);
    free(gapa_matrix[i]);
  }
  */
  free(gapa_mat);
  free(matrix);
  free(backtrace);

  fprintf(stderr,"Merge time %1.2f seconds\n",(clock()-start_time)/(float) CLOCKS_PER_SEC);
  return r;
}

static po_t *make_po(int seq_index) {
  po_t *po;
  int i;
  seq_t *seq;

  seq = sequences[seq_index];

  MA(po, sizeof(po_t));
  MA(po->nodes, sizeof(node_t)*seq->length);
  CA(po->out_edges, seq->length, sizeof(edge_t *));
  CA(po->in_edges, seq->length, sizeof(edge_t *));
  po->n_nodes = seq->length;
  po->n_seq = 0;
  po->seq = NULL;
  PUSH(po->seq, po->n_seq, sizeof(int));
  po->seq[po->n_seq++] = seq_index;

  for(i=0;i<seq->length;i++) {
    po->nodes[i].base = seq->seqstr[i];
    po->nodes[i].seq = NULL;
    po->nodes[i].seq_pos = NULL;

    PUSH(po->nodes[i].seq, 0, sizeof(int));
    po->nodes[i].seq[0] = seq_index;
    PUSH(po->nodes[i].seq_pos, 0, sizeof(int));
    po->nodes[i].seq_pos[0] = i;
    po->nodes[i].n_seq = 1;

    po->nodes[i].n_align = 0;
    po->nodes[i].aligned_nodes = NULL;

    if (i > 0) {
      po->in_edges[i] = malloc(sizeof(edge_t));
      po->in_edges[i]->node = i-1;
      po->in_edges[i]->next = NULL;
      po->out_edges[i-1] = malloc(sizeof(edge_t));
      po->out_edges[i-1]->node = i;
      po->out_edges[i-1]->next = NULL;
    }
  }  

  return po;
}

void display_po(po_t *po, FILE *f) {
  int i,j,k,l;
  int *used;
  int *seq_map;
  int *seq_rmap;
  int map_id;
  uchar **seq;
  node_t *node;
  int *highlight;

  seq_map = malloc(n_seq*sizeof(int));
  seq_rmap = malloc(po->n_seq*sizeof(int));
  for(i=0;i<n_seq;i++)
    seq_map[i] = -1;

  map_id = 0;
  seq = malloc(sizeof(uchar *)*po->n_seq);
  for(i=0;i<po->n_seq;i++) {
    seq[i] = calloc(po->n_nodes+1, sizeof(uchar));
    memset(seq[i], '-', po->n_nodes);
  }
  used = calloc(po->n_nodes, sizeof(int));
  highlight = calloc(po->n_nodes, sizeof(int));

  k = 0;
  for(i=0;i<po->n_nodes;i++) {
    if (used[i]) continue;
    node = po->nodes + i;
    used[i] = 1;
    for(j=0;j<node->n_seq;j++) {
      if (seq_map[node->seq[j]] == -1) {
	seq_rmap[map_id] = node->seq[j];
	seq_map[node->seq[j]] = map_id++;
      }
      if (sequences[node->seq[j]]->qual[node->seq_pos[j]] < 12) {
	seq[seq_map[node->seq[j]]][k] = node->base + 32;
      } else {
	seq[seq_map[node->seq[j]]][k] = node->base ;
      }
    }
    
    for(j=0;j<node->n_align;j++) {
      node_t *t_node;
      highlight[k] = 1;
      assert(used[node->aligned_nodes[j]] == 0);
      used[node->aligned_nodes[j]] = 1;
      t_node = po->nodes + node->aligned_nodes[j];
      for(l=0;l<t_node->n_seq;l++) {
	if (seq_map[t_node->seq[l]] == -1) {
	  assert(map_id < po->n_seq);
	  seq_rmap[map_id] = t_node->seq[l];
	  seq_map[t_node->seq[l]] = map_id++;
	}
	if (sequences[t_node->seq[l]]->qual[t_node->seq_pos[l]] < 12) {
	  seq[seq_map[t_node->seq[l]]][k] = t_node->base + 32;
	} else {
	  seq[seq_map[t_node->seq[l]]][k] = t_node->base;
	}	  
      }
    }
    
    k++;
  }

  /* Add null-string terminator */
  for(i=0;i<po->n_seq;i++) {
    seq[i][k] = 0;
  }

  i = 0;
  while(i<k) {
    // for(j=0;j<po->n_seq;j++) {
    // fprintf(f,"%-15.15s   %-60.60s\n",seq_names[seq_rmap[j]],seq[j] + i);
    //}
    for(j=0;j<po->n_seq;j++) {
      fprintf(f,"%-20.20s   ",sequences[seq_rmap[j]]->label);
      for(l=i;l<i+105;l++) {
	if (l>=k) break;
	if (highlight[l]) {
	  fprintf(f,"%c%c01;31m",27,91);
	  putc(seq[j][l], f);
	  fprintf(f,"%c%c00m",27,91);
	} else {
	  putc(seq[j][l], f);
	}
      }
      fprintf(f,"\n");
    }

    fprintf(f,"\n");
    i += 105;
  }


  for(i=0;i<po->n_seq;i++) {
    free(seq[i]);
  }
  free(seq);
  free(seq_map);
  free(seq_rmap);
  free(used);
  free(highlight);
  fflush(f);
}

static void po_dfs_visit(po_t *src, int *o, int *on, 
			 int *d, int *f, int *time, int n) {
  edge_t *edge;

  d[n] = (*time)++;
  edge = src->out_edges[n]; 
#if 0
  for(i=0;i<src->nodes[n].n_align;i++) {
    if (d[src->nodes[n].aligned_nodes[i]] == -1)
      po_dfs_visit(src, o, on, d, f, time, src->nodes[n].aligned_nodes[i]);
  }
#endif
  while(edge!=NULL) {
    if (d[edge->node] == -1) po_dfs_visit(src, o, on, d, f, time, edge->node);
    edge = edge->next;
  }
  f[n] = (*time)++;
  o[(*on)--] = n;
}

static void po_topological_sort(po_t *src) {
  int *d, *f, *o, *r, *m;
  int on, time, i, j;
  node_t *temp_nodes;
  edge_t **temp_inedge;
  edge_t **temp_outedge;

  d = malloc(src->n_nodes*sizeof(int));
  f = malloc(src->n_nodes*sizeof(int));
  o = calloc(src->n_nodes, sizeof(int));
  r = calloc(src->n_nodes, sizeof(int));
  m = calloc(src->n_nodes, sizeof(int));
  on = src->n_nodes - 1;
  time = 0;

  for(i=0;i<src->n_nodes;i++) {
    d[i] = f[i] = r[i] = -1;
  }
 
  for(i=0;i<src->n_nodes;i++) {
    if(d[i] == -1) po_dfs_visit(src, o, &on, d, f, &time, i);
  }
  assert(on == -1);
#if 0
  /* At this stage we have *a* topological sorting of the graph. Note
     carefully that many such orderings are possible. However, we would like
     the last base in the MSA to actually be the last node in the graph. Due
     to the un-ordered nature of the adjacency lists, this may not be what
     we get. So*/
  time = 0;
  for(i=0;i<src->n_nodes;i++) {
    if (r[o[i]] != -1) continue;
    m[time] = o[i];
    r[o[i]] = time++;
    for(j=0;j<src->nodes[o[i]].n_align;j++) {
      //assert( == -1);
      m[time] = src->nodes[o[i]].aligned_nodes[j];
      r[src->nodes[o[i]].aligned_nodes[j]] = time++;
    }
  }
  memcpy(o,m,sizeof(int)*src->n_nodes);
#endif
  /* Construct reverse mapping */
  for(i=0;i<src->n_nodes;i++)
    r[o[i]] = i;
  

  /* Rebuild the node array in topological sort order */
  temp_nodes = malloc(sizeof(node_t)*src->n_nodes);
  temp_inedge = malloc(sizeof(edge_t *)*src->n_nodes);
  temp_outedge = malloc(sizeof(edge_t *)*src->n_nodes);
  for(i=0;i<src->n_nodes;i++) {
    edge_t *edge;

    temp_nodes[i] = src->nodes[o[i]];
    temp_inedge[i] = src->in_edges[o[i]];
    temp_outedge[i] = src->out_edges[o[i]];
    for(j=0;j<temp_nodes[i].n_align;j++) {
      temp_nodes[i].aligned_nodes[j] = r[temp_nodes[i].aligned_nodes[j]];
    }
    edge = temp_inedge[i];
    while(edge!=NULL) {
      edge->node = r[edge->node];
      edge = edge->next;
    }
    edge = temp_outedge[i];
    while(edge!=NULL) {
      edge->node = r[edge->node];
      edge = edge->next;
    }
  }
  free(src->nodes);
  free(src->in_edges);
  free(src->out_edges);
  src->nodes = temp_nodes;
  src->in_edges = temp_inedge;
  src->out_edges = temp_outedge;
    

  free(d);
  free(f);
  free(o);
  free(r);
  free(m);
}

static int disjoint_members(node_t *a, node_t *b) {
  int i,j;

  for(i=0;i<a->n_seq;i++) {
    for(j=0;j<b->n_seq;j++) {
      if (a->seq[i] == b->seq[j]) return 0;
    }
  }

  return 1;
}

static void delete_edgelist(edge_t *edge) {
  edge_t *next;

  while(edge!=NULL) {
    next = edge->next;
    free(edge);
    edge = next;
  }

}

static void po_fuseall(po_t *src) {
  int i,j,dead_node,k;
  edge_t *edge, *t_edge, *newedge;
  node_t *node1, *node2;
  edge_t *incoming, *outgoing, *last;

  for(i=0;i<src->n_nodes;i++) {
    edge = src->out_edges[i];
    while(edge != NULL) {
      t_edge = edge->next;
      while(t_edge != NULL) {
	/* Look for another edge which leads to a node with the same base */
	if (src->nodes[t_edge->node].base == src->nodes[edge->node].base) {
	  node1 = src->nodes + edge->node;
	  node2 = src->nodes + t_edge->node;
	  assert(node1 != node2);

	  if (disjoint_members(node1, node2)) {
	    fprintf(stderr,"Attempting to fuse nodes\n");
	    /* t_edge is about to get wiped out, advance now and remember the
	       node # */
	    dead_node = t_edge->node;
	    t_edge = t_edge->next;

	  /* We are going to fuse node2 to node1 and unlink node2. First, 
	     append sequence membership to node1. */
	  for(j=0;j<node2->n_seq;j++) {
	    PUSH(node1->seq, node1->n_seq, sizeof(int));
	    node1->seq[node1->n_seq] = node2->seq[j];
	    PUSH(node1->seq_pos, node1->n_seq, sizeof(int));
	    node1->seq_pos[node1->n_seq] = node2->seq_pos[j];
	    node1->n_seq++;
	  }

	  /* Now, we redirect all incoming edges for node2 to node1 */
	  incoming = src->in_edges[dead_node];
	  while(incoming != NULL) {
	    /* Scan this prior node for an edge that already leads to the
	       fused node */
	    outgoing = src->out_edges[incoming->node];
	    while(outgoing!=NULL) {
	      if (outgoing->node == edge->node) break;
	      outgoing = outgoing->next;
	    }
	    /* If there isn't one, add one */
	    if (outgoing == NULL) {
	      MA(newedge, sizeof(edge_t));
	      newedge->next = src->out_edges[incoming->node];
	      newedge->node = edge->node;
	      src->out_edges[incoming->node] = newedge;
	    }
	    /* Now get rid of the one leading to the node we are fusing */
	    last = outgoing = src->out_edges[incoming->node];
	    while(outgoing!=NULL) {
	      if (outgoing->node == dead_node) {
		if (last == outgoing) {
		  src->out_edges[incoming->node] = outgoing->next;
		  free(last);
		} else {
		  last->next = outgoing->next;
		  free(outgoing);
		}
		break;
	      }
	      outgoing = outgoing->next;
	    }
	    assert(outgoing != NULL);
	    /* Now add the incoming edge to node1's incoming edges */
	    MA(newedge, sizeof(edge_t));
	    newedge->node = incoming->node;
	    newedge->next = src->in_edges[edge->node];
	    src->in_edges[edge->node] = newedge;
	    
	    incoming = incoming->next;
	  }

	  /* Redirect all outgoing edges from node2 to come from node1 */
	  outgoing = src->out_edges[dead_node];
	  while(outgoing != NULL) {
	    /* Scan this following node for an edge that already comes from
	       our fused node */
	    incoming = src->in_edges[outgoing->node];
	    while(incoming!=NULL) {
	      if (incoming->node == edge->node) break;
	      incoming = incoming->next;
	    }
	    if (incoming == NULL) {
	      MA(newedge, sizeof(edge_t));
	      newedge->next = src->in_edges[outgoing->node];
	      newedge->node = edge->node;
	      src->in_edges[outgoing->node] = newedge;
	    }
	    last = incoming = src->in_edges[outgoing->node];
	    while(incoming != NULL) {
	      if (incoming->node == dead_node) {
		if (last == incoming) {
		  src->in_edges[outgoing->node] = incoming->next;
		  free(last);
		} else {
		  last->next = incoming->next;
		  free(incoming);
		}
		break;
	      }
	      incoming = incoming->next;
	    }
	    assert(incoming!=NULL);

	    /* Finally, add the outgoing edge to node1's outgoing edges */
	    MA(newedge, sizeof(edge_t));
	    newedge->node = outgoing->node;
	    newedge->next = src->out_edges[edge->node];
	    src->out_edges[edge->node] = newedge;
	    
	    outgoing = outgoing->next;
	  }

	  /* Now deal with aligned nodes */
	  for(j=0;j<node2->n_align;j++) {
	    for(k=0;k<node1->n_align;k++) {
	      if (node1->aligned_nodes[k] == node2->aligned_nodes[j]) {
		fprintf(stderr,"Shared aligned node %d (%d) (%d)\n",
			node1->aligned_nodes[k],node1->n_align, node2->n_align);
		break;
	      }
	    }
	    if (k == node1->n_align) {
	      PUSH(node1->aligned_nodes, node1->n_align, sizeof(int));
	      node1->aligned_nodes[node1->n_align++] = node2->aligned_nodes[j];
	    }
	  }
	  delete_edgelist(src->in_edges[dead_node]);
	  delete_edgelist(src->out_edges[dead_node]);
	  src->in_edges[dead_node] = NULL;
	  src->out_edges[dead_node] = NULL;
	  src->nodes[dead_node].n_align = 0;
	  free(src->nodes[dead_node].aligned_nodes);
	  src->nodes[dead_node].aligned_nodes = NULL;
	  free(src->nodes[dead_node].seq);
	  free(src->nodes[dead_node].seq_pos);
	  src->nodes[dead_node].seq = NULL;
	  src->nodes[dead_node].seq_pos = NULL;
	  src->nodes[dead_node].n_seq = 0;

	  } else {
	    fprintf(stderr,"Nodes do not have disjoint sequence membership\n");
	    t_edge = t_edge->next;
	  }
	} else {
	  t_edge = t_edge->next;
	}
      }
      edge = edge->next;
    }
  }
}

static void fuse_identical_alignring_nodes(po_t *src) {
  int i,j,k;
  node_t *node, *node2;
  edge_t *edge;
  int next, list_size, do_fuse;

  int *fuse_map;
  int *map;
  int *used;
  int *list;
  
  node_t *temp_nodes;
  edge_t **temp_inedge;
  edge_t **temp_outedge;

  do_fuse = 0;
  MA(fuse_map, sizeof(int)*src->n_nodes);
  MA(map, sizeof(int)*src->n_nodes);
  for(i=0;i<src->n_nodes;i++)
    fuse_map[i] = -1;

  for(i=0;i<src->n_nodes;i++) {
    if (fuse_map[i] != -1) continue;
    if (src->nodes[i].n_align > 0) {
      node = src->nodes + i;
      for(j=0;j<node->n_align;j++) {
	if (src->nodes[node->aligned_nodes[j]].base == node->base) {
	  fuse_map[node->aligned_nodes[j]] = i;
	  fprintf(stderr,"Fusing nodes %d-%d\n",i,node->aligned_nodes[j]);
	  do_fuse = 1;
	}
      }
      for(j=0;j<node->n_align;j++) {
	if (fuse_map[node->aligned_nodes[j]] != -1) continue;
	node2 = src->nodes + node->aligned_nodes[j];
	for(k=j+1;k<node->n_align;k++) {
	  if (src->nodes[node->aligned_nodes[k]].base == node2->base) {
	    fuse_map[node->aligned_nodes[k]] = node->aligned_nodes[j];
	  fprintf(stderr,"Fusing nodes %d-%d\n",node->aligned_nodes[j],
		  node->aligned_nodes[k]);
	    do_fuse = 1;
	  }
	}
      }
    }
  }

  if (! do_fuse) {
    free(map);
    free(fuse_map);
    return;
  }

  next = 0;
  for(i=0;i<src->n_nodes;i++) {
    if (fuse_map[i]==-1) map[i] = next++;
    else 
      map[i] = map[fuse_map[i]];
  }

  MA(temp_nodes, sizeof(node_t)*next);
  MA(temp_inedge, sizeof(edge_t *)*next);
  MA(temp_outedge, sizeof(edge_t *)*next);

  for(i=0;i<src->n_nodes;i++) {
    if (fuse_map[i] == -1) {
      temp_nodes[map[i]] = src->nodes[i];
      temp_inedge[map[i]] = src->in_edges[i];
      temp_outedge[map[i]] = src->out_edges[i];
    }
  }

  /* Add sequences of fused nodes, prepend edge lists. We don't have to
     append align ring information since fused nodes are always in the
     same align ring. */
  CA(used, src->n_nodes, sizeof(int));
  CA(list, src->n_nodes, sizeof(int));
  for(i=0;i<src->n_nodes;i++) {
    if (fuse_map[i] == -1) continue;
    for(j=0;j<src->nodes[i].n_seq;j++) {
      node = temp_nodes + map[fuse_map[i]];
      PUSH(node->seq, node->n_seq, sizeof(int));
      PUSH(node->seq_pos, node->n_seq, sizeof(int));
      node->seq[node->n_seq] = src->nodes[i].seq[j];
      node->seq_pos[node->n_seq] = src->nodes[i].seq_pos[j];
      node->n_seq++;
    }
    edge = src->in_edges[i];
    while(edge != NULL && edge->next != NULL) edge = edge->next;
    if (edge!=NULL) {
      edge->next = temp_inedge[map[fuse_map[i]]];
      temp_inedge[map[fuse_map[i]]] = src->in_edges[i];
    }
    
    edge = src->out_edges[i];
    while(edge != NULL && edge->next != NULL) edge = edge->next;
    if (edge!=NULL) {
      edge->next = temp_outedge[map[fuse_map[i]]];
      temp_outedge[map[fuse_map[i]]] = src->out_edges[i];
    } 
  }    

  /* Now, walk through this mess and change all the referenced node indices
     to their mapped index. Note that we have to check all kinds of crap for
     multiple entries that may have occured due to fusing above */
  memset(used, 0x0, sizeof(int)*next);
  for(i=0;i<next;i++) {
    /* Avoid repeating work for nodes which were processed by previous node
       in the same align ring */
    if (used[i]) continue;
    node = temp_nodes + i;
    if (node->n_align > 0) {
      list_size = 1;
      list[0] = i;
      used[i] = 1;
      for(j=0;j<node->n_align;j++) {
	int n;
	n = map[node->aligned_nodes[j]];
	if (!used[n]) {
	  list[list_size++] = n;
	  used[n] = 1;
	}
      }
      for(j=0;j<list_size;j++) {
	node = temp_nodes + list[j];
	node->n_align = 0;
	free(node->aligned_nodes);
	node->aligned_nodes = NULL;
	for(k=0;k<list_size;k++) {
	  if (k == j) continue;
	  PUSH(node->aligned_nodes, node->n_align, sizeof(int));
	  node->aligned_nodes[node->n_align++] = list[k];
	}
      }
    }
  }

  memset(used,0x0,sizeof(int)*next);
  for(i=0;i<next;i++) {
    
    list_size = 0;
    edge = temp_inedge[i];
    while(edge!=NULL) {
      edge_t *t;

      if (! used[map[edge->node]]) {
	list[list_size++] = map[edge->node];
	used[map[edge->node]] = 1;
      }
      t = edge;
      edge = edge->next;
      free(t);
    } 
    temp_inedge[i] = NULL;
    for(j=0;j<list_size;j++) {
      MA(edge, sizeof(edge_t));
      edge->node = list[j];
      edge->next = temp_inedge[i];
      temp_inedge[i] = edge;
    }
    for(j=0;j<list_size;j++)
      used[list[j]] = 0;

    list_size = 0;
    edge = temp_outedge[i];
    while(edge!=NULL) {
      edge_t *t;

      if (! used[map[edge->node]]) {
	list[list_size++] = map[edge->node];
	used[map[edge->node]] = 1;
      }
      t = edge;
      edge = edge->next;
      free(t);
    }
    for(j=0;j<list_size;j++)
      used[list[j]] = 0;
    

    temp_outedge[i] = NULL;
    for(j=0;j<list_size;j++) {
      MA(edge, sizeof(edge_t));
      edge->node = list[j];
      edge->next = temp_outedge[i];
      temp_outedge[i] = edge;
    }
  }

  free(src->nodes);
  free(src->in_edges);
  free(src->out_edges);

  src->nodes = temp_nodes;
  src->in_edges = temp_inedge;
  src->out_edges = temp_outedge;
  src->n_nodes = next;

  free(list);
  free(used);
  free(fuse_map);
  free(map);
}



void scan_alignlists(po_t *po) {
  int *map;
  int i,j;

  map = calloc(po->n_nodes, sizeof(int));

  for(i=0;i<po->n_nodes;i++) {
    for(j=0;j<po->nodes[i].n_align;j++) {
      if (po->nodes[i].aligned_nodes[j] == i) {
	fprintf(stderr,"Node %d is aligned to itself\n",i);
      }
      map[po->nodes[i].aligned_nodes[j]]++;
    }
    for(j=0;j<po->n_nodes;j++) {
      if (map[j] > 1) {
	fprintf(stderr,"In node %d (base %c) %d is aligned %d times\n",
		i, po->nodes[i].base, j, map[j]);
      }
      map[j] = 0;
    }
  }
}


static void destroy_po(po_t *po) {
  int i;
  edge_t *edge;

  for(i=0;i<po->n_nodes;i++) {
    free(po->nodes[i].aligned_nodes);
    free(po->nodes[i].seq);
    free(po->nodes[i].seq_pos);
    edge = po->in_edges[i];
    while(edge!=NULL) {
      edge_t *t;
      t = edge;
      edge = edge->next;
      free(t);
    }
    edge = po->out_edges[i];
    while(edge!=NULL) {
      edge_t *t;
      t = edge;
      edge = edge->next;
      free(t);
    }
  }
  free(po->in_edges);
  free(po->out_edges);
  free(po->nodes);
  free(po->seq);
  free(po);
}

static int analyze_sgp_maxflow_coverage(po_t *po) {
  int i;
  int next_node, max, node_count;
  edge_t *edge;

  i = 0;
  node_count = 0;
  while(po->out_edges[i]!=NULL) {
    if (po->out_edges[i]->next != NULL) {

      next_node = po->out_edges[i]->node;
      edge = po->out_edges[i]->next;
      max = po->nodes[next_node].n_seq;
      while(edge!=NULL) {
	if (po->nodes[edge->node].n_seq > max) {
	  next_node = edge->node;
	  max = po->nodes[next_node].n_seq;
	}
	edge = edge->next;
      }
      i = next_node;
    } else {
      i = po->out_edges[i]->node;
    }
    node_count++;
  }

  fprintf(stderr,"Max-flow path covers %d nodes (%1.2f)\n",node_count,
	  node_count/(float) po->n_nodes * 100.0);
  
  return node_count;
}

inline int base_toint(uchar base) {

  switch(base) {
  case 'A':
    return 0;
  case 'C':
    return 1;
  case 'G':
    return 2;
  case 'T':
    return 3;
  }

  return -1;
}

inline uchar int_tobase(int i) {
  uchar base[4] = "ACGT";

  if (i <0 || i > 3) return 0;
  else return base[i];
}

#define PE(x) (pow(10,(-(x)/10)))

static int analyze_sgp_mismatches(po_t *po) {
  int i,j,k,l;
  int baseint;
  int rval;
  node_t *temp_node, *node;
  int *seqdiv;
  double pdb[4];
  double pbd[4];
  double pd, val;

  CA(seqdiv,n_seq,sizeof(int));

  for(i=0;i<po->n_nodes;i++) {
    if (po->nodes[i].n_align > 0) {
      int n_obs[4] = {0, 0, 0, 0};
      double obs[4] = { 1., 1., 1., 1.};
      int n;

      node = po->nodes + i;
      baseint = base_toint(node->base);
      if (baseint != -1) 
	for(j=0;j<node->n_seq;j++) {
	  n_obs[baseint]++;
	  val = PE(sequences[node->seq[j]]->qual[node->seq_pos[j]]);
	  if (val < obs[baseint]) obs[baseint] = val;
	}
      for(j=0;j<node->n_align;j++) {
	temp_node = po->nodes + node->aligned_nodes[j];
	baseint = base_toint(temp_node->base);
	if (baseint != -1) {
	  for(k=0;k<temp_node->n_seq;k++) {
	    n_obs[baseint]++;
	    val=PE(sequences[temp_node->seq[k]]->qual[temp_node->seq_pos[k]]);
	    if (val < obs[baseint]) obs[baseint] = val;
	  }
	}
      }

      pd = 0;
      for(j=0;j<4;j++) {
	pdb[j] = 1.0;
	for(k=0;k<4;k++) {
	  if (j == k) continue;
	  if (n_obs[k]>0) 
	    pdb[j] *= obs[k]*(0.333333333333333);
	}
	if (n_obs[j] > 0)
	  pdb[j] *= (1.0 - obs[j]);
	pd += pdb[j]*0.25;
      }

      for(j=0;j<4;j++) {
	pbd[j] = pdb[j]*0.25/pd;
	fprintf(stderr,"Base %c: \t%d Observations P = %1.7f\n",int_tobase(j),
		n_obs[j],pbd[j]);
	assert(n_obs[j] < po->n_seq);
      }

      /* Temporary code -- FIX ME

         This calculation repeated unnecessarily for other nodes in the
	 align ring */
      n = 0;
      for(j=0;j<4;j++) {
	if (pbd[j] > 0.25) n++;
      }
      if (n>1 && pbd[base_toint(node->base)] > 0.25) {
	for(j=0;j<node->n_seq;j++) {
	  seqdiv[node->seq[j]]++;
	}
      }
    }
  }

  rval = 0;
  for(i=0;i<po->n_seq;i++) {
    if ((seqdiv[po->seq[i]]/(float) sequences[po->seq[i]]->length)
	> 0.01)
      rval = 1;
    fprintf(stderr,"Sequence %s: %1.5f\n", sequences[po->seq[i]]->label,
	    seqdiv[po->seq[i]]/(float) sequences[po->seq[i]]->length);
  }
  free(seqdiv);
  return rval;
}

static int analyze_single_gene_profile(po_t *po) {

  if (analyze_sgp_mismatches(po)) return 0;
  return 1;
  //analyze_sgp_mismatches(po);
}

#if 0
static void analyze_single_gene_profile(po_t *po) {
  int i, j;
  int *used;
  int mismatch, corr_mismatch, max, n_unaligned, node_count;
  int aligned_start, aligned_end;
  
  mismatch = corr_mismatch = 0;
  CA(used, sizeof(int), po->n_nodes);
  for(i=0;i<po->n_nodes;i++) {
    if (used[i]) continue;
    used[i] = 1;
    if (po->nodes[i].n_align > 0) {
      max = po->nodes[i].n_seq;
      for(j=0;j<po->nodes[i].n_align;j++) {
	if (po->nodes[po->nodes[i].aligned_nodes[j]].n_seq > max) {
	  max = po->nodes[po->nodes[i].aligned_nodes[j]].n_seq;
	}
      }
      if (po->nodes[i].n_seq < max) { 
	if (po->nodes[i].n_seq > 1) {
	  corr_mismatch++;
	} else {
	  mismatch++;
	}
      }
      for(j=0;j<po->nodes[i].n_align;j++) {
	used[po->nodes[i].aligned_nodes[j]] = 1;
	if (po->nodes[po->nodes[i].aligned_nodes[j]].n_seq < max) {
	  if (po->nodes[po->nodes[i].aligned_nodes[j]].n_seq > 1) {
	    corr_mismatch++;
	  } else {
	    mismatch++;
	  }
	}
      }
    }
  }

  free(used);
  fprintf(stderr,"Mismatched = %d\tCorroborated Mismatch = %d\n",
	  mismatch, corr_mismatch);

  /* Scan for unmatched/unaligned nodes */
  aligned_start = 0;
  while(po->nodes[aligned_start].n_seq == 1 && 
	po->nodes[aligned_start].n_align == 0) aligned_start++;
  aligned_end = po->n_nodes-1;
  while(po->nodes[aligned_end].n_seq == 1 &&
	po->nodes[aligned_end].n_align == 0) aligned_end--;
  n_unaligned = 0;
  for(i=aligned_start;i<aligned_end;i++) {
    if (po->nodes[i].n_seq == 1 && po->nodes[i].n_align == 0) {
      n_unaligned++;
    }
  }
  fprintf(stderr,"Unaligned nodes %d\n", n_unaligned);

}
#endif

static void check_alignrings(po_t *po) {
  int al_nodes[1000];
  int n_align;
  int i,j;
  edge_t *edge;

  for(i=0;i<po->n_nodes;i++) {
    if (po->nodes[i].n_align > 0) {
      n_align = 0;
      for(j=0;j<po->nodes[i].n_align;j++) {
	al_nodes[n_align++] = po->nodes[i].aligned_nodes[j];
      }
      edge = po->in_edges[i];
      while(edge!=NULL) {
	for(j=0;j<n_align;j++) {
	  if (edge->node == al_nodes[j]) {
	    fprintf(stderr,"Edge within align ring\n");
	  }
	}
	edge = edge->next;
      }
      edge = po->out_edges[i];
      while(edge!=NULL) {
	for(j=0;j<n_align;j++) {
	  if (edge->node == al_nodes[j]) {
	    fprintf(stderr,"Edge within align ring\n");
	  }
	}
	edge = edge->next;
      }
    }
  }
}

void test_shit() {
  po_t **pos;
  int *tracking;
  int *po_size;
  int **po_members;
  int i,j,k;
  int contigs, singlets;
  po_t *a, *b, *r;
  int did_join;

  CA(pos, n_seq, sizeof(po_t *));
  CA(tracking, n_seq, sizeof(int));
  CA(po_members, n_seq, sizeof(int *));
  CA(po_size, n_seq, sizeof(int));

  for(i=0;i<n_seq;i++) {
    tracking[i] = -1;
  }
#if 1
  for(i=0;i<n_components;i++) {
    int fill;
    for(j=0;j<component_size[i];j++) {
      pos[components[i][j]] = make_po(components[i][j]);
    }
    while(component_size[i]>1) {
      fill = 0;
      for(j=0;j<component_size[i];j+=2) {
	clock_t stopwatch;
	if (j+1 == component_size[i]) break;
	stopwatch = clock();
	a = pos[components[i][j]];
	b = pos[components[i][j+1]];
	r = poa_align(a,b);
	fuse_identical_alignring_nodes(r);
	po_topological_sort(r);
	//display_po(r,  stderr);
	pos[components[i][j]] = NULL;
	pos[components[i][j+1]] = NULL;
	pos[components[i][fill++]] = r;
	destroy_po(a);
	destroy_po(b);
	fprintf(stderr,"Total time for poa_align call: %5.5f\n",
		(clock() - stopwatch)/(float) (CLOCKS_PER_SEC)); 
      }
      if (j<component_size[i]) {
	pos[components[i][fill++]] = pos[components[i][j]]; 
	pos[components[i][j]] = NULL;
      }
      component_size[i] = fill;
    }
  }
      
	
#else
  for(i=0;i<n_components;i++) {
    for(j=0;j<component_size[i];j++) {
      pos[components[i][j]] = make_po(components[i][j]);
      tracking[components[i][j]] = components[i][j];
      PUSH(po_members[components[i][j]],po_size[components[i][j]],sizeof(int));
      po_members[components[i][j]][po_size[components[i][j]]++] = 
	components[i][j];
    }

    for(j=0;j<component_size[i];j++) {
      int aint, bint, k;

      if (assembly_order[i][j].s1 == -1) {
	fprintf(stderr,"Root node %s\n",sequences[assembly_order[i][j].s2]->label);
	continue;
      }
      bint = tracking[assembly_order[i][j].s2];
      aint = tracking[assembly_order[i][j].s1];

      fprintf(stderr,"%s (%d) --> %s (%d) (%1.1f%% done)\n",
	      sequences[assembly_order[i][j].s1]->label,
	      po_size[tracking[aint]],
	      sequences[assembly_order[i][j].s2]->label,
	      po_size[tracking[bint]],
	      (j/(float) component_size[i])*100.0);
      a = pos[aint];
      b = pos[bint];
      r = poa_align(a, b);
      fuse_identical_alignring_nodes(r);
      po_topological_sort(r);
      display_po(r, stderr);

      pos[tracking[assembly_order[i][j].s1]] = r;
      pos[tracking[assembly_order[i][j].s2]] = NULL;
      destroy_po(a);
      destroy_po(b);
      for(k=0;k<po_size[bint];k++) {
	tracking[po_members[bint][k]] = aint;
	PUSH(po_members[aint], po_size[aint], sizeof(int));
	po_members[aint][po_size[aint]++] = po_members[bint][k];
      }
      free(po_members[bint]);
      po_members[bint] = NULL;
      po_size[bint] = 0;
    }
  }
#endif
  contigs = singlets = 0;
  for(i=0;i<n_seq;i++) {
    if (pos[i]!=NULL) {
      if (pos[i]->n_seq > 1) {
	fprintf(stdout,"Contig %d\n",contigs++);
	display_po(pos[i], stdout);
      } else {
	fprintf(stdout,"Singlet %d\n",singlets++);
	display_po(pos[i], stdout);
      }
    }
  }

}

static void print_edgelist(edge_t *edge) {

  while(edge!=NULL) {
    fprintf(stderr,"%d\t",edge->node);
    edge = edge->next;
  }
  fprintf(stderr,"\n");
}


#if 0
inline int check_membership(node_t *node, int seq_id) {
  int i;
  for(i=0;i<node->n_seq;i++)
    if (node->seq[i] == seq_id) return 1;

  return 0;
}

static void analyze_po(po_t *po) {

  for(i=0;i<po->n_seq;i++) {
    seq_id = po->seq[i];
    
    /* Find starting node */
    j = 0;
    while(j<po->n_nodes) {
      if (check_membership(po->nodes + j, seq_id)) break;
      j++;
    }
    assert(j<po->n_nodes);
    while(j<po->n_nodes) {
      edge = po->out_edges[j];
      while(edge!=NULL) {
	if (check_membership(po->nodes + edge->node, seq_id)) break;
	edge = edge->next;
      }
      if (edge == NULL) {
	fprintf(stderr,"Terminating node found\n");
	break;
      }
      
    
}
#endif


