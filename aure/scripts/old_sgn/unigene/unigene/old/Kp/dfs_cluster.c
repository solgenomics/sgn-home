#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <limits.h>
#include <assert.h>
#include <getopt.h>

#include "log_message.h"
#include "kp_types.h"

/* These variables are set by command line options. If there value is
   not NULL, the options they specify are enabled below. */
static uchar *chimera_file = NULL;
static uchar *database_name = NULL;
static int verbosity_level;

static int *chimeric = NULL;
static seqmeta_t *seqmeta = NULL;
static uchar *seqname_data = NULL;

static int n_seq;

static int **tree_edges;
static int *n_tedges;
static int **back_edges;
static int *n_bedges;

static int n_arti;
static int *arti_points;

typedef struct {
  int s2;
  int score;
} fasta_t;

static int *n_fasta;
static fasta_t **fasta_scores;

static int n_components;
static int *component_size;
static int **components;

static int scan_artipoint(int node, int l, int *level) {
  int i, child_level, min, atri_point;

  atri_point = 0;
  min = INT_MAX;

  assert(node < n_seq);
  level[node] = l;
  for(i=0;i<n_tedges[node];i++) {
    child_level = scan_artipoint(tree_edges[node][i], l+1, level);
    if (child_level == INT_MAX) {
      atri_point = 1;
    } else {
      if (child_level >= l) atri_point = 1;
      if (child_level < min) min = child_level;
    }
  }

  for(i=0;i<n_bedges[node];i++) {
    if (level[back_edges[node][i]] < min) min = level[back_edges[node][i]];
  }

  if (atri_point && l!=0) {
    arti_points[node] = 1;
    n_arti++;
  }

  return min;
}

static void scan_arti_points(int component) {
  int *level;

  if (component_size[component] <= 1) return;
  CA(level, n_seq, sizeof(int));
  if (n_tedges[components[component][0]] > 1) {
    arti_points[components[component][0]] = 1;
  }

  scan_artipoint(components[component][0], 0, level);
  free(level);
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



static void cc_dfs_visit(int *color, int **component, int *component_size, 
			 int n) {
  int i;

  color[n] = 1;

  PUSH(*component, *component_size, sizeof(int));
  (*component)[(*component_size)++] = n;
  assert(chimeric[n] == 0);
  for(i=0;i<n_fasta[n];i++) {
    if (chimeric[fasta_scores[n][i].s2]) continue;
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

static void connected_components() {
  int *color;
  int i;
  int clusters, singletons;

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
    if (chimeric[i]) continue;

    PUSH(components, n_components, sizeof(int *));
    PUSH(component_size, n_components, sizeof(int));
    components[n_components] = NULL;
    component_size[n_components] = 0;
    cc_dfs_visit(color, components + n_components, 
		 component_size + n_components, i);
    n_components++;
  }


  fprintf(stderr,"Found %d connected components\n", n_components);
  CA(arti_points, n_seq, sizeof(int));
  n_arti = 0;
  clusters = singletons = 0;

  for(i=0;i<n_components;i++) {
    int j;

    if (component_size[i] > 1) {

      fprintf(stdout,">Cluster %d (%d sequences)\n",clusters,
	      component_size[i]);

      scan_arti_points(i);
      if (database_name) {
	for(j=0;j<component_size[i];j++) {
	  fprintf(stdout,"%s ",seqname_data + 
		  seqmeta[components[i][j]].name_pos);
	}
      } else {
	for(j=0;j<component_size[i];j++) {
	  fprintf(stdout,"%d ",components[i][j]);
	}
      }
      fprintf(stdout,"\n");
      clusters++;
    } else {
      singletons++;
    }
  }

  fprintf(stdout,">Singletons (%d sequences)\n", singletons);
  if (database_name) {
    for(i=0;i<n_components;i++) {
      if (component_size[i] == 1) {
	fprintf(stdout,"%s ",seqname_data + 
		seqmeta[components[i][0]].name_pos);
      }
    }
  } else {
    for(i=0;i<n_components;i++) {
      if (component_size[i] == 1) {
	fprintf(stdout,"%d ",components[i][0]);
      }
    }
  }
  fprintf(stdout,"\n");

  free(color);
  free(back_edges);
  free(tree_edges);
  free(n_tedges);
  free(n_bedges);
  fprintf(stderr,"Clusters %d Singletons %d\n",clusters,singletons);
}

static void load_scores(FILE *f) {
  int i;
#if 0
  f = fopen(filename, "r");
  if (f == NULL) {
    logmsg(MSG_FATAL,"Failed opening file \"%s\" (%s)\n",filename,
	   strerror(errno));
  }
#endif  
  if (database_name == NULL) {
    fread(&n_seq, sizeof(uint), 1, f);
  } else {
    uint x;
    fread(&x, sizeof(uint), 1, f);
    if (x > n_seq) {
      logmsg(MSG_FATAL,"More sequences found in adjacency list input than in "
	     "sequence database index.\n");
    }
  }

  MA(n_fasta, sizeof(int)*n_seq);
  MA(fasta_scores, sizeof(fasta_t *)*n_seq);
  fread(n_fasta, sizeof(int), n_seq, f);
  for(i=0;i<n_seq;i++) {
    MA(fasta_scores[i], sizeof(fasta_t)*n_fasta[i]);
    fread(fasta_scores[i], sizeof(fasta_t), n_fasta[i], f);
  }

#if 0
  for(i=0;i<n_seq;i++) {
    fprintf(stderr,"For sequence %d loaded %d scores:\n",i,n_fasta[i]);
    for(j=0;j<n_fasta[i];j++) {
      fprintf(stderr,"%d\t%d\t",fasta_scores[i][j].s2,fasta_scores[i][j].score);
    }
    fprintf(stderr,"\n");

  }
#endif
}

static void usage(char *program_name) {

  fprintf(stderr,"\n\n%s:\n\n"
"Quick program to scan formatted sequence file against a pre-formatted \n"
"database of words (sub-sequence), to approximate alignment by linking \n"
"together consecutive sequences of matching words.\n"
"\n"
"Options:\n"
"--chimera=<chimera file> (-c) \n"
"    Filename of sequence ids (integers) which are (probably) chimeric. \n"
"    These sequences are excluded in clustering\n"
"--database=<basename> (-s) \n"
"    Basename of preformatted sequence 'database' from which homology reports\n"
"    are derived\n"
"--verbose=<integer> (-v)\n"
"    Verbosity level. 0 (normal) by default. Negative enables debugging messages\n"
"    Positive makes program quieter.\n"
"--help (-h)\n"
"    Prints this message.\n"
"\n"
"    Program expects binary format adjacency list on standard input, writes\n"
"    clustering output to standard output in text. If database is specified,\n"
,program_name);

}


static void parse_arguments(int argc, char *argv[]) {
  int option_index, commandline_error, rval;
  struct option longopts[] = {
    { "chimera", 1, NULL, 'c'},
    { "database", 1, NULL, 'd'},
    { "verbose", 1, NULL, 'v'},
    { "help", 1, NULL, 'h'},
    { NULL, 0, NULL, 0}
  };
  char *optstring = "c:d:v:h";

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
      database_name = strdup(optarg);
      break;
    case 'v':
      verbosity_level = atoi(optarg);
      break;
    case 'c':
      chimera_file = strdup(optarg);
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
  
  if (commandline_error) {
    logmsg(MSG_ERROR,"! Program halted due to command line option errors\n");
    usage(argv[0]);
    exit(-1);
  }

}

static void load_chimeras(uchar *chimera_file) {
  FILE *f;
  uchar inputline[1024];
  uint seq_id;
  int c;

  f = fopen(chimera_file, "r");
  if (f == NULL) {
    logmsg(MSG_FATAL,"Failed opening chimera file \"%s\" (%s)\n",
	   chimera_file, strerror(errno));
  }

  while((c = fgetc(f)) != EOF) {
    if (c == '>') {
      fgets(inputline, 1024, f);
      seq_id = atoi(inputline);
      if (seq_id < n_seq) {
	chimeric[seq_id] = 1;
      }
    }
  }

  fclose(f);
}

static void load_seqnames(uchar *database_name) {
  int l;
  uchar *temp;
  uint x;
  uint i;
  int total_size;
  FILE *f;

  l = strlen(database_name) + 6;
  temp = alloca(l);
  strcpy(temp, database_name);
  strcat(temp, ".ind");

  f = fopen(temp, "r");
  if (f == NULL) {
    logmsg(MSG_FATAL,"! Failed opening database index file \"%s\" (%s)\n",
	   temp, strerror(errno));
  }

  fread(&x, sizeof(uint), 1, f);
  if (x != INDFILE_MAGIC) {
    logmsg(MSG_FATAL,"! Database index file \"%s\" does not appear to be properly formatted\n",temp);
  }
  fread(&n_seq, sizeof(uint), 1, f);
  MA(seqmeta, sizeof(seqmeta_t)*n_seq);
  fread(seqmeta, sizeof(seqmeta_t), n_seq, f);
  total_size = 0;
  for(i=0;i<n_seq;i++) {
    total_size += seqmeta[i].name_length + 1;
  }
  MA(seqname_data, total_size);
  fread(seqname_data, total_size, 1, f);

  fclose(f);
}

int main(int argc, char *argv[]) {
  FILE *f;
  int i;

  configure_logmsg(MSG_DEBUG1);
  parse_arguments(argc, argv);

  n_seq = 0;
  if (database_name) load_seqnames(database_name);
  load_scores(stdin);
  CA(chimeric, n_seq, sizeof(int));
  if (chimera_file) load_chimeras(chimera_file);
  connected_components();

  f = fopen("articulations.txt","w");
  for(i=0;i<n_seq;i++) {
    if (arti_points[i]) {
      fprintf(f,"%d\n",i);
    }
  }
  fclose(f);

  return 0;
}
