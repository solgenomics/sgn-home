#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <math.h>

#include "kp_types.h"
#include "log_message.h"


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

void load_inputsequence(uchar *input_seqfile) {
  FILE *sf, *qf;
  uint line_no, input_length, name_start, name_end, i, j;
  int last_seqid;
  uint input_allocsize;
  uchar *inputline;
  seq_t *seq;
  
  n_seq = input_allocsize = 0;
  sequences = NULL;
  inputline = NULL;


  MA(seq, sizeof(seq_t));
  sf = openfile(input_seqfile, "r", "FASTA sequence file");

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
	     "name not found in header: %s\n",line_no, input_seqfile, 
	     inputline);
    }

    inputline[name_end] = 0;

    seq->label = strdup(inputline + name_start);
    seq->seqstr = NULL;
    seq->seq = NULL;

    input_length = read_fullline(&inputline, &input_allocsize, sf);
    line_no++;
    i = 0;
    while(input_length>0 && inputline[0]!='>') {
      j=0;
      while(j<input_length) {
	if (nucleotide(inputline[j])) {
	  PUSH(seq->seqstr, i, sizeof(uchar));
	  seq->seqstr[i++] = inputline[j];
	}
	j++;
      }
      input_length = read_fullline(&inputline, &input_allocsize, sf);
      line_no++;
    }
    seq->seqstr = realloc(seq->seqstr, (i+1)*sizeof(uchar));
    seq->seqstr[i] = 0;

    seq->length = i;
    seq->qual = NULL;

    PUSH(sequences, n_seq, sizeof(seq_t *));
    sequences[n_seq] = seq;
    n_seq++;
  }

  logmsg(MSG_INFO,"Loaded %d sequences from %s\n",n_seq, input_seqfile);
  fclose(sf);

  free(inputline);

}

void polya_truncate(void) {
  int i,j,k;
  int na, qual, pasta_qual;
  seq_t *seq;

  typedef struct {
    uint truncate_pos;
    uint a_length;
    uint a_qual;
    uint pasta_qual;
  } candidate_t;

  candidate_t *candidates;
  int cand, cand_allocate;

  logmsg(MSG_INFO,"Trimming PolyA and PolyT tails...\n");
  candidates = NULL;
  cand_allocate = 0;
  for(i=0;i<n_seq;i++) {
    j = 0;
    seq = sequences[i];

    cand = 0;
    while(j<seq->length) {
      if (seq->seqstr[j++] == 'A') {
	na = 1; qual = seq->qual[j-1];
	while(j<seq->length && seq->seqstr[j] == 'A') {
	  na++;
	  qual += seq->qual[j++];
	}
	if (na > 11) {
	  k = j;
	  pasta_qual = 0;
	  while(k<seq->length && (k-j)<na) {
	    pasta_qual += seq->qual[k++];
	  }
	  if (k>j && qual/na > (pasta_qual/(k-j))*1.5 && 
	      seq->length-j<seq->length/3) {	    
	    if (cand >= cand_allocate) {
	      cand_allocate += 10;
	      candidates = realloc(candidates, 
				   sizeof(candidate_t)*cand_allocate);
	      if (candidates == NULL) {
		logmsg(MSG_FATAL,"Failed allocating memory for ploy-A "
		       "truncation (%d bytes)\n",
		       sizeof(candidate_t)*cand_allocate);
	      }
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
	     "position %d\n", seq->label, candidates[select].truncate_pos);
      seq->seqstr = realloc(seq->seqstr, 
			    candidates[select].truncate_pos + 1);
      seq->seqstr[candidates[select].truncate_pos] = 0;
      seq->qual = 
	realloc(seq->qual, (candidates[select].truncate_pos + 1)*sizeof(int));
      seq->length = candidates[select].truncate_pos;
    }
  }

  for(i=0;i<n_seq;i++) {
    seq = sequences[i];
    j = seq->length;

    cand = 0;
    while(j>=0) {
      if (seq->seqstr[j--] == 'T') {
	na = 1; qual = seq->qual[j+1];
	while(j>=0 && seq->seqstr[j] == 'T') {
	  na++;
	  qual += seq->qual[j--];
	}
	if (na > 11) {
	  k = j;
	  pasta_qual = 0;
	  while(k>=0 && (j-k)<na) {
	    pasta_qual += seq->qual[k--];
	  }
	  if (k<j && qual/na > (pasta_qual/(j-k))*1.5 && j<seq->length/3) {
	    if (cand >= cand_allocate) {
	      cand_allocate += 10;
	      candidates = realloc(candidates, 
				   sizeof(candidate_t)*cand_allocate);
	      if (candidates == NULL) {
		logmsg(MSG_FATAL,"Failed allocating memory for ploy-A "
		       "truncation (%d bytes)\n",
		       sizeof(candidate_t)*cand_allocate);
	      }
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
      logmsg(MSG_DEBUG1,"Truncating leading pre-poly-T noise for sequence "
	     "%s at position %d\n", seq->label, 
	     candidates[select].truncate_pos);
      seq->length -= candidates[select].truncate_pos + 1;
      memmove(seq->seqstr, 
	      seq->seqstr + candidates[select].truncate_pos + 1,
	      seq->length);
      seq->seqstr = realloc(seq->seqstr, seq->length + 1);
      memmove(seq->qual, 
	      seq->qual + candidates[select].truncate_pos,
	      seq->length*sizeof(int));
      seq->qual = realloc(seq->qual, seq->length*(sizeof(int)));
    }
  }
  free(candidates);
  logmsg(MSG_INFO,"Finished trimming PolyA and PolyT tails.\n");
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

void compute_intsequence(void) {
  int i,j;
  seq_t *seq;
  double expected_error;

  for(i=0;i<n_seq;i++) {
    seq = sequences[i];
    if (seq->seq != NULL) {
      logmsg(MSG_FATAL,"Integer sequence already computed for sequence %s\n",
	     seq->label);
    }
    expected_error = 0.0;
    MA(seq->seq, sizeof(int)*seq->length);
    for(j=0;j<seq->length;j++) {
      seq->seq[j] = base_toint(seq->seqstr[j]);
      expected_error += pow(10.0, seq->qual[j]/-10.0);
    }
    seq->expected_error = expected_error;
  }
}

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
