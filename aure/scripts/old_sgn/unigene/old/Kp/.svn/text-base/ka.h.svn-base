#ifndef _KA_H
#define _KA_H
typedef unsigned char uchar;

typedef struct {
  uchar *label;
  uchar *seqstr;
  int *seq;
  int *qual;
  uint length;
  double expected_error;
} seq_t;

typedef struct {
  int s2;
  int score;
  double obs_diff;
} fasta_t;

typedef struct {
  int s1;
  int s2;
  int score;
} asm_order_t;

#define PUSH(a,l,t) \
if ((l) % 128 == 0)   { \
   (a) = realloc((a), (t)*((l) + 128)); \
   if ((a) == NULL) {                     \
      logmsg(MSG_FATAL,"Failed allocating memory at %s:%d (%d bytes)\n", \
             (t)*(l + 128),__FILE__,__LINE__);  \
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


#endif
