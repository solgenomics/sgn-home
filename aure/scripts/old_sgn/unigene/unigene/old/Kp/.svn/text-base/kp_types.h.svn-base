#ifndef _KP_TYPES_H
#define _KP_TYPES_H

typedef unsigned char uchar;

typedef struct {
  int name_length;
  int name_pos;
  int seq_length;
  uint seqstr_pos;
  uint seqbin_pos;
} seqmeta_t;

typedef struct {
  uint seq_id;
  uint seq_pos;
} word_t;

typedef struct {
  int n_words;
  uint start_pos;
} lookupmeta_t;

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

#define RA(p,n,s)     \
(p) = realloc(p, (n)*(s)); \
if ((p) == NULL) {    \
   logmsg(MSG_FATAL,"Failed allocating memory at %s:%d (%d bytes)\n", \
          __FILE__,__LINE__,(n)*(s)); \
}

#define INDFILE_MAGIC (0x10001217)
#define STRFILE_MAGIC (0x10001218)
#define BINFILE_MAGIC (0x10001219)
#define LOOKUP_MAGIC  (0x100013A1)

#endif
