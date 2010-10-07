#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

/* ---- This file's types and exports ---- */
#include "log_message.h"

static int verbosity_level;

static int logmessage_ready = 0;

/* Set variables which control where a message is copied. Negative value means
   don't modify the current setting. This function should be called at startup
   with 0,1,LOG_INFO, and switched to 1,0,LOG_WARNING when the program forks
   into background. Debugging builds however may do something different. */
int configure_logmsg(int verbosity_status) {
  
  if (logmessage_ready == 0) {
    logmessage_ready = 1;
  }

  verbosity_level = verbosity_status;

  return 0;
}

/* Write a message to the syslog or console, depending on its priority level
   and the messaging configuration via the variables above. */ 
void logmsg(int priority, char *s, ...) {
  va_list ap;

  if (logmessage_ready) {
    if (priority >= verbosity_level) {
      va_start(ap, s);
      vfprintf(stderr, s, ap);
    }
  } else {
    fprintf(stderr,"Call configure_logmsg() first!\n");
  }
  if (priority == MSG_FATAL) 
    exit(-1);
}
