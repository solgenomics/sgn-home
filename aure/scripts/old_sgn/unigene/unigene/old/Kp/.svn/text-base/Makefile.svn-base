COMM_OBJS=	log_message.o

LIBS= -lm
CFLAGS=-Wall -ggdb
#CFLAGS=-Wall -fomit-frame-pointer -funroll-loops -fexpensive-optimizations -O3 -march=pentiumpro #-pg

all: format_seqdata format_lookup scan_sequences dfs_cluster

%.o: %.c
	gcc -c $(CFLAGS) $<

scan_sequences: $(COMM_OBJS) scan_sequences.o
	gcc $(CFLAGS) -oscan_sequences scan_sequences.o $(COMM_OBJS) $(LIBS)

format_seqdata: $(COMM_OBJS) format_seqdata.o
	gcc $(CFLAGS) -oformat_seqdata format_seqdata.o $(COMM_OBJS) $(LIBS)

format_lookup: $(COMM_OBJS) format_lookup.o
	gcc $(CFLAGS) -oformat_lookup format_lookup.o $(COMM_OBJS) $(LIBS)

dfs_cluster: $(COMM_OBJS) dfs_cluster.o
	gcc $(CFLAGS) -odfs_cluster dfs_cluster.o $(COMM_OBJS) $(LIBS)

clean:
	rm -f *.o ka format_lookup format_seqdata scan_sequences
