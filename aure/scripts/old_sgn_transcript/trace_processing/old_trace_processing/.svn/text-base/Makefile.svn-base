# -*- Mode: Fundamental -*- # Makefile mode made my emacs hang...?
# This is the library processing aspect of the Unigene pipeline.

# I (Marty) inherited this from Koni, whose scripts were fairly hard
# to follow and mostly changed things in the database.  When we moved
# from MySQL to Postgres, it became necessary to rewrite those
# programs, and since the requirements of the programs change subtly
# and not-subtly from run to run, it seemed worthwhile to me to have
# the library processing take place entirely outside the database, and
# to load the results only when everything worked out.

# This is my second reimplementation of library processing; my first
# reimplementation involved an object-oriented facility for isolating
# the inputs and outputs of each stage of the process.  I realized
# only after I wrote it that the I/O managing objects were doing about
# half of what Unix make(1) does, but not as well.  ISTM now that make
# offers a reasonable amount of flexibility (e.g., you can pick up in
# the middle of the process by touching some files), and is a facility
# that one might want to get familiar with anyhow.  This Makefile has
# been kept as simple as possible, not taking advantage of any
# particularly fancy make features. 

# What you need:

# * a directory called "chromat" containing all your chromatograms,
#   /OR/ an empty directory called "chromat", an empty directory
#   called "phddir", a file called "rawseq.fasta", and a file called
#   "rawseq.qual" in case there are no chromatograms (e.g., if you are
#   processing sequence from Genbank).
# * a copy of UniVec in FASTA format
# * a directory containing our screening contaminants
# * an arabidopsis peptide blast database
# * an arabidopsis peptide self-blast in m8 format
# * (optionally) the vector your sequences were cloned with
# * (optionally) the 5' adaptor used in sequencing

# Here's the approximate command to run:

# make -f ../../Makefile UNIVEC=/data/local/unigene-pipeline/screening-contaminants/UniVec.seq PHRED_PARAMETER_FILE=/home/kreuter/sys/bin/phredpar.dat CONTAMINANTS_DIR=/data/local/unigene-pipeline/screening-contaminants VECTOR=../pcr4topo.fasta ATHDB=/data/local/unigene-pipeline/arabidopsis/ath1/ATH1_pep ADAPTOR=GTGGATCCAAGAATTCGGCACGAGG ATHSELFB=/path/to/arabidopsis-self-blast.m8

# General notes about make: dollar-signs in code fragments (e.g., in
# Perl) must be doubled, or else make will interpolate its own
# variables into the code.
#

# If you do need to start the processing from some mid-point (e.g., if
# there are no chromatograms, but only raw sequence FASTA), make sure
# to touch all the files that are produced by steps prior to the place
# where you want to start, so that make will consider those targets as
# satisfied.

#BLAST=/root/sgn-tools/bin/qsblast.pl # or blastall
BLAST=blastall
#UNIVEC=/home/kreuter/screening-contaminants/UniVec.seq
UNIVEC=/data/prod/public/screening/screening-contaminants/UniVec.seq
#CONTAMINANTS_DIR=/home/kreuter/screening-contaminants
CONTAMINANTS_DIR=/data/prod/public/screening/screening_contaminants
#ATHDB=/home/kreuter/arabidopsis/ath1/ATH1_pep
ATHDB=/data/prod/public/screening/arabidopsis/ath1/ATH_pep
#ATHSELFB=/home/kreuter/arabidopsis/ATH1-self-blast/ATH1_pep.blastp.m8
ATHSELFB=/data/prod/public/screening/arabidopsis/ATH1-self-blast/ATH1_pep.blastp.m8
VECTOR=
ADAPTOR=

# The targets in this file are written in reverse order for no
# particular reason.  Feel free to rearrange them.

all: chimera_screen_analyze.tab

chimera_screen_analyze.tab: PROG:='\
	  my %equivalence\;\
	  open my $$athselfbp, $$ARGV[0] or die "$$!"\;\
	  while (<$$athselfbp>) {\
	    chomp\;\
	    my @record = split\;\
	    $$id = shift @record\;\
	    @{$$equivalence{$$id}} = @record\;\
	  }\
	  my %blast\;\
	  keys %blast = 20000\;\
	  while (<STDIN>) {\
	    my ($$query, $$subject) = split\;\
	    $$blast{$$query} = $$subject unless (exists($$blast{$$query}))\;\
	  }\
	  my %clones_already_seen\;\
	  foreach my $$query (keys %blast) {\
	    my $$clone_id = $$query\;\
	    $$clone_id =~ s/:[35]//\;\
	    next if ($$clones_already_seen{$$clone_id})\;\
	    $$clones_already_seen{$$clone_id} = 1\;\
	    my $$partner = $$query\;\
	    ($$query =~ m/3$$/) ? $$partner =~ s/.$$/5/ : $$partner =~ s/.$$/3/\;\
	    if ($$blast{$$query} && $$blast{$$partner}) {\
	      if (($$blast{$$query} ne $$blast{$$partner}) &&\
	  	(!(grep { $$_ eq $$blast{$$query} }\
	  	   @{$$equivalence{$$blast{$$partner}}}))) {\
	        print "$$clone_id\t$$blast{$$query}\t$$blast{$$partner}\n"\;\
	      }\
	    }\
	  }'
chimera_screen_analyze.tab: chimera_screen.m8 athselfb-permuted.tab
	-rm chimera_screen_analyze.tab chimera_screen_analyze.tab.tmp
	perl -e $(PROG) athselfb-permuted.tab < chimera_screen.m8 > chimera_screen_analyze.tab.tmp
	mv chimera_screen_analyze.tab.tmp chimera_screen_analyze.tab

athselfb-permuted.tab: PROG:='my %equivalence\;\
	  keys %equivalence = 32000\;\
	  while (<STDIN>) {\
	    my ($$query, $$subject) = split\;\
	    push @{$$equivalence{$$query}}, $$subject \
	      unless (grep { $$_ eq $$subject } @{$$equivalence{$$query}})\;\
	    push @{$$equivalence{$$subject}}, $$query \
	      unless (grep { $$_ eq $$query } @{$$equivalence{$$subject}})\;\
	  }\
	  foreach my $$key (keys %equivalence) {\
	    printf "%s\n", join("\t", $$key, @{$$equivalence{$$key}})\
	  }'
athselfb-permuted.tab: $(ATHSELFB)
	-rm athselfb-permuted.tab athselfb-permuted.tab.tmp
	perl -e $(PROG) < $(ATHSELFB) > athselfb-permuted.tab.tmp
	mv athselfb-permuted.tab.tmp athselfb-permuted.tab

chimera_screen.m8: chimera_screen_input.fasta
	-rm chimera_screen.m8.tmp chimera_screen.m8
	$(BLAST) -i chimera_screen_input.fasta -d $(ATHDB) -p blastx -e 1e-6 -m 8 -o chimera_screen.m8.tmp -v 1 -b 1
	mv chimera_screen.m8.tmp chimera_screen.m8

chimera_screen_input.fasta: PROG='\
	  use Bio::SeqIO\;\
	  $$in=Bio::SeqIO->new(-file=>$$ARGV[0],-format=>"Fasta")\;\
	  while ($$s=$$in->next_seq) {\
	    $$length=length($$s->seq)\;\
	    if ($$length>300) {\
	      $$range=(($$length < 600) ? $$length/2 : 300)\;\
	      $$forw=substr($$s->seq, 0, $$range)\;\
	      $$revs=substr($$s->seq, -1*$$range)\;\
	      printf ">%s:5\n%s\n>%s:3\n%s\n", $$s->id, $$forw, $$s->id, $$revs\;\
	    }\
	  }'
chimera_screen_input.fasta: contaminant_screened.fasta
	-rm chimera_screen_input.fasta.tmp chimera_screen_input.fasta
	perl -e $(PROG) contaminant_screened.fasta > chimera_screen_input.fasta.tmp
	mv chimera_screen_input.fasta.tmp chimera_screen_input.fasta
	/root/sgn-tools/util/fastaFileCheck -m chimera_screen_input.fasta

contaminant_screened.fasta: PROG='\
	use strict\;\
	use Bio::SeqIO\;\
	my %bad\;\
	while(<STDIN>){ chomp\; $$bad{$$_}=1 }\
	my $$f=Bio::SeqIO->new(-file=>$$ARGV[0],-format=>"Fasta")\;\
	while (my $$s=$$f->next_seq){\
	    my $$fmt=">%s %d \n%s\n"\;\
	    printf $$fmt, $$s->id, length($$s->seq),$$s->seq\
	}'
contaminant_screened.fasta: contaminant_screen.m8 trimseq.fasta
	-rm contaminant_screened.fasta.tmp contaminant_screened.fasta
	cut -f1 contaminant_screen.m8 | perl -e $(PROG) trimseq.fasta > contaminant_screened.fasta.tmp
	mv contaminant_screened.fasta.tmp contaminant_screened.fasta
	/root/sgn-tools/util/fastaFileCheck -m contaminant_screened.fasta

contaminant_screen.m8: trimseq.fasta contaminant.db
	-rm contaminant_screen.m8.tmp contaminant_screen.m8
	$(BLAST) -i trimseq.fasta -d contaminant.db -p blastn -e 1e-6 -m 8 -o contaminant_screen.m8.tmp
	mv contaminant_screen.m8.tmp contaminant_screen.m8

contaminant.db: contaminant.fasta
	formatdb -p F -i contaminant.fasta -n contaminant.db -l contaminant.log
	touch contaminant.db

contaminant.fasta:
	-rm contaminant.fasta.tmp contaminant.fasta
	test -n "$(CONTAMINANTS_DIR)" || { echo "CONTAMINANTS_DIR is unset" && exit 1; }
	cat $(CONTAMINANTS_DIR)/*.seq > contaminant.fasta.tmp
	mv contaminant.fasta.tmp contaminant.fasta

qc.tab: trimseq.fasta
	sed -n '/^>/{s/^>//; s/[[:space:]]/	/g; p}' trimseq.fasta > qc.tab

# Note: Makefile syntax requires doubling dollar signs in commands.
trimseq.fasta: PROG='\
	use strict\;\
	use Bio::SeqIO\;\
	my $$in = Bio::SeqIO->new(-fh=>\*STDIN,-format=>"fasta")\;\
	my $$adaptor_end_length = 6\;\
	my $$adaptor_length = length ($$ARGV[0])\;\
	my $$adaptor = substr($$ARGV[0], -$$adaptor_end_length)\;\
	while (my $$s=$$in->next_seq){\
	  my $$l\;\
	  split"X+",$$s->seq\;\
	  map{ if (length($$_) > length($$l)){ $$l=$$_ } } @_\;\
	  if($$l){\
	    my $$pos = index($$l, $$adaptor)\;\
	    if ((0 < $$pos) && ($$pos < 2*$$adaptor_length)) {\
	      $$l = substr($$l, $$pos+$$adaptor_end_length)\;\
	    }\
	    my $$fmt=">%s %d %d\n%s\n"\;\
	    if (length($$l)) {\
	      printf $$fmt, $$s->id, index($$s->seq,$$l), length($$l),$$l\
	    }\
	  }\
	}'
trimseq.fasta: rawseq.fasta.screen
	-rm trimseq.fasta.tmp trimseq.fasta
	perl -e $(PROG) $(ADAPTOR) < rawseq.fasta.screen > trimseq.fasta.tmp
	mv trimseq.fasta.tmp trimseq.fasta
	/root/sgn-tools/util/fastaFileCheck -m trimseq.fasta

rawseq.fasta.screen: rawseq.fasta
	test -n "$(UNIVEC)" || { echo "UNIVEC is unset" && exit 1; }
	cross_match rawseq.fasta $(UNIVEC) $(VECTOR) -screen > cross_match.out 2>&1

# FIXME: maybe escape semicolons for make
raw.tab: phddir
	-rm raw.tab.tmp raw.tab
	find phddir -type f -exec awk 'NR==1{id=$$2} /BEGIN_DNA/{seq="";qual="";pos="";inseq=1;next} /END_DNA/{inseq=0} inseq ~ /1/{seq=seq $$1; if (qual=="") qual=$$2; else qual=qual " " $$2; if(pos=="") pos=$$3; else pos=pos " " $$3} END{fmt="%s\t%s\t%s\t%s\n";printf(fmt,id,seq,qual,pos);}' {} \; > raw.tab.tmp
	mv raw.tab.tmp raw.tab

# FIXME: maybe escape semicolons for make
rawseq.fasta rawseq.fasta.qual: phddir
	-rm rawseq.fasta.tmp rawseq.fasta.qual.tmp rawseq.fasta rawseq.fasta.qual
	find phddir -type f -exec awk -v seqf=rawseq.fasta.tmp -v qualf=rawseq.fasta.qual.tmp 'NR==1{id=$$2} /BEGIN_DNA/{seq="";qual="";inseq=1;next} /END_DNA/{inseq=0} inseq ~ /1/{seq=seq $$1; if (qual=="") qual=$$2; else qual=qual " " $$2} END{if (seq) {fmt=">%s\n%s\n"; printf(fmt, id, seq) >> seqf; printf(fmt, id, qual) >> qualf}}' {} \;
	mv rawseq.fasta.tmp rawseq.fasta
	mv rawseq.fasta.qual.tmp rawseq.fasta.qual

# Note: phred will quit unless PHRED_PARAMETER_FILE is set.
phddir: chromat
	mkdir phddir
	find chromat -maxdepth 1 -type f | xargs -n1 phred -trim_alt "" -trim_phd -trim_cutoff 0.03 -pd phddir
#	find chromat -maxdepth 1 -type f \( -name '*.scf' -o -name '*.ab?' \) | xargs -n1 phred -trim_alt "" -trim_phd -trim_cutoff 0.03 -pd phddir

.PHONY: clean
clean:
	-rm rawseq.fasta rawseq.fasta.qual rawseq.fasta.log rawseq.fasta.screen cross_match.out trimseq.fasta contaminant_screen.m8 contaminant.db.* contaminant.db contaminant.log contaminant.fasta contaminant_screened.fasta read-clone-length-direction.tab chimera_screen_input.fasta chimera_screen.m8 athselfb-permuted.tab  chimera_screen_analyze.tab
	-rm -r phddir

.PHONY: dbload
dbload: raw.tab contaminant-scores.tab qc.tab
	test -n "$(LIBRARY_NAME)" || { echo "LIBRARY_NAME is unset" && exit 1; }
	echo "\
	set search_path = sgn_dev;\
	begin;\
	create temporary table tmpest (read_id integer, read_name varchar,\
	       seq varchar, qscore varchar, call_positions varchar, \
	       status integer, flags integer);\
	create temporary table tmpflags (read_name varchar, flags integer);\
	create temporary table tmpchim (read_name varchar, seq1 varchar, seq2 varchar);--\
	create temporary table tmpqc (est_id integer, read_name varchar, hqi_start bigint, hqi_length bigint);--\
	\copy tmpest (read_name, seq, qscore, call_positions) from 'raw.tab'--\
	update tmpest set flags = 0;\
	\copy tmpflags (read_name, flags) from 'contaminant-scores.tab'--\
	\copy tmpchim  (read_name, seq1, seq2) from 'chimera_screen_analyze.tab'--\
	\copy tmpqc (read_name, hqi_start, hqi_length) from qc.tab--\
	update tmpest set flags = tmpflags.flags from tmpflags \
	       where tmpest.read_name = tmpflags.read_name;\
	update tmpest set flags = flags | cast(x'80' as integer) where\
	       read_name in (select read_name from tmpchim);\
	update tmpqc set read_name = substring (read_name from '(.*)\\\\.')\
	       where read_name ~ '.*\\\\.(scf|ab[i012])$$';\
	update tmpest set status = 0;\
	select count(*) from tmpest where flags = 0;\
	update tmpest set read_name = substring (read_name from '(.*)\\\\.')\
	       where read_name ~ '.*\\\\.(scf|ab[i012])$$';\
	update tmpest set read_id = seqread.read_id from seqread\
	       where tmpest.read_name=seqread.trace_name;\
	update tmpest set read_id = seqread.read_id from\
	       seqread, clone, library\
	       where tmpest.read_name=seqread.trace_name\
	       and seqread.clone_id = clone.clone_id\
	       and clone.library_id = library.library_id\
	       and library.library_shortname='$(LIBRARY_NAME)';\
	insert into sgn_dev.est (read_id, seq, qscore, call_positions, status,\
	       flags) select read_id, seq, qscore, call_positions, status,\
	       flags from tmpest;\
	update tmpqc set est_id = est.est_id from est, seqread \
	       where tmpqc.read_name=seqread.trace_name and est.read_id=seqread.read_id;\
	insert into sgn_dev.qc_report (est_id, hqi_start, hqi_length) \
	       select est_id, hqi_start, hqi_length from tmpqc;\
	select count(*) from est join seqread using (read_id) join clone \
	       using (clone_id) join library using (library_id) where \
	       library_shortname='$(LIBRARY_NAME)' and flags = 0;\
	commit;\
	" | sed 's/--/\n/g;s/;/;\n/g' | \
	psql --variable AUTOCOMMIT=off --variable ON_ERROR_STOP=t --echo-all --dbname $(DBNAME) --host $(DBHOST) --user $(DBUSER)
#	psql --variable AUTOCOMMIT=off --variable ON_ERROR_STOP=t --echo-all --dbname cxgn --host hyoscine.sgn.cornell.edu


contaminant-scores.tab: PROG='\
	  open C, "/home/kreuter/screening-contaminants/contaminants.tab"\
	    or die "$$!"\;\
	  %flags\; while (<C>) {chomp\;split\;$$flags{$$_[1]}=$$_[0]}\
	  while (<STDIN>) {\
	    chomp\;split\;\
	    $$fn=`grep \">$$_[1]\" /home/kreuter/screening-contaminants/*|sed \"s/:.*//\"|xargs -n1 basename`\;\
	    chomp ($$fn)\;\
	    $$seqflags{$$_[0]}=$$seqflags{$$_[0]}|hex($$flags{$$fn})\;\
	  }\
	  map{printf("%s\t%s\n",$$_, $$seqflags{$$_})}keys(%seqflags)\
	'
contaminant-scores.tab: contaminant_screen.m8
	-rm contaminant-scores.tab.tmp contaminant-scores.tab
	perl -e $(PROG)  < contaminant_screen.m8 > contaminant-scores.tab.tmp
	mv contaminant-scores.tab.tmp contaminant-scores.tab
