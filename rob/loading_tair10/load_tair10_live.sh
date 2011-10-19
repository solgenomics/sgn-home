set -e;
date;
gmod_bulk_load_gff3.pl  --gfffile TAIR10_GFF3_genes.gff --organism 'Arabidopsis thaliana' --fastafile  TAIR10_genome.fas --dbuser postgres --dbpass 'XXXXX' --dbhost db --dbname cxgn --noexon --remove_lock --recreate_cache;
date;
mx-run -I/home/rob/cxgn/bio-chado-loader/lib Bio::Chado::Loader::FASTA --db_dsn 'dbi:Pg:host=db;dbname=cxgn;user=postgres;password=XXXXXX' --type_name chromosome --organism_name 'Arabidopsis _ thaliana' TAIR10_genome.fas
date;
perl -MDBI -E "DBI->connect('dbi:Pg:host=db;dbname=cxgn;user=postgres;password=XXXXXX')->do(q|delete from public.feature where name like 'polypeptide-auto%' and organism_id = ( select organism_id from organism where species IN('thaliana','Arabidopsis thaliana'))|)"
date;
mx-run -I/home/rob/cxgn/bio-chado-loader/lib Bio::Chado::Loader::FASTA --db_dsn 'dbi:Pg:host=db;dbname=cxgn;user=postgres;password=XXXXXX' --type_regex 'type:(\S+)' --organism_name 'Arabidopsis _ thaliana' TAIR10_cdna_20101214.typed TAIR10_pep_20101214.typed