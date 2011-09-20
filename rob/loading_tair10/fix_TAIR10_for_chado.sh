perl -i -pE 's/^>/>TAIR10/; s/Chr(\d+)/sprintf("ch%02d",$1)/e' TAIR10_genome.fas
perl -i -pE 's/(^|ID=|Name=)Chr(\d+)/sprintf("${1}TAIR10ch%02d",$2)/eg' TAIR10_GFF3_genes.gff
perl -i -pE 's/ChrM\b/TAIR10mitochondria/g; s/ChrC\b/TAIR10chloroplast/g' TAIR10_GFF3_genes.gff
perl -i -pE 's/\tmRNA_TE_gene\t/\tmRNA\t/' TAIR10_GFF3_genes.gff
perl -i -pE 's/\tprotein\t/\tpolypeptide\t/' TAIR10_GFF3_genes.gff
perl -i -pE 's/Derives_from/transposable_element_name/ if /\ttransposable_element_gene\t/' TAIR10_GFF3_genes.gff
gff3_insert_sync_directives  -i TAIR10_GFF3_genes.gff
