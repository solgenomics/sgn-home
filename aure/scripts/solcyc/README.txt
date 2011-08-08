=============================================================================
== ANNOTATION PIPELINE FOR SolCyc GENES                                    ==
==                                                                         ==
== Created by Lukas Mueller                                                ==
== Updated by Aureliano Bombarely                                          ==
==                                                                         ==
== [2011-08-04]                                                            ==
=============================================================================


This pipeline has as goal the creation of a pathologic file with metabolic
relations between genes and pathways. It will use as pathway information
the data contained in AraCyc.

The annotation process can be divided into 3 steps:

1) Blast homology search with three datasets:
   + Arabidopsis pep.
   + Swissprot
   + GenBank NR

2) Blast result filtering and attachment of the gene annotations 
   (using blast_result_handle_tool.pl).

3) Creation of the pathologic file (using generate_pathologic_file.pl script)

Follow the PathwayTools instructions to create or modify a pathway database.


Please use perldoc or -h option with each script for more detail.
