-- -*- SQL -*-

SET SEARCH_PATH=sgn;
\set ON_ERROR_STOP on

\pset format unaligned
\pset tuples_only on
\timing

SELECT 'Dumping unigene build id ' || :ubid;  -- this ensures that :ubid is defined.

-- Make a table containing the sequence and quality scores for unigenes
-- in the current build.  Since this is a union, it could be done as two
-- queries, in case that's faster.
SELECT 'Accumulating sequence and quality score information.';
CREATE TEMPORARY TABLE unigene_sequence AS
       -- First select the sequences and scores for non-singletons
       SELECT unigene_id, seq, qscores AS qscore 
              FROM unigene 
	      JOIN unigene_consensi USING (consensi_id)
              WHERE unigene.consensi_id IS NOT NULL -- i.e., nr_members > 1
              AND unigene_build_id = :ubid
       -- and union it with trimmed sequences and scores for ESTs.
       -- (Trimming qscores is slow.  Maybe we ought to store offsets into
       -- qscore strings in qc_report?)
       UNION (SELECT unigene_id,
                     COALESCE(SUBSTRING(seq FROM CAST(hqi_start AS int) + 1
                                            FOR CAST(hqi_length AS int)), seq) AS seq,
                     COALESCE(array_to_string((string_to_array(qscore, ' '))
                                                [hqi_start+1:hqi_start+hqi_length], ' '),
                                              qscore) AS qscore
                      FROM unigene
                      JOIN unigene_member USING (unigene_id)
                      JOIN est USING (est_id)
                      JOIN qc_report USING (est_id)
                      WHERE consensi_id IS NULL -- (a.k.a. nr_members = 1)
                      AND unigene_build_id = :ubid);

-- Create a table containing the best hits for each unigene for each
-- blast_annotation the unigene has.
-- This is the slowest part of the program.
SELECT 'Accumulating blast annotations.';
CREATE TEMPORARY TABLE unigene_blast_annotations AS
       SELECT unigene_id, defline_id, blast_target_id, target_db_id,
              evalue, score
              FROM unigene u
              JOIN blast_annotations a ON (a.apply_id = u.unigene_id)
              JOIN blast_hits USING (blast_annotation_id) 
	      WHERE unigene_build_id = :ubid
	            AND blast_hit_id = (SELECT blast_hit_id
                                               FROM blast_hits h
	                                       WHERE h.blast_annotation_id =
	                                             a.blast_annotation_id
			                       AND score =
			                           (SELECT max (score)
                                                           FROM blast_hits h2
                                                           WHERE h2.blast_annotation_id =
                                                                 a.blast_annotation_id)
			                       LIMIT 1);

SELECT 'Formatting deflines for output.';
CREATE TEMPORARY TABLE unigene_defline AS
       SELECT unigene_id,
              '>SGN-U' || unigene_id || '\t' || g.comment || ' #' || build_nr ||
              ' [' || nr_members || ' ESTs aligned]' ||
              CASE WHEN arabidopsis_hit IS NOT NULL THEN arabidopsis_hit ELSE '' END ||
              CASE WHEN genbank_hit IS NOT NULL THEN genbank_hit ELSE '' END AS defline
              FROM groups g
              JOIN unigene_build ON (organism_group_id = group_id)
              JOIN unigene u USING (unigene_build_id) 
              LEFT JOIN (SELECT unigene_id, ' arabidopsis/peptide: ' || ad.defline ||
                                '(evalue: ' || evalue || ', score=' || score || ')'
                                AS arabidopsis_hit
                                FROM unigene_blast_annotations u1
                                LEFT JOIN blast_defline ad USING (defline_id)
                                WHERE u1.blast_target_id = 2)
                   AS aq USING (unigene_id)
              LEFT JOIN (SELECT unigene_id, ' genbank/nr: ' || gd.defline ||
                                '(evalue: ' || evalue || ', score=' || score || ')' 
                                AS genbank_hit
                                FROM unigene_blast_annotations u2
                                LEFT JOIN blast_defline gd USING (defline_id)
                                WHERE u2.blast_target_id = 1)
                   AS gq USING (unigene_id)
              JOIN unigene_sequence USING (unigene_id)
              WHERE unigene_build_id = :ubid;

CREATE TEMPORARY TABLE unigene_seq_dump AS
       SELECT defline, seq
              FROM unigene_defline
              JOIN unigene_sequence USING (unigene_id);

CREATE TEMPORARY TABLE unigene_qscore_dump AS
       SELECT defline, qscore
              FROM unigene_defline
              JOIN unigene_sequence USING (unigene_id);

SELECT 'Outputting to unigene-dump.fasta.';
\copy unigene_seq_dump to 'unigene-dump.fasta' with delimiter as '\n'
\! sed -i 's/\\t/        /' unigene-dump.fasta
SELECT 'Outputting to unigene-qscore-dump.fasta.';
\copy unigene_qscore_dump to 'unigene-qscore-dump.fasta' with delimiter as '\n'
\! sed -i 's/\\t/        /' unigene-qscore-dump.fasta
