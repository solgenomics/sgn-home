rollback;
begin;

create temporary table new_notes as
    select outer_f.feature_id as feature_id
           , ( select cvterm_id from cvterm where name = 'Note' limit 1 ) as type_id
           , array_to_string( array(
                select fp.value
                from featureprop fp
                join feature f using(feature_id)
                where fp.type_id IN( select cvterm_id from cvterm where name = 'Note' )
                  and f.feature_id = outer_f.feature_id
             ), ',' ) as new_note
    from feature outer_f
    where outer_f.organism_id = (select organism_id from organism where species = 'Solanum tuberosum')
      and 1 < ( select count(*)
                from featureprop
                where feature_id = outer_f.feature_id
                  and type_id IN( select cvterm_id from cvterm where name = 'Note' )
              );

select * from new_notes limit 10 offset 100;

delete from featureprop fp using new_notes n
  where fp.feature_id = n.feature_id
    and fp.type_id IN( select cvterm_id from cvterm where name = 'Note' );

insert into featureprop ( feature_id, type_id, value )
     select feature_id, type_id, new_note
       from new_notes;

commit;