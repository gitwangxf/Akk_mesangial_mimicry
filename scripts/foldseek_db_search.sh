#!/bin/bash

epi_pdb_dir=/path/to/Akk_epitopes_pdb
epi_db=/path/to/Akk_epitope_database
mes_pdb_dir=/path/to/mesangial_protein_pdb
mes_db=/path/to/mesangial_database
home_dir=/path/to/home_dir
tag=foldseek_result_tag


foldseek createdb $epi_pdb_dir $epi_db
foldseek createdb $mes_pdb_dir $mes_db
foldseek easy-search $epi_db $mes_db $home_dir/foldseek_result/$tag.e0.1.out tmp --threads 24 --num-iterations 3 --exhaustive-search -e 0.1 --max-seqs 5000 --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage"
