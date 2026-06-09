#!/bin/bash

# predict B-cell epitopes by discotope3
python /path/to/discotope3_web/discotope3/main.py --cpu_only --pdb_or_zip_file /path/to/input_protein_structure.pdb --out_dir /path/to/output_dir

# predict B-cell epitopes by BepiPred
python /path/to/bepipred3_CLI.py -i /path/to/input_protein_sequence.fa -o /path/to/output_filename -pred vt_pred

# predict B-cell epitopes by ElliPro
ellipro_jar="/path/to/ElliPro/ElliPro.jar"
java -jar "$ellipro_jar" \
        --input-file "/path/to/input_protein_structure.pdb" \
        --min-score 0.5 \
        --max-dist 6 \
        --output "/path/to/epitope_output_filename" \
        --table "/path/to/protrusion_output_filename"

# predict T-cell epitopes by netmhcii-ba/el
method=netmhciipan_ba # netmhciipan_el
allele="HLA-DRB1*01:01,HLA-DRB1*03:01,HLA-DRB1*04:01,HLA-DRB1*04:05,HLA-DRB1*07:01,HLA-DRB1*08:02,HLA-DRB1*09:01,HLA-DRB1*11:01,HLA-DRB1*12:01,HLA-DRB1*13:02,HLA-DRB1*15:01,HLA-DRB3*01:01,HLA-DRB3*02:02,HLA-DRB4*01:01,HLA-DRB5*01:01,HLA-DQA1*05:01/DQB1*02:01,HLA-DQA1*05:01/DQB1*03:01,HLA-DQA1*03:01/DQB1*03:02,HLA-DQA1*04:01/DQB1*04:02,HLA-DQA1*01:01/DQB1*05:01,HLA-DQA1*01:02/DQB1*06:02,HLA-DPA1*02:01/DPB1*01:01,HLA-DPA1*01:03/DPB1*02:01,HLA-DPA1*01:03/DPB1*04:01,HLA-DPA1*03:01/DPB1*04:02,HLA-DPA1*02:01/DPB1*05:01,HLA-DPA1*02:01/DPB1*14:01" # or other alleles
python /path/to/mhc_ii/mhc_II_binding.py $method $allele  /path/to/input_protein_sequence.fa > /path/to/output_netmhciipan_ba_result.csv
