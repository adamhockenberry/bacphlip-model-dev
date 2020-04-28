#!/bin/bash

hmmfile=../../Data/model_data/protein_domain_data/all_current_models.hmm

for faa_file in ../../Data/model_data/phage_data_nmicro2017/phage_fasta_files/*_6frame.faa; do
    outfile=${faa_file/fasta_files/hmmsearch_out}
    outfile=${outfile/_6frame.faa/.hmmsearch.out}
    echo $outfile
    hmmsearch $hmmfile $faa_file > $outfile 
    wait
done
