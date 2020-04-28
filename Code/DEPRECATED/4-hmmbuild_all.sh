#!/bin/bash

for infile in ../../Data/model_data/protein_domain_data/domain_alignments_and_hmms/*.afa; do
    outfile=${infile/.afa/.hmm}
    hmmbuild $outfile $infile
    wait
done

cat ../../Data/model_data/protein_domain_data/domain_alignments_and_hmms/*.hmm > ../../Data/model_data/protein_domain_data/all_current_models.hmm 
