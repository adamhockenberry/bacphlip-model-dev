#!/bin/bash

for infile in ./examples/*; do
    if [[ $infile == ./examples/genome_example.fasta ]]; then
        continue
    elif [[ $infile == ./examples/aa_example.faa ]]; then
        continue
    else
        rm -i $infile
    fi
done
