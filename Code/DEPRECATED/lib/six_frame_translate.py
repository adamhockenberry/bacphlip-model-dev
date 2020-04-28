from Bio import SeqIO
import argparse
import os

"""
Note that this code will fail on phages with non-standard genetic codes!
"""
############################################################################
###Constants
############################################################################
min_prot_length = 40
############################################################################
############################################################################

###Command line arguments    
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file",\
        required=True, help="Should be a valid path to a single genome (nucleotide) fasta file containing only 1 record/contig.")
parser.add_argument("-o", "--output_file",\
        required=True, help="Can be any valid path that does not currently exist.")
args = parser.parse_args()

###Basic error testing of input/output files
if not os.path.exists(args.input_file):
    print("Input file does not appear to exist. Exiting.")
    exit(1)
if os.path.exists(args.output_file):
    print("Specified output file ({}) appears to already exist. Remove before running program again. Exiting.".format(args.output_file))
    exit(1)

###More error testing to ensure that the input file contains data
nt_records = list(SeqIO.parse(args.input_file, 'fasta'))
if len(nt_records) == 1:
    nt_record = nt_records[0]
elif len(nt_records) == 0:
    print('Input fasta file appears to be empty. Exiting.')
    exit(1)
else:
    print('Input fasta file appears to contain more than one sequence record. Exiting.')
    exit(1)

###Run basic code
genome_id = nt_record.id
prots = []
for i in [0, 1, 2]:###Each of the three reading frames
    ###Translate the regular strand
    tempy = nt_record[i:].translate()
    ###Split according to stop codons
    seq = str(tempy.seq).split('*')
    ###Append the sequences if they are longer than the defined minimum
    for j in seq:
        if len(j) >= min_prot_length:
            prots.append(j)
    ###Repeat for the reverse complement
    tempy = nt_record.reverse_complement()[i:].translate()
    seq = str(tempy.seq).split('*')
    for j in seq:
        if len(j) >= min_prot_length:
            prots.append(j)

###Write out the results
with open(args.output_file, 'w') as outfile:
    for i, seq in enumerate(prots):
        outfile.write('>{}_{}\n{}\n'.format(genome_id, i, seq))

exit(0)
