import argparse
import pandas as pd
import os
from Bio import SearchIO

############################################################################
###Constants
############################################################################
total_protein_models = 371
############################################################################
############################################################################

###Command line arguments    
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file",\
        required=True, help="Should be a valid path to a hmmsearch output file.")
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

with open(args.input_file, 'r') as infile:
    results = list(SearchIO.parse(infile, 'hmmer3-text'))
    simple_res = []
    for i in results:
        if len(i.hits) > 0:
            simple_res.append((i.id, 1))
        else:
            simple_res.append((i.id, 0))
if len(simple_res) != total_protein_models:
    print('Appears to be an error, too many or too few results given expected number of protein models tested. Exiting')
    exit(1)

#single_df = pd.DataFrame(dict(simple_res), index=[index])

exit(0)
