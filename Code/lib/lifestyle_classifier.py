import argparse
import pandas as pd
import os
import joblib

############################################################################
###Constants
############################################################################
############################################################################
############################################################################

###Command line arguments    
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file",\
        required=True, help="Should be a valid path to a .tsv datatable.")
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



###Load classifier model
clf = joblib.load('../Data/rf_highMinAJH.joblib')
###Load dataset
single_df = pd.read_csv(args.input_file, sep='\t', index_col=0)
###Predict
class_probs = clf.predict_proba(single_df)
###Write output
with open(args.output_file, 'w') as outfile:
    outfile.write('{}\t{}\n'.format('Lytic', 'Temperate'))
    outfile.write('{}\t{}\n'.format(class_probs[0][0], class_probs[0][1]))
exit(0)
