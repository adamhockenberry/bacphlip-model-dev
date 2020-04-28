#!/bin/bash

###The default program path assumes hmmsearch and python are globally accessible
HMMSEARCH_PATH="hmmsearch"
PYTHON_PATH="python"

###Identify the base directory of this program to easily point towards library and data files
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
###Location of the hmm file to use for hmmsearch
HMM_DB=$DIR/lib/prot_models.hmm

###Help function
function usage() {
    echo "Welcome to BACPHLIP"
    echo ""
    echo -e "\t-h --help"
    echo -e "\t--environment=$HMMSEARCH_PATH"
    echo -e "\t--db-path=$DB_PATH"
    echo ""
}

###Parsing command line arguments
while [ "$1" != "" ]; do
    PARAM=`printf "%s\n" $1 | awk -F= '{print $1}'`
    VALUE=`printf "%s\n" $1 | sed 's/^[^=]*=//g'`
    if [[ $VALUE == $PARAM ]]; then
        shift
        VALUE=$1
    fi
    case $PARAM in
        -i | --genome_fasta)
            GENOME_FASTA=$VALUE
            ;;
        -j | --proteins_fasta)
            PROTEINS_FASTA=$VALUE
            ;;
        -h | --help)
            usage
            exit
            ;;
        --hmmsearch-path)
            HMMSEARCH_PATH=$VALUE
            ;;
        --python-path)
            PYTHON_PATH=$VALUE
            ;;
        *)
            echo "ERROR: unknown parameter \"$PARAM\". Program is exiting."
            usage
            exit 1
            ;;
    esac
    shift
done

###Checking that required arguments are satisfied
if [ -z $GENOME_FASTA ] && [ -z $PROTEINS_FASTA ] ; then
        echo "You must provide the path to a genome fasta file (or amino acid fasta file) as input with the -i (or -j) flag. Exiting."
    exit 1
fi
###And that both weren't provided, which would be pointless
if [ ! -z $GENOME_FASTA ] && [ ! -z $PROTEINS_FASTA ] ; then
        echo "You must provide the path to EITHER a genome fasta file or amino acid fasta file as input with the -i OR -j flag. Not both. Exiting."
    exit 1
fi

#######################################################################
#######################################################################
###Running the actual pipeline
#######################################################################
#######################################################################
echo ""
###Run 6-fold amino acid prediction if the input is a genome file
if [ ! -z $GENOME_FASTA ] ; then
    echo "###############################################################"
    echo "Running six-frame translation of genome input file."
    PROTEINS_FASTA="$GENOME_FASTA.6frame"
    $PYTHON_PATH $DIR/lib/six_frame_translate.py -i $GENOME_FASTA -o $PROTEINS_FASTA
    if [[ $? = 0 ]]; then
        echo "Program call appears to be successful."
    else
        echo "Something went wrong (error code $?). Exiting."
        exit 1
    fi
    echo "Results written to $PROTEINS_FASTA"
    echo ""
fi

###Run HMMSEARCH
echo "###############################################################"
echo "Running hmmsearch on amino acid fasta file located at $PROTEINS_FASTA"
HMMSEARCH_OUT="$PROTEINS_FASTA.hmmsearch"
hmmsearch $HMM_DB $PROTEINS_FASTA > $HMMSEARCH_OUT
if [[ $? = 0 ]]; then
    echo "Program call appears to be successful."
else
    echo "Something went wrong (error code $?). Exiting."
    exit 1
fi
echo "Results written to $HMMSEARCH_OUT"
echo ""

####Compile HMMsearch results.
echo "###############################################################"
echo "Parsing hmmsearch results located at $HMMSEARCH_OUT"
HMMSEARCH_DF="$HMMSEARCH_OUT.df"
$PYTHON_PATH $DIR/lib/process_hmmsearch.py -i $HMMSEARCH_OUT -o $HMMSEARCH_DF
if [[ $? = 0 ]]; then
    echo "Program call appears to be successful."
else
    echo "Something went wrong (error code $?). Exiting."
    exit 1
fi
echo "Results written to $HMMSEARCH_DF"
echo ""

###Run predictor
echo "Running classifier on $HMMSEARCH_DF"
PREDICTION_FILE="$PROTEINS_FASTA.bacphlip"
$PYTHON_PATH $DIR/lib/lifestyle_classifier.py -i $HMMSEARCH_DF -o $PREDICTION_FILE
if [[ $? = 0 ]]; then
    echo "Program call appears to be successful."
else
    echo "Something went wrong (error code $?). Exiting."
    exit 1
fi
echo ""
echo "Full program completed! Hopefully without any errors noted above!"
echo ""

