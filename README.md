# Compiling and/or updating your own model will require the steps listed here. 

## Project overview
This repository covers the designing and building of BACPHLIP (a bacteriophage lifestyle prediction tool). 

Users looking to classy the lifestyle of their favorite phage should not be here and should instead navigate directly to the separate BACPHLIP repository:
    xxx

The purpose of *this* repository is for full scientific transparency with regard to the associated BACPHLIP publication, which can be found at:
    xxx

Additionally, this repository can be cloned by others interested in building their own machine learning pipeline based off of this dataset and can also be used to update the protein domain database that underlies BACPHLIP for future releases (note that the updating is manual and BACPHLIP is fully distinct from this repository). 

## Usage notes
The pipeline consists of two initially distinct steps (1 and 2) that deal with compiling a set of phage and a set of protein domain data, respectively. Step 1 pertains to one particular phage training data set and will probably not need to be run ever again unless users simply want to replicate the results in their own hands using the associated data files (more on that in teh Requirements section). As a better / more expansive training set of phage data one day becomes available, it is likely that processing the data for this hypothetical set will require different code from what is used in Steps 1 and 2, but that code should essentially do the same thing and should not be too difficult to integrate within this framework.

By contrast, step 2 may be run to develop future releases of BACPHLIP since this step relies on the conserved domain database (which will be updated and grow over time). Thus, step 2 onwards may result in building an improved classifier that can serve as the basis for a new version of BACPHLIP even in the absence of more phage training data.

If you're not re-downloading any data files (this can take some time and we encourage you to follow the link to our archived data repository instead, see below), the whole pipeline can be run top-to-bottom on most common laptops (typing on a 2013 macbook pro) in an hour or two. 

## Requirements
Users will need several items if they wish to fully step through all of the associated code. In addition to python and several common packages that should be obvious within the notebooks, users will need to:

    1. Install HMMER3 (available from xxxx)
    2. Download the full associated dataset and unzip it within this repository (xxxx)
    3. A small (and optional) portion of the code also relies on fastANI (xxx) to cluster sequences. This is used for a more stringent test of model accuracy but is not strictly necessary

## Usage walkthrough
1. (ignore if only updating conserved domains) Step through `1-compile_phage_training_data.ipynb`. This should run with no errors and if you have downloaded/cloned the full repository this should do basically nothing except replace/re-write some files that were already there.

2. (optional: Download and replace the `cddid.tbl` file contained within `../Data/protein_domain_data/` to use an updated database. See link in the referenced notebook that follows.) Next, step through the notebook titled `2-compile_search_families.ipynb`. Depending on your goals here, you *may* want to ensure that you've deleted all existing `.afa` (and `.hmm` files, if they've been built) or alter the code accordingly to place an updated database in a new folder. In default formatting, the code will not re-download existing `.afa` files even though more recent updates to CDD may expand the size/quality of these files. I am of the opinion that this should be done more sporadically, on the order of every few months or yearly.

3. Now step through `3-hmmer_time.ipynb`. This is the part of the pipeline where you'll face a mountain of errors if you don't have HMMER3 installed. This will first create `.hmm` files for each of the protein domain `.afa` files. To check that the HMMER3 software is installed and accessible just type `hmmbuild` in your terminal. You should **not** see any variation of "command not found". The code assumes that the HMMER3 installation is accessible in your system path, but experienced useres may with to alter this code in order to point to a local install. This notebook then uses `hmmsearch` to create output files (containing information for all of the protein domains) for each of the viral sequences. This may take an hour or two to run. Note that I'm also running `hmmsearch` with defaults and it may be wise to tweak some sensitivity parameters at a future date/update. 

4. Now we're ready for the fun part and can step through the notebook `4-build-ml-classifier.ipynb`. This compiles full data tables from the `hmmsearch` output, splits data into training and testing sets, and ultimately fits (then saves) a Random Forest classifier.

5. To assess how well it all works, users can finally step through `5-assess-classifier.ipynb` (part of this code will require stepping through `cluster_seqs.ipynb` as well).
