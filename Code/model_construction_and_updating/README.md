# Compiling and/or updating your own model will require the steps listed here. 


Steps 1 and 2 pertain to one particular phage training data set and will probably not need to be run ever again. As a better / more expansive training set one day becomes available, it is likely that processing the data for this hypothetical set will require different code from what is used in Steps 1 and 2, but that code should essentially do the same thing.

If you're not re-downloading any data files (this can take some time), the whole pipeline can be run top-to-bottom on most common laptops (typing on a 2013 macbook pro) in an hour or two. Will re-iterate this later, but be sure that you have HMMER installed and accessible in the system path if you're planning on running the whole pipeline.  

1. (ignore for updating conserved domains) Step through `1-compile_phage_training_data.ipynb`. This should run with no errors and if you have downloaded/cloned the full repository this should do basically nothing except replace/re-write some files that were already there.

2. (ignore for updating conserved domains) In the terminal, run the command `python translate_6_frame.py`. This will take a bit of time, but you're actually just re-writing files that should again already exist for this particular phage dataset.

3. (optional) Download and replace the `cddid.tbl` file contained within `../../Data/model_data/protein_domain_data/` to use an updated database. See link in notebook below.

3. Next, step through the notebook titled `2-compile_search_families.ipynb`. Depending on your goals here, you *may* want to ensure that you've deleted all existing `.afa` (and `.hmm` files, if they've been built) or alter the code accordingly to place an updated database in a new folder. In default formatting, the code will not re-download existing `.afa` files even though more recent updates to CDD may expand the size/quality of these files.

4. In the terminal, run `bash hmmbuild_all.sh`. This will create `.hmm` files for each of the protein domain `.afa` files. It is likely that if you have any issues running through this code it will occur here. This portion of the pipeline relies on HMMER, so to check if you have this software installed and accessible just type `hmmbuild` in your terminal. You should **not** see any variation of "command not found". If you do see this message, you'll need to install HMMER (see: http://hmmer.org/, but note that if you are using anaconda you might want to instead visit: https://anaconda.org/bioconda/hmmer). The installation should be accessible in your system path, but you may alter `hmmbuild_all.sh` to point to a local install. Lastly, note that the final line of this command will concatenate all of the `.hmm` files into a single file. If you change any paths or want to change any file names, ensure that you do so here as well.

5. We'll return to HMMER so this next step will again require that it is properly installed. In the terminal run `bash hmmsearch_all.sh` to create `hmmsearch` output files (containing information for all of the protein domains) for each of the viral sequences. You may edit this file if you want to use a different `.hmm` file other than `all_current_models.hmm` or if you want to scan any other file full of protein sequences in the fasta file format. This may take an hour or two to run but since everything that we're working with is fairly small it did not seem worth trying to parallelize any computations. Note that I'm also running `hmmsearch` with defaults and it may be wise to tweak some sensitivity parameters. 

6. Now we're ready for the fun part and can step through the notebook ``.
