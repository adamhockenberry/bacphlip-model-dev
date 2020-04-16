# BacPhLiP - BacterioPhage Lifestyle Predictor

Welcome! The most important thing to do if you're new here is to check out our manuscript for details of the tool that you're looking at:


# Preparation

To use this tool you'll need it on your computer. For the latest version, directly download from github or go ahead and open a terminal on your computer, navigate to a directory that you'd like BacPhLiP to be installed in and type:

    `git clone ....`

For past versions, check out the releases page on the github. Download the selected release, move it somewhere logical on your computer, and unzip it. 

That's it! Installation complete. No really. BacPhLiP is actually just a bash script (`bacphlip.sh`) that strings together a few different sequential calls to various programs (most of which are provided in this repository [see `lib/`] and can be used individually by advanced users). As such, there is no need for any complicated installations (though, see the Warnings section because you might need to install some other tools first!). 

## Warnings / dependencies / disclaimers
Depending on how you tried to download/clone BacPhLiP you might have already encountered a problem. Which is to say, this is a terminal (also called command-line) based program meant for use on Linux/Unix machines. It may work on Windows terminals, but no promises are made. In fact, I'd appreciate any users getting in touch who know what (if any) changes I'd need to make in order to make the program accessible to the command-line Windows community. 

**Additionally, while this repository should provide **most** of what you need. It does have it's own set of requirements**

Most notably among these requirements is that BacPhLiP relies on the HMMER suite of tools to search a given genome for various protein domains. If you can't type "hmmsearch" into your terminal and see a program yelling at you for inputs, you're going to have a problem. But fear not! Go ahead and install HMMER via this link (I for one find it useful to use Anaconda to manage package installs and HMMER can be found on Anaconda and seemlessly installed, see this link). If you have (or have just gotten) the HMMER suite of tools but prefer to keep them referenced in a local directory (via their absolute path), this path can be provided as a command line argument (discussed later). This is also a good time to let you know that the tool was developed using HMMER version 3.3.x and no promises are made about compatibility with past or future versions of HMMER. 

Aside from `hmmsearch`, the other scripts that are necessary to analyze your genome rely on you having `python` (version 3.x) successfully installed (and accessible but again, the absolute path to your python install can be manually provided if you insist). In addition to so-called "base" `python`, BacPhLiP relies on a few relatively common packages. See the requirements.txt for a list of them all. If you don't have one of these requirements satisfied, you'll likely encounter errors. With the possible exception of `biopython`, if you use python for literally any kind of data analyis you'll likely have most of these packages installed already. And `biopython` is relatively wide-spread and easy to use as well, so make sure you have it (see this link). 

Also, I'm just putting this here to avoid any easy-to-catch errors, but if you're new to the terminal make sure that none of the files you want to run have spaces in their names. It'll *probably* wreak havoc and it's going to save you a lot of headaches in the future to get used to using under-scores "_" instead of spaces. 

# Usage

Finally, you've made it. Your best place to start is to check out the help page out by navigating to the BacPhLiP directory and typing (omit the `$ `, this just signifies that it's a terminal based command):

    `$ bash bacphlip.sh -h` (note that you might be able to get away with the shortcut `bacphlip.sh` depending on the locaiton of your bash environment)

and hitting "Enter". Ideally, you'll see a helpful menu of command line options.

After you check that out, the real test of your system and installation of the program dependencies will come by entering:

    `$ bash bacphlip.sh -i examples/genome_example.fasta`

You should see a few comments popping up as the program keeps you informed of its progress at the various steps of analysis. Those comments *should* make sense and not end in anything resembling an error message. You'll know that everything is successful if you go into the `examples/` folder and see a few different files:

    `genome_example.fasta.6frame` is a fasta file containing all the amino acid sequences above a certain length threshold in all 6 possible translational reading frames for the given genome input
    `genome_example.fasta.6frame.hmmsearch` is a large output file from the `hmmsearch` program testing your genome for the presence of various protein domains
    `genome_example.fasta.6frame.hmmsearch.df` simplifies those `hmmsearch` results into a tiny tab-separated file describing the presence (1) / absence (0) of the various domains
    `genome_example.fasta.6frame.bacphlip` is our final file, a simple tab-delimited text file showing the probabilities of the input genome as being either Lytic or Temperate. In the case of this example, it should be Lytic with overwhelming probability). Basically whichever category has a probability > 0.5 is the programs predictions for the lifestyle of your phage.

If you don't have all of those files (or they're empty, or some other problem...), it's likely that you actually had some errors that were reported and should guide you to which step of the pipe-line failed so that you can interrogate why. But it's also possible that I didn't catch some possible errors, and of course, keep me informed if this is the case so I can make the program better for you all.

# Final thoughts
And that's basically it! You should be able to run the `bacphlip.sh` command with any valid path to a genome fasta file as input (or alternatively a fasta file of the proteins, provided with the -j flag). As you can see from the example, the intermediate and final files of interest are automatically created in the same directory as the input file with the various names of the files appended to them. If you have a batch of files that you want to run BacPhLiP on, at the moment you'll have to just write your own bash `for` loop that runs `bacphlip.sh -i xyz` with the various xyz input files.

If you have any issues, questions, omments, or concerns feel free to file them in the github repository or to reach out to me directly.

--Adam J Hockenberry
