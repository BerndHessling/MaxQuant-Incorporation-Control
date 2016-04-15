# MaxQuant-Incorporation-Control
R-script to check incorporation level of heavy SILAC amino acids


## Script description

This R-script can be used to calculate and plot incorporation rates of heavy SILAC amino acids like any type of labelled Arg, Lys or Leu. It can be used after analyzing raw files by MaxQuant and takes the MaxQuant output file "peptides.txt" as input.
The result file is a pdf with plotted incorporation rates for the single amino acids and all identified peptides.


## Prerequisits

* A successful MaxQuant analysis must be performed in advance with the correct labelling (Multiplicity 2, specified light and heavy label).

*	[R](https://cran.r-project.org/bin/windows/base/) needs to be installed (package was developed for version 3.2.4).


## Running the script

The script can be started without editing. In a first pop-up window the path to the combined/txt must be defined, in a second one a project name can be given and in a third the label-type can be defined (currently Arg&Lys and Leu labelling are supported).
The resulting pdf file will be created in the same combined/txt folder.
