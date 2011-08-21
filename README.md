WigglePlotR
===========
A tool to create wiggle plots from RNA-Seq data.

Installation
------------
Download the WigglePlotR.R file, start R and then run:
	
	source("WigglePlotR.R")
	
Usage
-----
To create a plot, execute the following command:

	bedfile <- read.bedfile("path_to_bed_file.txt")
	WigglePlotR(ids, bamfiles, bedfile, total.reads=NULL, cex = 1, kernel.width = 50,
		exon.colors="black", intron.color="lightgray", bg.colors="transparent")

 * `ids` [character vector] -- IDs of the transcripts. Each ID must be present in the `bedfile`.

Acknowledgements
----------------
This script is based on the the wiggleplots.R script originally written by Adam Gower.

