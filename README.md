WigglePlotR
===========
A tool to create wiggle plots from RNA-Seq data.

Installation
------------
Download the WigglePlotR.R file, start R and then run:
	
	source("WigglePlotR.R")
	
Usage
-----
You first have to import a bedfile using the following command:

	bedfile <- ReadBedFile("path_to_bed_file.txt")

Then you can create a wiggle plot by running:

	WigglePlotR(ids, bamfiles, bedfile, total.reads=NULL, cex = 1, kernel.width = 0
		exon.colors="black", intron.color="lightgray", exon.events = TRUE)

 * `ids` [character vector] -- IDs of the transcripts. Each ID must be present in the `bedfile`.
 * `bamfiles` [charcter vector] -- Full path to each BAM file containing the aligned reads.
 * `bedfile` [data.frame] -- BED file imported with `ReadBedFile()`. Specifies the structure of transcripts.
 * `total.reads` [numeric vector] -- Optional, total number of in each BAM file. When specified, used to normalize the y-axis to RPKM values.
 * `cex` [number] -- Scaling factor for text and axes. Default is `cex = 1`
 * `kernel.width` [number] -- Width of the kernel used to smoothen the wiggle plots. Default value is `kernel.width = 1` which means no smoothing.
 * `exon.colors` [character vector] -- Color(s) for the reads overlapping exons.
 * `intron.colors` [character vector] -- Color(s) for the reads overlapping introns.
 * `exon.events` [boolean] -- If TRUE, first two transcripts will be colored differently from others and unique exons will be colored red.

Acknowledgements
----------------
This script is based on the the wiggleplots.R script originally written by Adam Gower.

