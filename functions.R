# Helper functions for the wiggleplots.R file
# 
# Author: alasoo
###############################################################################

read.bedfile <- function (bedfile) {
	# A utility function for parsing a BED annotation file
	# Adam Gower, 2010
	#
	# INPUT
	# bedfile     The name of a BED file to parse
	#
	# OUPUT
	# A data frame with columns labeled according to the conventions of the BED file format as outline by UCSC:
	# http://genome.ucsc.edu/FAQ/FAQformat.html#format1
	
	bed.colClasses=c("character", rep("integer",2), rep("character",3),
			rep("integer",2), "character", "integer", rep("character",2),
			"integer", "integer", "integer", "double");
	bed.colnames <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand",
			"thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts",
			"depth", "bases.at.depth", "feature.size", "percent.coverage");
	
	# Get the number of elements in the first line
	n <- ncol(read.table(bedfile, nrows = 1));
	
	if (n == 10) {
		# If the BED file only contains the first 6 elements, followed by 4 more, assume that it is coverageBed output
		bed <- read.table(bedfile, stringsAsFactors=FALSE, colClasses=bed.colClasses[c(1:6,13:16)]);
		# Label the columns accordingly
		colnames(bed) <- bed.colnames[c(1:6,13:16)];
	} else {
		# Otherwise, just read the file using the normal colClasses
		bed <- read.table(bedfile, stringsAsFactors=FALSE, colClasses=bed.colClasses[1:n]);
		# Label the columns accordingly
		colnames(bed) <- bed.colnames[1:n];
	}
	# Return the data frame
	return(bed);
}

DrawChevrons <- function(transcript){
	# Function for drawing chevrons on the wiggle plot that show transcription direction.
	#
	# INPUT:
	# transcript	Object of the class "Transcript"
	
	# OUTPUT
	# Add chevrons to the gene structure plot
	
	chevron.intervals <- 50;
	x <- seq(transcript@txStart, transcript@txEnd, length.out=chevron.intervals+1)[2:chevron.intervals];
	width <- diff(x)[1] * 0.25;
	height <- 0.2;
	direction <- c("-"=-1, "+"=1)[transcript@strand];
	chevron.starts <- pmin(x, x-width*direction);
	chevron.ends <- pmax(x, x-width*direction);
	# If there are introns, draw chevrons there; if not, draw them in white on the exon
	if (transcript@exonCount > 1) {
		x <- x[!mapply(function(starts, ends) any((starts < transcript@exonEnds) & (ends > transcript@exonStarts)), chevron.starts, chevron.ends)];
		chevron.col <- "black";
	} else {
		chevron.col <- "white";
	}
	# Draw the chevrons
	if (length(x) > 0){
		segments(x0=x, x1=x-width*direction, y0=0, y1=-height, col=chevron.col);
		segments(x0=x, x1=x-width*direction, y0=0, y1=height, col=chevron.col);
	}
}

DrawExonStructure <- function(transcript, color, start.pos){
	# Function to draw the exon structure for single transcripts
	
	# INPUT
	# transcript	Object of the class "Transcript"
	# color			color of the exons
	
	# OUTPUT
	# Add exons from one trancript to the gene structure plot
	
	# Draw a line across the transcription unit
	lines(x=c(transcript@txStart+1, transcript@txEnd), y=c(0,0));
	# Draw any exons before the coding sequence
	i <- 1;
	while (i < transcript@cdsStartExon) {
		rect(xleft=transcript@exonStarts[i], xright=transcript@exonEnds[i], ybottom=-0.5, ytop=0.5, col=color, border = NA);
		if (i < transcript@exonCount) i <- i + 1;
	}
	# Draw the exon that intersects with the start of the coding sequence
	rect(xleft=transcript@exonStarts[i], xright=transcript@cdsStart+1, ybottom=-0.5, ytop=0.5, col=color, border = NA);
	# This line checks to see if the CDS ends before the exon does (in case CDS is in one exon)
	rect(xleft=transcript@cdsStart+1, xright=min(transcript@cdsEnd, transcript@exonEnds[i]), ybottom=-0.75, ytop=0.75, col=color, border = NA);
	if (i < transcript@exonCount) i <- i + 1;
	# Draw the exons that overlap completely with the coding sequence
	while (i < transcript@cdsEndExon) {
		rect(xleft=transcript@exonStarts[i], xright=transcript@exonEnds[i], ybottom=-0.75, ytop=0.75, col=color, border = NA);
		if (i < transcript@exonCount) i <- i + 1;
	}
	# Draw the exon that intersects with the end of the coding sequence
	rect(xleft=max(transcript@cdsStart, transcript@exonStarts[i]), xright=transcript@cdsEnd, ybottom=-0.75, ytop=0.75, col=color, border = NA);
	rect(xleft=transcript@cdsEnd, xright=transcript@exonEnds[i], ybottom=-0.5, ytop=0.5, col=color, border = NA);
	if (i < transcript@exonCount) i <- i + 1;
	# Draw any exons after the coding sequence
	while (i <= transcript@exonCount) {
		rect(xleft=transcript@exonStarts[i], xright=transcript@exonEnds[i], ybottom=-0.5, ytop=0.5, col=color, border = NA);
		i <- i + 1;
	}

	#Add name
	text(x = start.pos, y = -1, transcript@name, cex = 2, pos = 4)
	
	#Draw chevrons
	DrawChevrons(transcript)
}

InitializeExonStructurePlot <- function(start.pos, end.pos){
	# Function to initialize exon structure plot
	# INPUT
	# start.pos		start position of the longest transcript
	# end.pos		end position of the longest transcript
	
	#par(mar=c(3,5,0.5,1));
	par(mar=c(0.5,5,0.5,1));
	plot(x=NULL, y=NULL, yaxt="n", xaxt="n", xlab=NA, ylab=NA, xlim=c(start.pos, end.pos), ylim=c(-1,1));
}



DrawUniqueExons <- function(primaryTranscript, otherTranscripts, bed){
	
	#require(plotrix)
	#Extract information about primary transcript
	primaryRecord <- subset(bed, name == primaryTranscript)
	primaryExonStarts <- as.integer(strsplit(primaryRecord$exonStarts, ",")[[1]])+primaryRecord$txStart + 1;
	primaryExonEnds <- as.integer(strsplit(primaryRecord$exonSizes, ",")[[1]])+primaryExonStarts;
	print(length(primaryExonStarts))
	print(length(primaryExonEnds))
	
	uniqueStartIndices <- c(1:length(primaryExonStarts))
	uniqueEndIndices <- c(1:length(primaryExonEnds))
	
	#Go through all other transcripts
	for (id in otherTranscripts){
		bed.record <- subset(bed, name == id)
		newStarts <- as.integer(strsplit(bed.record$exonStarts, ",")[[1]])+bed.record$txStart + 1;
		newEnds <- as.integer(strsplit(bed.record$exonSizes, ",")[[1]])+newStarts;
		
		newStartIndices <- c(1:length(primaryExonStarts))[primaryExonStarts %in% newStarts == FALSE]
		newEndIndices <- c(1:length(primaryExonEnds))[primaryExonEnds %in% newEnds == FALSE]
		
		uniqueStartIndices <- intersect(uniqueStartIndices, newStartIndices)
		uniqueEndIndices <- intersect(uniqueEndIndices, newEndIndices)
	}
	
	uniqueIndices = unique(c(uniqueStartIndices,uniqueEndIndices))
	exonStarts <- primaryExonStarts[uniqueIndices]
	exonEnds <- primaryExonEnds[uniqueIndices]
	
	print(exonStarts)
	print(exonEnds)	
	if (length(exonStarts) > 0){
		for (i in c(1:length(exonStarts))){
			rect(xleft=exonStarts[i], xright=exonEnds[i], ybottom=-1, ytop=-0.77, col="red", border = NA);
		}		
	}	
}


