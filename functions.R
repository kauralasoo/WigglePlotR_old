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

DrawChevrons <- function(start.pos, end.pos, exonStarts, exonEnds, exonCount, strand){
	# Function for drawing chevrons on the wiggle plot that show transcription direction.
	#
	# INPUT:
	# start.pos		start position of the transcript
	# end.pos		end position of the transcript
	# exonStarts	start coordinates of all exons
	# exonEnds		end coordinates of all exons
	# econCount		the total number of exons in transcript
	# strand		the strand of the transcript
	
	# OUTPUT
	# Add chevrons to the gene structure plot
	
	chevron.intervals <- 50;
	x <- seq(start.pos, end.pos, length.out=chevron.intervals+1)[2:chevron.intervals];
	width <- diff(x)[1] * 0.25;
	height <- 0.2;
	direction <- c("-"=-1, "+"=1)[strand];
	chevron.starts <- pmin(x, x-width*direction);
	chevron.ends <- pmax(x, x-width*direction);
	# If there are introns, draw chevrons there; if not, draw them in white on the exon
	if (exonCount > 1) {
		x <- x[!mapply(function(starts, ends) any((starts < exonEnds) & (ends > exonStarts)), chevron.starts, chevron.ends)];
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

DrawExonStructure <- function(bed.record, color, start.pos){
	# Function to draw the exon structure for single transcripts
	
	# INPUT
	# bed.record	the description of the transcript on bed format 
	# color			color of the exons
	
	# OUTPUT
	# Add exons from one trancript to the gene structure plot
	
	# Extract variables from bed record for convenience
	txStart <- bed.record$txStart;
	txEnd <- bed.record$txEnd;
	cdsStart <- bed.record$cdsStart;
	cdsEnd <- bed.record$cdsEnd;
	exonCount <- bed.record$exonCount;
	
	# Extract start positions of exons and CDS, converting from 0-based starts to 1-based starts
	exonStarts <- as.integer(strsplit(bed.record$exonStarts, ",")[[1]])+txStart+1;
	exonEnds <- as.integer(strsplit(bed.record$exonSizes, ",")[[1]])+exonStarts;
	cdsStartExon <- findInterval(cdsStart+1, exonStarts);
	cdsEndExon <- findInterval(cdsEnd, exonStarts);
	
	# Draw a line across the transcription unit
	lines(x=c(txStart+1, txEnd), y=c(0,0));
	# Draw any exons before the coding sequence
	i <- 1;
	while (i < cdsStartExon) {
		rect(xleft=exonStarts[i], xright=exonEnds[i], ybottom=-0.5, ytop=0.5, col=color, border = NA);
		if (i < exonCount) i <- i + 1;
	}
	# Draw the exon that intersects with the start of the coding sequence
	rect(xleft=exonStarts[i], xright=cdsStart+1, ybottom=-0.5, ytop=0.5, col=color, border = NA);
	# This line checks to see if the CDS ends before the exon does (in case CDS is in one exon)
	rect(xleft=cdsStart+1, xright=min(cdsEnd, exonEnds[i]), ybottom=-0.75, ytop=0.75, col=color, border = NA);
	if (i < exonCount) i <- i + 1;
	# Draw the exons that overlap completely with the coding sequence
	while (i < cdsEndExon) {
		rect(xleft=exonStarts[i], xright=exonEnds[i], ybottom=-0.75, ytop=0.75, col=color, border = NA);
		if (i < exonCount) i <- i + 1;
	}
	# Draw the exon that intersects with the end of the coding sequence
	rect(xleft=max(cdsStart, exonStarts[i]), xright=cdsEnd, ybottom=-0.75, ytop=0.75, col=color, border = NA);
	rect(xleft=cdsEnd, xright=exonEnds[i], ybottom=-0.5, ytop=0.5, col=color, border = NA);
	if (i < exonCount) i <- i + 1;
	# Draw any exons after the coding sequence
	while (i <= exonCount) {
		rect(xleft=exonStarts[i], xright=exonEnds[i], ybottom=-0.5, ytop=0.5, col=color, border = NA);
		i <- i + 1;
	}

	#Add name
	text(x = start.pos, y = -1, bed.record$name, cex = 2, pos = 4)
	
	#Draw chevrons
	DrawChevrons(txStart, txEnd, exonStarts, exonEnds, exonCount, bed.record$strand)
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

CreatePileups <- function(start.pos, end.pos, chrom, bamfiles, total.reads){
	#Function to retrieve reads from BAM files and create pileups.
	
	# INPUT
	# start.pos		start position of the longest transcript
	# end.pos		end position of the longest transcript
	# chrom			chromosome where the transcript is located
	# bamfiles		list of BAM files to be analyzed
	# total.reads	vector of total number of reads per each BAM file
	
	# OUTPUT
	# list(reads = reads, counts = counts) - List of two lists, showing the reads and read counts from each BAM files.
	
	n <- length(bamfiles);
	reads <- list();
	counts <- list();
	for (i in 1:n) {
		# Prepare a parameter object to look for any alignments falling within the given range
		param <- ScanBamParam(which=GRanges(seqnames=Rle(chrom),
						ranges=IRanges(start.pos, end.pos)),
				what=c("qwidth","pos"));
		# Retrieve the positions of all alignments falling within the given range
		reads[[i]] <- scanBam(bamfiles[i], param=param)[[1]];
		# Perform a pileup of the alignments
		if (length(reads[[i]]$pos)) {
			coordinates <- reads[[i]]$pos - reads[[i]]$pos[1] + 1;
			counts[[i]] <- rep(as.integer(0), length=max(coordinates)+33);
			for (j in 1:length(coordinates)) {
				start <- coordinates[j];
				end <- start + reads[[i]]$qwidth[j] - 1;
				counts[[i]][start:end] <- counts[[i]][start:end] + 1;
			}
		} else {
			counts[[i]] <- integer(0);
		}
	}
	
	# If total reads were requested, normalize all counts to million reads per sample
	if (!is.null(total.reads)) {
		for (i in 1:n) {
			counts[[i]] <- counts[[i]] / (total.reads[i] / 1E6);
		}
	}
	
	return(list(reads = reads, counts = counts))
}

DrawWigglePlots <- function(start.pos, end.pos, exonStarts, exonEnds, pileups, total.reads, sample.names, bg.colors, intron.color, exon.colors){
	
	reads = pileups$reads
	counts = pileups$counts
	n <- length(counts)
	exonPositions <- unlist(mapply(seq, exonStarts, exonEnds, SIMPLIFY=FALSE));
	
	
	# Get the maximum across all counts and all samples
	min.counts <- Inf;
	max.counts <- 0;
	for (i in 1:n) {
		nonzero <- which(counts[[i]] > 0);
		first.read <- reads[[i]]$pos[1];
		nonzero.exonic <- intersect(nonzero, exonPositions-first.read+1);
		max.counts <- max(max.counts, max(counts[[i]][nonzero.exonic]));
		min.counts <- min(min.counts, min(counts[[i]][nonzero]));
	}
	
	par(mar=c(0.5,5,0.5,1));
	for (i in 1:n) {
		# Set up each plot, scaled to the maximum count across all samples
		#ylab <- sprintf("%s%s", sample.names[i], ifelse(!is.null(total.reads), "\nRPKM", "\nreads"));
		ylab <- sprintf("%s", sample.names[i]);
		plot(x=NA, y=NA, xaxt="n", xlab=NA, ylab=ylab,
				cex.axis=2, cex.lab = 2, xlim=c(start.pos, end.pos), ylim=c(0, max.counts));
		# Draw a background over just the plotting area
		rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=par("usr")[3], ytop=par("usr")[4], col=bg.colors[i]);
		# If there are counts, plot them as individual line segments
		nonzero <- which(counts[[i]] > 0);
		first.read <- reads[[i]]$pos[1];
		nonzero.exonic <- intersect(nonzero, exonPositions-first.read+1);
		nonzero.intronic <- setdiff(nonzero, exonPositions-first.read+1);
		if (length(nonzero.intronic)) {
			segments(x0=first.read+nonzero.intronic-1, x1=first.read+nonzero.intronic-1,
					y0=0, y1=counts[[i]][nonzero.intronic],
					col=intron.color);
		}
		if (length(nonzero.exonic)) {
			segments(x0=first.read+nonzero.exonic-1, x1=first.read+nonzero.exonic-1,
					y0=0, y1=counts[[i]][nonzero.exonic],
					col=exon.colors[i]);
		}
	}
}

CreateFullGeneModel <- function(transcripts, bed){
	# Function to create the full gene model from the transcripts
	
	# INPUT
	# transcripts		vector of all transcript IDs
	# bed				BED file 
	
	# OUTPUT
	# list(start.pos = start.pos, end.pos = end.pos, exonStarts = exonStarts, exonEnds = exonEnds, exonCount = length(exonStarts))
	# start.pos		minimal start position of all transcripts
	# end.pos		maximal end position of all transcripts
	# exonStarts	start coordinates of exons
	# exonEnds		end cooridnates of exons
	# exonCount		total number of exons
	
	startPositions = list()
	endPositions = list()
	exonStarts = c()
	exonEnds = c()
	
	i = 0
	for (id in transcripts){
		i = i + 1
		bed.record <- subset(bed, name == id)
		
		if (nrow(bed.record) > 1){ #If multiple rows, choose first.
			print(paste("ERROR: Transcript", id, "has more than one record in BED file. Selecting first.", sep = " "))
			bed.record = bed.record[1,]	
		}
		
		startPositions = append(startPositions, bed.record$txStart + 1)
		endPositions = append(endPositions, bed.record$txEnd)
		
		if (i == 1){
			exonStarts <- as.integer(strsplit(bed.record$exonStarts, ",")[[1]])+bed.record$txStart + 1;
			exonEnds <- as.integer(strsplit(bed.record$exonSizes, ",")[[1]])+exonStarts;
		}
		else{
			newStarts <- as.integer(strsplit(bed.record$exonStarts, ",")[[1]])+bed.record$txStart + 1;
			newEnds <- as.integer(strsplit(bed.record$exonSizes, ",")[[1]])+newStarts;
			newStartIndices <- c(1:length(newStarts))[newStarts %in% exonStarts == FALSE]
			newEndIndices <- c(1:length(newStarts))[newStarts %in% exonStarts == FALSE]
			newExons = unique(c(newStartIndices, newEndIndices))
			exonStarts = c(exonStarts, newStarts[newExons])
			exonEnds = c(exonEnds, newEnds[newExons])	
		}		
	}
	start.pos = min(unlist(startPositions))
	end.pos = max(unlist(endPositions))
	
	return(list(start.pos = start.pos, end.pos = end.pos, exonStarts = exonStarts, exonEnds = exonEnds, exonCount = length(exonStarts)))
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


