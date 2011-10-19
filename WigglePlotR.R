
WigglePlotR <- function(ids, bamfiles, bedfile, total.reads=NULL, cex = 1, kernel.width = 1,
		exon.colors=rgb(99/255,99/255,99/255), intron.color=rgb(189/255,189/255,189/255), exon.events = FALSE) {
	# Function for creating "wiggle plots" from BAM alignment files and BED annotation files
	# Originally written by Adam Gower, modified by Kaur Alasoo.
	#
	# INPUT
	# ids             A character vector of IDs that correspond to features in the BED file, e.g., c("NM_001234")
	# bamfiles        A character vector of BAM file names, one for each alignment that will be used to generate a wiggle plot
	# bedfile         The BED file as data.frame that will be used to determine the chromosomal coordinates of the features in 'ids'
	# cex             How many times character size should be increased? Usful when creating jpg images.
	# kernel.width    If > 1, kernel smoothing of the read profiles will be used.
	# total.reads     An integer vector of the total number of reads for each sample; if supplied, used to scale the y-axes of plots
	# exon.colors     A character vector of color names with which to draw base positions that overlap with exonic regions
	# intron.color    A color name with which to draw base positions that overlap with intronic regions
	# exon.events     If TRUE, first two transcripts will be colored differently from others and unique exons will be colored red.
	#
	# OUTPUT
	# Draws a single pane for each feature in 'ids', containing a single wiggle plot for each alignment,
	#     as well as a representation of the structure of the gene at the bottom of the page
	
	# Need the Rsamtools package to read BAM files
	require(Rsamtools)
	
	#define colors
	transcript.colors = c(rgb(55/255,126/255,184/255), rgb(127/255,205/255,187/255)) #Blue and green
	unique.exon.color = rgb(228/255, 26/255, 28/255) #Red
	
	if (!is.null(names(bamfiles))) {
		# Get the sample names from the names of the filename vector, if available
		sample.names <- names(bamfiles);
	} else {
		# Otherwise, extract the names of the files after the path and before the .bam extension
		sample.names <- sub("(^.*/)*(.*)\\.bam$","\\2", bamfiles);
	}
	
	#Read all transcripts from the BED file
	transcript.list = ReadTranscriptsFromBed(ids, bedfile)	
	#Create full gene model (contains all exons)
	full.gene.model = CreateFullGeneModel(transcript.list)
	
	#Retrieve alignments and create pileups
	pileups = CreatePileups(full.gene.model, bamfiles, total.reads, kernel.width)
	
	# Create the plot layout
	n <- length(bamfiles)
	m <- length(transcript.list)
	layout(matrix(1:(n+m),n+m,1), heights = c(rep(4,n),rep(1,m)))
	par(mar=c(2,50,2,1), bg="transparent")
	
	#Draw wiggle plots
	if (length(exon.colors) == 1) exon.colors <- rep(exon.colors, n)
	DrawWigglePlots(full.gene.model, pileups, total.reads, sample.names,
			intron.color, exon.colors, cex = cex)
	
	#Draw exon structures
	first <- TRUE
	i = 0
	for (transcript in transcript.list){
		i = i + 1
		#Initialize plot and draw exon structure
		par(mar=c(0.5,5,0,1))
		plot(x=NULL, y=NULL, yaxt="n", xaxt="n", xlab=NA, ylab=NA, 
				xlim=c(full.gene.model@txStart, full.gene.model@txEnd), 
				ylim=c(-1,1), frame.plot = FALSE)
		
		#Choose the correct color for transcripts
		if (exon.events){
			if (i < 3){ color = transcript.colors[1]} else { color = transcript.colors[2] }
		}
		else{
			color = transcript.colors[1]
		}
		
		DrawExonStructure(transcript, color, full.gene.model, cex)
		
		#Mark unique exons on the first transcript
		if (first){
			DrawUniqueExons(transcript.list, unique.exon.color, exon.events)
			first <- FALSE
		}
	}
}

#### Helper functions and data structures ####

setClass("Transcript", representation(
				#Class to store transcript information
				name = "character",
				txStart = "numeric",
				txEnd = "numeric",
				cdsStart = "numeric",
				cdsEnd = "numeric",
				exonCount = "numeric",
				exonStarts = "numeric",
				exonEnds = "numeric",
				cdsStartExon = "numeric",
				cdsEndExon = "numeric",
				strand = "character",
				chrom = "character")
)

#Define a constructor for the class that loads tanscript information from a BED file
setGeneric("Transcript", function(bed.record) standardGeneric("Transcript"))
setMethod("Transcript", signature(bed.record = "data.frame"), function(bed.record){
			#Calculate values based on the bed record
			exonStarts <- as.integer(strsplit(bed.record$exonStarts, ",")[[1]])+bed.record$txStart+1
			exonEnds <- as.integer(strsplit(bed.record$exonSizes, ",")[[1]])+exonStarts
			cdsStartExon <- findInterval(bed.record$cdsStart+1, exonStarts)
			cdsEndExon <- findInterval(bed.record$cdsEnd, exonStarts)

			transcript = new("Transcript")
			transcript@name = bed.record$name
			transcript@strand = bed.record$strand
			transcript@txStart = bed.record$txStart + 1
			transcript@txEnd = bed.record$txEnd
			transcript@cdsStart = bed.record$cdsStart + 1
			transcript@cdsEnd = bed.record$cdsEnd
			transcript@exonStarts = exonStarts
			transcript@exonEnds = exonEnds
			transcript@cdsStartExon = cdsStartExon
			transcript@cdsEndExon = cdsEndExon
			transcript@exonCount = bed.record$exonCount
			transcript@chrom = bed.record$chrom
			return(transcript)
		})

ReadBedFile <- function (bedfile) {
	# A utility function for parsing a BED annotation file
	# Adam Gower, 2010
	#
	# INPUT
	# bedfile     The name of a BED file to parse
	#
	# OUTPUT
	# A data frame with columns labeled according to the conventions of the BED file format as outline by UCSC:
	# http://genome.ucsc.edu/FAQ/FAQformat.html#format1
	
	bed.colClasses=c("character", rep("integer",2), rep("character",3),
			rep("integer",2), "character", "integer", rep("character",2),
			"integer", "integer", "integer", "double");
	bed.colnames <- c("chrom", "txStart", "txEnd", "name", "score", "strand",
			"cdsStart", "cdsEnd", "itemRgb", "exonCount", "exonSizes", "exonStarts",
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
	
	# If this is a 'short' BED file with only 6 columns, use the data to add extra columns to the full BED file specification 
	if (ncol(bed) == 6) {
		bed <- cbind(bed, data.frame(thickStart=bed$chromStart, thickEnd=bed$chromEnd, itemRgb=0,
						blockCount=1, blockSizes=sprintf("%d,", bed$chromEnd-bed$chromStart), blockStarts="0,",
						stringsAsFactors=FALSE));
	}
	
	# Keep only the rows with the standard chromosome IDs and remove duplicate rows.
	bed <- subset(bed, chrom %in% sprintf("chr%s", c(as.character(1:22),"X","Y","M")))
	bed <- unique(bed)
	
	# Return the data frame
	return(bed)
}

ReadTranscriptsFromBed <- function(transcript.ids, bed.file){
	# Read transcripts from BED file and convert them to a list of Transcript objects.
	# 
	# INPUT
	# transcript.ids	vector of transcript ids
	# bed.file			BED file
	#
	# OUTPUT
	#transcript.list	list of transcript objects
	
	transcript.list = list()
	i = 1
	for (transcript.id in transcript.ids){
		bed.record = subset(bed.file, name == transcript.id)
		#Skip transcript if it's not found in the BED file.
		if (nrow(bed.record) == 0){
			print(paste("ERROR: Transcript", transcript.id , "has no record in the BED file. Skipping.", sep = " "))
			next
		}
		#If there are multiple rows (same transcript on different strands), then take the first one.
		if (nrow(bed.record) > 1){
			print(paste("ERROR: Transcript", transcript.id, "has more than one record in BED file. Selecting first.", sep = " "))
			bed.record = bed.record[1,]	
		}
		transcript = Transcript(bed.record)
		transcript.list[[i]] = transcript
		i = i + 1
	}
	return(transcript.list)
}

CreateFullGeneModel <- function(transcript.list){
	# Function to create the full gene model from the transcripts.
	
	# INPUT
	# transcript.list	list of all transcript objects.
	
	# OUTPUT
	# full.gene.model	Object of class "Transcript" specifying the full gene model.
	
	#Define variables
	startPositions = list()
	endPositions = list()
	exonStarts = c()
	exonEnds = c()
	strand = NULL
	chrom = NULL
	
	i = 1
	for (transcript in transcript.list){
		startPositions = append(startPositions, transcript@txStart)
		endPositions = append(endPositions, transcript@txEnd)
		
		if (i == 1){  #Copy most of the infromation from the first transcript
			exonStarts <- transcript@exonStarts
			exonEnds <- transcript@exonEnds
			strand <- transcript@strand
			chrom <- transcript@chrom
			i = 2
		}
		else{  #Add unique exons from the other transcripts
			newStarts <- transcript@exonStarts
			newEnds <- transcript@exonEnds
			newStartIndices <- c(1:length(newStarts))[newStarts %in% exonStarts == FALSE]
			newEndIndices <- c(1:length(newStarts))[newStarts %in% exonStarts == FALSE]
			newExons = unique(c(newStartIndices, newEndIndices))
			exonStarts = c(exonStarts, newStarts[newExons])
			exonEnds = c(exonEnds, newEnds[newExons])	
		}		
	}
	
	txStart = min(unlist(startPositions))
	txEnd = max(unlist(endPositions))
	exonCount = length(exonStarts)
	
	full.gene.model = new("Transcript", txStart = txStart, txEnd = txEnd,
			exonStarts = exonStarts, exonEnds = exonEnds, exonCount = exonCount,
			strand = strand, chrom = chrom, name = "All exons")
	return(full.gene.model)
}

CreatePileups <- function(full.gene.model, bamfiles, total.reads, kernel.width = 1){
	#Function to retrieve reads from BAM files and create pileups.
	
	# INPUT
	# full.gene.model	Transcript object specifying the full gene model.
	# bamfiles		list of BAM files to be analyzed
	# total.reads	vector of total number of reads per each BAM file
	
	# OUTPUT
	# list(reads = reads, counts = counts) - List of two lists, showing the reads
	#										and read counts from each BAM files.
	
	start.pos = full.gene.model@txStart
	end.pos = full.gene.model@txEnd
	chrom = full.gene.model@chrom
	
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
	
	# If given a kernel width, use a kernel to smooth the raw counts values
	if (kernel.width > 1){
		smooth.counts = list()
		for (i in 1:n){
			counts.vector = counts[[i]]
			counts.vector[is.na(counts.vector)] = 0
			k1 = kernel("daniell", 50)
			counts.vector = kernapply(counts.vector, k1)
			smooth.counts[[i]] = counts.vector
		}
		counts = smooth.counts
	}
	
	return(list(reads = reads, counts = counts))
}

DrawWigglePlots <- function(full.gene.model, pileups, total.reads, sample.names, 
							intron.color, exon.colors, cex){
	
	reads = pileups$reads
	counts = pileups$counts
	n <- length(counts)

	exonPositions <- unlist(mapply(seq, full.gene.model@exonStarts, 
					full.gene.model@exonEnds, SIMPLIFY=FALSE));
	
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
		ylab <- sprintf("%s", sample.names[i]);
		plot(x=NA, y=NA, xaxt="n", xlab=NA, ylab=ylab, cex.axis = cex, cex.lab = cex, 
				xlim=c(full.gene.model@txStart, full.gene.model@txEnd),
				ylim=c(0, max.counts), frame.plot = FALSE);
		# If there are counts, plot them as individual line segments
		nonzero <- which(counts[[i]] > 0);
		first.read <- reads[[i]]$pos[1];
		nonzero.exonic <- intersect(nonzero, exonPositions-first.read+1);
		nonzero.intronic <- setdiff(nonzero, exonPositions-first.read+1);
		if (length(nonzero.intronic)) {
			segments(x0=first.read+nonzero.intronic-1, x1=first.read+nonzero.intronic-1,
					y0=0, y1=counts[[i]][nonzero.intronic], col=intron.color)
		}
		if (length(nonzero.exonic)) {
			segments(x0=first.read+nonzero.exonic-1, x1=first.read+nonzero.exonic-1,
					y0=0, y1=counts[[i]][nonzero.exonic], col=exon.colors[i])
		}
	}
}

DrawUniqueExons <- function(transcript.list, color, exon.events){
	#Plot a red rectangle under unique exons of the first transcript or draw the unique
	# exon in red (exon events).
	#
	# INPUT
	# transcript.list	list of all Transcript objects
	
	#Extract information about primary transcript
	primaryExonStarts <- transcript.list[[1]]@exonStarts
	primaryExonEnds <- transcript.list[[1]]@exonEnds
	
	uniqueStartIndices <- c(1:length(primaryExonStarts))
	uniqueEndIndices <- c(1:length(primaryExonEnds))
	
	#Go through all other transcripts
	for (i in 2:length(transcript.list)){
		newStarts <- transcript.list[[i]]@exonStarts
		newEnds <- transcript.list[[i]]@exonEnds
				
		newStartIndices <- c(1:length(primaryExonStarts))[primaryExonStarts %in% newStarts == FALSE]
		newEndIndices <- c(1:length(primaryExonEnds))[primaryExonEnds %in% newEnds == FALSE]
		
		uniqueStartIndices <- intersect(uniqueStartIndices, newStartIndices)
		uniqueEndIndices <- intersect(uniqueEndIndices, newEndIndices)
	}
	
	uniqueIndices = unique(c(uniqueStartIndices,uniqueEndIndices))
	exonStarts <- primaryExonStarts[uniqueIndices]
	exonEnds <- primaryExonEnds[uniqueIndices]
	
	#Mark unique exons with rectangles
	if (length(exonStarts) > 0){
		for (i in c(1:length(exonStarts))){
			if (exon.events){
				rect(xleft=exonStarts[i], xright=exonEnds[i], ybottom=-0.6, ytop=1, col=color, border = NA)
			}
			else{
				rect(xleft=exonStarts[i], xright=exonEnds[i], ybottom=-0.8, ytop=-0.6, col=color, border = NA)
			}
		}		
	}	
}

DrawExonStructure <- function(transcript, color, full.gene.model, cex){
	# Function to draw the exon structure for single transcripts
	
	# INPUT
	# transcript	Object of the class "Transcript"
	# color			color of the exons
	
	# OUTPUT
	# Add exons from one trancript to the gene structure plot
	
	# Draw a line across the transcription unit
	lines(x=c(transcript@txStart+1, transcript@txEnd), y=c(0.2,0.2), col = color);
	
	#Define bottom and top coordinates
	exon.bottom = -0.4
	exon.top = 0.8
	cds.bottom = -0.6
	cds.top = 1
	
	#Iterate over all exons
	for (i in c(1:transcript@exonCount)){
		#Draw exons before CDS
		if (i < transcript@cdsStartExon){
			rect(xleft=transcript@exonStarts[i], xright=transcript@exonEnds[i], 
					ybottom = exon.bottom, ytop = exon.top, col=color, border = NA)
		}
		#Draw exons that intersect with CDS start
		else if (i == transcript@cdsStartExon){
			rect(xleft=transcript@exonStarts[i], xright=transcript@cdsStart+1,
					ybottom = exon.bottom, ytop = exon.top, col=color, border = NA)
			rect(xleft=transcript@cdsStart+1, xright=min(transcript@cdsEnd, transcript@exonEnds[i]), 
					ybottom = cds.bottom, ytop = cds.top, col=color, border = NA)
		}
		#Draw exons that overlap with CDS
		else if (i < transcript@cdsEndExon){
			rect(xleft=transcript@exonStarts[i], xright=transcript@exonEnds[i], 
					ybottom = cds.bottom, ytop = cds.top, col=color, border = NA)
		}
		#Draw exons that intersect with CDS end
		else if (i == transcript@cdsEndExon){
			rect(xleft=max(transcript@cdsStart, transcript@exonStarts[i]), xright=transcript@cdsEnd, 
					ybottom = cds.bottom, ytop = cds.top, col=color, border = NA)
			rect(xleft=transcript@cdsEnd, xright=transcript@exonEnds[i], 
					ybottom= exon.bottom, ytop = exon.top, col=color, border = NA)
		}
		#Draw exons after the CDS
		else{
			rect(xleft=transcript@exonStarts[i], xright=transcript@exonEnds[i], 
					ybottom = exon.bottom, ytop = exon.top, col=color, border = NA);
		}	
	}
	
	#Add name
	text(x = full.gene.model@txStart, y = -0.9, transcript@name, cex = cex, pos = 4,
			col = color)
	
	#Calculate the fraction of chevrons needed
	transcript.length = transcript@txEnd - transcript@txStart
	gene.length = full.gene.model@txEnd - full.gene.model@txStart
	fraction = transcript.length/gene.length
	chevron.intervals <- round(50*fraction)
	
	#Draw chevrons
	x <- seq(transcript@txStart, transcript@txEnd, length.out=chevron.intervals+1)[2:chevron.intervals];
	width <- diff(x)[1] * 0.25;
	height <- 0.2;
	direction <- c("-"=-1, "+"=1)[transcript@strand];
	chevron.starts <- pmin(x, x-width*direction);
	chevron.ends <- pmax(x, x-width*direction);
	# If there are introns, draw chevrons there; if not, draw them in white on the exon
	if (transcript@exonCount > 1) {
		x <- x[!mapply(function(starts, ends) any((starts < transcript@exonEnds) & (ends > transcript@exonStarts)), chevron.starts, chevron.ends)];
		chevron.col <- color;
	} else {
		chevron.col <- "white";
	}
	# Draw the chevrons
	if (length(x) > 0){
		segments(x0=x, x1=x-width*direction, y0=+0.2, y1=-height+0.2, col=chevron.col);
		segments(x0=x, x1=x-width*direction, y0=0+0.2, y1=height+0.2, col=chevron.col);
	}
}

Test = function() {
	#Test that the plot works
	ids = c('NM_044472','NM_001039802','NM_001791')
	#ids = c('chr2:204299601:204299957:+@chr2:204302616:204302740:+@chr2:204307752:204310898:+.A',
	#		'chr2:204299601:204299957:+@chr2:204302616:204302740:+@chr2:204307752:204310898:+.B',
	#		'NM_006139')
	bamfiles = c("/Users/alasoo/projects/Tripathi/alignments/ThP/ThP.shrimp.final.sorted.bam", 
			"/Users/alasoo/projects/Tripathi/alignments/Th0/Th0.shrimp.final.sorted.bam")
	names(bamfiles) = c("ThP", "Th0")
	bedfile1 = ReadBedFile("/Users/alasoo/projects/Tripathi/annotations/bed/SE.hg18.bed")
	bedfile = ReadBedFile("/Users/alasoo/projects/Tripathi/annotations/refGene/refGene.hg18.270711.bed")
	bedfile = rbind(bedfile, bedfile1)
	WigglePlotR(ids, bamfiles, bedfile, cex = 1, kernel.width = 50)
}

Test()


