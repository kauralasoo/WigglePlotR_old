#Load functions from functions file
source('/Users/alasoo/workspace/WigglePlotR/functions.R')

#Define a class to store transcript infromation
setClass("Transcript", representation(
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

ReadTranscriptsFromBed <- function(transcript.ids, bed.file){
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

CreatePileups <- function(full.gene.model, bamfiles, total.reads){
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
	
	return(list(reads = reads, counts = counts))
}

DrawWigglePlots <- function(full.gene.model, pileups, total.reads, sample.names, bg.colors, intron.color, exon.colors){
	
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
		#ylab <- sprintf("%s%s", sample.names[i], ifelse(!is.null(total.reads), "\nRPKM", "\nreads"));
		ylab <- sprintf("%s", sample.names[i]);
		plot(x=NA, y=NA, xaxt="n", xlab=NA, ylab=ylab, cex.axis=2, cex.lab = 2, 
				xlim=c(full.gene.model@txStart, full.gene.model@txEnd), ylim=c(0, max.counts));
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

DrawUniqueExons <- function(transcript.list){
	#Plot a red rectangle under unique exons of the first transcript.
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
			rect(xleft=exonStarts[i], xright=exonEnds[i], ybottom=-1, ytop=-0.77, col="red", border = NA);
		}		
	}	
}

wiggleplots <- function(primaryTranscript, otherTranscripts, bamfiles, bedfile, total.reads=NULL,
		exon.colors="black", intron.color="lightgray", bg.colors="transparent") {
	
	# Function for creating "wiggle plots" from BAM alignment files and BED annotation files
	# Adam Gower, 2010
	#
	# INPUT
	# ids             A character vector of IDs that correspond to features in the BED file, e.g., "NM_001234"
	# bamfiles        A character vector of BAM file names, one for each alignment that will be used to generate a wiggle plot
	# bedfile         The name of a BED file that will be used to determine the chromosomal coordinates of the features in 'ids'
	# total.reads     An integer vector of the total number of reads for each sample; if supplied, used to scale the y-axes of plots
	# exon.colors     A character vector of color names with which to draw base positions that overlap with exonic regions
	# intron.color    A color name with which to draw base positions that overlap with intronic regions
	# bg.colors       A character vector of color names with which to shade the background of each wiggle plot
	#
	# OUTPUT
	# Draws a single page for each feature in 'ids', containing a single wiggle plot for each alignment,
	#     as well as a representation of the structure of the gene at the bottom of the page
	
	# Need the Rsamtools package to read BAM files
	require(Rsamtools);
	
	n <- length(bamfiles);
	m <- length(otherTranscripts)
	
	if (!is.null(names(bamfiles))) {
		# Get the sample names from the names of the filename vector, if available
		sample.names <- names(bamfiles);
	} else {
		# Otherwise, extract the names of the files after the path and before the .bam extension
		sample.names <- sub("(^.*/)*(.*)\\.bam$","\\2", bamfiles);
	}
	
	#if (!is.null(names(ids))) {
	#	# Get the ID names from the names of the ID vector, if available
	#	id.names <- names(ids);
	#} else {
	#	# Otherwise, use the ID itself
	#	id.names <- ids;
	#}
	
	# Parse the BED file
	#cat(sprintf("Parsing BED file %s.\n", bedfile));
	#bed <- read.bedfile(bedfile);
	bed = bedfile
	# If this is a 'short' BED file with only 6 columns, use the data to add extra columns to the full BED file specification 
	if (ncol(bed) == 6) {
		bed <- cbind(bed, data.frame(thickStart=bed$chromStart, thickEnd=bed$chromEnd, itemRgb=0,
						blockCount=1, blockSizes=sprintf("%d,", bed$chromEnd-bed$chromStart), blockStarts="0,",
						stringsAsFactors=FALSE));
	}
	# Change some of the column names
	colnames(bed) <- sub("chromStart","txStart", colnames(bed));
	colnames(bed) <- sub("chromEnd","txEnd", colnames(bed));
	colnames(bed) <- gsub("thick","cds", colnames(bed));
	colnames(bed) <- gsub("block", "exon", colnames(bed));
	# Keep only the rows with the standard chromosome IDs (no "random" or haplotype-specific chromosomes)
	bed <- subset(bed, chrom %in% sprintf("chr%s", c(as.character(1:22),"X","Y","M")));
	# Remove duplicate rows
	bed <- unique(bed);
	
	if (length(exon.colors) == 1) exon.colors <- rep(exon.colors, n);
	if (length(bg.colors) == 1) bg.colors <- rep(bg.colors, n);
	
	#Read all transcripts from the BED file
	transcript.list = ReadTranscriptsFromBed(c(primaryTranscript, otherTranscripts), bed)	
	#Create full gene model (contains all exons)
	full.gene.model = CreateFullGeneModel(transcript.list)
	
	### Retrieve alignments and create pileups ###
	cat("Retrieving alignments and creating pileups.\n");
	pileups = CreatePileups(full.gene.model, bamfiles, total.reads)
	
	# Create the plot layout
	layout(matrix(1:(n+m+1),n+m+1,1), heights = c(rep(2,n),rep(1,m)));
	par(mar=c(2,50,2,1), bg="transparent");
	
	##### DRAW WIGGLE PLOTS #####
	DrawWigglePlots(full.gene.model, pileups, total.reads, sample.names, bg.colors, intron.color, exon.colors)
	
	##### ADD EXON STRUCTURES #####
	i <- 1
	for (transcript in transcript.list){
		#Initialize plot and draw exon structure
		par(mar=c(0.5,5,0.5,1))
		plot(x=NULL, y=NULL, yaxt="n", xaxt="n", xlab=NA, ylab=NA, 
			xlim=c(full.gene.model@txStart, full.gene.model@txEnd), ylim=c(-1,1))
		DrawExonStructure(transcript, "black", full.gene.model@txStart)
		
		#Mark unique exons on the first transcript
		if (i == 1){
			DrawUniqueExons(transcript.list)
			i <- 0
		}
	}
}   

primary_tx = 'NM_044472'
other_tx = c('NM_001039802','NM_001791')
#primary_tx = 'chr4:152243458:152243681:+@chr4:152244779:152244880:+@chr4:152245074:152245254:+.A'
#other_tx = c('chr4:152243458:152243681:+@chr4:152244779:152244880:+@chr4:152245074:152245254:+.B')
bamfiles = c("/Users/alasoo/projects/Tripathi/alignments/ThP/ThP.shrimp.final.sorted.bam", 
		"/Users/alasoo/projects/Tripathi/alignments/Th0/Th0.shrimp.final.sorted.bam")
names(bamfiles) = c("ThP", "Th0")
bedfile = read.bedfile("/Users/alasoo/projects/Tripathi/annotations/refGene/refGene.hg18.270711.bed")
record = wiggleplots(primary_tx, other_tx, bamfiles, bedfile)


