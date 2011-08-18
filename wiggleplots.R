#Load functions from functions file
source('/Users/alasoo/workspace/WigglePlotR/functions.R')

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
	
	# Obtain the record from the BED file corresponding to the primary transcript
	primaryBed <- subset(bed, name == primaryTranscript)
	if (nrow(primaryBed) == 0){
		print(paste("ERROR: Primary transcript", primaryTranscript, "has no record in the BED file. Skipping.", sep = " "))
		return()
	}
	else if (nrow(primaryBed) > 1){
		print(paste("ERROR: Primary transcript", primaryTranscript, "has more than one record in BED file. Selecting first.", sep = " "))
		primaryBed = primaryBed[1,]	
	}
	chrom <- primaryBed$chrom
	strand <- primaryBed$strand
	
	#Check if other transcripts acutally exist in the BED file:
	for (transcript_id in otherTranscripts){
		bed.record = subset(bed, name == transcript_id)
		if (nrow(bed.record) == 0){
			otherTranscripts = otherTranscripts[-which(otherTranscripts == transcript_id)] #Remove if not in bed file
		}			
	}
	
	### EXTRACT PROPTERTIES OF THE LONGEST TRANSCRIPT ###
	fullGeneModel <- CreateFullGeneModel(c(primaryTranscript,otherTranscripts), bed)
	start.pos <- fullGeneModel$start.pos
	end.pos <- fullGeneModel$end.pos
	exonCount <- fullGeneModel$exonCount
	exonStarts <- fullGeneModel$exonStarts
	exonEnds <- fullGeneModel$exonEnds
	
	### Retrieve alignments and create pileups ###
	cat("Retrieving alignments and creating pileups.\n");
	pileups = CreatePileups(start.pos, end.pos, chrom, bamfiles, total.reads)
	
	# Layout has n+2 rows: one for each sample, plus 1 title row and 1 annotation row
	layout(matrix(1:(n+m+1),n+m+1,1), heights = c(rep(2,n),rep(1,m)));
	par(mar=c(2,50,2,1), bg="transparent");
	# Write the title bar
	#frame();
	#with(primaryBed, {
	#			text(x=0.5, y=0.5,
	#					#sprintf("%s\n(%s:%s-%s)", primaryTranscript, chrom, format(start.pos, big.mark=","), format(end.pos, big.mark=",")), cex=2);
	#					sprintf("%s", primaryTranscript, format(end.pos, big.mark=",")), cex=1);
	#		});
	
	##### DRAW WIGGLE PLOTS #####
	DrawWigglePlots(start.pos, end.pos, exonStarts, exonEnds, pileups, total.reads, sample.names, bg.colors, intron.color, exon.colors)
	
	##### DRAW EXON STRUCTURE OF TRANSCRIPT #####
	InitializeExonStructurePlot(start.pos, end.pos)
	
	#Draw exon structures of the primary transcript
	DrawExonStructure(primaryBed, "black", start.pos)
	#Draw unique exons of the primary transcript
	DrawUniqueExons(primaryTranscript, otherTranscripts, bed)
	
	# Draw chevrons to show direction of transcription
	#DrawChevrons(start.pos, end.pos, exonStarts, exonEnds, exonCount, strand)
	
	#Draw exon structures of other transcripts
	for (id in otherTranscripts){
		InitializeExonStructurePlot(start.pos, end.pos)
		secondaryBed = subset(bed, name == id)
		DrawExonStructure(secondaryBed, "black", start.pos)
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
wiggleplots(primary_tx, other_tx, bamfiles, bedfile)


