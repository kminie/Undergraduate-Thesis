## annotateChIPpeaks.R
## Version 10.2, Wednesday, December 10th, 2014
## Created By: Tyler Kolisnik, and Mark Bieda.

## DESCRIPTION:
# annotatePeaks takes in a set of peak information in .bed format
# TSS information is obtained from the specified database
# This TSS information is used to determine which genes each peak is associated with
# See the bioconductor website for lists of available TxDb and org.??.eg.db databases

# Run Speed: This program usually takes less than 30 seconds to execute. 

## DIRECTIONS:
# Change the appropriate parameters.
# Run the program.
# Output is two .txt files, one that contains all the Data with the file name format: runName_upstreamRangeup_downstreamRangedown_allData.txt
# and one with only the gene symbols with the file name format: runName_upstreamRangeup_downstreamRangedown_geneSymbols.txt

## EXAMPLE INPUT:
# inputFile <- "/home/tylerk/massInputData/GATA1_peaks.bed" 
# outputDirectory <- "/home/tylerk/massTestResults/peakAnnotation" 
# runName <- "massTest1" 
# transcriptionDB <- "TxDb.Hsapiens.UCSC.hg19.knownGene"
# annotationDB <- "org.Hs.eg.db"
# upstreamRange <- 5000 
# downstreamRange <- 2000 

## See these examples when using different species or genome assembly:
# All annotation packages can be found at the bioconductor website: http://www.bioconductor.org/packages/3.0/data/annotation/

## Example Use of Different TxDB. For Mmusculus assembly mm10:
# Visit Bioconductor website link above, find correct TxDB package in list, click link and 
# Install by running in R:
# source("http://bioconductor.org/biocLite.R")
# biocLite("TxDb.Mmusculus.UCSC.mm10.ensGene")
# Then be sure to change transcriptDB parameter to "TxDb.Mmusculus.UCSC.mm10.ensGene", example:
# transcriptDB <- "TxDb.Mmusculus.UCSC.mm10.ensGene" 

## You will also have to install the correct organism annotation Db. Example is for Mouse (mus musculus).
# Install by running in R:
# source("http://bioconductor.org/biocLite.R")
# biocLite("org.Mm.eg.db")
# Then be sure to change the annotationDB parameter to "org.Mm.eg.db", example:
# annotationDB <- "org.Mm.eg.db"

##########################   SET VARIABLES HERE   ####################################

## The input file and output directory must be complete file paths with no trailing /
inputFile <- "/home/tylerk/AFX1/ctTest/gluco-chip/diffPeaks/MAnorm_uniquePeaks2_cut.txt"
outputDirectory <- "/home/tylerk/AFX1/ctTest/gluco-chip/diffPeaks" # The directory where output files will be located
runName <- "uniqueToUntreated" # runName will be used as a prefix in the output file names
transcriptionDB <- "TxDb.Mmusculus.UCSC.mm9.knownGene" # Specify the correct TxDb package for the species being studied.
annotationDB <- "org.Mm.eg.db" # Specify the correct orgDB annotation database package for the species being studied.
upstreamRange <- 5000 # Numeric maximum distance upstream of TSS to look for a peak
downstreamRange <- 2000 # Numeric maximum distance downstream of TSS to look for a peak

######################################################################################

library(transcriptionDB, character.only=TRUE)
assign("txdb", get(transcriptionDB))
library(annotationDB, character.only=TRUE)
assign("annotdb", get(annotationDB))
library(GenomicRanges)
library(rtracklayer)

peaks <- import.bed(con=inputFile, asRangedData=F)
start(peaks) <- start(peaks)-1

#function for extending the TSS to the ranges specified by upstreamRange and downstreamRange
extendTSS <- function(x){
	if(x[5]=="-"){
		x[length(x)+1] <- as.integer(x[7])
		x[length(x)+1] <- as.integer(x[7])-downstreamRange
		x[length(x)+1] <- as.integer(x[7])+upstreamRange
	}
	else{
		x[length(x)+1] <- as.integer(x[6])
		x[length(x)+1] <- as.integer(x[6])-upstreamRange
		x[length(x)+1] <- as.integer(x[6])+downstreamRange
	}
	return(x)
}

## get information from the TxDb and extend the TSSs
keys <- keys(txdb, keytype="TXID")
columns=c("TXSTART", "TXEND","TXSTRAND", "TXCHROM", "GENEID", "TXNAME")
TSSdata <- select(txdb, keys=keys, columns=columns, keytype="TXID")
TSSdata <- data.frame(t(apply(TSSdata, 1, extendTSS)))
TSSdata$V8 <- as.numeric(levels(TSSdata$V8))[TSSdata$V8]
TSSdata$V9 <- as.numeric(levels(TSSdata$V9))[TSSdata$V9]
TSSdata$V10 <- as.numeric(levels(TSSdata$V10))[TSSdata$V10]

## create Granges object for the TSS ranges
TSSranges <- GRanges(seqnames = Rle(TSSdata$TXCHROM),
ranges = IRanges(as.numeric(TSSdata$V9), as.numeric(TSSdata$V10)),
strand = Rle(TSSdata$TXSTRAND), entrezID=TSSdata$GENEID, UCSCID=TSSdata$TXNAME, TSS=TSSdata$V8)

## find overlaps between peak data and TSS ranges
peakData <- data.frame(as.character(seqnames(peaks)), ranges(peaks), peaks@elementMetadata) 
TSSdata <- data.frame(TSSranges@elementMetadata)
overlap <- findOverlaps(peaks,TSSranges)

## extract overlapping peaks and associated data
outputData <- data.frame(as.character(seqnames(peaks))[queryHits(overlap)], ranges(peaks)[queryHits(overlap)], 
peaks@elementMetadata[queryHits(overlap),],ranges(TSSranges)[subjectHits(overlap)], TSSranges@elementMetadata[subjectHits(overlap),], 
strand(TSSranges)[subjectHits(overlap)])
colnames(outputData)[1] <- "chromosome"
colnames(outputData)[13] <- "strand"
colnames(outputData)[which(colnames(outputData)=="start")] <- "peak start"
colnames(outputData)[which(colnames(outputData)=="end")] <- "peak end"
colnames(outputData)[which(colnames(outputData)=="name")] <- "peak name"
colnames(outputData)[which(colnames(outputData)=="UCSCID")] <- "UCSC_ID"
colnames(outputData)[which(colnames(outputData)=="TXSTART")] <- "TSS"
colnames(outputData)[which(colnames(outputData)=="GENEID")] <- "entrez ID"
colnames(outputData)[which(colnames(outputData)=="start.1")] <- "TSS_region_start"
colnames(outputData)[which(colnames(outputData)=="end.1")] <- "TSS_region_end"

## remove extra data fields
drops <- c("width", "TXCHROM", "width.1")
outputData <- outputData[,!(names(outputData) %in% drops)]

## get additional data from the annotationDB
keys <- keys(annotdb, keytype="ENTREZID")
columns=c("SYMBOL","GENENAME")
annotData <- select(annotdb, keys=keys, columns=columns, keytype="ENTREZID")

## add new annotation data to the previous data list
colnames(annotData)[which(colnames(annotData)=="ENTREZID")] <- "entrezID"
colnames(annotData)[which(colnames(annotData)=="SYMBOL")] <- "Gene_Symbol"
colnames(annotData)[which(colnames(annotData)=="GENENAME")] <- "Full Gene Name"
outputData <- merge(outputData, annotData, by="entrezID")
outputData <- outputData[c(2:6, 10, 7:8, 11, 9, 1, 12:13)]

## write data tables to file
symbolList <- unique(outputData$Gene_Symbol)
outputFile <- paste(outputDirectory, "/", runName, "_", upstreamRange, "up_", downstreamRange, "down_allData.txt", sep="")
symbolFile <- paste(outputDirectory, "/", runName, "_", upstreamRange, "up_", downstreamRange, "down_geneSymbols.txt", sep="")
write.table(outputData, file=outputFile,append=FALSE, quote=FALSE, sep="\t", eol = "\n", row.name=FALSE, col.name=TRUE)
write.table(symbolList, file=symbolFile,append=FALSE, quote=FALSE, sep="\t", eol = "\n", row.name=FALSE, col.name=FALSE)
