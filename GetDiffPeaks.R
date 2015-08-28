##GetDiffPeaks.R
##Created by Tyler Kolisnik

#This is a program to get differential peaks from two ChIP-seq datasets.
#Input the read file and peak file from each of two datasets. The number 1 corresponds to the first dataset, 
# and the number 2 corresponds to the second dataset. These must match up.
#BP shift values can be found in the MACS peak file (*_peaks.xls) after "# d =".
#This progam is usually followed by running the annotatePeaks.R Rscript for each of the output files:
#"MAnorm_uniquePeaks1_cut.txt"
#"MAnorm_uniquePeaks2_cut.txt"
# in order to annotate the peaks to genes. 


#This program must be run in the command line, as it requires access to the path, also 
#Bedtools must be installed.
#The files MAnorm.sh, MAnorm.R, 
################# Begin Parameter Block ############
workingDirectory<- "/home/tylerk/AFX1_Ready/ChIP-GE"
chipReadData1<- "/home/tylerk/AFX1_Ready/ChIP-GE/mcf7Reads.bed"
chipPeakData1<- "/home/tylerk/AFX1_Ready/ChIP-GE/mcf7Peaks.bed"
chipReadData2<- "/home/tylerk/AFX1_Ready/ChIP-GE/osteoReads.bed"
chipPeakData2<- "/home/tylerk/AFX1_Ready/ChIP-GE/osteoPeaks.bed"
bpShift1 <- "150"
bpShift2 <- "150"
################## End Parameter Block #############


system(paste("./MAnorm_Input_Column_Adjust.sh", chipPeakData1,chipPeakData2, chipReadData1,chipReadData2, sep=" "), wait=TRUE)
newchipReadData1<- paste(workingDirectory,"read1.bed",sep="/")
newchipPeakData1<- paste(workingDirectory,"peak1.bed",sep="/")
newchipReadData2<- paste(workingDirectory,"read2.bed",sep="/")
newchipPeakData2<- paste(workingDirectory,"peak2.bed",sep="/")
system(paste("./MAnorm_Edited.sh", chipPeakData1,chipPeakData2, chipReadData1,chipReadData2,bpShift1, bpShift2, sep=" "), wait=TRUE)

system(paste("./MAnorm_output_formatting.sh ",workingDirectory,"/MAnorm_result.xls", sep=""), wait=TRUE)