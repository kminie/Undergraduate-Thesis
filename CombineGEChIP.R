## CombineGEChIP.R
## Version 16, March 2nd, 2015
## Created By: Tyler Kolisnik and Mark Bieda.

## Description: This is a program to analyze gene expression microarray data and compare it to ChIP-seq data from two samples and generate meaningful
## graphical displays of their correlations.
## This program is used to determine the pharmacogenomic applications of transcription-factor altering drugs (i.e. histone deacetylase inhibitors) 
## or to compare gene expression levels between cell types. 
## It allows the user to determine what genes the transcription factor is present at uniquely in one input, and look at 
## the gene expression levels of those specific genes, allowing for example, the determination of the effect of a drug on these genes.


## Dependencies: 
library(simpleaffy) # Required for data normalization.
library(marray) # Required for Heatmap Generation.
library(gplots) # Required for Heatmap Generation.
library(limma) # Required for DEG (Differentially Expressed Gene) Selection.
library(hgu133plus2.db) # Must be changed if using a different array other than Affymetrix hgu133plus2.0 array. Other points in code must be changed as well.
library(RColorBrewer) # Required for Heatmap Generation.
library(org.Hs.eg.db) # Organism Database ID (used by GOstats), must be changed if using an organism other than human (default).
library(GOstats) # Required for Gene Ontology Analysis.
library(gage) # Required for Pathway Analysis.
library(pathview) # Required for Pathway Analysis.
library(ggplot2)
library(plotrix)



##################################################################### Begin Parameter Block #####################################################################

chipInputData <- "/home/tylerk/AFX1/ctTest/mcf7-chip/nohomerresults/encs5r000dwc_5000up_2000down_allData.txt" #mcf7-cell DATA
chipInputData2 <- "/home/tylerk/AFX1/ctTest/osteo-chip/results/osteo_5000up_2000down_allData.txt" #osteo-cell DATA
chipInputData_lostPeaks <- "/home/tylerk/AFX1/ctTest/mcf7-osteo/lostinosteo_5000up_2000down_allData.txt"  #Peaks present in one but not other.

#For gene expression data, they must both be from the same array type.
#inputDirectory must contain .CEL files corresponding to the same celltype/experiment type as in chipInputData
#inputDirectory 2 must contain .CEL files corresponding to the same celltype/experiment type as in chipInputData2
inputDirectory <- "/home/tylerk/AFX1/ctTest/mcf7-expr" # Path of directory where input is located, must contain .CEL files and a covdesc file, no trailing /.
inputDirectory2 <- "/home/tylerk/AFX1/ctTest/osteo-expr" 
outputDirectory <- "/home/tylerk/AFX1/ctTest/mcf7-osteo" # Path of directory where output will be located, no trailing /.
baseCondition <- "MCF7"
comparator <- "Osteo"
isItAnExonArray <- "yes"
exonArrayPackageName <- "pd.huex.1.0.st.v2" 
dbpackagename <- "huex10sttranscriptcluster.db"
histOrTFName <- "H3K27ac" #The name of the modification or transcription factor you are analyzing, so the graph can be labeled appropriately
graphOutputDirectory <- "/home/tylerk/AFX1/ctTest/mcf7-osteo"



##################################################################### End Parameter Block #######################################################################

##################################################################### Begin Coding Block #####################################################################
######################### Begin Data Entry and Manipulation Block ############################
#Load libraries specific to array type
library(dbpackagename, character.only=TRUE)
library(exonArrayPackageName, character.only=TRUE)
assign("dbpackagename", get(dbpackagename))

options(max.print=100)
setwd(inputDirectory)

################## PROCESS CEL FILES IN FIRST INPUT DIRECTORY ##########################
if(isItAnExonArray == "yes"){
############If an exon array ##########
print("Processing as Exon Array")
rawData <- read.celfiles(list.celfiles(),pkgname=exonArrayPackageName)
dataSet <- rma(rawData, target="core")
exprset <- exprs(dataSet)
exprMeans <- rowMeans(exprset) 
exprMatrix <- as.matrix(exprMeans)
probeIDs <- rownames(exprMatrix)
linkprobes <- select(dbpackagename, keys= probeIDs, columns = "ENTREZID", keytype = "PROBEID")
exprsetlinkedtogenes <-merge(exprMatrix, linkprobes, by.x=0, by.y="PROBEID")
############End  an exon array ##########
}else{
#############If an affy array ###########
print("Processing as Affy Array")
setwd(inputDirectory)
celfiles <- list.celfiles()
dataSet <- justRMA(filenames=celfiles)
exprset <- exprs(dataSet)
exprMeans <- rowMeans(exprset) 
exprMatrix <- as.matrix(exprMeans)
probeIDs <- rownames(exprMatrix)
linkprobes <- select(dbpackagename, keys= probeIDs, columns = "ENTREZID", keytype = "PROBEID")
exprsetlinkedtogenes <-merge(exprMatrix, linkprobes, by.x=0, by.y="PROBEID")
#############End if an affy array ###########
}
entrezID_column_position <- grep("ENTREZID", names(exprsetlinkedtogenes)) #This gets the position of the column with EntrezIDs, It needs to be done this way as $ is not valid for atomic vectors.
exprsetlinkedtogenes <- exprsetlinkedtogenes[,c(entrezID_column_position, (1:ncol(exprsetlinkedtogenes))[-entrezID_column_position])]
exprsetlinkedtogenes[2]<-NULL
expmatrix<-data.matrix(exprsetlinkedtogenes)
rownames(expmatrix)<-expmatrix[,1]
microarrayData_1 <- expmatrix

################## END PROCESS CEL FILES IN FIRST INPUT DIRECTORY ##########################

################## PROCESS CEL FILES IN SECOND INPUT DIRECTORY ##########################
setwd(inputDirectory2)
if(isItAnExonArray == "yes"){
  ############If an exon array ##########
  print("Processing as Exon Array")
  rawData2 <- read.celfiles(list.celfiles(),pkgname=exonArrayPackageName)
  dataSet2 <- rma(rawData2, target="core")
  exprset2 <- exprs(dataSet2)
  exprMeans2 <- rowMeans(exprset2) 
  exprMatrix2 <- as.matrix(exprMeans2)
  probeIDs2 <- rownames(exprMatrix2)
  linkprobes2 <- select(dbpackagename, keys= probeIDs2, columns = "ENTREZID", keytype = "PROBEID")
  exprsetlinkedtogenes2 <-merge(exprMatrix2, linkprobes2, by.x=0, by.y="PROBEID")
  ############If an exon array ##########
}else{
  #############If an affy array ###########
  print("Processing as Affy Array")
  setwd(inputDirectory2)
  celfiles2 <- list.celfiles()
  dataSet2 <- justRMA(filenames=celfiles2)
  exprset2 <- exprs(dataSet2)
  exprMeans2 <- rowMeans(exprset2) 
  exprMatrix2 <- as.matrix(exprMeans2)
  probeIDs2 <- rownames(exprMatrix2)
  linkprobes2 <- select(dbpackagename, keys= probeIDs2, columns = "ENTREZID", keytype = "PROBEID")
  exprsetlinkedtogenes2 <-merge(exprMatrix2, linkprobes2, by.x=0, by.y="PROBEID")
  #############End if an affy array ###########
}
entrezID_column_position2 <- grep("ENTREZID", names(exprsetlinkedtogenes2)) #This gets the position of the column with EntrezIDs, It needs to be done this way as $ is not valid for atomic vectors.
exprsetlinkedtogenes2 <- exprsetlinkedtogenes2[,c(entrezID_column_position2, (1:ncol(exprsetlinkedtogenes2))[-entrezID_column_position2])]
exprsetlinkedtogenes2[2]<-NULL
expmatrix2<-data.matrix(exprsetlinkedtogenes2)
rownames(expmatrix2)<-expmatrix2[,1]
microarrayData_2 <- expmatrix2
################## END PROCESS CEL FILES IN SECOND INPUT DIRECTORY ##########################


##Get gene list from First ChIP-data
firstChIPData<- read.delim(chipInputData)
ChIPgeneList<- firstChIPData$entrezID

##Get gene list from second ChIP-seq Data
secondChIPData<- read.delim(chipInputData2)
ChIPgeneList2<- secondChIPData$entrezID

##get list of 'lost' genes from ChIP-seq Data
fullChIPData_lostPeaks<- read.delim(chipInputData_lostPeaks)
ChIPgeneList_lostPeaks<- fullChIPData_lostPeaks$entrezID

## This outputs a list of EntrezIDs where the histone or transcription factor signal was Present in the first set of ChIP data
IDlistP <- microarrayData_1[microarrayData_1[,1] %in% ChIPgeneList,] #Subset the microarray data based on the chip entrezids
IDlistP <- IDlistP[!is.na(IDlistP[,1]),] # delete all NA rows
IDlistP <- subset(IDlistP, !duplicated(IDlistP[,1])) # removes all duplicates

## This outputs a list of EntrezIDs where the histone or transcription factor signal was Present in the second set of ChIP data
IDlistP_2 <- microarrayData_2[microarrayData_2[,1] %in% ChIPgeneList2,] #Subset the microarray data based on the chip entrezids
IDlistP_2 <- IDlistP_2[!is.na(IDlistP_2[,1]),] # delete all NA rows
IDlistP_2 <- subset(IDlistP_2, !duplicated(IDlistP_2[,1])) # removes all duplicates

## This outputs a list of EntrezIDs where the histone or transcription factor signal was Not Present in the first set of ChIP data
IDlistN <- microarrayData_1[!microarrayData_1[,1] %in% ChIPgeneList,]
IDlistN <- IDlistN[!is.na(IDlistN[,1]),] # deletes all NA rows
IDlistN <- subset(IDlistN, !duplicated(IDlistN[,1])) # removes all duplicates

## This outputs a list of EntrezIDs where the histone or transcription factor signal was Not Present in the first set of ChIP data
IDlistN_2 <- microarrayData_2[!microarrayData_2[,1] %in% ChIPgeneList2,]
IDlistN_2 <- IDlistN_2[!is.na(IDlistN_2[,1]),] # deletes all NA rows
IDlistN_2 <- subset(IDlistN_2, !duplicated(IDlistN_2[,1])) # removes all duplicates

## This outputs a list of EntrezIDs where the histone or transcription factor signal was 'lost' in the first set of data
IDlistL <- microarrayData_1[microarrayData_1[,1] %in% ChIPgeneList_lostPeaks,] #Subset the microarray data based on the chip entrezids
IDlistL <- IDlistL[!is.na(IDlistL[,1]),] # delete all NA rows
IDlistL <- subset(IDlistL, !duplicated(IDlistL[,1])) # removes all duplicates

## This outputs a list of EntrezIDs where the histone or transcription factor signal was 'lost' in the second set of data
IDlistL_2 <- microarrayData_2[microarrayData_2[,1] %in% ChIPgeneList_lostPeaks,] #Subset the microarray data based on the chip entrezids
IDlistL_2 <- IDlistL_2[!is.na(IDlistL_2[,1]),] # delete all NA rows
IDlistL_2 <- subset(IDlistL_2, !duplicated(IDlistL_2[,1])) # removes all duplicates

## T Test:
ttest <- t.test(IDlistL[,2],IDlistL_2[,2],paired=TRUE)
pval <- ttest$p.value
write.list(ttest,file="ttestdata.txt")
##First Graph: a bar graph Gene Expression Levels of The "lost" genes in both sets of data:
setwd(graphOutputDirectory)
graphData1 <- IDlistL[,2]
graphData2 <- IDlistL_2[,2]
names(graphData1) <- NULL
names(graphData2) <- NULL
## Create the data for the legend
legendName1 <- paste(baseCondition,"-",comparator,"+ genes in",baseCondition)
legendName2 <- paste(baseCondition,"-",comparator,"+ genes in",comparator)
legend.data <- c(legendName1,legendName2)


sd1 <- sd(graphData1)
sd2 <- sd(graphData2)
sdList <- matrix(c(sd1,sd2))


mean1 <- mean(graphData1)
mean2 <- mean(graphData2)
meanList <- c(mean1, mean2)
names(meanList) <- legend.data
barplotData <- as.data.frame(meanList)

##Calculate p value to display
if(pval<0.01){ 
  legendpval <- "< 0.01"
} else if(pval<0.05){
  legendpval <- "< 0.05"
} else if(pval==0.05){
  legendpval <- "= 0.05"
} else if(pval>0.05){
  legendpval <- "> 0.05"
}
#Bar Plot is Plotted Here:

pdf(paste(histOrTFName,baseCondition,"_",comparator,"barplot.pdf", sep=""),height=9.5)
limits <- aes(ymax = meanList + sdList, ymin=meanList - sdList)
dodge <- position_dodge(width=0.5)
PNbarplot <- qplot(rownames(barplotData),barplotData[,1],data=barplotData, geom="bar",stat="identity", ylab= "Avg Gene Expression Levels", xlab= paste("Differential H3K27ac Presence"), ylim=c(0,10), fill=barplotData[,1], guides(fill=FALSE))
PNbarplot + theme(legend.position="none") #removes legend
PNbarplot + geom_path(x=c(1.25,1.75),y=c(7.9,7.9)) + 
  annotate("text",x=2,y=10,label=paste("** = p value", legendpval)) +
  annotate("text", x=1.5, y=8, label="**") +
  geom_errorbar(limits, position=dodge, width=0.1) +
  theme(legend.position="none")
dev.off()

########################################### END FIRST GRAPH ##################################################################

############ Create Histogram #################


pdf(paste(histOrTFName,baseCondition,"_",comparator,"_histogram.pdf", sep=""),height=6,width=4)
histP<-hist(rbind(graphData1),plot=FALSE,breaks="Sturges")$counts
histN<-hist(rbind(graphData2),plot=FALSE,breaks="Sturges")$counts


histPadj <- hist(rbind(graphData1),plot=FALSE,breaks=seq(0.5,15,1))$counts/sum(hist(rbind(graphData1),plot=FALSE,breaks=1:15)$counts)
histNadj <- hist(rbind(graphData2),plot=FALSE,breaks=seq(0.5,15,1))$counts/sum(hist(rbind(graphData2),plot=FALSE,breaks=1:15)$counts)

barp(rbind(histPadj,histNadj),staxx=F,col=c("red","gray"), xlab=("Gene Expression Levels"), ylab=("Relative Frequency")) 
legend('topright', legend=legend.data, col=c("red", "gray"), pch=15, bty="n", cex=0.5)
dev.off()
############ End Create Histogram #################

##############  Create Bar Plot   ###################

sdP <- sd(graphDataP)
sdN <- sd(graphDataN)
sdL <- sd(graphDataL)
sdOtherL <- sd(cat[,1])
sdList <- matrix(c(sdL,sdOtherL))


meanP <- mean(graphDataP)
meanN <- mean(graphDataN)
meanL <- mean(graphDataL)
meanOtherL <- mean(cat[,1])
meanList <- c(meanL, meanOtherL)
names(meanList) <- legend.data
barplotData <- as.data.frame(meanList)

##Calculate p value to display
if(ttestdata$p.value<0.01){ 
  legendpval <- "< 0.01"
} else if(ttestdata$p.value<0.05){
  legendpval <- "< 0.05"
} else if(ttestdata$p.value==0.05){
  legendpval <- "= 0.05"
} else if(ttestdata$p.value>0.05){
  legendpval <- "> 0.05"
}
#Bar Plot is Plotted Here:

pdf(paste(histOrTFName,"tylers_barplot.pdf", sep=""),height=9.5)
limits <- aes(ymax = meanList + sdList, ymin=meanList - sdList)
dodge <- position_dodge(width=0.5)
PNbarplot <- qplot(rownames(barplotData),barplotData[,1],data=barplotData, geom="bar",stat="identity", ylab= "Avg Gene Expression Levels", xlab= paste("Differential H3K27ac Presence"), ylim=c(0,10), fill=barplotData[,1], guides(fill=FALSE))
PNbarplot + theme(legend.position="none") #removes legend
PNbarplot + geom_path(x=c(1.25,1.75),y=c(7.9,7.9)) + 
annotate("text",x=2,y=10,label=paste("** = p value", legendpval)) +
annotate("text", x=1.5, y=8, label="**") +
geom_errorbar(limits, position=dodge, width=0.1) +
theme(legend.position="none")
dev.off()

##############  End Create Bar Plot   ###################

################################### End Graphing Block ##################################

##################################################################### End Coding Block #####################################################################