# Producing a clustered heatmap based on the presence of the secreted genes
# Initially made by Moez Valliani on May 24, 2017 (for CAZymes)

# Please run all the code present in here.
# For changing what data you want displayed, please search for: #! CHANGE DATA TO WHAT IS NEEDED

# TROUBLESHOOTING:
#   If the kclust() says the data is not a numeric dataset, set the argument 'forceNumeric' to TRUE

##### INITIAL/IMPORTANT DATA INFORMATION #####
setwd("SET_DIRECTORY")
if(!exists("kclust", mode="function")) source("../_R-functions.R")


#### Analyze Data ####

#testData<-read.delim("FILENAMW.tsv", row.names = 1) # The data is already log10 transformed with an log10(e-value=0) set to -180 (evalue = 1e-180)
#colnames(testData); colnames(testData)<-renameSH(colnames(testData)); dim(testData); colnames(testData); 

testData<-as.data.frame(rbind(c(1,2,0),c(0,1,0),c(1,1,0),c(0,1,1),c(1,1,0),c(1,0,0),c(0,0,0),c(1,1,1),c(0,1,1))); testData; colSums(testData)
numIsolate<-dim(testData)[2]
dim(analyze(Dataset = testData, whiteBelowAbsOne = FALSE, colorMode=TRUE, transformData=abs, rowLabels=paste(rownames(testData), sep="|"), showLegend=FALSE))[1] #, DEBUG = 1)
as.data.frame(rev(table(rowSums(sapply(testData,ifelse,1,0))))); legend("bottomleft", legend = numIsolate:1, col = rev(rainbow(numIsolate+1)[-(numIsolate+1)]), lty= 1, lwd = 2, cex=.8)
