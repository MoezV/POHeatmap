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

#testData.data<-read.delim("FILENAMW.tsv", row.names = 1) # The data is already log10 transformed with an log10(e-value=0) set to -180 (evalue = 1e-180)
#colnames(testData.data); colnames(testData.data)<-renameSH(colnames(testData.data)); dim(testData.data); colnames(testData.data); 

testData.data<-as.data.frame(rbind(c(0,1,0),c(1,1,0),c(0,1,1),c(1,0,1),c(0,0,0),c(1,1,0))); testData.data
numIsolate<-dim(testData.data)[2]
dim(analyze(Dataset = testData.data, Series = "", colorMode=TRUE, transformData=abs, rowLabels=paste(rownames(testData.data), sep="|"), showLegend=FALSE))[1] #, DEBUG = 1)
as.data.frame(rev(table(rowSums(sapply(testData.data,ifelse,1,0))))); legend("bottomleft", legend = numIsolate:1, col = rev(rainbow(numIsolate+1)[-(numIsolate+1)]), lty= 1, lwd = 2, cex=.8)
