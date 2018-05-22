# Several functions used in my heatmap analyses
# Moez Valliani, 2018/05/02 (combining functions from old scripts in one location)
#
# To include in an R project, run:
# if(!exists("kclust", mode="function")) source("C:/Users/Moez V/Desktop/Hsiang/R/_R-functions.R")  # Obtained from https://stackoverflow.com/questions/6456501/how-to-include-source-r-script-in-other-scripts


# Heatmap plotting area
library(gplots)
library(RColorBrewer)
##### DEFINING IMPORTANT FUNCTIONS #####
{
  # This section have various functions to produce the result we want without keep repeating it
  # Not all functions here may be useful to you do do not be afraid to remove some - LEAVE kclust() IN THOUGH
  
  renameSH <- function (isoList){
    if(!exists("convertNameTable")){convertNameTable<-read.delim("../renameSH.txt", row.names=1, header = F, colClasses = "character");}
    isoList<-tolower(isoList)
    for(iso in isoList){
      isoList <- replace(isoList, match(iso,isoList), convertNameTable[iso,1])
    }
    return(isoList)
  }
  
  
  # Clustered heatmap (k = number of subclusters). NoKlabel means for the names returned
  # For more information, please run: ?heatmap.2  (or help("heatmap.2"))
  # colRowFormat=function(z){rainbow(z)}){
  
  
  # Arguments (Default value if present):
  #   Data     = The matrix data you want the clustered heatmap for. Expects a dataset with named columns & rows
  #   Title    = Define if you want a title in the top center (great way to remind yourself which set you are looking at)
  #   K        = The number of clusters you want to be displayed SEPARATELY! 1 by default (shows overall cluster)
  #   noKlabel =  If TRUE, it does not append "K=#" to the label
  #   forceNumeric        = Set to TRUE (or 1 ... or anything except "", 0 and null) if your data is stating it cannot be displayed as it does not contain a numeric value.
  #   minColRng/maxColRng = Force a minimum/maximum value for the colour range 
  #   colRowData          = An array of the colours to be displayed with the rownames() set as the same as the data provided
  
  kclust <- function(Data, axis="column", orderPresence=TRUE, colorMode=FALSE, reorderFunction=ReorderFunction, 
                     Title="", K=1, noSplitK=TRUE, noKlabel=TRUE, forceNumeric=FALSE, colScale=NULL,
                     minColRng=FALSE, maxColRng=FALSE, colRowData=NULL, label_at_top=TRUE, whiteBelowAbsOne=TRUE
  ){
    if(forceNumeric){
      rName<-rownames(Data)
      Data<-sapply(Data,as.numeric)
      rownames(Data)<-rName
    }
    
    origData<-Data
    if(orderPresence==TRUE){
            Data[Data>1]<-1;
      ##Data[Data<1 && Data>-1]<-0;
      #      Data[Data<-1]<-(-0.5);
      reOrder<-order(rowSums(Data), decreasing = TRUE)
      colRowData<-colRowData[reOrder]
#      colRowData<-NULL
      Data<-origData[reOrder,]
    }
    #    reOrder<-order(rowSums(Data), decreasing = TRUE); Data<-Data[reOrder,]; colRowData<-colRowData[reOrder]
    # rowSums(Data*2^(1:dim(Data)[2])), 
    #reOrder<-order(colSums(Data), decreasing = TRUE); Data<-Data[,reOrder];
    
    hc.rows <- hclust(dist(Data)) # transpose the matrix and cluster columns
    hc.cols <- hclust(dist(t(Data)))  # draw heatmap for first cluster
    
    # Setting up the colour range with data set to 1 at the most
    
    if(colorMode==FALSE){
      # BW/presence mode
      Colpal<-c("white","black")
      ColRange <--1:1
    }else{
      IncrRng<-1.0 # Double the range
      totData=max(abs(floor(min(Data))), ceiling(max(Data)))*2  # dim(Data)[1] # max(Data)
      
      if(colorMode=="rainbow"){
        minN=ifelse(minColRng==FALSE,min(Data)-1,minColRng)
        maxN=ifelse(maxColRng==FALSE,max(Data),maxColRng)
        ColRange <- seq(minN,maxN, length.out=totData/2+2) # Note that totData above was doubled; here, we do not want that otherwise the scaling is off
        # warning(paste0("max(", maxN, ") - min(", minN, ") = ", maxN-minN, "; rngN = ", rngN, " ; totData = ", totData))
        #ColRange <- seq(ifelse(minColRng==FALSE,min(Data)-1,minColRng),ifelse(maxColRng==FALSE,max(Data),maxColRng), length.out=totData+2)
        
        Colpal <- c("white","black",colorRampPalette(rev(c("red","yellow","green","aquamarine","blue","purple")))(n = totData/2-1))
        if(totData==3){Colpal <- c("white","black","purple","green")}
      }else if(whiteBelowAbsOne==TRUE){
        Colpal <- c(colorRampPalette(rev(brewer.pal(11, "RdBu"))[1:5])(n = totData*IncrRng/2),
                    rep.int("white",2),
                    rev(colorRampPalette(brewer.pal(11, "RdBu")[1:5])(n=totData*IncrRng/2)))
        ColRange <- c(seq(ifelse(minColRng==FALSE,-1*(totData/2),minColRng),-1,1/IncrRng),
                      c(-0.99, 0, 0.99),
                      seq(1,ifelse(maxColRng==FALSE,totData/2,maxColRng),1/IncrRng))
      }else{
        Colpal <- c(colorRampPalette(rev(brewer.pal(11, "RdBu"))[1:5])(n = totData*IncrRng/2),
                    rep.int("white",2),
                    rev(colorRampPalette(brewer.pal(11, "RdBu")[1:5])(n=totData*IncrRng/2)))
        ColRange <- c(seq(ifelse(minColRng==FALSE,-1*(totData/2),minColRng),-1,1/IncrRng),
                      c(-0.01, 0, 0.01),
                      seq(1,ifelse(maxColRng==FALSE,totData/2,maxColRng),1/IncrRng))
      }

      
    }
    if(0||exists("DEBUG") && DEBUG){warning(paste0("Colpal: ", length(Colpal), " /// ColRange: ", length(ColRange)))}
    if(0||exists("DEBUG") && DEBUG){warning(paste0(" ", (Colpal)));warning(paste0(" ", (ColRange)))}
    if(!is.null(colScale)){Colpal <- colScale}
    
    # Splitting dendrogram in to define # of clusters defined by k (i.e. k-clusters [ = 1])
    cutData<-cutree(hc.rows,k=K)
    hm2<-c()
    for(N in 1:K){
      newTitle<-Title
      if(K>1){
        newTitle<-paste0(Title,"/nK=",N)
        subData<-Data[cutData==N,]
      }else{subData<-Data; newTitle<-Title}
      
      if(is.null(colRowData)){colRowData<-rep.int("#FFFFFF",dim(subData)[1]); names(colRowData)<-rownames(subData); colorProvided<-FALSE}
      
      hm<-heatmap.2(subData, main=newTitle, dendrogram=c(ifelse(axis=="column","column", ifelse(axis=="row","row",ifelse(axis=="none","none","both")))),
                    Rowv=(axis!="column"), Colv=(axis!="row"), revC=FALSE, 
                    RowSideColors=colRowData, col=Colpal, breaks=ColRange, density.info="none", 
                    key = TRUE, key.title="No. occurences", trace="none", symm=F,symkey=F,symbreaks=T, 
                    scale="none", na.rm=FALSE, cexCol=1, cexRow=1, reorderfun=function(d, w) reorder(d, w, agglo.FUN = reorderFunction))
      
      hm2<-c(hm2,if(noKlabel){rev(colnames(hm$carpet))}else{paste0("K=",N," | ",rev(colnames(hm$carpet)))})
      if(label_at_top){
        # Found information from https://stackoverflow.com/questions/17117753/r-how-to-build-angled-column-headings-above-columns-in-heatmap-2-pass-text-pl
        # Got locations of a pixel box and the avg point for n = 22 plotted heatmap
        # for small screen, set xrange max to 0.965. The value below is for the max x for the fullscreen on my laptop
        newPos.adj<-(0.005);newPos.xrange<-c(0.235,1.0152); newPos.n<-dim(subData)[2]; newPos.box<-(newPos.xrange[2]-newPos.xrange[1])/
          newPos.n; newPos.mid<-newPos.box/2; newPos.rng<-seq(newPos.xrange[1]+newPos.mid-newPos.adj, newPos.xrange[2]-newPos.mid-newPos.adj, newPos.box)
          text(labels(hm$carpet[,1]),x=newPos.rng, y=rep(0.81, newPos.n), xpd=TRUE, adj = 0, srt=90, cex = 0.80)
      }
      if(N<K){x=readline(prompt = paste0("Paused at k=",N,". Press <Enter> to continue..."))}
      
      
    }
    return(hm2)
  }
  
geomean <- function(x,MARGIN=2,noDim=FALSE,DEBUG=FALSE){
  # Obtaining a (modified) geometric mean [transformation using ln aka log_e]
  # x = the data set for the values to be obtained through a geometric mean (nth root of a product of the values)
  # MARGIN
  #    0 = Obtain a single value (return the geomean of the entire dataset)
  #    1 = Only apply to the rows for the data (across, per gene)
  #    2 = Only apply to the columns for the data (top down, per library)
  # na.rm = set to TRUE if you want NA values to be removed
  
  # This is actually a modified geometric mean. As log(0) will result in -Inf,
  # data points with 0s are ignored than accounted at the very end when the
  # second mean is obtained.
  
  # Geomean = antilog [ mean (log values of non-0 values) ]
  # Mod Geomean = mean of Geomean + the excluded 0 values for each time frame
  
  ## Old script:
  ##  Data<-log(x)
  ##  logx<-sapply(x,function(z)ifelse(z==0,log_of_0,log(z)))
  ##  if(MARGIN==0){return(exp(mean(logx, na.rm=na.rm)))}
  ##  else{return(apply(logx,MARGIN,function(z)exp(mean(z, na.rm=na.rm))))}
  x2<-log(x);x2<-sapply(x2, function(z)ifelse(is.infinite(z),NA,z))
  if(MARGIN==0){x3<-mean(x2,na.rm=TRUE)}else{x3<-apply(x2,2,mean,na.rm=TRUE)}
  x4<-exp(x3)
  x5<-x4
  if(MARGIN==0){
    lenMax<-length(x)
    n<-lenMax
    x5<-sum(x5,na.rm=T)*(sum(!is.na(x2)))/n
  }
  else{
    lenMax<-dim(x)[2]
    n<-dim(x)[1]
    for(j in 1:lenMax){
      nn<-sum(sapply(x2[,j],is.na))}
    if(nn>0){x5[j]<-(x5[j]*(n-nn)/nn)}#mean(c(x5[j], rep.int(0, n)))}
  }
  if(DEBUG){return(paste0("x: ",x,"    x2: ",x2,"    x3: ",x3,"    x4: ",x4,"    x5: ",x5))}
  return(x5)
}

plotSlope<-function(x, xlabel="X", ylabel="Y",title="",subtitle="",txtX=2,txtY=10){
  # Shows the plot with the slope and line equation displayed. As per http://stackoverflow.com/questions/13114539/how-to-add-rmse-slope-intercept-r2-to-r-plot
  fit<-lm(x)
  rmse <- round(sqrt(mean(resid(fit)^2)), 2)
  coefs <- coef(fit)
  b0 <- round(coefs[1], 2)
  b1 <- round(coefs[2],2)
  r2 <- round(summary(fit)$r.squared, 2)
  eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~ 
                  r^2 == .(r2) * "," ~~ RMSE == .(rmse))
  plot(x,main=title,sub=subtitle,xlab=xlabel,ylab=ylabel)
  abline(fit)
  text(txtX, txtY, eqn, pos = 4)
  
}

logTransform <- function(x,Transformed_0_value_for_log=-Inf,MARGIN=2, Base=logBase, applyInt=FALSE){
  # A simple function to transform the dataset given. Any value with 0 with be replaced with given whereas non-0 will be transformed using log
  if(MARGIN==0){y<-sapply(x,function(z) ifelse(z==0,Transformed_0_value_for_log,log(z, base=Base)))}
  else{y<-apply(x,MARGIN,function(z) ifelse(z==0,Transformed_0_value_for_log,log(z, base=Base)))}
  if(applyInt==TRUE){y<-as.integer(y)}
  return(y)
}
  
  ##### DATA ENTRY #####
  #! CHANGE DATA TO WHAT IS NEEDED
  # Produce core results based on the geometric mean e-values (the query isolates were used as the isolate grouping and the subject hits were subject to the geometric mean)
  
  
  # A function to make everything more efficient.
  analyze<-function(Dataset, DataCols=1:numIsolate,Series="", SeriesCol=min(numIsolate+1, dim(Dataset)[1]), colScale=NULL,
                    colorMode=TRUE, rowLabels=NULL, transformData=NULL, showLegend=TRUE, Title="", 
                    clustAxis="column", orderPresence=TRUE, ReorderFunction=sum, Ksplit=1, whiteBelowAbsOne=TRUE,
                    ForceNumeric=FALSE, DEBUG=FALSE ){
    # Series    = what subset of data is wanted (e.g. AA, CBM, etc. Leave blank for all)
    # SeriesCol = the column in the provided Dataset which refers to the Series data
    # DataCols  = the columns where the data we want viewed
    # Dataset   = the expected data
    data<-Dataset[,DataCols]; data.title<-"All"
    #  data<-data[order(rowSums(data), decreasing = TRUE),]
    
    if(!is.null(rowLabels)){rownames(data)<-rowLabels}
    if(Series!=""){
      data<-data[grep(paste0("^",Series), Dataset[,SeriesCol]),]
      data.title <- Series
    }
    nIso <- dim(data)[2]
    
    if(!is.null(transformData)){rowLabels<-rownames(data); data<-sapply(data,transformData); rownames(data)<-rowLabels}
    #print(head(data));stop("Debug")
    if(DEBUG){warning("DEBUG mode");return(data);}
    
    # Setting up the row colours wanted to show
    rowValues<-rowSums(data > 0)
    rowColors<-rainbow(nIso+1) 
    colRowData<-rowColors[rowValues]
    names(colRowData)<-ifelse(is.null(rowLabels),rownames(data),rowLabels)
    if(DEBUG){warning("DEBUG mode");return(colRowData);}
    
    
    if(Title!=""){data.title<-Title}
    write.table(
      kclust(data, axis = clustAxis, forceNumeric = ForceNumeric, Title = data.title, orderPresence=orderPresence,
             colorMode = colorMode, colRowData=colRowData, colScale=colScale, K = Ksplit,
             whiteBelowAbsOne=whiteBelowAbsOne, reorderFunction=ReorderFunction),
      "x.tsv",sep='\t',quote=F, append=F, col.names=T)
    
    # A legend for the provided row colours
    if(showLegend==TRUE){
      legend("bottomleft",
             legend = 1:numIsolate,
             col = rowColors, 
             lty= 1, lwd = 2, cex=.8
      )
    }
    return(data)
  }
  
}