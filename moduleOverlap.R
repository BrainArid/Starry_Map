#countsMatrix - an integer matrix of counts representing pairwise intersection between module sets
#lengths1 - vector of the lengths of the modules found in the rows of countMatrix
#lengths2 - vector of the lengths of the modules found in the cols of countMatrix
#subuniverseSize - size of the gene universe
starryMutualInformation <- function(countsMatrix, lengths1, lengths2, subuniverseSize)
{
   #Calculate symmetric uncertainty mutual information- p(x,y)log2(p(x,y)/(p(x)p(y)))
    I<-0;
    Hi<-0;
    Hj<-0;
    rS<-rowSums(countsMatrix);
    cS<-colSums(countsMatrix);
    
    #calculate the entropy Hi and Hj
    for(i in 1:length(lengths1))
    {
      p<-rS[i]/subuniverseSize;
      Hi<-Hi-(p*log2(p));
    }
    for(j in 1:length(lengths2))
    {
      p<-cS[j]/subuniverseSize;
      Hj<-Hj-(p*log2(p));
    }
    for(i in 1:length(lengths1))
      for(j in 1:length(lengths2))
      {
        if(countsMatrix[i,j]>0)
          I <- I + countsMatrix[i,j]* log2(subuniverseSize*countsMatrix[i,j]/(lengths1[i]*lengths2[j]));
      }
  return(2*I/subuniverseSize/Hi/Hj);
}

#input is expected to be module per row otherwise specify clusts#Columned=TRUE
#clustsFiles[1/2] - files containing cluster elements organized into either rows or columns
#clusts[1/2]Columned - indicates if clustsFiles[1/2] is organized in columns
#clusts[1/2]Names - indicates if clustsFiles[1/2] contains header with cluster/module names
#outDir - output directory address
#threshold - Fisher's Exact test significance threshold
#clustHeatMap - legacy to be removed
#fisherExact - indicates if FET is performed
#histogram - legacy to be removed
#geneUnivers3[1/2] - a string vector of genes that define the gene universe
#sort - when boolean this indicates if the resulting starry map is ordered by module size, 
#       otherwise, a numeric vector to define the module orderings
starryMap <- function(clustsFile1, clusts1Columned=FALSE, clusts1Names=FALSE, 
                          clustsFile2, clusts2Columned=FALSE, clusts2Names=FALSE, 
                          outDir, threshold, clustHeatMap=FALSE, fisherExact=FALSE, 
                          histogram=FALSE, geneUniverse1, geneUniverse2, sort1=TRUE, sort2=TRUE)
{ 
  output<-list();
  
  print("opening files");
  clusts1 <- read.table(file=clustsFile1,sep = "\t",stringsAsFactors = FALSE,fill = TRUE)
  clusts2 <- read.table(file=clustsFile2,sep = "\t",stringsAsFactors = FALSE,fill = TRUE)
  print("files opened");
  
  #rotate if necessary
  if(clusts1Columned)
    clusts1<-t(clusts1)
  if(clusts2Columned)
    clusts2<-t(clusts2)
  #remove first column if module names found
  if(clusts1Names)
  {
    row.names(clusts1)<-factor(clusts1[,1])
    clusts1<- clusts1[,-1] 
  }else
  {
    n<- paste0("module", row.names(clusts1))
    row.names(clusts1)<- factor(n,levels=n)
    rm(n);
  }
  if(clusts2Names)
  {
    row.names(clusts2)<-factor(clusts2[,1])
    clusts2<- clusts2[,-1]
  } else
  {
    n<- paste0("module", row.names(clusts2))
    row.names(clusts2)<- factor(n,levels=n)
    rm(n);
  }
  #remove any columns containing NAs
  clusts1 <- clusts1[,colSums(is.na(clusts1))==0];
  clusts2 <- clusts2[,colSums(is.na(clusts2))==0];
  
  #get intersection of universes
  #browser();
  geneUniverseIntersect<-intersect(x=geneUniverse1,y=geneUniverse2);
  totalGenes<-length(geneUniverseIntersect);
  subUniverse<-union(clusts1,clusts2);
  subUniverseSize<-length(subUniverse);
  
  #intersect each clust with universe Intersection
  lengths1<-vector(length = dim(clusts1)[1])
  for(i in 1:dim(clusts1)[1])
  {
    temp <- intersect(x=clusts1[i,][clusts1[i,]!=""],
                             y=geneUniverseIntersect);
    clusts1[i,]<-"";
    clusts1[i,1:length(temp)]<-temp;
    lengths1[i]<-length(clusts1[i,][clusts1[i,]!=""]);
  }
  
  lengths2<-vector(length = dim(clusts2)[1])
  for(i in 1:dim(clusts2)[1])
  {
    temp <- intersect(x=clusts2[i,][clusts2[i,]!=""],
                      y=geneUniverseIntersect);
    clusts2[i,]<-"";
    clusts2[i,1:length(temp)]<-temp;
    lengths2[i]<-length(clusts2[i,][clusts2[i,]!=""]);
  }
  remove(temp);
  
  #sort clusters' modules by module length
  if(class(sort1)=="numeric" || sort1!=FALSE)
  {
    if(class(sort1)=="numeric" || class(sort1)=="integer" || class(sort1)=="integer"){
      order1<- sort1;
    }
    else
    {
      order1<- order(lengths1,decreasing=TRUE);
    }
    clusts1 <- clusts1[order1, ];
    lengths1<- lengths1[order1];
    output$order1 <- order1;
  }
  if(class(sort2)=="numeric" || class(sort2)=="integer" || sort2!=FALSE)
  {
    if(class(sort2)=="numeric" || class(sort2)=="integer"){
      order2<- sort2;
    }
    else
    {
      order2<- order(lengths2,decreasing=TRUE);
    }
    clusts2 <- clusts2[order2, ];
    lengths2<- lengths2[order2];
    output$order2 <- order2;
  }
  output$lengths1<-lengths1;
  output$lengths2<-lengths2;
  countsMat <- matrix(nrow = dim(clusts1)[1], ncol= dim(clusts2)[1]);
  row.names(countsMat)<-row.names(clusts1)
  colnames(countsMat)<-row.names(clusts2)
  countsMat[,]<-0;
  rMat <- matrix(nrow = dim(clusts1)[1], ncol= dim(clusts2)[1]);
  row.names(rMat)<-row.names(clusts1)
  colnames(rMat)<-row.names(clusts2)
  rMat[,]<-0;
  gMat <- matrix(nrow = dim(clusts1)[1], ncol= dim(clusts2)[1]);
  row.names(gMat)<-row.names(clusts1)
  colnames(gMat)<-row.names(clusts2)
  gMat[,]<-0;
  bMat <- matrix(nrow = dim(clusts1)[1], ncol= dim(clusts2)[1]);
  row.names(bMat)<-row.names(clusts1)
  colnames(bMat)<-row.names(clusts2)
  bMat[,]<-0;
  if(fisherExact)
  {
    feMat <- matrix(nrow = dim(clusts1)[1], ncol= dim(clusts2)[1],dimnames=list(row.names(clusts1),row.names(clusts2)));
    row.names(feMat)<-row.names(clusts1)
    colnames(feMat)<-row.names(clusts2)
    feMat[,]<-0;
    feMat_est<-feMat;
  }
  
  print("Processing Overlaps...")
  for(i in 1:dim(clusts1)[1])
  {
    for(j in 1:dim(clusts2)[1])
    {
      print(paste0("Cross comparing row ", i, " of ", dim(clusts1)[1]," and column ", j, " of ", dim(clusts2)[1]));
      countsMat[i,j] <- length(intersect(clusts1[i,][clusts1[i,]!=""],clusts2[j,][clusts2[j,]!=""]))
      r<- lengths1[i];
      g<- countsMat[i,j];
      b<- lengths2[j];
      x<- r + b - g;
      rMat[i,j] <- r/x;
      gMat[i,j] <- g/x;
      bMat[i,j] <- b/x;
      if(fisherExact)
      {
        ft<-fisher.test(matrix(data=c(g,r-g,b-g,totalGenes-(r+b-g)),nrow=2,ncol=2))
        feMat[i,j] <- ft$p.value
        feMat_est[i,j] <- ft$estimate
      }
    }
  }
  
  output$intCounts <- countsMat;
  output$starryMI <- starryMutualInformation(countsMatrix=countsMat,lengths1=lengths1,lengths2=lengths2,subuniverseSize=subUniverseSize);
  
  #visualize
  print("Drawing pretty pictures...")
  library("ggplot2")
  library("reshape")
  library("plyr")
  library("scales")

  rMat.m <- melt(rMat)
  rMat.m$X1 <- factor(rMat.m$X1, levels = row.names(clusts1))
  rMat.m$X2 <- factor(rMat.m$X2, levels = row.names(clusts2))

  gMat.m <- melt(gMat)
  gMat.m$X1 <- factor(gMat.m$X1, levels = row.names(clusts1))
  gMat.m$X2 <- factor(gMat.m$X2, levels = row.names(clusts2))

  #if(clustHeatMap)
  #{
  #  p <- ggplot(gMat.m, aes(y=X1, x=X2))+ 
  #    geom_tile(aes(y=X1,x=X2,fill = rgb(red = 0,green=value,blue=0)), colour = "black")+scale_fill_identity()+
  #    ylab(label=basename(clustsFile1))+
  #    xlab(label=basename(clustsFile2))+
  #    theme(axis.text.x=element_text(angle=90));
  #  ggsave(filename = paste0(outDir, basename(clustsFile1), "_VS_", basename(clustsFile2), "_INTERSECT.png"),plot = p,width=3,height=3)
  #}
   if(fisherExact)
   {
      feMat.m <- melt(round(feMat,5))
      feMat.m$X1 <- factor(feMat.m$X1, levels = row.names(clusts1))
      feMat.m$X2 <- factor(feMat.m$X2, levels = row.names(clusts2))
      
      feMat_est.m <- melt(round(feMat_est,5))
      feMat_est.m$X1 <- factor(feMat_est.m$X1, levels = row.names(clusts1))
      feMat_est.m$X2 <- factor(feMat_est.m$X2, levels = row.names(clusts2))
#       p <- ggplot(feMat.m, aes(y=X1, x=X2))+
#         geom_tile(aes(y=X1,x=X2,fill = rgb(red = 0,green=1-value,blue=0)), colour = "black")+scale_fill_identity()+
#         ylab(label=basename(clustsFile1))+
#         xlab(label=basename(clustsFile2))+
#         theme(axis.text.x=element_text(angle=90));
#       ggsave(filename = paste0(outDir, basename(clustsFile1), "_VS_", basename(clustsFile2), "_FISHER.png"),plot = p,width=3,height=3)
   }
  
  bMat.m <- melt(bMat)
  bMat.m$X1 <- factor(bMat.m$X1, levels = row.names(clusts1))
  bMat.m$X2 <- factor(bMat.m$X2, levels = row.names(clusts2))
  maxSize<-1;
  colMat.m <- data.frame(X1=rMat.m$X1,X2=rMat.m$X2,r=rMat.m$value,g=gMat.m$value,b=bMat.m$value,fe=feMat.m$value,feset=feMat_est.m$value,feSize=(feMat.m$value<threshold)*maxSize/2,feEst=feMat_est.m$value>1)
  colMat.m$X1 <- factor(colMat.m$X1, levels = row.names(clusts1))
  colMat.m$X2 <- factor(colMat.m$X2, levels = row.names(clusts2))
  output$colMat <- colMat.m;
  if(clustHeatMap)
  {
    p <- ggplot(colMat.m, aes(y=X1, x=X2)) +
      geom_tile(aes(y=X1,x=X2,fill = rgb(red = r,green=g,blue=b),width=maxSize,height=maxSize), colour = "black")+
      scale_fill_identity()+
      ylab(label=basename(clustsFile1))+
      xlab(label=basename(clustsFile2))+
      theme(axis.text=element_text(size=9),axis.text.x=element_text(angle=-90,hjust=0,vjust=0.5));
    if(fisherExact)
    {
      p <- p + geom_tile(aes(y=X1,x=X2,fill = rgb(red = 0,green=feEst,blue=0),width=feSize,height=feSize), colour = "black");
    }
    ggsave(filename = paste0(outDir, basename(clustsFile1), "_VS_", basename(clustsFile2), "_COMPOSITE.png"),plot = p,width=7,height=7)
    output$fig <- p;
  }

#   if(histogram)
#   {
#     
#     p <- ggplot(colMat.m[colMat.m$g>0.25,], aes(x=g))+ geom_histogram(binwidth=0.02)+
#       xlab(label="Intersection Size over Union Size")+
#       ylab(label="module count");
#     ggsave(filename = paste0(outDir, basename(clustsFile1), "_VS_", basename(clustsFile2), "_HIST.png"),plot = p,width=3,height=3)
#     
#   }

  #list white tiles:
  if(dim(colMat.m[colMat.m$fe<threshold,])[1]>0)
  {
    print(colMat.m[colMat.m$fe<threshold,]);
    
    x <- colMat.m[colMat.m$fe<threshold,]$X1
    y <- colMat.m[colMat.m$fe<threshold,]$X2
    
    commonModules <- list();
    commonModules$profiles <- list();
    commonModules$clusts1 <- list();
    commonModules$clusts2 <- list();
    
    relativeRank <- function(L1, L2)
    {      
      l1 <- length(L1);
      l2 <- length(L2);
      m1<-match(L1,L2);
      m2<-match(L2,L1);
      
      return(list(r1=m1-seq(from = 1,to = l1,1),r2=m2-seq(from = 1,to = l2,1)));
    }
    
    print(paste0("Selecting modules for output that meet overlap threshold of ", threshold))
    for(i in 1:length(x))
    {
      print(paste0("Getting overlap ", i, " of ", length(x)) )
      commonModules$profiles[i] <- list(colMat.m[colMat.m$X1==x[i] & colMat.m$X2==y[i],c("X1","X2","r","g","b", "fe")]);
      #if(fisherExact)
      #{
      #  commonModules$profiles[[i]]$fe <- feMat.m[feMat.m$X1==x[i] & feMat.m$X2==y[i],c("value")]
      #}
      commonModules$clusts1[i] <- list(clusts1[x[i],clusts1[x[i],]!=""]);
      commonModules$clusts2[i] <- list(clusts2[y[i],clusts2[y[i],]!=""]);
    }
    
    #and thier profiles:
    gSortedOrder <- order(matrix(data=unlist(commonModules$profiles),nrow=length(commonModules$profiles),ncol=6,byrow=TRUE)[,6],decreasing=TRUE)
    commonModules <- list(profiles=commonModules$profiles[gSortedOrder],
                          clusts1=commonModules$clusts1[gSortedOrder],
                          clusts2=commonModules$clusts2[gSortedOrder]);
    output$commonModules<-commonModules;
  }
  
  #temp code for Tara
  #browser();
  #options(scipen=999)
  #df <- melt(data = output$intCounts);
  #names(df)<- c("A Module Label","B Module Label", "Overlap");
  #df[,"A Module Size"] <- melt(matrix(data = output$lengths1,nrow = dim(output$intCounts)[1],ncol = dim(output$intCounts)[2]))[,3];
  #df[,"B Module Size"] <- melt(matrix(data = output$lengths2,nrow = dim(output$intCounts)[1],ncol = dim(output$intCounts)[2],byrow = TRUE))[,3];
  #df[,"FE pValue"] <- melt(round(feMat,8))[,3];
  #df[,"Size Ratio"] <- apply(df,1,function(i)
  #  {
  #    return(round(min(as.numeric(i["A Module Size"]),
  #               as.numeric(i["B Module Size"])) / 
  #             max(
  #               as.numeric(i["A Module Size"]),
  #               as.numeric(i["B Module Size"])),8));
  #  })
  #write.csv(df[,c(1,2,4,5,7,3,6)],file = paste0(basename(clustsFile1),"_VS_",basename(clustsFile2),".csv"),quote=FALSE,row.names = FALSE);
  return(output);
}

#input is expected to be module per row otherwise specify clusts#Columned=TRUE
#inputs are analogous to those in staryMap
starryMapList <- function(clustsFiles, clustsColumned, clustsNames, 
                      outDir, threshold, clustHeatMap=FALSE, fisherExact=FALSE,
                      histogram=FALSE, geneUniverses, sorts)
{
  #all self-comparisons first
  outputs<-matrix(list(),nrow = length(clustsFiles), ncol=length(clustsFiles))
  for(i in 1:length(clustsFiles))
  {
     outputs[[i,i]]<-starryMap(clustsFile1 = clustsFiles[i],clusts1Columned = clustsColumned[i],clusts1Names = clustsNames[i],
                               clustsFile2 = clustsFiles[i],clusts2Columned = clustsColumned[i],clusts2Names = clustsNames[i],
                               outDir, threshold, clustHeatMap, fisherExact,
                               histogram, geneUniverses[i], geneUniverses[i], sort[i], sort[i])
  }
  
  #all other pair-wise comparisons
  for(i in 1:length(clustsFiles))
  {
    for(j in i+1:length(clustsFiles))
    {
      outputs[[i,j]]<-starryMap(clustsFile1 = clustsFiles[i],clusts1Columned = clustsColumned[i],clusts1Names = clustsNames[i],
                                clustsFile2 = clustsFiles[j],clusts2Columned = clustsColumned[j],clusts2Names = clustsNames[j],
                                outDir, threshold, clustHeatMap, fisherExact,
                                histogram, geneUniverses[i], geneUniverses[j], sort[i], sort[j])
    }
  }
}

#prints starryMap results in text form
print.StarryMap <- function(starryMap)
{
  
}