#run moduleOverlap
source("moduleOverlap.R")

print("Reading in command line arguments.");
args <- commandArgs(trailingOnly = TRUE);
print(paste0("commandArgs: ",args));

if(length(args) > 0)
{
  #Parse arguments (we expec the form --argName=argValue)
  parseArgs <- function (x) 
  {
    s<- unlist(strsplit(sub("^--","",x), "="));
    return(list(V1=s[1],V2=paste(s[-1],collapse = "=")))
  }
  argsDF <- as.data.frame(do.call("rbind", lapply(X = args,FUN = parseArgs)));
  args <- as.list(argsDF$V2)
  names(args) <- argsDF$V1
  rm(argsDF)
} else
{
  args <- list();
}

print(paste0("commandArgs: ",args));

#initialize arguments if 
initializeBooleanArg <- function(arg, default){
  if(is.null(arg))
  {
    arg <- default;
  } else if(is.character(arg))
  {
    arg <- as.logical(arg);
  }
  return(arg);
}

initializeStringArg <- function(arg, default){
  if(is.null(arg))
  {
    arg <- default;
  } else if(!is.character(arg))
  {
    arg <- as.character(arg);
  }
  return(arg);
}

initializeFloatArg <- function(arg, default){
  if(is.null(arg))
  {
    arg <- default;
  } else if(!is.numeric(arg))
  {
    arg <- as.numeric(arg);
  }
  return(arg);
}

initializeIntArg <- function(arg, default){
  if(is.null(arg))
  {
    arg <- default;
  } else if(!is.integer(arg))
  {
    arg <- as.integer(arg);
  }
  return(arg);
}

args$dir1 <- initializeStringArg(arg=args$dir1, default="Test Sets/");
args$dir2 <- initializeStringArg(arg=args$dir2, default="Test Sets/");
args$clustsFile1 <- initializeStringArg(arg=args$clustsFile1, default="GSE57872StarryMap.txt");
args$clusts1Columned <- initializeBooleanArg(arg=args$clusts1Columned, default=TRUE);
args$clusts1Names <- initializeBooleanArg(arg=args$clusts1Names, default=TRUE);
args$clustsFile2 <- initializeStringArg(arg=args$clustsFile2, default="GSE48865StarryMap.txt");
args$clusts2Columned <- initializeBooleanArg(arg=args$clusts2Columned, default=TRUE);
args$clusts2Names <- initializeBooleanArg(arg=args$clusts2Names, default=TRUE);
args$outDir <- initializeStringArg(arg=args$outDir, default="Figures");
args$threshold <- initializeFloatArg(arg=args$threshold, default=0.001);
args$clustHeatMap <- initializeBooleanArg(arg=args$clustHeatMap, default=TRUE);
args$fisherExact <- initializeBooleanArg(arg=args$fisherExact, default=TRUE);
args$permutationTest <- initializeBooleanArg(arg=args$permutationTest, default=FALSE);
args$histogram <- initializeBooleanArg(arg=args$histogram, default=TRUE);
args$geneUniverseFile1 <- initializeStringArg(arg=args$geneUniverseFile1, default="Test Sets/geneOrder.txt");
args$geneUniverseFile2 <- initializeStringArg(arg=args$geneUniverseFile2, default="Test Sets/GSE48865_geneOrder.txt");

if(DEBUG)
{
  dir1<-args$dir1
  dir2<-args$dir2
  clustsFile1<-paste0(args$dir1, args$clustsFile1)
  clusts1Columned<-args$clusts1Columned
  clusts1Names<-args$clusts1Names
  clustsFile2<-paste0(args$dir2, args$clustsFile2)
  clusts2Columned<-args$clusts2Columned
  clusts2Names<-args$clusts2Names
  outDir<-args$outDir
  threshold<-args$threshold
  permutationTest<-args$permutationTest
  clustHeatMap<-args$clustHeatMap
  fisherExact<-args$fisherExact
  histogram<-args$histogram
}

#compute single-cell vs bulk sample gene-level comparison figure
geneUniverse1 <- as.character(read.table(file=args$geneUniverseFile1,skip=1,header=FALSE)[,2]);
geneUniverse2 <- as.character(read.table(file=args$geneUniverseFile2,skip=1,header=FALSE)[,2]);

output <-moduleOverlap(paste0(args$dir1, args$clustsFile1),args$clusts1Columned, 
                       args$clusts1Names, paste0(args$dir2, args$clustsFile2), args$clusts2Columned, 
                       args$clusts2Names, args$outDir, args$threshold, args$clustHeatMap, 
                       args$fisherExact, args$permutationTest, args$histogram, 
                       geneUniverse1=geneUniverse1, geneUniverse2=geneUniverse2,
                       sort1=c(20,3,12,13,15,10,7,11,1,2,5,6,9,16,4,8,19,21,14,17,18),
                       sort2=c(24,37,34,1,22,23,46,26,4,39,43,10,21,18,38,11,17,25,2,3,13,29,40,5,7,47,6,16,19,20,27,30,41,44,8,9,14,15,28,31,32,33,35,36,42,45,48,12));
p3<-output$fig;
output <-moduleOverlap(paste0(args$dir1, args$clustsFile1),args$clusts1Columned, 
                       args$clusts1Names, paste0(args$dir1, args$clustsFile1),args$clusts1Columned, 
                       args$clusts1Names, args$outDir, args$threshold, args$clustHeatMap, 
                       args$fisherExact, args$permutationTest, args$histogram, 
                       geneUniverse1=geneUniverse1, geneUniverse2=geneUniverse1);
p1<-output$fig;
output <-moduleOverlap(paste0(args$dir2, args$clustsFile2),args$clusts2Columned, 
                       args$clusts2Names, paste0(args$dir2, args$clustsFile2), args$clusts2Columned, 
                       args$clusts2Names, args$outDir, args$threshold, args$clustHeatMap, 
                       args$fisherExact, args$permutationTest, args$histogram, 
                       geneUniverse1=geneUniverse2, geneUniverse2=geneUniverse2);
p2<-output$fig;

#compute pairwise single-cell CODENSE output comparisons 
moduleFiles1 <- list.files(pattern="G=5_E=3_D=0\\.[4-8]_Q=0\\.4_B=0\\.4_S=80_C=0\\.6\\.modulesFO\\.hgnc\\.david\\.module$",path=args$dir1)
moduleFiles2 <- moduleFiles1;#list.files(pattern="G=5_E=1_D=0\\.[4-8]_Q=0\\.4_B=0\\.4_S=80_C=0\\.6\\.modulesFO\\.hgnc\\.david\\.module$",path=args$dir2)

stats <- matrix(data=0,nrow=length(moduleFiles1),ncol=length(moduleFiles2));
orders<-list(d4=c(16,3,12,15,13,10,7,1,9,6,2,11,5,8,14,4,17,18,19),
             d5=c(20,3,12,22,13,15,10,7,1,11,17,6,9,2,5,23,4,8,18,21,14,16,19),
             d6=c(20,3,12,13,15,10,7,11,1,2,5,6,9,16,4,8,19,21,14,17,18),
             d7=c(16,3,12,14,20,17,22,10,18,7,11,5,6,9,1,2,4,15,8,13,19,21),
             d8=c(18,3,14,16,11,13,15,20,24,9,10,22,5,7,8,21,23,2,1,4,6,12,17,19)
)
i<-1;
j<-1;
for(moduleFile1 in moduleFiles1)
{
  for(moduleFile2 in moduleFiles2)
  {
    output <-moduleOverlap(paste0(args$dir1, moduleFile1), args$clusts1Columned, 
                           args$clusts1Names, paste0(args$dir1, moduleFile2), args$clusts1Columned, 
                           args$clusts1Names, args$outDir, args$threshold, args$clustHeatMap, 
                           args$fisherExact, args$permutationTest, args$histogram,
                           geneUniverse1=geneUniverse1, geneUniverse2=geneUniverse1,
                           sort1=orders[[i]],sort2=orders[[j]]);
    stats[i,j]<-output$stat;
    j<-j+1;
  }
  j<-1;
  i<-i+1;
}

#...


fileName<-paste0(args$outDir, basename(args$clustsFile1), "_VS_", basename(args$clustsFile2), "_MATCHING_MODULES.csv");
if(is.null(commonModules))
{
  write("No consensus.",fileName);
} else 
{
  unlistAndWrite <- function(x, index, file, append=TRUE, ncolumns=1000, sep=";")
  {
    write(unlist(x[[index]]),file,ncolumns=ncolumns,append = TRUE, sep=sep);
  }
  
  #if(args$enrich)
  #{
  #  GOCharts1 <- commonModules$clusts1GO
  #  GOCharts2 <- commonModules$clusts2GO
  #  commonModules$clusts1GO <- NULL;
  #  commonModules$clusts2GO <- NULL;
  #}
  
  write(paste0(basename(args$clustsFile1), ";", basename(args$clustsFile2)),fileName,append=FALSE);#FALSE to clear file contents
  write(paste0("number of matching modules:,",length(commonModules[[1]])),fileName,append=TRUE);
  write("X1 index;X2 index;r;g;b;fe",fileName,append=TRUE);
  
  for(i in 1:length(commonModules[[1]]))
  {
    lapply(commonModules, unlistAndWrite, i, fileName);
    write(x = paste("\n",args$clustsFile1,";;;;;;;;;;",args$clustsFile2,"\n"),fileName,append = TRUE);
    #       if(args$enrich)
    #       {
    #         write(x = "BIOLOGICAL PROCESS",fileName,append = TRUE);
    #         write.table(x=cbind(GOCharts1[[i]]$BP,x=rep(x="",times = length(GOCharts1[[i]]$BP[,1])),GOCharts2[[i]]$BP),file = fileName,append = TRUE,sep = ";",quote = FALSE,row.names = FALSE,col.names = TRUE);
    #         write(x = "MOLECULAR FUNCTION",fileName,append = TRUE);
    #         write.table(x=cbind(GOCharts1[[i]]$MF,x=rep(x="",times = length(GOCharts1[[i]]$MF[,1])),GOCharts2[[i]]$MF),file = fileName,append = TRUE,sep = ";",quote = FALSE,row.names = FALSE,col.names = TRUE);
    #         write(x = "CELL CYCLE",fileName,append = TRUE);
    #         write.table(x=cbind(GOCharts1[[i]]$CC,x=rep(x="",times = length(GOCharts1[[i]]$CC[,1])),GOCharts2[[i]]$CC),file = fileName,append = TRUE,sep = ";",quote = FALSE,row.names = FALSE,col.names = TRUE);
    #         write(x = "\n",fileName,append = TRUE);
    #       }
  }
}
