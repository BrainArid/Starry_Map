source("StarryModuleSet.R");

StarryMap <- setRefClass(Class = "StarryMap", 
                      fields = list(
                        #public
                        moduleSet1="StarryModuleSet",
                        moduleSet2="StarryModuleSet",
                        interMat="matrix", #matrix of intersection counts
                        
                        #private
                        .totalElements="integer",
                        .interMat="matrix",
                        .exclusMat1="matrix",
                        .exclusMat2="matrix",
                        .unionMat="matrix",
                        .fePMat="matrix",
                        .feEstMat="matrix"
                      ),
                      
                      methods = list(
                        initialize=function(...)
                        {
                          .self$initFields(...);
                          .self$.recalcTotalElements();
                          return(.self);
                        },
                        show=function()
                        {
                          cat("An object of class ", class(.self), "\n");
                          .self$moduleSet1$show();
                          .self$moduleSet2$show();
                          invisible(NULL);
                        },
                        
                        #contains side effects
                        intersect.elements=function(otherElements)
                        {
                          .self$moduleSet1$intersect.elements(otherElements);
                          .self$moduleSet2$intersect.elements(otherElements);
                        },
                        
                        #contains side effects
                        intersect.starryModule=function(otherStarryModule)
                        {
                          .self$moduleSet1$intersect.starryModule(otherStarryModule);
                          .self$moduleSet2$intersect.starryModule(otherStarryModule);
                        },
                        
                        #contains side effects
                        .recalcTotalElements=function()
                        {
                          .self$.totalElements <- length(union(.self$moduleSet1$universe$elements, .self$moduleSet1$universe$elements));
                        },
                        
                        #reduces all elements to the intersection of module set universes
                        reduceElements=function()
                        {
                          #copy module set 1's universe
                          interUniverse<- .self$moduleSet1$universe$copy();
                          #intersect with module set 2's universe
                          interUniverse$intersect.starryModule(.self$moduleSet2$universe);
                          #intersect moduleSets
                          .self$moduleSet1$intersect.starryModule(interUniverse);
                          .self$moduleSet2$intersect.starryModule(interUniverse);
                          #recalculate total elements 
                          .self$.totalElements <- interUniverse$length();
                        },
                        
                        .calcPairwiseProfile=function()
                        {
                          mods1 <- .self$moduleSet1;
                          mods2 <- .self$moduleSet2;
                          print("Processing Overlaps...")
                          
                          .self$.interMat <- sapply(X=1:mods2$length(), 
                                                 FUN = function(j)
                                                 { sapply(X=1:mods1$length(),
                                                          FUN = function(i)
                                                          {
                                                            length(intersect(x = mods1$modules[[i]]$elements,
                                                                      y = mods2$modules[[j]]$elements));
                                                          })
                                                 });
                          .self$.exclusMat1 <- rep(mods1$lengths(use.names = FALSE),times=mods2$length()) - .self$.interMat;
                          .self$.exclusMat2 <- rep(mods2$lengths(use.names = FALSE),each=mods1$length()) - .self$.interMat;
                          .self$.unionMat <- .self$.exclusMat1 + .self$.exclusMat2 + .self$.interMat;
                        },
                        
                        .calcFisherExact=function()
                        {
                          mods1 <- .self$moduleSet1;
                          mods2 <- .self$moduleSet2;
                          .self$.fePMat<- matrix(data=0,nrow=mods1$length(),ncol=mods2$length());
                          .self$.feEstMat<- matrix(data=0,nrow=mods1$length(),ncol=mods2$length());
                          
                          sapply(X=1:mods1$length(), 
                                 FUN = function(i)
                                 { sapply(X=1:mods2$length(),
                                          FUN = function(j)
                                          {
                                            ft<-fisher.test(matrix(data=c(.self$.interMat[i,j],
                                                                          .self$.exclusMat1[i,j],
                                                                          .self$.exclusMat2[i,j], 
                                                                          .self$.totalElements-.self$.unionMat[i,j]),nrow=2,ncol=2));
                                            .self$.fePMat[i,j] <- ft$p.value;
                                            .self$.feEstMat[i,j] <- ft$estimate;
                                          })
                                 });
                        },
                        
                        .calcFDRAdjustment=function()
                        {
                          
                        },
                        
                        writeToFile(file, onlyTopTriangle=FALSE)
                        {
                          write.table();
                        }
                      )
);

testSMSet1 <-StarryModuleSet$new();
testSMSet1$initialize.byFiles("tests/GSE48865StarryMap.txt", "tests/GSE48865_geneOrder.txt",named=TRUE,byColumn=TRUE);
testSMSet2 <-StarryModuleSet$new();
testSMSet2$initialize.byFiles("tests/GSE57872StarryMap.txt", "tests/geneOrder.txt",named=TRUE,byColumn=TRUE);

testStarryMap<- StarryMap$new(moduleSet1=testSMSet1, moduleSet2=testSMSet2);
testInterSM <- testStarryMap$copy();
testInterSM$.calcPairwiseProfile();
testInterSM$.calcFisherExact();
