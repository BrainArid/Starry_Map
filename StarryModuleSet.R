source("StarryModule.R");

#StarryModuleSet Class - a set of StarryModules and the superset of elements comprising the collective element universe.
StarryModuleSet <- setRefClass(Class ="StarryModuleSet",
                            fields = list(
                              #public
                              modules="list",
                              universe="StarryModule"
                              
                              #private
                              #.lengths="integer"
                            ),
                            
                            methods = list(
                              
                              initialize.byFiles=function(moduleSetFile,universeFile,named=TRUE,byColumn=TRUE)
                              {
                                mods <- read.table(file=moduleSetFile, sep="\t", stringsAsFactors=FALSE, fill=TRUE);
                                
                                #rotate if necessary
                                if(byColumn)
                                  mods<-t(mods)
                                #remove first column if module names found
                                if(named)
                                {
                                  row.names(mods)<-factor(mods[,1])
                                  mods<- mods[,-1] 
                                }else
                                {
                                  n<- paste0("module", 1:dim(mods)[1])
                                  row.names(mods)<- factor(n,levels=n)
                                  rm(n);
                                }
                                #remove any columns containing NAs
                                mods <- mods[,colSums(is.na(mods))==0];
                                
                                .self$modules <- apply(cbind(as.character(row.names(mods)),mods), 1,
                                                         function(x, name){
                                                           name <- x[1];
                                                           elements <- x[-1];
                                                           elements <- elements[elements!=""];
                                                           StarryModule$new(name=name,elements=elements);
                                                         })
                                
                                #read in universe
                                univ <- read.table(file=universeFile, skip=0,header=FALSE,stringsAsFactors = FALSE)[,1];
                                .self$universe <- StarryModule$new(elements=univ);
                              },
                              
                              show=function()
                              {
                                cat("An object of class ", class(.self), "\n",
                                    "\tUniverse of length ", .self$universe$length(),": ", .self$universe$head());
                                if(.self$universe$length()>6) cat("...");
                                cat("\n\tContaining ",base::length(.self$modules), " modules:\n");
                                lapply(.self$modules,function(module){module$show(); cat("\n");})
                                invisible(NULL);
                              },
                              
                              length=function()
                              {
                                return(base::length(.self$modules));
                              },
                              
                              lengths=function(use.names=TRUE)
                              {
                                lengs <- lapply(.self$modules,function(module){module$length();});
                                if(use.names){
                                  return(lengs);
                                } else {
                                  return(as.vector(unlist(lengs)));
                                }
                              },
                              
                              #contains side effects
                              intersect.elements=function(otherElements)
                              {
                                #intersect all modules
                                lapply(.self$modules,function(module){module$intersect.elements(otherElements)});
                                #intersect universe
                                .self$universe$intersect.elements(otherElements);
                              }, 
                              
                              #contains side effects
                              intersect.starryModule=function(otherStarryModule)
                              {
                                #intersect all modules
                                lapply(.self$modules,function(module){module$intersect.starryModule(otherStarryModule)});
                                #intersect universe
                                .self$universe$intersect.starryModule(otherStarryModule);
                              }
                            )
);


testSM_1 <- StarryModule$new(name="Cell Cycle", elements=c("A","B","C","D","E","F","G"));
testSM_2 <- StarryModule$new(name="Metabolism", elements=c("F","G","H","I","J"));
testU <- StarryModule$new(name="Universe", elements=union(testSM_1$elements, testSM_2$elements));
testSMSet <- StarryModuleSet$new(modules=c(testSM_1,testSM_2),universe=testU);

testSMSet <-StarryModuleSet$new();
testSMSet$initialize.byFiles("tests/GSE48865StarryMap.txt", "tests/GSE48865_geneOrder.txt",named=TRUE,byColumn=TRUE);
