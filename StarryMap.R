StarryMap <- setClass(Class = "StarryMap", 
                      
                      slots = c(moduleSet1="StarryModuleSet",
                                moduleSet2="StarryModuleSet"));
setMethod("show", "StarryMap",
          definition = function(object){
            cat("An object of class ", class(object), "\n");
            show(object@moduleSet1);
            show(object@moduleSet2);
            invisible(NULL);
          })
   
   
#    "clusts1Name","clustsclustsFile1",clusts1Columned=FALSE, clusts1Names=FALSE, 
#                                         clustsFile2, clusts2Columned=FALSE, clusts2Names=FALSE, 
#                                         outDir, threshold, clustHeatMap=FALSE, fisherExact=FALSE, 
#                                         histogram=FALSE, geneUniverse1, geneUniverse2, sort1=TRUE, sort2=TRUE)
#            )

testSM_1 <- StarryModule(name="Cell Cycle", elements=c("A","B","C","D","E","F","G"));
testSM_2 <- StarryModule(name="Metabolism", elements=c("F","G","H","I","J"));
testU <- StarryModule(name="Universe", elements=union(testSM_1@elements, testSM_2@elements));
testSMSet <- StarryModuleSet(modules=c(testSM_1,testSM_2),universe=testU);
testStarryMap<- StarryMap(moduleSet1=testSMSet, moduleSet2=testSMSet);
