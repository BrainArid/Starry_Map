#StarryModuleSet Class - a set of StarryModules and the superset of elements comprising the collective element universe.
StarryModuleSet <- setClass(Class = 
                              "StarryModuleSet",
                            slots = c(
                              modules="list",
                              universe="StarryModule"
                            )
);

setMethod("show", "StarryModuleSet",
          definition = function(object){
            cat("An object of class ", class(object), "\n",
                "\tUniverse of length ", length(object@universe),": ", head(object@universe));
            if(length(object@universe)>6) cat("...");
            cat("\n\tContaining ",length(object@modules), " modules:\n");
            show(object@modules);
            invisible(NULL);
          })

#testSM_1 <- StarryModule(name="Cell Cycle", elements=c("A","B","C","D","E","F","G"));
#testSM_2 <- StarryModule(name="Metabolism", elements=c("F","G","H","I","J"));
#testU <- StarryModule(name="Universe", elements=union(testSM_1@elements, testSM_2@elements));
#testSMSet <- StarryModuleSet(modules=c(testSM_1,testSM_2),universe=testU);
