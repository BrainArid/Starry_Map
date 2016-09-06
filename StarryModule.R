#StarryModule Class - contains a single set of unique elements.
StarryModule <- setRefClass(Class = "StarryModule",
                            
                            fields = list(
                              
                              #public
                              name="character",
                              elements="character"
                              
                              #private
                            ),
                            
                            methods = list(
                              show=function(){
                                cat("An object of class ", class(.self), "\n",
                                    "\tName: ", .self$name, "\n",
                                    "\tContains ", base::length(.self$elements),
                                    " elements: \n\t\t", utils::head(.self$elements));
                                if(base::length(.self$elements) > 6) 
                                  cat("...")
                                invisible(NULL);
                              },
                              
                              length=function(){
                                return(base::length(.self$elements))
                              },
                              head=function(){
                                return(utils::head(.self$elements))
                              },
                              tail=function(){
                                return(utils::tail(.self$elements))
                              },
                              #contains side effects
                              intersect.elements=function(otherElements){
                                .self$elements <- base::intersect(.self$elements,otherElements)
                              },
                              #contains side effects
                              intersect.starryModule=function(otherStarryModule)
                              {
                                intersect.elements(otherStarryModule$elements);
                              }
                            )
);

testSM1 <- StarryModule$new(name="Cell Cycle", elements=c("A","B","C","D","E","F","G"));
testSM2 <- StarryModule$new(name="Metabolism", elements=c("A","B","G","H","I"));
