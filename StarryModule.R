#StarryModule Class - contains a single set of unique elements.
StarryModule <- setClass(Class = 
                           "StarryModule",
                         slots = c(
                           name="character",
                           elements="character"
                         )
);
#StarryModule overridden methods
setMethod("show", "StarryModule",
          function(object){
            cat("An object of class ", class(object), "\n",
                "\tName: ", object@name, "\n",
                "\tContains ", length(object@elements),
                " elements: \n\t\t", head(object@elements));
            if(length(object@elements) > 6) 
              cat("...")
            invisible(NULL);
          })
setMethod("length","StarryModule",
          function(x){
            return(length(x@elements))
          })
setMethod("head","StarryModule",
          function(x){
            return(head(x@elements))
          })
setMethod("tail","StarryModule",
          function(x){
            return(tail(x@elements))
          })

#testSM <- StarryModule(name="Cell Cycle", elements=c("A","B","C","D","E","F","G"));