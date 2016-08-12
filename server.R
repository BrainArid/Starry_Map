
args <- list();

#initialize arguments if 
initializeBooleanArg <- function(arg, default){
  if(is.null(arg))
  {
    arg <- default;
  } else if(!is.logical(arg))
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

args$dir1 <- initializeStringArg(arg=args$dir1, default="Test Sets\\");
args$dir2 <- initializeStringArg(arg=args$dir2, default="Test Sets\\");
args$clustsFile1 <- initializeStringArg(arg=args$clustsFile1, default="GSE57872StarryMap.txt");
args$clusts1Columned <- initializeBooleanArg(arg=args$clusts1Columned, default=TRUE);
args$clusts1Names <- initializeBooleanArg(arg=args$clusts1Names, default=TRUE);
args$clustsFile2 <- initializeStringArg(arg=args$clustsFile2, default="GSE48865StarryMap.txt");
args$clusts2Columned <- initializeBooleanArg(arg=args$clusts2Columned, default=TRUE);
args$clusts2Names <- initializeBooleanArg(arg=args$clusts2Names, default=TRUE);
args$outDir <- initializeStringArg(arg=args$outDir, default="Figures");
args$threshold <- initializeFloatArg(arg=args$threshold, default=0.001);
args$clustHeatMap <- initializeBooleanArg(arg=args$clustHeatMap, default=FALSE);
args$fisherExact <- initializeBooleanArg(arg=args$fisherExact, default=TRUE);
args$permutationTest <- initializeBooleanArg(arg=args$permutationTest, default=FALSE);
args$histogram <- initializeBooleanArg(arg=args$histogram, default=TRUE);
args$geneUniverseFile1 <- initializeStringArg(arg=args$geneUniverseFile1, default="Test Sets\\geneOrder.txt");
args$geneUniverseFile2 <- initializeStringArg(arg=args$geneUniverseFile2, default="Test Sets\\GSE48865_geneOrder.txt");

NumberModSets <- 3;
filesU <- reactiveValues();
filesS <- reactiveValues();
fileS_IsCol <- reactiveValues();
geneUniverses <- reactiveValues();
starryMaps <- list();

server <- function(input, output, session){
  
  #NumberModSets<<-reactive({
    #if(!(is.null(isolate(input$inputTabs))) && isolate(input$inputTabs) == "+")
    #{
    #   isolate(NumberModSets()) + 1;
    #}
    #else
    #{
    #  isolate(NumberModSets());
    #}
  #});
  output$fileUAllInput <- renderUI({
    if(input$fileU_IsAll)
      {fileInput("fileUAll", label = h5(paste0("Element Universe (U)")));}
  })
  
  output$fileS_IsColAllInput <- renderUI({
    if(input$fileS_Col_IsAll)
    {checkboxInput("fileS_IsColAll", label = h5(paste0("Columned?")),value = TRUE);}
  })
  
  
  output$ModuleSetInputTabs <- renderUI({
    
    #if(!(is.null(input$inputTabs)) && input$inputTabs == "+")
    # NumberModSets <<- NumberModSets + 1;
    
    Tabs <- lapply(1:NumberModSets, 
                   function(x){
                     uElement<-NULL;
                     colElement<-NULL;
                     if(!input$fileU_IsAll)
                     {
                       uElement<-fileInput(paste0("fileU",x), label = h5(paste0("Element Universe ",x," (U",x,")")));
                     }
                     if(input$fileS_Col_All=="3")
                     {
                       colElement<-checkboxInput(paste0("fileS",x,"_IsCol"),"Columned",TRUE);
                     }
                     tabPanel(title = paste0("Module Set",x),
                              value = x,
                              fileInput(paste0("fileS",x), label = h5(paste0("Module Set ",x," (S",x,")"))),
                              uElement,
                              colElement); 
                   });
    #if(!is.null(input$inputTabs) && input$inputTabs == "+")
    #  updateTabsetPanel(session,inputId = "inputTabs",selected = paste0("Module Set",NumberModSets));
    
    Tabs[[NumberModSets+1]] <- tabPanel(title = "+", value="+");
    do.call(what=tabsetPanel, c(Tabs,id="inputTabs"));
  });
  
#   Is <- reactive({
#     l<-as.list(1:NumberModSets);
#     names(l)<-lapply(l,as.character);
#     return(l);
#   })
  for(I in 1:NumberModSets)
  {
    local({
      i<-I;
      filesS[[as.character(i)]]<-reactive({
        if(input$useExample)
        {
          if(i=="1"){return("Test Sets/GSE48865StarryMap.txt")}
          else{return("Test Sets/GSE57872StarryMap.txt")}
        }
        else
        {
          return(input[[paste0("fileS",i)]]$datapath);
        }
      })
      
      filesU[[as.character(i)]]<-reactive({
        if(input$useExample)
        {
          if(i=="1"){return("Test Sets/GSE48865_geneOrder.txt")}
          else{return("Test Sets/geneOrder.txt")}
        }
        else if(input$fileU_IsAll)
        {
          return(input[["fileUAll"]]$datapath);
        }
        else
        {
          return(input[[paste0("fileU",i)]]$datapath);
        }
      })
      
      fileS_IsCol[[as.character(i)]]<-reactive({
        if(input$useExample)
        {
          if(i=="1"){return(TRUE)}
          else{return(TRUE)}
        }
        else if(!(input$fileS_Col_All=="3"))
        {
          if(input$fileS_Col_All=="1"){return(TRUE);}
          else(return(FALSE));
        }
        else
        {
          return(input[[paste0("fileS",i,"_IsCol")]]);
        }
      })
    })
  }
   Is<-as.list(1:NumberModSets);
   names(Is)<-lapply(Is,as.character)
#   filesS <- lapply(Is, function(I){
#     i<-as.character(I);
#     reactive({
#       if(input$useExample)
#       {
#         if(i=="1"){return("Test Sets/GSE48865StarryMap.txt")}
#         else{return("Test Sets/GSE57872StarryMap.txt")}
#       }
#       else
#       {
#         return(input[[paste0("fileS",i)]]$datapath);
#       }
#     })
#   })
#   #convert to reactiveValues object
#   filesS <- do.call(reactiveValues, filesS);
  
  #filesU[["1"]] <- reactive({"Test Sets/GSE48865_geneOrder.txt"});
  #filesS[["1"]] <- reactive({"Test Sets/GSE48865StarryMap.txt"});
  #filesU[["2"]] <- reactive({"Test Sets/geneOrder.txt"});
  #filesS[["2"]] <- reactive({"Test Sets/GSE57872StarryMap.txt"});
  #fileS_IsCol[["1"]] <- reactive({TRUE});
  #fileS_IsCol[["2"]] <- reactive({TRUE});
  
  #Is <- 1:NumberModSets
  geneUniverses <- lapply(Is, function(I){
    i<-as.character(I);
    return(reactive({
      g<-read.table(file=filesU[[i]](),skip=0,header=FALSE);
      return(apply(g,1,as.character));
    }))
  })
  
  #Is <- 1:NumberModSets;
  starryMaps <- lapply(Is, function(I){
    Js <- I:NumberModSets;
    return(lapply(Js,function(J){
      i<-force(as.character(I));
      j<-force(as.character(J));
      gU1<-geneUniverses[[I]];
      gU2<-geneUniverses[[J]];
      #k<-J+(I-1)*NumberModSets-(I*(I-1)/2);
      return(reactive({starryMap(filesS[[i]](),fileS_IsCol[[i]](), 
                                 args$clusts1Names, filesS[[j]](), fileS_IsCol[[j]](), 
                                 args$clusts1Names, args$outDir, args$threshold, args$clustHeatMap, 
                                 args$fisherExact, args$permutationTest, args$histogram, 
                                 geneUniverse1=gU1(),
                                 geneUniverse2=gU2());#need to add sorts 
      }))
    }))
  })
  
  output$MI <-renderText({
    return("Mutual Information: Here");
  });
  
  pStarry<-NULL;
  output$plot <-renderPlot({
    N <- NumberModSets;
    
    ps<-matrix(data = list(), nrow = N, ncol = N);
    starryTheme<-theme(axis.text=element_blank(),axis.text.x=element_text(angle=-90,hjust=0,vjust=0.5), axis.title=element_blank());

    #starryTheme<-theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=0.5), axis.title=element_blank());
    maxSize<-1;
    
    for(j in 1:N)
      for(i in 1:N)
      {
        inx<-(j-1)*N + i;
        if((N-i+1)<(j))#y=-1x+NumberModSets
          next;
        
        k<-i+(j-1)*N-(j*(j-1)/2);
        colMat<-starryMaps[[i]][[j]]()$colMat;
        
        #reorder
        l1<-sapply(starryMaps[[i]][[1]]()$order1,function(x){which(starryMaps[[i]][[j]]()$order1==x)});
        l2<-sapply(starryMaps[[j-1 + i-1 %% N + 1]][[1]]()$order1,function(x){which(starryMaps[[i]][[j]]()$order2==x)});
        colMat$X1<-factor(colMat$X1,levels=levels(colMat$X1)[l1])
        colMat$X2<-factor(colMat$X2,levels=levels(colMat$X2)[l2]);
        
        p<-ggplot(colMat, aes(y=X2, x=X1)) +
          geom_tile(aes(y=X2,x=X1,fill=hsv(r,b,g),width=maxSize,height=maxSize), colour = "black") +
          scale_fill_identity() +
          scale_x_discrete(labels=function(x){return(strtrim(x, 15))}) +
          scale_y_discrete(labels=function(x){return(strtrim(x, 15))}) +
          starryTheme;
        
        if(args$fisherExact)
          p <- p + geom_tile(aes(y=X2,x=X1,fill = rgb(red = 0,green=feEst,blue=0),width=feSize*maxSize,height=feSize*maxSize), colour = "black");
        #i j > prod> draw order > ps i j
        #1 1 > 1 1 > 7 > 1 3
        #1 2 > 1 2 > 8 > 2 3
        #1 3 > 1 3 > 9 > 3 3
        #2 1 > 2 2 > 5 > 2 2 
        #2 2 > 2 3 > 6 > 3 2
        #2 3 > null > 4 > 1 2
        #3 1 > 3 3 > 3 > 3 1
        #3 2 > null > 1 > 1 1
        #3 3 > null > 2 > 2 1
        ps[[j-1 + i-1 %% N + 1, N - i + 1]] <- p;
      }
    
    pStarry<-multiplot(plotlist=ps, cols=N);
    png(filename = "Starryplot.png",width=300,height=300,units = "px",res=300)
    print(pStarry);
    dev.off();
    print(pStarry); 
  });  
  
  starryMI<-reactive({
    N <- NumberModSets;
    m<-matrix(0,nrow = N,ncol=N);
    for(j in 1:N)
      for(i in 1:N)
      {
        if((N-i+1)<(j))#y=-1x+NumberModSets
          next;
        m[N - i + 1, N-(j-1 + i-1 %% N + 1)+1]<-starryMaps[[i]][[j]]()$starryMI
      }
    return(m);
  });
  
  pMI<-NULL;
  output$MIPlot <- renderPlot({
 
    starryMI.m<-melt(starryMI());
    starryMI.m<-starryMI.m[which(starryMI.m$X1>=starryMI.m$X2),]#filter top-triangle
    names(starryMI.m)<-c("X1","X2","MI")
    p<-ggplot(starryMI.m, aes(x=X1,y=X2)) + geom_tile(aes(x=X1,y=X2,fill=MI)) +
      geom_text(aes(label=round(as.numeric(MI),3)),colour="white",size=8) +
      theme(axis.text=element_blank(), axis.title=element_blank());
    pMI<-p;
    ggplot2::ggsave("MIplot.png",pMI,width=3,height=2.5);
    print(pMI);
  });
  
  output$stabPlot <-renderPlot({
    N <- NumberModSets;
    
    stabVec<-vector(length=N-1);
    for(i in 1:(N-1))
    {
      stabVec[i]<-starryMI()[i+1,i];
    }
    stabVec.m<-data.frame(X=1:(N-1),MI=stabVec);
    p<-ggplot(stabVec.m, aes(x=X,y=MI)) + geom_line() +
      theme(axis.text=element_blank(), axis.title=element_blank());
    print(p);
  })

  output$constellation<-renderPlot({
    s<-starryMI();
    s[upper.tri(s)] = t(s)[upper.tri(s)];
    g<-graph.adjacency(s,mode="undirected",weighted=TRUE);
    p<-plot.igraph(g,layout=layout.fruchterman.reingold(g, weights=E(g)$weight));
    print(p);
  })
  
  #to download output 
  output$saveStarry <- downloadHandler(
    filename = function() {svDialogs::dlgSave(title="Save To");},
    content = function(file) {
      ggsave(pStarry, file)
    }
  )
  
  output$saveMI <- downloadHandler(
    filename = function() { file.choose(new = TRUE) },
    content = function(file) {
      ggsave(pMI, file)
    }
  )
} 