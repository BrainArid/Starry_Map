library(shiny);
library(shinydashboard)
source("moduleOverlap.R")
source("multiplot.R")
library(ggplot2)
library("svDialogs")
library(igraph);

header <- dashboardHeader(title = "Starry Map")

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Input",tabName="inputTabView",icon=icon("dashboard")),
    menuItem("Vizualize",tabName="visualTabView",icon=icon("th"))
  ));
inputTab<-tabItem(tabName="inputTabView",
                  fluidRow(box(width=2,checkboxInput("useExample", label = "Use Example?", value = TRUE)),
                  box(width=2,checkboxInput("fileU_IsAll","Use same U for all S files?",TRUE), uiOutput("fileUAllInput")),
                  box(width=2,radioButtons("fileS_Col_All","S files organized by:",choices = c("Columns"="1","Rows"="2","Different per file"="3")))),
                  fluidRow(box(width=6,numericInput("NumberModSets",label="Number of Module Sets",value=2,step = 1))),
                  fluidRow(box(width=6,uiOutput("ModuleSetInputTabs"),background = "navy"))
);

visualTab<-tabItem(tabName="visualTabView",
                   div(style = 'overflow-y: scroll', 
                    box(title = "Starry Map",plotOutput("plot"), downloadButton('saveStarry', 'Download')),
                    box(title = "Set-wise Mutual Information",plotOutput("MIPlot"), downloadButton('saveMI', 'Download')),
                    box(title = "Stability",plotOutput("stabPlot")),
                    box(title = "Constellation",plotOutput("constellation"))
                    )
                   );

body <- dashboardBody(tabItems(inputTab,visualTab));

ui<-dashboardPage(header, sidebar, body);
