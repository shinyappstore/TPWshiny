# TPWshiny: an interactive RShiny app to explore the therapeutic response of NCI-60 cell lines
# version: 1.1 (release date: July 2023)
#
# PLEASE CITE:
# Title: TPWshiny: an interactive Shiny app to explore the therapeutic response of NCI-60 cell lines
# For authors and journal information, please refer to the BRP website: https://brb.nci.nih.gov/TPWshiny/


####
## DISCLAIMER:
####
# All users must agree to the following conditions:
# - Please cite the published manuscript in all studies using TPWshiny
# - The package is offered without support and the user will not hold the National Cancer Institute liable for any damages resulting from the use of TPWshiny.
# - This resource is intended for purely research purposes. It should not be used for emergencies or medical or professional advice.
# - THE SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT ARE DISCLAIMED. 
#  IN NO EVENT SHALL THE NATIONAL CANCER INSTITUTE(THE PROVIDER), THE NATIONAL INSTITUTES OF HEALTH, THE U.S. GOVERNMENT OR THE INDIVIDUAL DEVELOPERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
#  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
#  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  
####
## Changes:
## v1.1 - removed unsupported MRAN link for missing packages. Missing packages managed using the groundhog package.
## v1.1 - modified "Single Gene" tab to display "Baseline (no drug)" in the Parameters/Concentration drop-down
## 
####
## COPYRIGHT and LICENSE:
####
## TPWshiny is provided under Open Source Software Distribution Agreement. 
## Please refer to the BRP website to see the full terms of the License https://brb.nci.nih.gov/TPWshiny/ 
## For technical issues, please contact TPWshiny Support at ncitpwsupport@mail.nih.gov
####


#############
# 1 IMPORTS #
#############

# required packages
packages <- c("DT",
              "shiny",
              "ggplot2",
              "data.table",
              "shinyWidgets",
              "heatmaply",
              "plotly",
              "gridExtra",
              "UpSetR",
              "gplots",
              "dendextend",
              "devtools",
              "venn")


# find missing packages
missingPackages <- packages[!(packages %in% installed.packages()[,"Package"])]

# install missing packages
if(length(missingPackages)) {
  library("groundhog")
  groundhog.library(missingPackages, date = "2023-06-01")
}
 

# import packages
lapply(packages, library, character.only = TRUE)


server_name = "tpwb.nci.nih.gov"

ping <- function(x){
	mtry <- try(read.table(paste0("https://",x,"/GeneExpressionNCI60/rshiny/cellLines.txt")), 
				silent = TRUE)
	if (class(mtry) != "try-error") {
	  TRUE
	} else {
	  FALSE
	}
}


check_web_connection <- ping(server_name)


if(check_web_connection==FALSE) {
  ui<-shinyUI(fluidPage(
    
    # Application title
    titlePanel("TPWshiny"),
    
    # Sidebar with a slider input for number of observations
    sidebarLayout(
      sidebarPanel(
        br(),
        hr(style ="border: 1px solid red;"),
        p("ERROR CONNECTING TO THE SERVER!", style="padding: 15px; margin: 0px;border-style:solid; border-color: #e9e9e9; color:red; font-weight: bold;"),
        hr( style ="border: 1px solid red;"),
        br()
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        br(),
        hr(),
        p("Problem loading data from server. Please check your internet connection." ),
        br(),
        p("If the problem persists, contact the NCI-TPW support team at ncitpwsupport@mail.nih.gov"),
        hr(),
        br()
      )
    )
  ))
  
  server <- shinyServer(function(input, output, session) {
    # IMPORTANT!
    # this is needed to terminate the R process when the
    # shiny app session ends. Otherwise, you end up with a zombie process
    session$onSessionEnded(function() {
      stopApp()
    })
  })
 
} else {
  
###############
# 2 CONSTANTS #
###############
if(check_web_connection==TRUE) {
  source_url(paste0("https://",server_name,"/Scripts/heatmap3.R"))
zz <- file(".TPWshiny_logfile.txt", open = "wt")
sink(zz, type = "message")

list.as.matrix <- function(x, byrow=FALSE, filler=NA){
  maxlen <- max(sapply(x,length))
  xm <- sapply(x, 
               function(xs){fillen <- maxlen-length(xs)
               if (fillen>0) {c(xs,rep(filler,fillen))} else xs
               })
  if (byrow)return(t(xm)) else return(xm)
}


drug_names <- c("azacytidine", "bortezomib", "cisplatin", "dasatinib", "doxorubicin", "erlotinib", "geldanamycin", 
                "gemcitabine","lapatinib", "paclitaxel","sirolimus", "sorafenib", "sunitinib","topotecan", "vorinostat")


source_directory=paste0("https://",server_name,"/GeneExpressionNCI60/rshiny/venn/");
correlation_24_low <- read.table( paste0(source_directory,"GI50_top100genes_24_low_correlation.txt"), header=TRUE, as.is = TRUE)
correlation_24_high <- read.table(paste0(source_directory,"GI50_top100genes_24_high_correlation.txt"), header=TRUE, as.is = TRUE)
correlation_6_low <- read.table(paste0(source_directory,"GI50_top100genes_6_low_correlation.txt"), header=TRUE, as.is = TRUE)
correlation_6_high <- read.table(paste0(source_directory,"GI50_top100genes_6_high_correlation.txt"), header=TRUE, as.is = TRUE)
correlation_2_low <- read.table(paste0(source_directory,"GI50_top100genes_2_low_correlation.txt"), header=TRUE, as.is = TRUE)
correlation_2_high <- read.table(paste0(source_directory,"GI50_top100genes_2_high_correlation.txt"), header=TRUE, as.is = TRUE)

correlation = c()
correlation[["2 h"]][["low"]] = correlation_2_low
correlation[["2 h"]][["high"]] = correlation_2_high
correlation[["6 h"]][["low"]] = correlation_6_low
correlation[["6 h"]][["high"]] = correlation_6_high
correlation[["24 h"]][["low"]] = correlation_24_low
correlation[["24 h"]][["high"]] = correlation_24_high


time_points <- c("2 h", "6 h", "24 h")
conc_levels <- c("low", "high")

negative_correlation = c()
positive_correlation = c()
correlation_all = c()

for(which_drug in drug_names){
  for(timepoint in time_points){
    for(conclevel in conc_levels){
      correlation_current = correlation[[timepoint]][[conclevel]]
      index_drug = which(colnames(correlation_current)== which_drug)
      negative_correlation[[timepoint]][[conclevel]][[which_drug]] = correlation_current[correlation_current[,index_drug+1]<0,index_drug]
      positive_correlation[[timepoint]][[conclevel]][[which_drug]] = correlation_current[correlation_current[,index_drug+1]>0,index_drug]
      correlation_all[[timepoint]][[conclevel]][[which_drug]] = correlation_current[,index_drug]
    }
  }
}


# font for graphs
fonts <- list(sans = "Helvetica")

source_directory_main=paste0("https://",server_name,"/GeneExpressionNCI60/rshiny/");

# gene names
geneNames <- sort(scan(paste0(source_directory_main,"genes.txt"), character()))
geneNames_reduced <- c("A1CF", "RPL6", "CTDNEP1","RPL10A", "CBX3", "RPL17", "PSMB2", "KHDRBS1", "ERH", "SRSF9", "ABCF1")

# pathway names
pathwayNames <- sort(scan(paste0(source_directory_main,"pathways.txt"), character()))

# receptor group names
recepNames <- scan(paste0(source_directory_main,"receptors.txt"), character(), sep = "\n")

# transcriptional factor group names
transfactorNames <- sort(scan(paste0(source_directory_main,"transfactors.txt"), character(), sep = "\n"))

# cell line names
cellLineNames <- sort(scan(paste0(source_directory_main,"cellLines.txt"), character(), sep = "\n"))

# drug names
drugNames <- c("Azacytidine",
               "Bortezomib",
               "Cisplatin",
               "Dasatinib",
               "Doxorubicin",
               "Erlotinib",
               "Geldanamycin",
               "Gemcitabine",
               "Lapatinib",
               "Paclitaxel",
               "Sirolimus",
               "Sorafenib",
               "Sunitinib",
               "Topotecan",
               "Vorinostat")

# tissue names
tissueNames <- c("Breast",
                 "CNS",
                 "Colon",
                 "Leukemia",
                 "Lung",
                 "Melanoma",
                 "Ovarian",
                 "Prostate",
                 "Renal")

tissueColors <- c("BREAST" = "#f0251d",
                    "CNS" = "#b4e752",
                    "COLON" = "#f9a416",
                    "LEUKEMIA" = "#37df35",
                    "LUNG" = "#ea0ea6",
                    "MELANOMA"= "#5aedbc",
                    "OVARIAN"= "#1110c5",
                    "PROSTATE" = "#a409e9",
                    "RENAL" = "#3fbbf3")

# stub for GeneCards search
geneStub <- "https://www.genecards.org/cgi-bin/carddisp.pl?gene="

# color schemes
colorschemes <- c("Blue-White-Red", 
                  "Blue-Black-Red",
                  "Blue-White-Orange",
                  "Blue-Black-Orange",
                  "Green-White-Red",
                  "Green-Black-Red",
                  "Cyan-Navy")

TF_with_no_heatmap <-
  c("There is only one target gene DISC1 for transcription factor ATF5. No heatmap will be generated.",
    "There is only one target gene FOS for transcription factor ELK3. No heatmap will be generated.",
    "There is only one target gene FOS for transcription factor ELK4. No heatmap will be generated.",
    "There is only one target gene CCL5 for transcription factor FOSB. No heatmap will be generated.",
    "There is only one target gene TH for transcription factor HOXA4. No heatmap will be generated.",
    "There is only one target gene TGM1 for transcription factor HOXA7. No heatmap will be generated.",
    "No target gene has been found for transcription factor HOXB1. No heatmap will be generated.",
    "No target gene has been found for transcription factor HOXB2. No heatmap will be generated.",
    "There is only one target gene NKX2-1 for transcription factor HOXB3. No heatmap will be generated.",
    "No target gene has been found for transcription factor HOXB8. No heatmap will be generated.",
    "No target gene has been found for transcription factor HOXB9. No heatmap will be generated.",
    "There is only one target gene LCT for transcription factor HOXC11. No heatmap will be generated.",
    "There is only one target gene HOXD9 for transcription factor HOXD10. No heatmap will be generated.",
    "There is only one target gene ADRA1B for transcription factor NFIX. No heatmap will be generated.",
    "There is only one target gene INS for transcription factor PAX4. No heatmap will be generated.",
    "There is only one target gene KDM5B for transcription factor PAX9. No heatmap will be generated.",
    "No target gene has been found for transcription factor POU2F3. No heatmap will be generated.",
    "There is only one target gene CDK6 for transcription factor SMAD6. No heatmap will be generated.",
    "There is only one target gene DRD1 for transcription factor TFAP2B. No heatmap will be generated.",
    "There is only one target gene ALDH1A1 for transcription factor TLX1. No heatmap will be generated.")


names(TF_with_no_heatmap) <-
  c("ATF5",
    "ELK3",
    "ELK4",
    "FOSB",
    "HOXA4",
    "HOXA7",
    "HOXB1",
    "HOXB2",
    "HOXB3",
    "HOXB8",
    "HOXB9",
    "HOXC11",
    "HOXD10",
    "NFIX",
    "PAX4",
    "PAX9",
    "POU2F3",
    "SMAD6",
    "TFAP2B",
    "TLX1")

}
  
####################
# 3 USER INTERFACE #
####################


# navigation bar
ui <- navbarPage(
  id = "mainNavPage",
  #title
  tags$head(
    tags$style(
      HTML(".navbar-default .navbar-brand {margin: 15px;}
           .navbar > .container .navbar-brand, .navbar > .container-fluid .navbar-brand { margin-left: 0px;}")
      )
    ),
              
  title = "TPWshiny",
  selected =  p("Single Gene", style="padding: 15px; margin: 0px;border-style:solid; border-color: #e9e9e9;"),
  
  #############
  # 3.1 GENES #
  #############
  
  # tab title   
  tabPanel( p("Single Gene", style="padding: 15px; margin: 0px;border-style:solid; border-color: #e9e9e9;"),
          
     fluidRow(
       
       ##################
       # 3.1.1 LEFT BAR #
       ##################
        
       # begin left bar
       column(2, 
              
              ######################
              # 3.1.1.1 DATA INPUT #
              ######################
              
              wellPanel(
                
                # panel header
                h4("Dataset"),
                
                # input gene
                selectInput(inputId = "inGene", 
                            label =  HTML("Gene <br/>(type to trigger autofill)"), 
                            choices = geneNames_reduced),
                
            
                # link to GeneCards URL
                uiOutput("geneURL"),
                br(),
                # input drug
                selectInput(inputId = "inDrug", 
                            label = "Drug", 
                            
                            # does NOT include baseline - baseline is added as a conc measure option
                            choices = c(drugNames))),
              
              
              
              ######################
              # 3.1.1.2 PARAMETERS #
              ######################
              
              
              
              conditionalPanel(
                
                # if scatterplot, barplot, or non-baseline timeplot,
                # show the parameters panel
                condition = "!(input.gene_activePlot == 'Data')",
                
                wellPanel(
                  
                  # panel header
                  h4("Parameters"),
                  
                  
                  
                  # if scatterplot or barplot,
                  # show the tissue checkboxes
                  conditionalPanel(
                    condition = "!(input.gene_activePlot == 'Data')",
                    
                    # formatting for checkboxes
                    tags$head(
                      tags$style(
                        HTML(
                          ".checkbox-inline { 
                              margin-left: 0px;
                              margin-right: 10px;
                            }
                            .checkbox-inline+.checkbox-inline {
                              margin-left: 0px;
                              margin-right: 10px;
                            }"))),
                    
                    # tissue input
                    checkboxGroupInput(inputId = "inTissue",
                                       label = "Tissue",
                                       choices = tissueNames, 
                                       selected = tissueNames,
                                       inline = FALSE),
                    
                    # select all tissues
                    actionLink("selectall", "- Select All"),
                    br(),
                    # deselect all tissues
                    actionLink("deselectall", "- Clear selection")
                  ),
                  
                  br(),br(),
                  
                  # if there's a drug selected,
                  # display dosage level
                  conditionalPanel(
                    condition = "!(input.inDrug == 'baseline')",
                    
                    # input dosages
                    selectInput(inputId = "inDosage",
                                label = "Concentration",
                                choices = c("Low", "High", "Baseline (no drug)"))),
                  
                  
                  
                  # if scatterplot or barplot,
                  # display time
                  conditionalPanel(
                    condition = "input.gene_activePlot == 'Scatterplot' | 
                   input.gene_activePlot == 'Barplot'",
                    
                    # input time
                    selectInput(inputId = "inTime",
                                label = "Time",
                                choices = c("2 hrs", "6 hrs", "24 hrs"))),
                  
                  
                  
                  # if scatterplot,
                  # display resistance measure
                  conditionalPanel(
                    condition = "input.gene_activePlot == 'Scatterplot'",
                    
                    # input measure
                    selectInput(inputId = "inMeasure",
                                label = "Resistance Measure", 
                                choices = c("GI50" = "GI50",
                                            "Doubling Time" = "DoubleTime",
                                            "Multidrug Resistance" = "MultiRes"),
                                selected = "GI50"))))),
       
       
       
       ####################
       # 3.1.2 MAIN PANEL #
       ####################
       
       
       
       # begin main panel
       column(8,
              
              # formatting headers
              tags$head(
                tags$style(type='text/css', 
                           ".nav-tabs {font-size: 16px} ")),          
              tabsetPanel(
                id = "gene_activePlot",
                type = "tabs",
                
                
                
                # display scatterplot
                tabPanel("Scatterplot",
                         br(),
                         plotlyOutput("scatterplot", 
                                      width = "100%")),
                
                
                
                
                # display barplot
                tabPanel("Barplot",
                         br(),
                         plotlyOutput("bar", 
                                      width = "100%")),
                
                
                
                # display timeplot
                tabPanel("Timeplot",
                         br(),
                         plotlyOutput("timeplot", 
                                      width = "100%")),
                
                
                
                
                # display data table
                tabPanel("Data",
                         DT::dataTableOutput("table"),
                         style = "height:700px; overflow-y: scroll;overflow-x: scroll;"))),
       
       
       
       ###################
       # 3.1.3 RIGHT BAR #
       ###################
       
       
       # begin right bar
       column(2,
              
              
              
              ####################
              # 3.1.3.1 DOWNLOAD #
              ####################
              
              conditionalPanel(
                condition = "input.inTissue.length >= 1 || input.gene_activePlot == 'Data'",
              
              wellPanel(
                
                # panel header
                h4("Download"),
                
                
                # input download type
                selectInput(inputId = "gene_downloadType",
                            label = "Target",
                            choices = c("Full Dataset" = "full", 
                                        "Graphed Dataset" = "graphed")),
                
                # download button
                downloadButton("gene_download", 'Download'))),
              
              
              
              #######################
              # 3.1.3.2 CORRELATION #
              #######################
              
              
              
              # if scatterplot,
              # display correlations
              conditionalPanel(
                condition = "input.gene_activePlot == 'Scatterplot' && input.inTissue.length >= 1",
                
                
                wellPanel(
                  
                  
                  # panel header
                  h4("Correlation"),
                  
                  # display pearson's correlation
                  h5("Pearson's:"),
                  h5(textOutput("gene_pearsons")),
                  
                  # display spearman's correlation
                  h5("Spearman's"),
                  h5(textOutput("gene_spearmans")),
                  
                  #display kendall correlation
                  h5("Kendall's"),
                  h5(textOutput("gene_kendall"))),
                
                
                
                ######################
                # 3.1.3.3 REGRESSION #
                ######################
                
                
                
                wellPanel(
                  
                  # panel header
                  h4("Regression"),
                  
                  # input regression type
                  selectInput(inputId = "gene_regType",
                              label = "Regression Type",
                              choices = c("None" = "none",
                                          "Linear" = "linear",
                                          "Polynomial" = "poly",
                                          "Loess" = "loess")),
                  
                  # if a regresion method is selected,
                  # allow toggling of the confidence interval
                  conditionalPanel(
                    condition = "!(input.gene_regType == 'none')",
                    
                    # input confidence interval
                    sliderInput(inputId = "gene_confInt",
                                label = "Confidence Interval",
                                min = 0,
                                max = 0.99,
                                value = 0.95,
                                step = 0.01)))),
              
              
              
              ###################
              # 3.1.3.4 SORTING #
              ###################
              
              # if barplot,
              # enable sorting mechanisms
              conditionalPanel(
                condition = "input.gene_activePlot == 'Barplot' && input.inTissue.length >= 1",
                
                wellPanel(
                  
                  # panel header
                  h4("Sorting"),
                  
                   
                  selectInput(inputId = "gene_barSort",
                              label = "Sort By",
                              choices = c("None" = "None",
                                          "Gene Expression" = "Gene.Expression",
                                          "GI50" = "GI50",
                                          "Doubling Time" = "DoubleTime",
                                          "Multidrug Resistance" = "MultiRes")),
                  
                  conditionalPanel(
                    condition = "!(input.gene_barSort == 'Tissue.Origin')",
                    
                    selectInput(inputId = "gene_barOrder",
                                label = "Order",
                                choices = c("Ascending" = "up",
                                            "Descending" = "down"))),
                  
                  checkboxInput(inputId = "gene_bargroup",
                                label = "Group by Tissue",
                                value = TRUE)
                  )),
              
              
              
              #####################
              # 3.1.3.5 QUARTILES #
              #####################      
              
              
              # if timeplot,
              # show quartile option
              conditionalPanel(
                condition = "input.gene_activePlot == 'Timeplot' && input.inTissue.length >= 1",
                
                wellPanel(
                  
                  # panel header
                  h4("Quartiles"),
                  
                  # formatting for checkboxes
                  tags$head(
                    tags$style(
                      HTML(
                        ".checkbox-inline { 
                margin-left: 0px;
                margin-right: 10px;
              }
              .checkbox-inline+.checkbox-inline {
                margin-left: 0px;
                margin-right: 10px;
              }"))),
                  
                  # quartile measure
                  selectInput(inputId = "gene_quartMeasure",
                              label = "Segment By", 
                              choices = c("None" = "none",
                                          "GI50" = "GI50",
                                          "Doubling Time" = "DoubleTime",
                                          "Multidrug Resistance" = "MultiRes"),
                              selected = "none"),
                  
                  
                  # if a segmenting method is selected,
                  # show quartiles
                  conditionalPanel(
                    condition = "!(input.gene_quartMeasure == 'none')",
                    
                    # quartiles
                    checkboxGroupInput(inputId = "gene_quartiles",
                                       label = "Quartiles",
                                       choices = c("First" = "first",
                                                   "Second" = "second",
                                                   "Third" = "third",
                                                   "Fourth" = "fourth"), 
                                       selected = NULL,
                                       inline = TRUE))))))),


  
  ################
  # 3.2 HEATMAPS #
  ################
  
  
  
  # tab title
  tabPanel(p("Geneset Patterns", style="padding: 15px; margin: 0px;border-style:solid; border-color: #e9e9e9;"),
           
     fluidRow(
       
       
       
       ##################
       # 3.2.1 LEFT BAR #
       ##################
       
       
       
       # set sidebar
       column(2, 
              
              
              
              ######################
              # 3.2.1.1 DATA INPUT #
              ######################             
              
              
              
              wellPanel(
                
                # panel title
                h4("Dataset"),
                
                # input grouping option
                selectInput(inputId = "heat_groupType",
                            label = "Gene Set",
                            choices = c("Pathways" = "pathway",
                                        "Receptors" = "recep",
                                        "Transcription Factor" = "tf",
                                        "User-Defined" = "manual")),
                
                conditionalPanel(condition = "input.heat_groupType == 'tf'",
                selectInput(inputId = "heat_grouptfType",
                            label = "Parameter",
                            choices = c("Target genes of a specific TF" = "transfactor",
                                        "Drug effects on all TF" = "b",
                                        "Drug effects on specific TF" = "c")),
                # if selected transcription factor, show transcription factor options
                conditionalPanel(condition = "input.heat_grouptfType == 'transfactor' || input.heat_grouptfType == 'c'",
                
                # input transfactors
                selectInput(inputId = "heat_inTransFactor",
                label = "Transcription Factor",
                choices = transfactorNames))),
                
                
                # if selected pathway, show pathway options
                conditionalPanel(condition = "input.heat_groupType == 'pathway'",
                
                  # input pathway
                  selectInput(inputId = "heat_inPath", 
                              label = "Pathway", 
                              choices = pathwayNames)),
                
                # if selected receptors, show receptor options
                conditionalPanel(condition = "input.heat_groupType == 'recep'",
                                 
                   # input receptor
                   selectInput(inputId = "heat_inRecep", 
                               label = "Receptor", 
                               choices = recepNames)),
                
                
                conditionalPanel(condition = "input.heat_groupType == 'manual'",
                   
                   # input transfactors
                   selectizeInput(inputId = "heat_inGenes",
                                 label = "Genes (Paste Comma-separated values)",
                                 choices = geneNames_reduced,
                                 multiple = TRUE,
                                 selected = NULL,
                                 options = list(
                                   splitOn = I("(function() { return /[,; ]/; })()"),
                                   create = I("function(input, callback){
                                                return {
                                                value: input,
                                                label: input,
                                                text: input };}")))),
                               
                                 
                conditionalPanel(condition = "input.heat_groupType == 'manual' || input.heat_groupType == 'recep' || input.heat_groupType == 'pathway' || input.heat_grouptfType == 'transfactor' || input.heat_grouptfType == 'b'",
                # input drug
                selectInput(inputId = "heat_inDrug", 
                            label = "Drug", 
                            choices = drugNames)),
                
                
                # input dosages
                selectInput(inputId = "heat_inDose",
                            label = "Concentration",
                            choices = c("Low" = "L",
                                        "High" = "H")),
                
                # input time
                selectInput(inputId = "heat_inTime",
                            label = "Time",
                            choices = c("2 hrs" = "2hr",
                                        "6 hrs" = "6hr",
                                        "24 hrs" = "24hr")))),
             
             
             
       ####################
       # 3.2.2 MAIN PANEL #
       ####################
       
       
       # start main panel
       column(8,
              wellPanel(style = "background: white; border: 0px solid white; -webkit-box-shadow: inset 0 0px 0px rgba(0,0,0,0);box-shadow: inset 0 0px 0px rgba(0,0,0,0);",
                        uiOutput("genesNotFound"),
                        tags$br(),
                      # plot heatmap
                      plotlyOutput("heatmap",
                                    height = 750
                                 ))
              
        ),
       
       
       
       ###################
       # 3.2.3 RIGHT BAR #
       ###################
       
       
       
       # start right bar
       column(2,
              
              
              ####################
              # 3.2.3.1 PLOTTING #
              ####################
              
              
              wellPanel(
                
                # panel title
                h4("Plotting"),
                
                # select clusteirng option
                selectInput(inputId = "heat_cluster", 
                            label = "Clustering", 
                            choices = c("Both" = "both",
                                        "None" = "none",
                                        "Row" = "row",
                                        "Column" = "column")),
                
                # select seriation
                ##selectInput(inputId = "heat_seriation", 
                            ##label = "Seriation", 
                            ##choices = c("Mean" = "mean",
                                        ##"None" = "none",
                                        ##"Optimal Leaf Ordering" = "OLO",
                                        ##"Gruvaeus and Wainer" = "GW")),
              
                # input color scheme
                selectInput(inputId = "heat_inColor",
                            label = "Colorscheme",
                            choices = colorschemes)),
              
              
              
              #####################
              # 3.2.3.2 DOWNLOAD #
              #####################
              
              
              
              wellPanel(
                
                
                # panel title
                h4("Download"),
                
                # download button
                downloadButton("heat_download", 'Download Dataset'))))),
  
  
  #####################
  # 3.4 OVERLAY DRUGS #
  #####################
  
  
  
  # tab title
  tabPanel(p("Time profiles (by drug)", style="padding: 15px; margin: 0px;border-style:solid; border-color: #e9e9e9;"),
           
     fluidRow(
       
       
       
       ##################
       # 3.4.1 LEFT BAR #
       ##################
       
       
       
       # set left bar
       column(2,
              
              
              
          ######################
          # 3.4.1.1 DATA INPUT #
          ######################
          
          
          
          wellPanel(
            
            # panel title
            h4("Dataset"),
            
            # input gene
            selectInput(inputId = "drugs_inGene", 
                        label = HTML("Gene <br/>(type to trigger autofill)"), 
                        choices = geneNames_reduced),
            
            # link to GeneCards URL
            uiOutput("drugs_geneURL"),
            br(),
            
            # input drugs
            selectizeInput(inputId = "drugs_inDrugs",
                           label = "Drugs",
                           choices = drugNames,
                           multiple = TRUE,
                           selected = drugNames[1])),
          
          
          
          ######################
          # 3.4.1.2 PARAMETERS #
          ######################
          
          
          
          wellPanel(
            
            # panle title
            h4("Parameters"),
            
            
            conditionalPanel(
              
              condition = "input.drugs_activePlot == 'Single Tissue'",
            
            # input tissue
            selectInput(inputId = "drugs_inTissue",
                        label = "Tissue",
                        choices = tissueNames,
                        selected = tissueNames[1])),
            
            # input dosages
            selectInput(inputId = "drugs_inDosage",
                        label = "Concentration",
                        choices = c("Low", "High")))),
       
       
       
       ####################
       # 3.4.2 MAIN PANEL #
       ####################
       
       
       
       column(8,
              
          # formatting headers
          tags$head(
            tags$style(type='text/css', 
                       ".nav-tabs {font-size: 16px} ")),          
          tabsetPanel(
            id = "drugs_activePlot",
            type = "tabs",
            
            # display timeplot
            tabPanel("Single Tissue",
                     br(),    
          
              # output timeplot
              plotlyOutput("drugs_timePlot",
                           width = "100%")),
            
            
            tabPanel("Stratify By Tissue",
               
              br(),
              
              plotOutput("drugs_gridPlot")),
  
          
            # display data table
            tabPanel("Data",
                     DT::dataTableOutput("drugs_table"),
                     style = "height:700px; overflow-y: scroll;overflow-x: scroll;"))),
       
       
       
       ###################
       # 3.4.3 RIGHT BAR #
       ###################
       
       
       # set right bar
       column(2,
              
        ####################
        # 3.4.3.1 DOWNLOAD #
        ####################
              
              
        wellPanel(
          
          #panel title
          h4("Download"),
          
          # download button
          downloadButton("downloadDrugs", 'Download Dataset'))))),

  
  
  tabPanel(p("Common genes for sensitivity/resistance", style="padding: 15px; margin: 0px;border-style:solid; border-color: #e9e9e9;"), 
           fluidRow(
             
             ##################
             # 3.4.1 LEFT BAR #
             ##################
             
             
             # set left bar
             column(2,
                    ######################
                    # 3.4.1.1 DATA INPUT #
                    ######################
                    wellPanel(
                      p("Expression of genes most correlated with log(GI50)"),
                      p("(uses Pearson correlation)"),
                      br(),
                      selectInput('time_point', 'Time point', c("2 hrs" = "2 h", "6 hrs" = "6 h", "24 hrs" = "24 h"), selected="24 h", selectize=TRUE),
                      br(),
                      selectInput('drug_conc', 'Drug concentration', c ("High" = "high", "Low" = "low"), selected="low", selectize=TRUE),
                      br(),
                      checkboxGroupInput("drug_variable", "Drugs:",
                                         choices = list("Azacytidine"="Azacytidine", 
                                                        "Bortezomib"="Bortezomib", 
                                                        "Cisplatin"="Cisplatin", 
                                                        "Dasatinib"="Dasatinib", 
                                                        "Doxorubicin"="Doxorubicin", 
                                                        "Erlotinib"="Erlotinib", 
                                                        "Geldanamycin"="Geldanamycin", 
                                                        "Gemcitabine"="Gemcitabine",
                                                        "Lapatinib"="Lapatinib", 
                                                        "Paclitaxel"="Paclitaxel",
                                                        "Sirolimus"="Sirolimus", 
                                                        "Sorafenib"="Sorafenib", 
                                                        "Sunitinib"="Sunitinib",
                                                        "Topotecan"="Topotecan", 
                                                        "Vorinostat"="Vorinostat"),
                                         selected = c("Cisplatin","Topotecan", "Doxorubicin", "Gemcitabine"))
                      
                      )),
             column(8,
                    
                    # formatting headers
                    tags$head(
                      tags$style(type='text/css', 
                                 "#venn ~ .tab-content {padding-left: 20px;     
                                    border: 1px solid #ddd;
                                    border-top-color: white;
                                    border-right-color: white;}
                                 .nav-tabs {font-size: 16px}")),          
                    tabsetPanel(
                      id = "venn",
                      type = "tabs",
                      
                      tabPanel("Negatively correlated genes",
                               tags$br(),
                               uiOutput("upset_title"),
                               fluidRow(
                                 column(12, plotOutput(outputId = "upset1",  height = "90%", width="95%"))
                               ),
                               tags$br(),
                               tags$hr(),
                               uiOutput("venn_title"),
                               fluidRow(
                                 column(12,  plotOutput(outputId = "venn1", height = "500px", width="95%"))
                               ),
                               tags$hr(),
                               tags$p("Negatively correlated genes lists"),
                               fluidRow(
                                 column(12, DT::dataTableOutput("negative_table"))
                               ) 
                      ),
                      tabPanel("Positively correlated genes",
                               tags$br(),
                               uiOutput("upset_title2"),
                               fluidRow(
                                 column(12, plotOutput(outputId = "upset2",  height = "90%", width="95%"))
                               ),
                               tags$hr(),
                               uiOutput("venn_title2"),
                               fluidRow(
                                 column(12,  plotOutput(outputId = "venn2", height = "500px", width="95%"))
                               ),
                               tags$br(),
                               tags$hr(),
                               tags$p("Positively correlated genes lists"),
                               fluidRow(
                                 column(12, DT::dataTableOutput("positive_table"))
                               )
                      )
                    )
             )
             
             )#close fluidRow
           ),#close tabPanel venn
  tabPanel(p("Help...", style="padding: 15px; margin: 0px;border-style:solid; border-color: #e9e9e9;"), 
           fluidRow(
             # set left bar
             column(3,
                    wellPanel(
                      p("If you have any questions about TPWshiny please contact the Support Team.")
                    )
                    ),
             column(8,
                    fluidRow(  tags$br(),
                               p("Find more information about TPWShiny and Documentation at the following link"),
                               a("User Manual and more information", href = "https://brb.nci.nih.gov/TPWshiny/"),
                               tags$br(),
                               tags$br(),
                               p("If you have additional  questions or concerns email the Support Team at: ncitpwsupport@mail.nih.gov"),
                               tags$br(),
                               p("To facilitate debugging any issue that may be caused by the use of R packages different than the ones used and tested by our Team, please include the information below in your support request email."),
                               
                               uiOutput("packageVersion")
                      )
                    )
             
           )#close fluidRow
  )#close tabPanel about       
         
)#close main  navbarPage


############
# 4 SERVER #
############



server <- function(input, output, session) {
  # IMPORTANT!
  # this is needed to terminate the R process when the
  # shiny app session ends. Otherwise, you end up with a zombie process
  session$onSessionEnded(function() {
    stopApp()
  })
  
  
  options(warn = -1) 
  updateSelectizeInput(session, 'inGene', choices = geneNames, server = TRUE)
  updateSelectizeInput(session, 'drugs_inGene', choices = geneNames, server = TRUE)
  updateSelectizeInput(session, 'heat_inGenes', choices = geneNames, server = TRUE)
  
  ##############
  # 4.1 VALUES #
  ##############
  
  
  
  # reactive values
  values <- reactiveValues(
   
    react_val_genesNotFound = "",
    # database
    data = read.csv(paste0("https://",server_name,"/GeneExpressionNCI60/rshiny/singlegene/", geneNames[1], "_", tolower(drugNames[1]), "_formatted.csv"), 
                    header = TRUE, 
                    sep = ","),
    
    baselinedata = read.csv(paste0("https://",server_name,"/GeneExpressionNCI60/rshiny/singlegene/", geneNames[1], "_", "baseline", "_formatted.csv"), 
                    header = TRUE, 
                    sep = ","),
    
    # data to plot
    displaydata = read.csv(paste0("https://",server_name,"/GeneExpressionNCI60/rshiny/singlegene/", geneNames[1], "_", tolower(drugNames[1]), "_formatted.csv"),
                           header = TRUE, 
                           sep = ","),
     
    # heatmap database
    heatData = read.csv(file = paste0("https://",server_name,"/GeneExpressionNCI60/rshiny/heatmap/6hr_L_", pathwayNames[1], "_", tolower(drugNames[1]), "_formatted.csv", sep = ""), 
                        header = TRUE, 
                        sep = ","),
    
    # drugs database
    drugsData = read.csv(file = paste0("https://",server_name,"/GeneExpressionNCI60/rshiny/singlegene/", geneNames[1], "_", tolower(drugNames[1]), "_formatted.csv"), 
                         header = TRUE,
                         sep = ","),
    
    drugsFullData = read.csv(file = paste0("https://",server_name,"/GeneExpressionNCI60/rshiny/singlegene/", geneNames[1], "_", tolower(drugNames[1]), "_formatted.csv"),
                             header = TRUE,
                             sep = ","))
  
  
  
  #################
  # 4.2 LISTENERS #
  #################
  
  
  
  #####################
  # 4.2.1 SINGLE GENE #
  #####################
  
  
  
  #  update single gene database
  observeEvent({
    
    # listen for gene and drug updates
    input$inGene
    input$inDrug},{
      
      # update gene-drug path
      #fileName <- paste("data/", input$inGene, "_", input$inDrug, "_formatted.csv", sep = "") 
      
      source_directory=paste0("https://",server_name,"/GeneExpressionNCI60/rshiny/singlegene/");
      fileName_server <- paste0(source_directory, input$inGene, "_", tolower(input$inDrug), "_formatted.csv", sep = "") 
      fileName = fileName_server
      # update database
      values$data <- read.csv(file = fileName, header = TRUE, sep = ",")
      values$displaydata <- values$data
     
      filenName_baselinedata  <- paste0(source_directory, input$inGene, "_", "baseline", "_formatted.csv", sep = "") 
      values$baselinedata <- read.csv(file = filenName_baselinedata, header = TRUE, sep = ",")
      
      # update gene and drug
      values$gene <- input$inGene
      values$drug <- input$inDrug})
  

  
  # update single gene display data
  observeEvent({
    
    # listen for any parameter change
    input$inTissue
    input$inTime
    input$inDosage
    input$inMeasure
    values$displaydata},{
      
      # grab gene database
      data <- data.frame(values$data)
      
      # update column header
      setnames(data, 
               old = c(input$inMeasure), 
               new = c("Resistance.Measure"), 
               skip_absent = TRUE)
      
      # catch no selected tissues
      if(is.null(input$inTissue)){
        tissue <- tissueNames
      }else{
        tissue <- input$inTissue
      }
      
      # filter by tissue origin
      display <- subset(data, Tissue.Origin%in%toupper(tissue))
      
      # filter by time
      display <- subset(display, Time == input$inTime)
      
      # filter by dose
      highDose <- data[1, 'Dose']
      lowDose <- data[nrow(data),'Dose']
      
     
      if(input$inDosage == "Baseline (no drug)"){ 
        basedata = values$baselinedata
        currentdata = values$data
        for(cl in cellLineNames) {
           basedata[ basedata$"Cell.Line" == cl, "GI50"] = currentdata[currentdata$"Cell.Line" == cl,"GI50"][1]
         }
       display <- subset(basedata, Tissue.Origin%in%toupper(tissue))
        # update column header
        setnames(display, 
                 old = c(input$inMeasure), 
                 new = c("Resistance.Measure"), 
                 skip_absent = TRUE)
        
        display <- subset(display, Time == input$inTime)
        display <- subset(display, Dose == 0)
      }
      else {
        # find low and high doses
        if(input$inDosage == "Low"){
          targetDose <- lowDose
        }else{
          targetDose <- highDose
        }
        display <- subset(display, Dose == targetDose)
      }
      
      # store data
      values$displaydata <- display})
  
  
  # tissue selector
  observeEvent(
    
    # listen for select all click
    input$selectall,{
      
        # check all tissue types
        updateCheckboxGroupInput(session, inputId = "inTissue", label = "Tissue",
                                 choices = tissueNames, selected = tissueNames,
                                 inline = FALSE)
      }
    
    )
  
  # tissue selector
  observeEvent(
    
    # listen for select all click
    input$deselectall,{
      
       updateCheckboxGroupInput(session, inputId = "inTissue", label = "Tissue",
                                 choices = tissueNames, selected = NULL,
                                 inline = FALSE)
        
    
    }
    
  )
  
  
  
  
  ##################
  # 4.2.2 HEATMAPS #
  ##################
  
  
  
  # update heatmap database
  observeEvent({
    
    
    # listen for path and drug updates
    input$heat_groupType
    input$heat_grouptfType
    input$heat_inPath
    input$heat_inRecep
    input$heat_inTransFactor
    input$heat_inGenes
    input$heat_inDrug
    input$heat_inDose
    input$heat_inTime},{
      
      values$react_val_genesNotFound = ""
      
      if(input$heat_groupType == 'manual'){
       
        data <- read.table(text = "",
                           colClasses = c("character", rep("numeric", 60)),
                           col.names =  c("Gene", cellLineNames),
                           header = TRUE)
        
        genes <- input$heat_inGenes
        
        
        
        for(gene in genes){
          
           if(!(gene %in% geneNames)) { 
             values$react_val_genesNotFound <- paste(values$react_val_genesNotFound, gene)
             next; 
          }
          
          # get gene-drug path
          source_directory=paste0("https://",server_name,"/GeneExpressionNCI60/rshiny/singlegene/");
          fileName_server <- paste0(source_directory, gene, "_", tolower(input$heat_inDrug), "_formatted.csv", sep = "") 
          fileName = fileName_server
            # grab gene data
            geneData <- read.csv(file = fileName, header = TRUE, sep = ",")
             
            # find dose values
            highDose <- geneData[1, 'Dose']
            lowDose <- geneData[nrow(geneData),'Dose']
          # find target dose
          if(input$heat_inDose == "L"){
            targetDose <- lowDose
          }else{
            targetDose <- highDose}
          
          # filter by dose
          geneData <- subset(geneData, Dose == targetDose)
          
          # get time
          time <- input$heat_inTime
          
          # extract number
          time <- gsub("([0-9]+).*$", "\\1", input$heat_inTime)
          
          # add back space
          time <- paste(time, "hrs", sep = " ")
          
          # filter by time
          geneData <- subset(geneData, Time == time)
          
          geneData <- geneData[complete.cases(geneData$Gene.Expression), ]
         
          
          # get cell lines
          cellLines <- as.vector(geneData$Cell.Line)
          
          # get expression data
          expressions <- as.vector(geneData$Gene.Expression)
          
          geneData <- read.table(text = "",
                                 colClasses = c("character", rep("numeric", length(cellLines))),
                                 col.names = c("Gene", cellLines),
                                 header = TRUE)
          
          geneData[1, 1] <- c(gene)
          geneData[1, -1] <- expressions
          
         
          data <- rbind(geneData, data[, names(geneData)])
          
        }
        
        
        data[!sapply(data, function(x) all(is.na(x)))]
        old_names = colnames(data)
        new_names =  gsub("\\.", "-", old_names)
        new_names =  replace(new_names, new_names=="X786-0", "786-0")
        
        colnames(data) = new_names
        
      
        values$heatData <- data
        
        if(input$heat_groupType == 'manual'){
          # download heatmap data
          output$heat_download <- downloadHandler(
            
            # set file name
            filename = function(){
              "dataset" = paste("Heatmap_", input$heat_inTime, "_", input$heat_inDose, "_", input$heat_inDrug, "_UserDefinedGeneSet", '.csv', sep='')},
            
            # write heatmap data to csv
            content = function(file){
              "dataset" = write.csv(values$heatData, file)})}
        
      }
      if((input$heat_groupType == 'pathway')||(input$heat_groupType == 'recep')){
      
      group <- switch(input$heat_groupType,
                      "pathway" = input$heat_inPath,
                      "recep" = paste("Receptor", match(input$heat_inRecep, recepNames), sep = ""))
      
      # grab path
      source_directory=paste0("https://",server_name,"/GeneExpressionNCI60/rshiny/heatmap/");
      fileName_server <- paste0(source_directory, input$heat_inTime,"_", input$heat_inDose, "_", group, "_",  tolower(input$heat_inDrug), "_formatted.csv") 
      fileName = fileName_server
  
  
      # open file
      data <- read.csv(file = fileName, header = TRUE, sep = ",",check.names=FALSE)
      # remove GI50 row
      data <- data[-1,]
      data <- data[, colSums(is.na(data)) == 0]
      
      # store data
      values$heatData <- data
      
      if(input$heat_groupType == 'pathway'){
      # download heatmap data
      output$heat_download <- downloadHandler(
        
        # set file name
        filename = function(){
          "dataset" = paste("Heatmap_", input$heat_inTime, "_", input$heat_inDose, "_", input$heat_inDrug, "_", input$heat_inPath, '.csv', sep='')},
        
        # write heatmap data to csv
        content = function(file){
          "dataset" = write.csv(values$heatData, file)})}
      
      
      if(input$heat_groupType == 'recep') {
        # download heatmap data
        output$heat_download <- downloadHandler(
          
          # set file name
          filename = function(){
            "dataset" = paste("Heatmap_", input$heat_inTime, "_", input$heat_inDose, "_", input$heat_inDrug, "_", input$heat_inRecep, '.csv', sep='')},
          
          # write heatmap data to csv
          content = function(file){
            "dataset" = write.csv(values$heatData, file)})}
      }
      
      if((input$heat_groupType == 'tf')&&(input$heat_grouptfType == 'transfactor')){
      
        group <- switch(input$heat_grouptfType,"transfactor" = input$heat_inTransFactor)
        
        source_directory=paste0("https://",server_name,"/GeneExpressionNCI60/rshiny/heatmap/");
        fileName_server <- paste0(source_directory, input$heat_inTime,"_", input$heat_inDose, "_", group, "_",  tolower(input$heat_inDrug), "_formatted.csv") 
        fileName = fileName_server
      
       tryCatch({
        data <- read.csv(file = fileName, header = TRUE, sep = ",",check.names=FALSE)
        data <- data[-1,]
        data <- data[, colSums(is.na(data)) == 0]
        },warning=function(cond) { message(cond)},    error = function(e){  message(e)} )
        
        
        values$heatData <- data
        
        # download heatmap data
        output$heat_download <- downloadHandler(
          
          # set file name
          filename = function(){
            "dataset" = paste("Heatmap_", input$heat_inTime, "_", input$heat_inDose, "_", input$heat_inTransFactor, "_", input$heat_inDrug, '.csv', sep='')},
          
          # write heatmap data to csv
          content = function(file){
            "dataset" = write.csv(values$heatData, file)})
      }
      
      if((input$heat_groupType == 'tf')&&(input$heat_grouptfType == 'b')){
        
        source_directory=paste0("https://",server_name,"/GeneExpressionNCI60/rshiny/heatmap/");
        fileName_server <- paste0(source_directory, input$heat_inTime,"_", input$heat_inDose, "_", "alltf", "_",  tolower(input$heat_inDrug), "_formatted.csv") 
        fileName = fileName_server
        
          data <- read.csv(file = fileName, header = TRUE, sep = ",",check.names=FALSE)
          data <- data[-1,]
          data <- data[, colSums(is.na(data)) == 0]
        
        values$heatData <- data
       
        # download heatmap data
        output$heat_download <- downloadHandler(
          
          # set file name
          filename = function(){
            "dataset" = paste("Heatmap_", input$heat_inTime, "_", input$heat_inDose, "_", "alltf_", input$heat_inDrug, '.csv', sep='')},
          
          # write heatmap data to csv
          content = function(file){
            "dataset" = write.csv(values$heatData, file)})
      }
      
      if((input$heat_groupType == 'tf')&&(input$heat_grouptfType == 'c')){
        
        source_directory=paste0("https://",server_name,"/GeneExpressionNCI60/rshiny/heatmap/");
        fileName_server <- paste0(source_directory, input$heat_inTime,"_", input$heat_inDose, "_", input$heat_inTransFactor, "_", "alldrug", "_formatted.csv") 
        fileName = fileName_server
          data <- read.csv(file = fileName, header = TRUE, sep = ",",check.names=FALSE)
        
        values$heatData <- data
        
        # download heatmap data
        output$heat_download <- downloadHandler(
          
          # set file name
          filename = function(){
            "dataset" = paste("Heatmap_", input$heat_inTime, "_", input$heat_inDose, "_", input$heat_inTransFactor, "_alldrug", '.csv', sep='')},
          
          # write heatmap data to csv
          content = function(file){
            "dataset" = write.csv(values$heatData, file)})
      }
      })
  
  
  ######################
  # 4.2.3 DRUG OVERLAY #
  ######################
  
  # update drug overlay database
  observeEvent({
    
    # listen for drug and gene changes
    input$drugs_inGene
    input$drugs_inDrugs
    input$drugs_inTissue
    input$drugs_inDosage},{
      
      data <- data.frame("Cell Line" = character(), 
                         "Tissue Origin" = character(), 
                         "GI50" = numeric(),
                         "DoubleTime" = numeric(), 
                         "MultiRes" = numeric(), 
                         "Dose" = numeric(), 
                         "Time" = character(), 
                         "Gene Expression" = numeric(),
                         "Drug" = character())
      
      drugs <- input$drugs_inDrugs
      
      if(is.null(drugs))
        drugs <- drugNames
      
      # collect files
      for(drug in drugs){
        
        # get filePath
        source_directory=paste0("https://",server_name,"/GeneExpressionNCI60/rshiny/singlegene/");
        fileName_server <- paste0(source_directory, input$drugs_inGene, "_", tolower(drug), "_formatted.csv", sep = "") 
        fileName = fileName_server
        
        # get drug data
        drugData <- read.csv(file = fileName, header = TRUE, sep = ",")
        
        # calculate target dose
        highDose <- drugData[1, 'Dose']
        lowDose <- drugData[nrow(drugData),'Dose']
         
        if(input$drugs_inDosage == "Low"){
          targetDose <- lowDose
        }else{
          targetDose <- highDose}
        
        # filter by tissue origin
        drugData <- subset(drugData, Dose == targetDose)
        
        
        # add drug name column
        drugData$Drug = rep(drug, nrow(drugData))
        
        # merge data
        data <-rbind(data, drugData)
         
        }
      
      values$drugsData <- data})
  
  
  # update full overlay database
  observeEvent({
    
    # listen for drug and gene changes
    input$drugs_inGene
    input$drugs_inDrugs
    input$drugs_inTissue
    input$drugs_inDosage},{
      
      data <- data.frame("Cell Line" = character(), 
                         "Tissue Origin" = character(), 
                         "GI50" = numeric(),
                         "DoubleTime" = numeric(), 
                         "MultiRes" = numeric(), 
                         "Dose" = numeric(), 
                         "Time" = character(), 
                         "Gene Expression" = numeric(),
                         "Drug" = character())
      
      drugs <- input$drugs_inDrugs
      
      if(is.null(drugs))
        drugs <- drugNames
      
      # collect files
      for(drug in drugs){
        
        # get filePath
        source_directory=paste0("https://",server_name,"/GeneExpressionNCI60/rshiny/singlegene/");
        fileName_server <- paste0(source_directory, input$drugs_inGene, "_", tolower(drug), "_formatted.csv", sep = "") 
        fileName = fileName_server
        
        # get drug data
          drugData <- read.csv(file = fileName, header = TRUE, sep = ",")
        
           # add drug name column
          drugData$Drug = rep(drug, nrow(drugData))
          
          # merge data
          data <-rbind(data, drugData)
        }
      
      values$drugsFullData <- data})
  
  
  #############
  # 4.3 TABLE #
  #############
  
  output$table <- DT::renderDataTable(
    
    # output datatable
    DT::datatable({ basedata = values$baselinedata
                    currentdata = values$data
                    for(cl in cellLineNames) {
                      basedata[ basedata$"Cell.Line" == cl, "GI50"] = currentdata[currentdata$"Cell.Line" == cl,"GI50"][1]
                    }
                    rbind(values$data, basedata)
                  }, 
                  colnames=c("Cell Line",
                             "Tissue Origin",
                             "log(GI50)",
                             "Doubling Time",
                             "Multidrug Resistance",
                             "Concentration",
                             "Time",
                             "Gene Expression"),
                  options = list(paging = FALSE, searching = FALSE),
    rownames = FALSE) %>%
    formatRound(columns=c('Gene.Expression'), digits=3))
  
  
  output$drugs_table <- DT::renderDataTable(
    
    # output datatable
    DT::datatable(values$drugsFullData,  
                  colnames=c("Cell Line",
                             "Tissue Origin",
                             "log(GI50)",
                             "Doubling Time",
                             "Multidrug Resistance",
                             "Concentration",
                             "Time",
                             "Gene Expression",
                             "Drugs"),
                  options = list(paging = FALSE, searching = FALSE),
    rownames = FALSE)  %>%
      formatRound(columns=c('Gene.Expression'), digits=3))
  
  
  
  ###################
  # 4.4 SCATTERPLOT #
  ###################
  
  # display scatterplot
  output$scatterplot <- renderPlotly({
    if(length(input$inTissue) < 1) {
      stop("\nPlease select at least one tissue type from the list in the side panel.")
    }
    
    # grab data
    data <- values$displaydata
    
    # calculate target dose
    highDose <- values$data[1, 'Dose']
    lowDose <- values$data[nrow(values$data),'Dose']
    if(input$inDosage == "Low"){
      targetDose <- lowDose
    }else{
      targetDose <- highDose}
    
    
    # get proper ylabel
    # get proper xlabel
    #if(input$inDrug == 'baseline'){
    if(input$inDosage == "Baseline (no drug)"){
      ylabel <- "Gene Expression (log2)"
      xlabel <- switch(input$inMeasure, "GI50" = paste0("log(GI50) ",values$drug), 
                       "DoubleTime" = "Doubling Time", 
                       "MultiRes" = "Multidrug Resistance")
    }else{
      ylabel <- "Gene Expression (log2, relative to baseline)"
      xlabel <- switch(input$inMeasure, "GI50" = "log(GI50)", 
                     "DoubleTime" = "Doubling Time", 
                     "MultiRes" = "Multidrug Resistance")
    }

    # formatting plot
    #if(input$inDrug == 'baseline'){
    if(input$inDosage == "Baseline (no drug)"){
      title <- paste0("Gene expression of ", values$gene, " at baseline (", input$inTime, ")" , " vs. ", xlabel)
    }else{
      title <- paste0("Gene expression of ", values$gene, " following ", input$inTime, " exposure to ", targetDose, " nM of ", values$drug, " vs. ", xlabel)
    }
    
    data_rounded = data
    data_rounded$Gene.Expression = round(data_rounded$Gene.Expression,3)
    ggformatted <- ggplot(data_rounded,
                          aes(x = Resistance.Measure,
                              y = Gene.Expression))+ 
      labs(title = title, x = xlabel, y = ylabel)+ 
      theme(plot.title=element_text(size=12, hjust = 0.5),
            axis.text.x=element_text(size=10), 
            axis.text.y=element_text(size=10),
            axis.title.x=element_text(size=12),
            axis.title.y=element_text(size=12),
            legend.text=element_text(size= 10),
            legend.title=element_text(size= 10),
            plot.margin=unit(c(1.0,0,0,0),"cm")) + 
      scale_colour_manual(values = tissueColors)+
      theme(legend.title = element_blank())
    
    # coloring tissue origins
    plot <- ggformatted + geom_point(
      aes(color = Tissue.Origin,
          Cell.Line = Cell.Line), 
      size = 2)
    
    # linear regression line
    if(input$gene_regType == "linear"){
      plot <- plot + stat_smooth(method = "lm", 
                                 formula = y ~ x,
                                 level = input$gene_confInt)}
    
    # polynomial regression line
    if(input$gene_regType == "poly"){
      plot <- plot + stat_smooth(method = "lm", 
                                 formula = y ~ poly(x, 2),
                                 level = input$gene_confInt)}
    
    # loess regression
    if(input$gene_regType == "loess"){
      plot <- plot + geom_smooth(method = "loess", 
                                 level = input$gene_confInt)}
    
    # to plotly
    plot <- ggplotly(plot, 
                     height = 700,
                     tooltip = c("Cell.Line",
                                 "x",
                                 "y"))%>% 
      layout(legend = list(x = 1, y = 1,title=list(text=' Tissue of origin<br>'))) }) 
  
  
  
  ################
  # 4.5 BAR PLOT #
  ################
  
  # display barplot
  output$bar <- renderPlotly({
  
    if(length(input$inTissue) < 1) {
      stop("\nPlease select at least one tissue type from the list in the side panel.")
    }
    
    # grab data
    data <- values$displaydata
    
    # calculate target dose
    highDose <- values$data[1, 'Dose']
    lowDose <- values$data[nrow(values$data),'Dose']
    if(input$inDosage == "Low"){
      targetDose <- lowDose
    }else{
      targetDose <- highDose}
    
    # get sorting measure
    sortby <- input$gene_barSort 
    
    # account for renamining
    if(input$gene_barSort == input$inMeasure){
      sortby <- "Resistance.Measure"}
    
    # sorting
    if(!(sortby == "None")){
      data <- switch(input$gene_barOrder,
                     "up" = data[order(data[[sortby]]), ],
                     "down" = data[rev(order(data[[sortby]])), ])}
    
    if(input$gene_bargroup){
      data <- data[order(data[["Tissue.Origin"]]),]
    }
    
        
    
    # get title
    #if(input$inDrug == "baseline"){
    if(input$inDosage == "Baseline (no drug)"){
        title <- paste0("Gene expression of ", values$gene, " at baseline (", input$inTime,")")
    }else{
      title <- paste0("Gene expression of ", values$gene, " following ",input$inTime," exposure to ", targetDose, " nM of ", values$drug)
    }
    
    # get proper ylabel
    #if(input$inDrug == 'baseline'){
    if(input$inDosage == "Baseline (no drug)"){
        ylabel <- "Gene Expression (log2)"
    }else{
      ylabel <- "Gene Expression (log2, relative to baseline)"
    }
    
    data_rounded = data
    data_rounded$Gene.Expression = round(data_rounded$Gene.Expression,3)
    # set theme
    ggformatted <- ggplot(data_rounded, aes(x = Cell.Line, 
                                    y = Gene.Expression, 
                                    fill = Tissue.Origin,
                                    text = paste('Cell Line:', Cell.Line,
                                                 '<br>Gene Exp: ', Gene.Expression,
                                                 '<br>Tissue Origin: ', Tissue.Origin))) + 
      
      
      # labels
      labs(title = title, x = "Cell Line", y = ylabel, fill = "Tissue Origin") + 
      
      # text size
      theme(plot.title=element_text(size=12, hjust = 0.5),
            axis.text.x=element_text(size=10), 
            axis.text.y=element_text(size=9),
            axis.title.x=element_text(size=12),
            axis.title.y=element_text(size=12),
            legend.text=element_text(size= 10),
            legend.title=element_text(size= 10),
            panel.grid.major.y = element_blank(),
            plot.margin=unit(c(1.5,0,0,0),"cm")) + 
       scale_fill_manual(values = tissueColors)+
      
      # sort by cell line and reverse the order
      scale_x_discrete(limits=rev(data$Cell.Line))+
      theme(legend.title = element_blank())
    
    # add bars
    plot <- ggformatted + geom_bar(stat = "identity") +
      # flip coordinates
      coord_flip()
    
    # show plot
    p <- ggplotly(plot, width = 800, height = 700, tooltip = c("text")) %>% 
      layout(legend = list(x = 1, y = 1,title=list(text=' Tissue of origin<br>'))) 
    
    })
  
  #################
  # 4.6 TIME PLOT #
  #################
  
  output$timeplot <- renderPlotly({
    
    if(length(input$inTissue) < 1) {
      stop("\nPlease select at least one tissue type from the list in the side panel.")
    }
    
    data <- subset(values$data, Tissue.Origin%in%toupper(input$inTissue))
    # recompute data to include all times
    highDose <- data[1, 'Dose']
    lowDose <- data[nrow(data),'Dose']
   
    if(input$inDosage == "Baseline (no drug)"){ 
      data <- subset(values$baselinedata, Tissue.Origin%in%toupper(input$inTissue))
    }
    else {
      if(input$inDosage == "Low"){
        targetDose <- lowDose
      }else{
        targetDose <- highDose}
      
      data <- subset(data, Dose == targetDose)
    }
    
    # if quartiles are selected for segmentation
    if(!(input$gene_quartMeasure == "none") &&
       !is.null(input$gene_quartiles) && 
       !(length(input$gene_quartiles) == 4)){
      
      # calculate quartiles
      bottom <- as.double(quantile(data[[input$gene_quartMeasure]], 0.25, na.rm = TRUE))
      middle <- as.double(quantile(data[[input$gene_quartMeasure]], 0.50, na.rm = TRUE))
      top <- as.double(quantile(data[[input$gene_quartMeasure]], 0.75, na.rm = TRUE))
      
      # filter first quartile
      if(!("first" %in% input$gene_quartiles)){
        data <- data[data[[input$gene_quartMeasure]] >= bottom, ]}
      
      # filter second quartile
      if(!("second" %in% input$gene_quartiles)){
        data <- data[data[[input$gene_quartMeasure]] >= middle | 
                       data[[input$gene_quartMeasure]] < bottom, ]}
      
      # filter third quartile
      if(!("third" %in% input$gene_quartiles)){
        data <- data[data[[input$gene_quartMeasure]] >= top | 
                       data[[input$gene_quartMeasure]] < middle, ]}
      
      # filter fourth quartile
      if(!("fourth" %in% input$gene_quartiles)){
        data <- data[data[[input$gene_quartMeasure]] < top, ]}}
    
    
    # get title
    #if(input$inDrug == "baseline"){
    if(input$inDosage == "Baseline (no drug)"){ 
      title <- paste0("Gene expression of ", values$gene, " at baseline")
    }else{
      title <- paste0("Gene expression of ", values$gene, " following exposure to ", targetDose, " nM of ", values$drug)}
    
    # format plot    
    ggformatted <- ggplot(data) +
      
      # labels
      labs(title = title, x = "Time", y = "Gene Expression (log2)") + 
      
      # text size
      theme(plot.title=element_text(size=12, hjust = 0.5),
            axis.text.x=element_text(size=10), 
            axis.text.y=element_text(size=10),
            axis.title.x=element_text(size=12),
            axis.title.y=element_text(size=12),
            legend.text=element_text(size= 10),
            legend.title=element_text(size= 10),
            plot.margin=unit(c(1.5,0,0,0),"cm")) + 
      
      # coloring
      scale_color_manual(name = "Tissue of origin", values= tissueColors) +
    
      # force scale order
      scale_x_discrete(limits=c("2 hrs","6 hrs","24 hrs"))+
      theme(legend.title = element_blank())
    
    # plot lines
    plot <- ggformatted + geom_line(
      aes(x = Time, 
          y = Gene.Expression, 
          color = Tissue.Origin,
          group = Cell.Line), size = 1)
    
    # to plotly
    ggplotly(plot, 
             height = 700)%>% 
      layout(legend = list(x = 1, y = 1,title=list(text=' Tissue of origin<br>')))
    })
  
  
  
  #################
  # 4.7 DOWNLOADS #
  #################
  
  
  # download all
  output$gene_download <- downloadHandler(
    
    # set file name
    filename = function(){
      
      switch(input$gene_downloadType,
             "full" = paste(input$inGene, "_", input$inDrug, '.csv', sep=''),
             "graphed" = paste(input$inGene, "_", ifelse(input$inDosage == "Baseline (no drug)", 'GeneExprNoDrug', input$inDrug), "_", input$inTime, "_", ifelse(input$inDosage == "Baseline (no drug)", paste0('GI50',input$inDrug), input$inDosage), '.csv', sep=''))},
    
    
    
     # write displaydata to csv
    content = function(file) {
      basedata = values$baselinedata
      currentdata = values$data
      for(cl in cellLineNames) {
        basedata[ basedata$"Cell.Line" == cl, "GI50"] = currentdata[currentdata$"Cell.Line" == cl,"GI50"][1]
      }
      
      switch(input$gene_downloadType, 
             "full" = write.csv(rbind(values$data,basedata), file),
             "graphed" = write.csv(values$displaydata, file))})
  
  
  
  # download drug data
  output$downloadDrugs <- downloadHandler(
    
    # set file name
    filename = function(){
      paste(input$drugs_inGene, "_MultiDrug_", input$drugs_inTissue, "_",input$drugs_inDosage, '.csv', sep='')},
    
    # write drug data to csv
    content = function(file) {
      write.csv(values$drugsFullData, file)})
  
  
  
  ###################
  # 4.8 CORRELATION #
  ###################
  
  # pearsons correlation
  output$gene_pearsons <- renderText({
    if(length(input$inTissue) < 1) {    stop()    }
    
    # load x and y
    x <- values$displaydata$Resistance.Measure
    y <- values$displaydata$Gene.Expression
    
    # conduct correlation
    correlation <- cor.test(x, y, use="complete.obs", method = "pearson")
    
    # format
    paste( round(correlation$estimate, 4), 
           " (p = ", 
           round(correlation$p.value, 4),
           ")", sep = "")})
  
  # spearmans correlation
  output$gene_spearmans <- renderText({
    if(length(input$inTissue) < 1) {    stop()    }
    
    # load x and y
    x <- values$displaydata$Resistance.Measure
    y <- values$displaydata$Gene.Expression
    
    # conduct correlation
    correlation <- cor.test(x, y, use="complete.obs", method = "spearman")
    
    # format
    paste( round(correlation$estimate, 4), 
           " (p = ", 
           round(correlation$p.value, 4),
           ")", sep = "")})
  
  # kendall correlation
  output$gene_kendall <- renderText({
    if(length(input$inTissue) < 1) {    stop()    }
    
    # load x and y
    x <- values$displaydata$Resistance.Measure
    y <- values$displaydata$Gene.Expression
    
    # conduct correlation
    correlation <- cor.test(x, y, use="complete.obs", method = "kendall")
    
    # format
    paste( round(correlation$estimate, 4), 
           " (p = ", 
           round(correlation$p.value, 4),
           ")", sep = "")})
  
  
  
  
  ################
  # 4.9 HEATMAPS #
  ################
  
  output$heatmap <- renderPlotly({
    
    # load data
    data <- values$heatData
    
    
    if((input$heat_groupType == 'tf')&&(input$heat_grouptfType == 'transfactor')) {
      if(input$heat_inTransFactor  %in% names(TF_with_no_heatmap)) {
          stop(paste0(TF_with_no_heatmap[input$heat_inTransFactor]))
      }
    }
    
    if(class(data) =='function') {
      stop(error_message_serverfile)  
    } 
    
  
    if(nrow(data) < 2) {
     stop("\nPlease select at least two genes from the text-box in the side panel.")
    }
    
    rownames <- data[, 1]
    
    data <- data[,-1]
    
    row.names(data) <- rownames
    mat_data = as.matrix(data)
    
    
    # convert colorscheme to list
    colorscheme <- switch(input$heat_inColor, 
                          "Blue-White-Red" = colorpanel(75, "blue", "white", "red"),
                          "Blue-Black-Red" = colorpanel(75, "blue", "black", "red"),
                          "Blue-White-Orange" = colorpanel(75, "blue", "white", "orange"),
                          "Blue-Black-Orange" = colorpanel(75, "blue", "black", "orange"),
                          "Green-White-Red" = colorpanel(75, "green", "white", "red"),
                          "Green-Black-Red" = colorpanel(75, "green", "black", "red"),
                          "Cyan-Navy" = colorpanel(75, "cyan", "navy"))
    
    # get dose level
    dose <- switch(input$heat_inDose,
                   "L" = "low",
                   "H" = "high")
    
    # get dose level
    labeltime <- switch(input$heat_inTime,
                   "2hr" = "2 hrs",
                   "6hr" = "6 hrs",
                   "24hr" = "24 hrs")
  
    # get title
    title <- switch(input$heat_groupType,
                   "pathway" = paste("Expression of genes in", toupper(substr(input$heat_inPath, 3, nchar(input$heat_inPath) - 7)), "pathway \nfollowing",labeltime,"exposure to\n", dose, "concentration of", input$heat_inDrug, sep = " "),
                   "recep" = paste("Expression of genes in the receptor group\n", input$heat_inRecep, "following",labeltime, "exposure to\n", dose, "concentration of", input$heat_inDrug, sep = " "),
                   "manual" = paste("Expression of selected genes following",labeltime,"exposure to\n", dose, "concentration of", input$heat_inDrug, sep = " "),
    
    
    title <- switch(input$heat_grouptfType,
                    "transfactor" = paste("Expression of genes targeted by transcription factor",input$heat_inTransFactor, "\nfollowing", labeltime,"exposure to",dose,"concentration of",input$heat_inDrug, sep = " "),
                    "b" = paste("Gene expression of all transcription factors \nfollowing",labeltime, "exposure to", dose, "concentration of", input$heat_inDrug, sep = " "),
                    "c" = paste("Gene expression of transcription factor", input$heat_inTransFactor,"\nfollowing",labeltime, "exposure to", dose, "concentration of different drugs", sep = " ")))
    
    # determine clustering options
    clusterRow <- input$heat_cluster == "row" || input$heat_cluster == "both"
    clusterCol <- input$heat_cluster == "column" || input$heat_cluster == "both"
    
    
    total_row <- data[, 1]
    # heatmap.2
  
    jpeg(file=".heatmap_tmp.jpg")
    original_heatmap <- heatmap.3(mat_data, Rowv = clusterRow, Colv = clusterCol)
    
    
    
    dev.off()
    limit = max(abs(data),na.rm=TRUE)
    
 
    # If column have value in between 0 to 74, the heatmap will be small in size
    
    if(length(total_row) < 75){
      if(input$heat_cluster == 'both'){
    heatmaply(mat_data,
              main = title,
                trace = "none",
                Rowv = rev(original_heatmap$rowDendrogram),
                Colv = rev(original_heatmap$colDendrogram),
                dendrogram = input$heat_cluster,
                colors = colorscheme,limits = c(-limit, limit),
                row_dend_left = TRUE,
                subplot_heights=c(0.1, 0.70),subplot_widths=c(0.8, 0.2),column_text_angle = 90,fontsize_row = 12, fontsize_col = 12,
                plot_method = "plotly",colorbar_len = 0.15,colorbar_yanchor='top', colorbar_xpos=1.2, colorbar_ypos=0.90,key.title = "Color Key")%>%
          layout(height=1000)}
      
      else if(input$heat_cluster == 'row'){
        heatmaply(mat_data,
                  main = title,
                  trace = "none",margins = c(50, 50, 150, 40),
                  Rowv = rev(original_heatmap$rowDendrogram),
                  dendrogram = input$heat_cluster,
                  colors = colorscheme,limits = c(-limit, limit),
                  row_dend_left = TRUE,
                  branches_lwd=0.15,
                  subplot_widths=c(0.8, 0.2),column_text_angle = 90,fontsize_row = 12, fontsize_col = 12,
                  plot_method = "plotly",colorbar_len = 0.15,colorbar_yanchor='top', colorbar_xpos=1.2, colorbar_ypos=0.95,key.title = "Color Key")%>%
          layout(height=1000)}
      
      
      else if(input$heat_cluster == 'column'){
        heatmaply(mat_data,
                  main = title,
                  trace = "none",margins = c(50, 50, 150, 40),
                  Colv = rev(original_heatmap$colDendrogram),
                  dendrogram = input$heat_cluster,
                  colors = colorscheme,limits = c(-limit, limit),
                  subplot_heights=c(0.06, 0.93),column_text_angle = 90,fontsize_row = 12, fontsize_col = 12,
                  plot_method = "plotly",colorbar_len = 0.15,colorbar_yanchor='top', colorbar_xpos=1.3, colorbar_ypos=0.95,key.title = "Color Key")%>%
          layout(height=1000)}
      
      else if(input$heat_cluster == 'none'){
        heatmaply(mat_data,
                  main = title,
                  trace = "none",margins = c(50, 50, 150, 40),
                  dendrogram = input$heat_cluster,
                  colors = colorscheme,limits = c(-limit, limit),
                  column_text_angle = 90,fontsize_row = 12, fontsize_col = 12,
                  plot_method = "plotly",colorbar_len = 0.15,colorbar_yanchor='top', colorbar_xpos=1.3, colorbar_ypos=0.95,key.title = "Color Key")%>%
          layout(height=1000)}}
    
    
    # If column have value in between 75 to 149, the heatmap will be medium in size
    
    else if((length(total_row) >= 75)&&(length(total_row) < 150)){
      # display d3 heatmap
      if(input$heat_cluster == 'both'){
        heatmaply(mat_data,
                  main = title,
                  trace = "none",
                  Rowv = rev(original_heatmap$rowDendrogram),
                  Colv = rev(original_heatmap$colDendrogram),
                  dendrogram = input$heat_cluster,
                  colors = colorscheme,limits = c(-limit, limit),
                  row_dend_left = TRUE,
                  subplot_heights=c(0.1, 0.90),subplot_widths=c(0.8, 0.2),column_text_angle = 90,fontsize_row = 12, fontsize_col = 12,
                  plot_method = "plotly",colorbar_len = 0.07,colorbar_yanchor='top', colorbar_xpos=1.2, colorbar_ypos=0.90,key.title = "Color Key")%>%
          layout(height=2500)}
      
      else if(input$heat_cluster == 'row'){
        heatmaply(mat_data,
                  main = title,
                  trace = "none",margins = c(50, 50, 150, 40),
                  Rowv = rev(original_heatmap$rowDendrogram),
                  dendrogram = input$heat_cluster,
                  colors = colorscheme,limits = c(-limit, limit),
                  row_dend_left = TRUE,
                  subplot_widths=c(0.8, 0.2),column_text_angle = 90,fontsize_row = 12, fontsize_col = 12,
                  plot_method = "plotly",colorbar_len = 0.07,colorbar_yanchor='top', colorbar_xpos=1.0, colorbar_ypos=0.95,key.title = "Color Key")%>%
          layout(height=2500)}
      
      
      else if(input$heat_cluster == 'column'){
        heatmaply(mat_data,
                  main = title,
                  trace = "none",
                  Colv = rev(original_heatmap$colDendrogram),
                  dendrogram = input$heat_cluster,
                  colors = colorscheme,limits = c(-limit, limit),
                  subplot_heights=c(0.06, 0.93),column_text_angle = 90,fontsize_row = 12, fontsize_col = 12,
                  plot_method = "plotly",colorbar_len = 0.07,colorbar_yanchor='top', colorbar_xpos=1.3, colorbar_ypos=0.90,key.title = "Color Key")%>%
          layout(height=2500)}
      
      else if(input$heat_cluster == 'none'){
        heatmaply(mat_data,
                  main = title,
                  trace = "none",margins = c(50, 50, 150, 40),
                  dendrogram = input$heat_cluster,
                  colors = colorscheme,limits = c(-limit, limit),
                  column_text_angle = 90,fontsize_row = 12, fontsize_col = 12,
                  plot_method = "plotly",colorbar_len = 0.07,colorbar_yanchor='top', colorbar_xpos=1.3, colorbar_ypos=0.95,key.title = "Color Key")%>%
          layout(height=2500)}}
    
    
    # If column have value in between 150 to 299, the heatmap will be large in size
    
    else if((length(total_row) >= 150)&&(length(total_row) < 300)){
      if(input$heat_cluster == 'both'){
        heatmaply(mat_data,
                  main = title,
                  trace = "none",
                  Rowv = rev(original_heatmap$rowDendrogram),
                  Colv = rev(original_heatmap$colDendrogram),
                  dendrogram = input$heat_cluster,
                  colors = colorscheme,limits = c(-limit, limit),
                  row_dend_left = TRUE,
                  subplot_heights=c(0.05, 0.90),subplot_widths=c(0.8, 0.2),column_text_angle = 90,fontsize_row = 12, fontsize_col = 12,
                  plot_method = "plotly",colorbar_len = 0.05,colorbar_yanchor='top', colorbar_xpos=1.2, colorbar_ypos=0.95,key.title = "Color Key")%>%
          layout(height=3500)}
      
      else if(input$heat_cluster == 'row'){
        heatmaply(mat_data,
                  main = title,
                  trace = "none",margins = c(50, 50, 150, 40),
                  Rowv = rev(original_heatmap$rowDendrogram),
                  dendrogram = input$heat_cluster,
                  colors = colorscheme,limits = c(-limit, limit),
                  row_dend_left = TRUE,
                  subplot_widths=c(0.8, 0.2),column_text_angle = 90,fontsize_row = 12, fontsize_col = 12,
                  plot_method = "plotly",colorbar_len = 0.05,colorbar_yanchor='top', colorbar_xpos=1.2, colorbar_ypos=0.99,key.title = "Color Key")%>%
          layout(height=3500)}
      
      
      else if(input$heat_cluster == 'column'){
        heatmaply(mat_data,
                  main = title,
                  trace = "none",
                  Colv = rev(original_heatmap$colDendrogram),
                  dendrogram = input$heat_cluster,
                  colors = colorscheme,limits = c(-limit, limit),
                  subplot_heights=c(0.05, 0.90),column_text_angle = 90,fontsize_row = 12, fontsize_col = 12,
                  plot_method = "plotly",colorbar_len = 0.05,colorbar_yanchor='top', colorbar_xpos=1.2, colorbar_ypos=0.95,key.title = "Color Key")%>%
          layout(height=3500)}
      
      else if(input$heat_cluster == 'none'){
        heatmaply(mat_data,
                  main = title,
                  trace = "none",margins = c(50, 50, 150, 40),
                  dendrogram = input$heat_cluster,
                  colors = colorscheme,limits = c(-limit, limit),
                  column_text_angle = 90,fontsize_row = 12, fontsize_col = 12,
                  plot_method = "plotly",colorbar_len = 0.07,colorbar_yanchor='top', colorbar_xpos=1.2, colorbar_ypos=0.99,key.title = "Color Key")%>%
          layout(height=2500)}}
    
    
    # If column have value more than >=300, the heatmap will be very large in size
    
    else if(length(total_row) >= 300){
      if(input$heat_cluster == 'both'){
        heatmaply(mat_data,
                  main = title,
                  trace = "none",
              Rowv = rev(original_heatmap$rowDendrogram),
              Colv =  rev(original_heatmap$colDendrogram),
                  dendrogram = input$heat_cluster,
                  colors = colorscheme,limits = c(-limit, limit),
              row_dend_left = TRUE,
                  subplot_heights=c(0.05, 0.90),subplot_widths=c(0.8, 0.2),column_text_angle = 90,fontsize_row = 12, fontsize_col = 12,
                  plot_method = "plotly",colorbar_len = 0.04,colorbar_yanchor='top', colorbar_xpos=10, colorbar_ypos=0.95,key.title = "Color Key")%>%
          layout(height=4000)}
      
      else if(input$heat_cluster == 'row'){
        heatmaply(mat_data,
                  main = title,
                  trace = "none",margins = c(50, 50, 150, 40),
                  Colv = rev(original_heatmap$colDendrogram),
                  dendrogram = input$heat_cluster,
                  colors = colorscheme,limits = c(-limit, limit),
                  row_dend_left = TRUE,
                  subplot_widths=c(0.8, 0.2),column_text_angle = 90,fontsize_row = 12, fontsize_col = 12,
                  plot_method = "plotly",colorbar_len = 0.04,colorbar_yanchor='top', colorbar_xpos=1.2, colorbar_ypos=0.95,key.title = "Color Key")%>%
          layout(height=4000)}
      
      
      else if(input$heat_cluster == 'column'){
        heatmaply(mat_data,
                  main = title,
                  trace = "none",
                  Colv = rev(original_heatmap$colDendrogram),
                  dendrogram = input$heat_cluster,
                  colors = colorscheme,limits = c(-limit, limit),
                  subplot_heights=c(0.06, 0.93),column_text_angle = 90,fontsize_row = 12, fontsize_col = 12,
                  plot_method = "plotly",colorbar_len = 0.04,colorbar_yanchor='top', colorbar_xpos=1.2, colorbar_ypos=0.90,key.title = "Color Key")%>%
          layout(height=4000, width=1100)}
      
      else if(input$heat_cluster == 'none'){
        heatmaply(mat_data,
                  main = title,
                  trace = "none",margins = c(50, 50, 150, 40),
                  dendrogram = input$heat_cluster,
                  colors = colorscheme,limits = c(-limit, limit),
                  column_text_angle = 90,fontsize_row = 12, fontsize_col = 12,
                  plot_method = "plotly",colorbar_len = 0.04,colorbar_yanchor='top', colorbar_xpos=1.3, colorbar_ypos=0.99,key.title = "Color Key")%>%
          layout(height=4000)}}

    })
  
  

  ###########################
  # 4.10 OVERLAYED TIMEPLOT #
  ###########################
  
  
  # timeplot from input
  output$drugs_timePlot <- renderPlotly({
    
    # get drugs data
    data <- values$drugsData
    
   
    # filter by tissue origin
    data <- subset(data, Tissue.Origin == toupper(input$drugs_inTissue))
    
    # get dose level
    dosetitle <- switch(input$drugs_inDosage,
                   "Low" = "low",
                   "High" = "high")
    
    # plot title
    title <- paste("Gene expression of", input$drugs_inGene, "following exposure to", dosetitle,"concentration \nof the selected drugs on", input$drugs_inTissue, "cell lines", sep = " ")
    
    # format plot    
    ggformatted <- ggplot(data) +
      
      # labels
      labs(title = title, x = "Time", y = "Gene Expression (log2)") + 
      
      # text size
      theme(plot.title=element_text(size=12, hjust = 0.5),
            axis.text.x=element_text(size=10), 
            axis.text.y=element_text(size=10),
            axis.title.x=element_text(size=12),
            axis.title.y=element_text(size=12),
            legend.text=element_text(size= 10),
            legend.title=element_text(size= 12),
            plot.margin=unit(c(1.5,0,0,0),"cm")) + 
      
      # coloring
      scale_color_discrete(name="Drug") + 
      
      # force scale order
      scale_x_discrete(limits=c("2 hrs","6 hrs","24 hrs"))
    
    # plot lines
    plot <- ggformatted + geom_line(
      aes(x = Time, 
          y = Gene.Expression, 
          color = Drug,
          group = Cell.Line), size = 1)
    
    ggplotly(plot,
             height = 700)})
  
  
  
  # generate timeplots for specified tissue
  overlaygraph <- function(tissue, legend){
    
    # get drugs data
    data <- values$drugsData
    
    # filter by tissue origin
    data <- subset(data, Tissue.Origin == toupper(tissue))
    
    # plot title
    title <- paste(tissue, "cell lines", sep = " ")
    
    # format plot    
    ggformatted <- ggplot(data) +
      
      # labels
      labs(title = title, x = "Time", y = "Gene Expression (log2)") + 
      
      # text size
      theme(plot.title=element_text(size=12, hjust = 0.5),
            axis.text.x=element_text(size=8), 
            axis.text.y=element_text(size=8),
            axis.title.x=element_text(size=10),
            axis.title.y=element_text(size=10),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            legend.position = legend) + 
      
      # coloring
      scale_color_discrete(name= paste0("Drug (",input$drugs_inDosage,")")) + 
      
      # force scale order
      scale_x_discrete(limits=c("2 hrs","6 hrs","24 hrs"))
    
    # plot lines
    plot <- ggformatted + geom_line(
      aes(x = Time, 
          y = Gene.Expression,
          color = Drug,
          group = interaction(Cell.Line, Drug)), size = 1)
    
    plot}
  
  
  output$drugs_gridPlot <- renderPlot({
    
    # generate graphs
    p1 <- overlaygraph("Breast", "none")
    p2 <- overlaygraph("CNS", "none")
    p3 <- overlaygraph("Colon", "none")
    p4 <- overlaygraph("Leukemia", "none")
    p5 <- overlaygraph("Lung", "none")
    p6 <- overlaygraph("Melanoma", "none")
    p7 <- overlaygraph("Prostate", "none")
    p8 <- overlaygraph("Ovarian", "none")
    p9 <- overlaygraph("Renal", "none")
    
    # grab legend
    tmp <- ggplot_gtable(ggplot_build(overlaygraph("Breast", "right")))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    
    # arrange graphs
    graphs <- arrangeGrob(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3)
    
    # add legend
    grid.arrange(legend, graphs, ncol = 2, widths = c(2, 10))}, 
      height = 700)
  
  
  #############
  # 4.11 MISC #
  #############
  
  
  # display time and dosage selection for upset/venn titles
  output$upset_title <- renderUI({tagList(p(paste0("UpSet plot (time point: ",input$time_point,", drug concentration: ", input$drug_conc,")")))})
  output$venn_title <- renderUI({tagList(p(paste0("Venn Diagram (time point: ",input$time_point,", drug concentration: ", input$drug_conc,")")))})
  output$upset_title2 <- renderUI({tagList(p(paste0("UpSet plot (time point: ",input$time_point,", drug concentration: ", input$drug_conc,")")))})
  output$venn_title2 <- renderUI({tagList(p(paste0("Venn Diagram (time point: ",input$time_point,", drug concentration: ", input$drug_conc,")")))})
  
  output$genesNotFound <- renderUI(
    if(nchar(values$react_val_genesNotFound)>0) {
      {tagList(p(paste0("Query genes not in NCI-TPW: ",values$react_val_genesNotFound)))}
    } else {
      return(NULL)
    }
  )   
  
  output$packageVersion <- renderUI({
    x <- list(tags$p(paste0("DT: tested on v0.10, currently using v",packageVersion("DT"))),
              tags$p(paste0("shiny: tested on v1.4.0, currently using v",packageVersion("shiny"))),
              tags$p(paste0("ggplot2: tested on v3.2.1, currently using v",packageVersion("ggplot2"))),
              tags$p(paste0("data.table: tested on v1.12.6, currently using v",packageVersion("data.table"))),
              tags$p(paste0("shinyWidgets: tested on v0.5.0, currently using v",packageVersion("shinyWidgets"))),
              tags$p(paste0("heatmaply: tested on v1.0.0, currently using v",packageVersion("heatmaply"))),
              tags$p(paste0("plotly: tested on v4.9.1, currently using v",packageVersion("plotly"))),
              tags$p(paste0("gridExtra: tested on v2.3, currently using v",packageVersion("gridExtra"))),
              tags$p(paste0("UpSetR: tested on v1.4.0, currently using v",packageVersion("UpSetR"))),
              tags$p(paste0("gplots: tested on v3.0.1.2, currently using v",packageVersion("gplots"))),
              tags$p(paste0("dendextend: tested on v1.13.2, currently using v",packageVersion("dendextend"))),
              tags$p(paste0("devtools: tested on v2.2.1, currently using v",packageVersion("devtools"))),
              tags$p(paste0("venn: tested on v1.7, currently using v",packageVersion("venn")))
          )
      tagList(x)
  })
  
  # get gene URL
  output$geneURL <- renderUI({
    
    # define link
    url <- a("GeneCard information", href = paste(geneStub, input$inGene, sep= ""))
    
    # return hyperlink
    tagList(url)})
  
  # get gene URL for drugs overlay
  output$drugs_geneURL <- renderUI({
    
    # define link
    url <- a("GeneCard information", href = paste(geneStub, input$inGene, sep= ""))
    
    # return hyperlink
    tagList(url)})
  
  
  output$negative_table = DT::renderDataTable({
    
    timepoint = input$time_point
    conclevel = input$drug_conc
    
    negative_list_selected = c()
    for(dr in input$drug_variable){
      drugname = tolower(dr)
      negative_list_selected[[dr]] = negative_correlation[[timepoint]][[conclevel]][[drugname]]
    }
    
    if(length(input$drug_variable) >= 2) {
      ItemsList <- gplots::venn(negative_list_selected, show.plot = FALSE)
      genes_intersections <- list.as.matrix(attributes(ItemsList)$intersections)
      resorted_colnames = colnames(genes_intersections)
      genes_intersections <- genes_intersections[ ,order(resorted_colnames) ]
    } else {
      if(length(input$drug_variable) < 1) {
        stop("\nPlease select at least one drug from the list in the side panel.")
      }
      
      genes_intersections = negative_list_selected[[input$drug_variable]]
      genes_intersections = as.matrix(genes_intersections)
      colnames(genes_intersections) = input$drug_variable
    }
    
    datatable(genes_intersections, 
              class = 'cell-border stripe',
              extensions = list("Buttons" = NULL,
                                "ColReorder"=NULL),
              selection = list(target = "column"), 
              options = list(dom = 'Bt',
                             buttons = list(list(extend = 'csv', 
                                                 filename= paste0('TPW_logGI50_negcorrelatedgenes','_',timepoint,'_',conclevel)),
                                            list(extend = 'excel', 
                                                 filename= paste0('TPW_logGI50_negcorrelatedgenes','_',timepoint,'_',conclevel))),
                             colReorder = TRUE,
                             ordering = FALSE, 
                             searching = FALSE, 
                             pageLength = 100)
    )
    
  })
  
  
  output$positive_table = DT::renderDataTable({
    
    timepoint = input$time_point
    conclevel = input$drug_conc
    
    positive_list_selected = c()
    for(dr in input$drug_variable){
      drugname = tolower(dr)
      positive_list_selected[[dr]] = positive_correlation[[timepoint]][[conclevel]][[drugname]]
    }
    
    if(length(input$drug_variable) >= 2) {
      ItemsList <- gplots::venn(positive_list_selected, show.plot = FALSE)
      
      genes_intersections <- list.as.matrix(attributes(ItemsList)$intersections)
      
      resorted_colnames = colnames(genes_intersections)
      genes_intersections <- genes_intersections[ ,order(resorted_colnames) ]
    } else {
      if(length(input$drug_variable) < 1) {
        stop("\nPlease select at least one drug from the list in the side panel.")
      }
      genes_intersections = positive_list_selected[[input$drug_variable]]
      genes_intersections = as.matrix(genes_intersections)
      colnames(genes_intersections) = input$drug_variable
    }
    
    datatable(genes_intersections, 
              class = 'cell-border stripe',
              extensions = list("Buttons" = NULL,
                                "ColReorder"=NULL),
              selection = list(target = "column"), 
              options = list(dom = 'Bt',
                             buttons = list(list(extend = 'csv', 
                                                 filename= paste0('TPW_logGI50_poscorrelatedgenes','_',timepoint,'_',conclevel)),
                                            list(extend = 'excel', 
                                                 filename= paste0('TPW_logGI50_poscorrelatedgenes','_',timepoint,'_',conclevel))),
                             colReorder = TRUE,
                             ordering = FALSE, 
                             searching = FALSE, 
                             pageLength = 100)
    )
    
  })
  
  
  
  output$venn1 <- renderPlot({
    
    if(length(input$drug_variable) < 1) {
      stop("\nPlease select at least one drug from the list in the side panel.")
    }
    
    timepoint = input$time_point
    conclevel = input$drug_conc
    
    negative_list_selected = c()
    for(dr in input$drug_variable){
      drugname = tolower(dr)
      negative_list_selected[[dr]] = negative_correlation[[timepoint]][[conclevel]][[drugname]]
    }
    
    if(length(input$drug_variable) > 7) {
      venn::venn(8, ilab=TRUE, zcolor = "style")
    } else {
      venn::venn(negative_list_selected, borders=FALSE, ilabels=TRUE, zcolor = "style", ilcs = 1.2, sncs=1.1, opacity = 0.3, cexil = 1.5, cexsn = 1.3)
      zeroset <- matrix(1000*c(0,1,1,0,0,0,0,1,1,0), ncol = 2)
      lines(zeroset, col='white', lwd=5)
    }
    
    
  })
  
  
  output$venn2 <- renderPlot({
    
    if(length(input$drug_variable) < 1) {
      stop("\nPlease select at least one drug from the list in the side panel.")
    }

    timepoint = input$time_point
    conclevel = input$drug_conc
    
    positive_list_selected = c()
    for(dr in input$drug_variable){
      drugname = tolower(dr)
      positive_list_selected[[dr]] = positive_correlation[[timepoint]][[conclevel]][[drugname]]
    }
    
    
    if(length(input$drug_variable) > 7) {
      venn::venn(8, ilab=TRUE, zcolor = "style")
    } else {
      venn::venn(positive_list_selected, ilabels=TRUE, zcolor = "style", ilcs = 1.2, sncs=1.1, opacity = 0.3, cexil = 1.5, cexsn = 1.3)
      zeroset <- matrix(1000*c(0,1,1,0,0,0,0,1,1,0), ncol = 2)
      lines(zeroset, col='white', lwd=5)
    }
    
    
  })
  
  
  
  output$upset1 <- renderPlot({
    
    if(length(input$drug_variable) < 2) {
      stop("\nPlease select at least two drugs from the list in the side panel.")
    }
    
    timepoint = input$time_point
    conclevel = input$drug_conc
    
    selected_list = c()
    for(dr in input$drug_variable){
      drugname = tolower(dr)
      selected_list[[dr]] = correlation_all[[timepoint]][[conclevel]][[drugname]]
    }
    
    selected_names = input$drug_variable;
    
    all_genes = Reduce(union, selected_list)
    
    fig1_combined=matrix(ncol=length(selected_list),nrow=length(all_genes))
    
    rownames(fig1_combined) =  all_genes
    colnames(fig1_combined) =  selected_names
    
    
    negative_list_selected = c()
    for(dr in input$drug_variable){
      drugname = tolower(dr)
      negative_list_selected[[dr]] = negative_correlation[[timepoint]][[conclevel]][[drugname]]
    }
    
    fig1_combined[,]=0
    for(i in 1:length(selected_list)) {
      fig1_combined[, i] <- all_genes %in% unlist(negative_list_selected[i])
    }
    
    fig1_dataframe = data.frame(fig1_combined[,])
    
    
    upset(fig1_dataframe, sets = colnames(fig1_combined), 
          point.size = 3.0,  line.size = 0.9,
          mainbar.y.label = "", 
          sets.x.label = "",
          sets.bar.color = "black",
          main.bar.color = "grey30",
          order.by = "freq", 
          text.scale = 2, mb.ratio = c(0.55, 0.45),
          empty.intersections = NULL)
    
  },height = 450)
  
  
  
  
  output$upset2 <- renderPlot({
    
    
    if(length(input$drug_variable) < 2) {
      stop("\nPlease select at least two drugs from the list in the side panel.")
    }
    
    timepoint = input$time_point
    conclevel = input$drug_conc
    
    selected_list = c()
    for(dr in input$drug_variable){
      drugname = tolower(dr)
      selected_list[[dr]] = correlation_all[[timepoint]][[conclevel]][[drugname]]
    }
    
    
    selected_names = input$drug_variable;
    
    all_genes = Reduce(union, selected_list)
    
    fig1_combined=matrix(ncol=length(selected_list),nrow=length(all_genes))
    
    rownames(fig1_combined) =  all_genes
    colnames(fig1_combined) =  selected_names
    
    
    positive_list_selected = c()
    for(dr in input$drug_variable){
      drugname = tolower(dr)
      positive_list_selected[[dr]] = positive_correlation[[timepoint]][[conclevel]][[drugname]]
    }
    
    fig1_combined[,]=0
    for(i in 1:length(selected_list)) {
      fig1_combined[, i] <- all_genes %in% unlist(positive_list_selected[i])
    }
    
    fig1_dataframe = data.frame(fig1_combined[,])
    
    
    upset(fig1_dataframe, sets = colnames(fig1_combined), 
          point.size = 3.0,  line.size = 0.9,
          mainbar.y.label = "", 
          sets.x.label = "",
          sets.bar.color = "black",
          main.bar.color = "grey30",
          order.by = "freq", 
          text.scale = 2, mb.ratio = c(0.55, 0.45),
          empty.intersections = NULL)
    
  },height = 450)
  
  
} #close of server






############
# 5 LAUNCH #
############




}


# launch app
shinyApp(ui = ui, server = server, options = list(port=7990))
