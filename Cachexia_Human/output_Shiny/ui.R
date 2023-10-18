library(shiny) 
library(shinyhelper) 
library(data.table) 
library(Matrix) 
library(DT) 
library(magrittr) 
sc1conf = readRDS("sc1conf.rds")
sc1def  = readRDS("sc1def.rds")



### Start server code 
shinyUI(fluidPage( 
### HTML formatting of error messages 
 
tags$head(tags$style(HTML(".shiny-output-error-validation {color: red; font-weight: bold;}"))), 
list(tags$style(HTML(".navbar-default .navbar-nav { font-weight: bold; font-size: 16px; }"))), 
 
   
### Page title 
titlePanel("Pryce etal SingleCell Human Cachexia Positive"),  
navbarPage( 
  NULL,  
 ### Tab1.a1: cellInfo vs geneExpr on dimRed 
  tabPanel( 
    HTML("CellInfo vs GeneExpr"), 
    h4("Cell information vs gene expression on reduced dimensions"), 
    "In this tab, users can visualise both cell information and gene ",  
    "expression side-by-side on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("sc1a1drX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                           selected = sc1def$dimred[1]), 
            selectInput("sc1a1drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                        selected = sc1def$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("sc1a1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc1a1togL % 2 == 1", 
          selectInput("sc1a1sub1", "Cell information to subset:", 
                      choices = sc1conf[grp == TRUE]$UI, 
                      selected = sc1def$grp1), 
          uiOutput("sc1a1sub1.ui"), 
          actionButton("sc1a1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc1a1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("sc1a1tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc1a1tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("sc1a1siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("sc1a1psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("sc1a1fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("sc1a1asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("sc1a1txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information"), 
        fluidRow( 
          column( 
            6, selectInput("sc1a1inp1", "Cell information:", 
                           choices = sc1conf$UI, 
                           selected = sc1def$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc1a1tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc1a1tog1 % 2 == 1", 
              radioButtons("sc1a1col1", "Colour (Continuous data):", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("sc1a1ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("sc1a1lab1", "Show cell info labels", value = TRUE), 
              checkboxInput("sc1a1legend", "Show legend", value = FALSE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc1a1oup1.ui"))), 
        downloadButton("sc1a1oup1.pdf", "Download PDF"), 
        downloadButton("sc1a1oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc1a1oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc1a1oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)), br(), 
        actionButton("sc1a1tog9", "Toggle to show cell numbers / statistics"), 
        conditionalPanel( 
          condition = "input.sc1a1tog9 % 2 == 1", 
          h4("Cell numbers / statistics"), 
          radioButtons("sc1a1splt", "Split continuous cell info into:", 
                       choices = c("Quartile", "Decile"), 
                       selected = "Decile", inline = TRUE), 
          dataTableOutput("sc1a1.dt") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression"), 
        fluidRow( 
          column( 
            6, selectInput("sc1a1inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc1a1tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc1a1tog2 % 2 == 1", 
              radioButtons("sc1a1col2", "Colour:", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("sc1a1ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ) , 
        fluidRow(column(12, uiOutput("sc1a1oup2.ui"))), 
        downloadButton("sc1a1oup2.pdf", "Download PDF"), 
        downloadButton("sc1a1oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc1a1oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc1a1oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
 
 
  ### Tab1.a2: cellInfo vs cellInfo on dimRed 
  tabPanel( 
    HTML("CellInfo vs CellInfo"), 
    h4("Cell information vs cell information on dimension reduction"), 
    "In this tab, users can visualise two cell informations side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("sc1a2drX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                           selected = sc1def$dimred[1]), 
            selectInput("sc1a2drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                        selected = sc1def$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("sc1a2togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc1a2togL % 2 == 1", 
          selectInput("sc1a2sub1", "Cell information to subset:", 
                      choices = sc1conf[grp == TRUE]$UI, 
                      selected = sc1def$grp1), 
          uiOutput("sc1a2sub1.ui"), 
          actionButton("sc1a2sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc1a2sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("sc1a2tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc1a2tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("sc1a2siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("sc1a2psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("sc1a2fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("sc1a2asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("sc1a2txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information 1"), 
        fluidRow( 
          column( 
            6, selectInput("sc1a2inp1", "Cell information:", 
                           choices = sc1conf$UI, 
                           selected = sc1def$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc1a2tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc1a2tog1 % 2 == 1", 
              radioButtons("sc1a2col1", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("sc1a2ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("sc1a2lab1", "Show cell info labels", value = TRUE), 
              checkboxInput("sc1a2legend1", "Show legend", value = FALSE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc1a2oup1.ui"))), 
        downloadButton("sc1a2oup1.pdf", "Download PDF"), 
        downloadButton("sc1a2oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc1a2oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc1a2oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Cell information 2"), 
        fluidRow( 
          column( 
            6, selectInput("sc1a2inp2", "Cell information:", 
                           choices = sc1conf$UI, 
                           selected = sc1def$meta2) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc1a2tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc1a2tog2 % 2 == 1", 
              radioButtons("sc1a2col2", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("sc1a2ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("sc1a2lab2", "Show cell info labels", value = TRUE), 
              checkboxInput("sc1a2legend2", "Show legend", value = FALSE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc1a2oup2.ui"))), 
        downloadButton("sc1a2oup2.pdf", "Download PDF"), 
        downloadButton("sc1a2oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc1a2oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc1a2oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
 
 
  ### Tab1.a3: geneExpr vs geneExpr on dimRed 
  tabPanel( 
    HTML("GeneExpr vs GeneExpr"), 
    h4("Gene expression vs gene expression on dimension reduction"), 
    "In this tab, users can visualise two gene expressions side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("sc1a3drX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                           selected = sc1def$dimred[1]), 
            selectInput("sc1a3drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                        selected = sc1def$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("sc1a3togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc1a3togL % 2 == 1", 
          selectInput("sc1a3sub1", "Cell information to subset:", 
                      choices = sc1conf[grp == TRUE]$UI, 
                      selected = sc1def$grp1), 
          uiOutput("sc1a3sub1.ui"), 
          actionButton("sc1a3sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc1a3sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("sc1a3tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc1a3tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("sc1a3siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("sc1a3psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("sc1a3fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("sc1a3asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("sc1a3txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Gene expression 1"), 
        fluidRow( 
          column( 
            6, selectInput("sc1a3inp1", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc1a3tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc1a3tog1 % 2 == 1", 
              radioButtons("sc1a3col1", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("sc1a3ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc1a3oup1.ui"))), 
        downloadButton("sc1a3oup1.pdf", "Download PDF"), 
        downloadButton("sc1a3oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc1a3oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc1a3oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression 2"), 
        fluidRow( 
          column( 
            6, selectInput("sc1a3inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc1a3tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc1a3tog2 % 2 == 1", 
              radioButtons("sc1a3col2", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("sc1a3ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc1a3oup2.ui"))), 
        downloadButton("sc1a3oup2.pdf", "Download PDF"), 
        downloadButton("sc1a3oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc1a3oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc1a3oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
 ### Tab1.b2: Gene coexpression plot 
 tabPanel( 
    HTML("Gene Dual Coex"), 
    h4("Coexpression of two genes on reduced dimensions"), 
    "In this tab, users can visualise the coexpression of two genes ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
    column( 
      3, h4("Dimension Reduction"), 
          selectInput("sc1b2drX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                      selected = sc1def$dimred[1]), 
          selectInput("sc1b2drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                      selected = sc1def$dimred[2]), 
          actionButton("sc1b2tog1", "Toggle plot controls"), 
          conditionalPanel( 
            condition = "input.sc1b2tog1 % 2 == 1", 
            radioButtons("sc1b2col1", "Colour:", 
                         choices = c("Red (Gene1); Blue (Gene2)", 
                                     "Orange (Gene1); Blue (Gene2)", 
                                     "Red (Gene1); Green (Gene2)", 
                                     "Green (Gene1); Blue (Gene2)"), 
                         selected = "Red (Gene1); Blue (Gene2)"), 
            radioButtons("sc1b2ord1", "Plot order:", 
                         choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                         selected = "Max-1st", inline = TRUE) 
    ), 
    br(),br(), 
    actionButton("sc1b2togL", "Toggle to subset cells"), 
    conditionalPanel( 
      condition = "input.sc1b2togL % 2 == 1", 
      selectInput("sc1b2sub1", "Cell information to subset:", 
                  choices = sc1conf[grp == TRUE]$UI, 
                  selected = sc1def$grp1), 
      uiOutput("sc1b2sub1.ui"), 
      actionButton("sc1b2sub1all", "Select all groups", class = "btn btn-primary"), 
      actionButton("sc1b2sub1non", "Deselect all groups", class = "btn btn-primary") 
    ), 
    br(),br(), 
    actionButton("sc1b2tog0", "Toggle graphics controls"), 
    conditionalPanel( 
      condition = "input.sc1b2tog0 % 2 == 1", 
      sliderInput("sc1b2siz", "Point size:", 
                  min = 0, max = 4, value = 1.25, step = 0.25), 
      radioButtons("sc1b2psz", "Plot size:", 
                   choices = c("Small", "Medium", "Large"), 
                   selected = "Medium", inline = TRUE), 
      radioButtons("sc1b2fsz", "Font size:", 
                   choices = c("Small", "Medium", "Large"), 
                   selected = "Small", inline = TRUE), 
      radioButtons("sc1b2asp", "Aspect ratio:", 
                   choices = c("Square", "Fixed", "Free"), 
                   selected = "Square", inline = TRUE), 
      checkboxInput("sc1b2txt", "Show axis text", value = FALSE) 
    ) 
  ), # End of column (6 space) 
  
  column( 
    4,style="border-left: 2px solid black", h4("Gene Expression 1"), 
    selectInput("sc1b2inp1", "Gene 1:", choices=NULL) %>%  
      helper(type = "inline", size = "m", fade = TRUE, 
             title = "Gene expression to colour cells by", 
             content = c("Select gene to colour cells by gene expression", 
                         paste0("- Gene expression are coloured in a ", 
                                "White-Red colour scheme which can be ", 
                                "changed in the plot controls"))), 
    selectInput("sc1b2inp2", "Gene 2:", choices=NULL) %>% 
      helper(type = "inline", size = "m", fade = TRUE, 
             title = "Gene expression to colour cells by", 
             content = c("Select gene to colour cells by gene expression", 
                         paste0("- Gene expression are coloured in a ", 
                                "White-Blue colour scheme which can be ", 
                                "changed in the plot controls"))), 
    #style="border-right: 2px solid black", 
    uiOutput("sc1b2oup1.ui"), 
    downloadButton("sc1b2oup1.pdf", "Download PDF"), 
    downloadButton("sc1b2oup1.png", "Download PNG"), br(), 
    div(style="display:inline-block", 
        numericInput("sc1b2oup1.h", "PDF / PNG height:", width = "138px", 
                     min = 4, max = 20, value = 8, step = 0.5)), 
    div(style="display:inline-block", 
    numericInput("sc1b2oup1.w", "PDF / PNG width:", width = "138px", 
                 min = 4, max = 20, value = 10, step = 0.5)), 
    uiOutput("sc1b2oup2.ui"), 
    downloadButton("sc1b2oup2.pdf", "Download PDF"), 
    downloadButton("sc1b2oup2.png", "Download PNG"), 
    br(), h4("Cell numbers"), 
    dataTableOutput("sc1b2.dt") 
  ), # End of column (6 space) 
  
  column( 
    4, style="border-left: 2px solid black", h4("Gene Expression 2"), 
    selectInput("sc1b2inp3", "Gene 1:", choices=NULL) %>% 
      helper(type = "inline", size = "m", fade = TRUE, 
             title = "Gene expression to colour cells by", 
             content = c("Select gene to colour cells by gene expression", 
                         paste0("- Gene expression are coloured in a ", 
                         "White-Red colour scheme which can be ", 
                         "changed in the plot controls"))), 
    selectInput("sc1b2inp4", "Gene 2:", choices=NULL) %>% 
      helper(type = "inline", size = "m", fade = TRUE, 
             title = "Gene expression to colour cells by", 
             content = c("Select gene to colour cells by gene expression", 
                         paste0("- Gene expression are coloured in a ", 
                                "White-Blue colour scheme which can be ", 
                                "changed in the plot controls"))), 
    #style="border-right: 2px solid black", 
    uiOutput("sc1b2oup3.ui"), 
    downloadButton("sc1b2oup3.pdf", "Download PDF"), 
    downloadButton("sc1b2oup3.png", "Download PNG"), br(), 
    div(style="display:inline-block", 
        numericInput("sc1b2oup3.h", "PDF / PNG height:", width = "138px", 
                     min = 4, max = 20, value = 8, step = 0.5)), 
    div(style="display:inline-block", 
    numericInput("sc1b2oup3.w", "PDF / PNG width:", width = "138px", 
                 min = 4, max = 20, value = 10, step = 0.5)), 
    uiOutput("sc1b2oup4.ui"), 
    downloadButton("sc1b2oup4.pdf", "Download PDF"), 
    downloadButton("sc1b2oup4.png", "Download PNG"), 
    br(), h4("Cell numbers"), 
    dataTableOutput("sc1b2_2.dt") 
  ) 
 ) 
 ),     # End of tab (2 space) 
 
 
 
### Tab1.a1b: SPLIT CELL INFO 
 tabPanel( 
   HTML("Split Dataset - Cell Info"), 
   h4("Cell information vs gene expression on reduced dimensions"), 
      "In this tab, users can visualise both cell information and gene ", 
      "expression side-by-side on low-dimensional representions.", 
   br(),br(), 
   fluidRow( 
     column( 
       3, h4("Dimension Reduction"), 
       fluidRow( 
         column( 
           12, 
           selectInput("sc1a1bdrX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                           selected = sc1def$dimred[1]), 
           selectInput("sc1a1bdrY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                       selected = sc1def$dimred[2])) 
         ) 
       ), # End of column (6 space) 
       column( 
         3, 
         actionButton("sc1a1btogL", "Toggle to subset cells"), 
         conditionalPanel( 
           condition = "input.sc1a1btogL % 2 == 1", 
           selectInput("sc1a1bsub1", "Cell information to subset:", 
                       choices = sc1conf[grp == TRUE]$UI, 
                       selected = sc1def$grp1), 
           uiOutput("sc1a1bsub1.ui"), 
           actionButton("sc1a1bsub1all", "Select all groups", class = "btn btn-primary"), 
           actionButton("sc1a1bsub1non", "Deselect all groups", class = "btn btn-primary") 
         ) 
       ), # End of column (6 space) 
       column( 
       6, 
       actionButton("sc1a1btog0", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.sc1a1btog0 % 2 == 1", 
         fluidRow( 
           column( 
             6, 
             sliderInput("sc1a1bsiz", "Point size:", 
                         min = 0, max = 4, value = 1.25, step = 0.25), 
             radioButtons("sc1a1bpsz", "Plot size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE), 
             radioButtons("sc1a1bfsz", "Font size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE) 
           ), 
           column( 
             6, 
             radioButtons("sc1a1basp", "Aspect ratio:", 
                          choices = c("Square", "Fixed", "Free"), 
                          selected = "Square", inline = TRUE), 
             checkboxInput("sc1a1btxt", "Show axis text", value = FALSE) 
           ) 
         ) 
       ) 
     )  # End of column (6 space) 
   ),   # End of fluidRow (4 space) 
   fluidRow( 
     column( 
       6, style="border-right: 2px solid black", 
       h4("Split dataset"), 
       fluidRow( 
         column( 
           12, 
           selectInput("sc1a1bsplit1", "Split dataset:", 
                       choices = sc1conf[sc1conf$split==TRUE]$UI), 
           uiOutput("sc1a1bsplit_select1_out1"), 
           uiOutput("sc1a1bsplit_select2_out1") 
         ) 
       ) 
     ) 
   ), 
   fluidRow( 
     column( 
       6, style="border-right: 2px solid black", 
       h4("Cell information"), 
       fluidRow( 
         column( 
           6, 
           selectInput("sc1a1binp1", "Cell information:", 
                       choices = sc1conf$UI, 
                       selected = sc1def$meta1) %>% 
                       helper(type = "inline", size = "m", fade = TRUE, 
                              title = "Cell information to colour cells by", 
                              content = c("Select cell information to colour cells", 
                                          "- Categorical covariates have a fixed colour palette", 
                                          paste0("- Continuous covariates are coloured in a ", 
                                          "Blue-Yellow-Red colour scheme, which can be ", 
                                          "changed in the plot controls"))) 
           ), 
           column( 
           6, 
           actionButton("sc1a1btog1", "Toggle plot controls"), 
           conditionalPanel( 
             condition = "input.sc1a1btog1 % 2 == 1", 
             radioButtons("sc1a1bcol1", "Colour (Continuous data):", 
                          choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                          selected = "Blue-Yellow-Red"), 
             radioButtons("sc1a1bord1", "Plot order:", 
                          choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                          selected = "Original", inline = TRUE), 
             checkboxInput("sc1a1blab1", "Show cell info labels", value = TRUE), 
             checkboxInput("sc1a1blegend", "Show legend", value = FALSE) 
           ) 
         ) 
       ), 
       fluidRow( 
         column( 
           12, 
           uiOutput("sc1a1boup1.ui") 
         ) 
       ), 
       downloadButton("sc1a1boup1.pdf", "Download PDF"), 
       downloadButton("sc1a1boup1.png", "Download PNG"), br(), 
       div(style="display:inline-block", 
         numericInput("sc1a1boup1.h", "PDF / PNG height:", width = "138px", 
                      min = 4, max = 20, value = 6, step = 0.5)), 
       div(style="display:inline-block", 
         numericInput("sc1a1boup1.w", "PDF / PNG width:", width = "138px", 
                      min = 4, max = 20, value = 8, step = 0.5)), br(), 
       actionButton("sc1a1btog9", "Toggle to show cell numbers / statistics"), 
#       conditionalPanel( 
#         condition = "input.sc1a1btog9 % 2 == 1", 
#         h4("Cell numbers / statistics"), 
#         radioButtons("sc1a1bsplt", "Split continuous cell info into:", 
#                      choices = c("Quartile", "Decile"), 
#                      selected = "Decile", inline = TRUE), 
#         dataTableOutput("sc1a1b.dt") 
#       ) 
     ), # End of column (6 space) 
     column( 
       6, 
       br(),br(), br(),br(),br(), br(), 
# fluidRow( 
#   column( 
#     6, selectInput("sc1a1binp1", "Cell information:", 
#                    choices = sc1conf$UI, 
#                    selected = sc1def$meta1) %>% 
#       helper(type = "inline", size = "m", fade = TRUE, 
#              title = "Cell information to colour cells by", 
#              content = c("Select cell information to colour cells", 
#                          "- Categorical covariates have a fixed colour palette", 
#                          paste0("- Continuous covariates are coloured in a ", 
#                                 "Blue-Yellow-Red colour scheme, which can be ", 
#                                 "changed in the plot controls"))) 
#   ), 
#   column( 
#     6, actionButton("sc1a1btog1", "Toggle plot controls"), 
#     conditionalPanel( 
#       condition = "input.sc1a1btog1 % 2 == 1", 
#       radioButtons("sc1a1bcol1", "Colour (Continuous data):", 
#                    choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
#                    selected = "Blue-Yellow-Red"), 
#       radioButtons("sc1a1bord1", "Plot order:", 
#                    choices = c("Max-1st", "Min-1st", "Original", "Random"), 
#                    selected = "Original", inline = TRUE), 
#       checkboxInput("sc1a1blab1", "Show cell info labels", value = TRUE) 
#     ) 
#   ) 
# ), 
# fluidRow( 
#   column( 
#     6, selectInput("sc1a1binp2", "Gene name:", choices=NULL) %>% 
#       helper(type = "inline", size = "m", fade = TRUE, 
#              title = "Gene expression to colour cells by", 
#              content = c("Select gene to colour cells by gene expression", 
#                          paste0("- Gene expression are coloured in a ", 
#                                 "White-Red colour scheme which can be ", 
#                                 "changed in the plot controls"))) 
#   ), 
#   column( 
#     6, actionButton("sc1a1btog2", "Toggle plot controls"), 
#     conditionalPanel( 
#       condition = "input.sc1a1tbog2 % 2 == 1", 
#       radioButtons("sc1a1cbol2", "Colour:", 
#                    choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
#                    selected = "White-Red"), 
#       radioButtons("sc1a1bord2", "Plot order:", 
#                    choices = c("Max-1st", "Min-1st", "Original", "Random"), 
#                    selected = "Max-1st", inline = TRUE) 
#     ) 
#   ) 
# ) , 
       fluidRow( 
         column( 
           12, 
           uiOutput("sc1a1boup2.ui") 
         ) 
       ), 
       downloadButton("sc1a1boup2.pdf", "Download PDF"), 
       downloadButton("sc1a1boup2.png", "Download PNG"), br(), 
       div(style="display:inline-block", 
         numericInput("sc1a1boup2.h", "PDF / PNG height:", width = "138px", 
                      min = 4, max = 20, value = 6, step = 0.5)), 
       div(style="display:inline-block", 
         numericInput("sc1a1boup2.w", "PDF / PNG width:", width = "138px", 
                      min = 4, max = 20, value = 8, step = 0.5)) 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
### Tab1.a1c: geneExpr vs geneExpr on dimRed 
 tabPanel( 
   HTML("Split dataset - Gene Expression"), 
   h4("Gene expression vs gene expression on dimension reduction"), 
      "In this tab, users can visualise two gene expressions side-by-side ", 
      "on low-dimensional representions.", 
   br(),br(), 
 
 
 
   fluidRow( 
     column( 
       3, 
       h4("Dimension Reduction"), 
       fluidRow( 
         column( 
           12, 
           selectInput("sc1a1cdrX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                       selected = sc1def$dimred[1]), 
           selectInput("sc1a1cdrY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                       selected = sc1def$dimred[2])) 
         ) 
       ), # End of column (6 space) 
       column( 
         3, 
         actionButton("sc1a1ctogL", "Toggle to subset cells"), 
         conditionalPanel( 
           condition = "input.sc1a1ctogL % 2 == 1", 
           selectInput("sc1a1csub1", "Cell information to subset:", 
                       choices = sc1conf[grp == TRUE]$UI, 
                       selected = sc1def$grp1), 
           uiOutput("sc1a1csub1.ui"), 
           actionButton("sc1a1csub1all", "Select all groups", class = "btn btn-primary"), 
           actionButton("sc1a1csub1non", "Deselect all groups", class = "btn btn-primary") 
         ) 
       ), # End of column (6 space) 
       column( 
         6, 
         actionButton("sc1a1ctog0", "Toggle graphics controls"), 
         conditionalPanel( 
           condition = "input.sc1a1ctog0 % 2 == 1", 
           fluidRow( 
             column( 
               6, 
               sliderInput("sc1a1csiz", "Point size:", 
                           min = 0, max = 4, value = 1.25, step = 0.25), 
               radioButtons("sc1a1cpsz", "Plot size:", 
                            choices = c("Small", "Medium", "Large"), 
               selected = "Medium", inline = TRUE), 
               radioButtons("sc1a1cfsz", "Font size:", 
               choices = c("Small", "Medium", "Large"), 
               selected = "Medium", inline = TRUE) 
             ), 
             column( 
               6, 
               radioButtons("sc1a1casp", "Aspect ratio:", 
                            choices = c("Square", "Fixed", "Free"), 
                            selected = "Square", inline = TRUE), 
               checkboxInput("sc1a1ctxt", "Show axis text", value = FALSE) 
             ) 
           ) 
         ) 
       )  # End of column (6 space) 
     ),   # End of fluidRow (4 space) 
 
 
 
 
     fluidRow( 
       column( 
         6, 
         #style="border-right: 2px solid black", 
         h4("Split dataset"), 
         fluidRow( 
           column( 
             12, 
             selectInput("sc1a1csplit1", "Split dataset:", 
               choices = sc1conf[sc1conf$split==TRUE]$UI), 
             uiOutput("sc1a1csplit_select1_out1"), 
             uiOutput("sc1a1csplit_select2_out1") 
           ) 
         ) 
       ) 
     ), 
     fluidRow( #within column 6 - gene expression 1 columns 
       column( 
         6, 
         #style="border-right: 2px solid black", 
         h4("Gene expression 1"), 
         fluidRow( 
           column( 
             6, 
             selectInput("sc1a1cinp1", "Gene name:", choices=NULL) %>% 
               helper(type = "inline", size = "m", fade = TRUE, 
                      title = "Gene expression to colour cells by", 
                      content = c("Select gene to colour cells by gene expression", 
                      paste0("- Gene expression are coloured in a ", 
                             "White-Red colour scheme which can be ", 
                             "changed in the plot controls"))) 
           ), 
           column( 
             6, 
               actionButton("sc1a1ctog1", "Toggle plot controls"), 
               conditionalPanel( 
                 condition = "input.sc1a1ctog1 % 2 == 1", 
                 radioButtons("sc1a1ccol1", "Colour:", 
                              choices = c("White-Red", "Blue-Yellow-Red", 
                                          "Yellow-Green-Purple"), 
                              selected = "White-Red"), 
                 radioButtons("sc1a1cord1", "Plot order:", 
                              choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                              selected = "Max-1st", inline = TRUE) 
               ) 
             ) 
           ), 
           fluidRow( 
             column( 
               12, 
               uiOutput("sc1a1coup1.ui") 
             ) 
           ), 
           downloadButton("sc1a1coup1.pdf", "Download PDF"), 
           downloadButton("sc1a1coup1.png", "Download PNG"), br(), 
           div(style="display:inline-block", 
             numericInput("sc1a1coup1.h", "PDF / PNG height:", width = "138px", 
                          min = 4, max = 20, value = 6, step = 0.5)), 
           div(style="display:inline-block", 
             numericInput("sc1a1coup1.w", "PDF / PNG width:", width = "138px", 
                          min = 4, max = 20, value = 8, step = 0.5)) 
         ), # End of column (6 space) 
         column( 
           6, 
           br(), br(), br(), br(), 
           fluidRow( 
           # column( 
           #   6, selectInput("sc1a1cinp1", "Gene name:", choices=NULL) %>% 
           #     helper(type = "inline", size = "m", fade = TRUE, 
           #            title = "Gene expression to colour cells by", 
           #            content = c("Select gene to colour cells by gene expression", 
           #                        paste0("- Gene expression are coloured in a ", 
           #                               "White-Red colour scheme which can be ", 
           #                               "changed in the plot controls"))) 
           # ), 
             column( 
               6, 
               actionButton("sc1a1ctog2", "Toggle plot controls"), 
               conditionalPanel( 
                 condition = "input.sc1a1ctog2 % 2 == 1", 
                 radioButtons("sc1a1ccol2", "Colour:", 
                              choices = c("White-Red", "Blue-Yellow-Red", 
                                          "Yellow-Green-Purple"), 
                              selected = "White-Red"), 
                 radioButtons("sc1a1cord2", "Plot order:", 
                              choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                              selected = "Max-1st", inline = TRUE) 
               ) 
             ) 
           ), 
           fluidRow( 
             column( 
               12, 
               uiOutput("sc1a1coup2.ui") 
             ) 
           ), 
           downloadButton("sc1a1coup2.pdf", "Download PDF"), 
           downloadButton("sc1a1coup2.png", "Download PNG"), br(), 
           div(style="display:inline-block", 
             numericInput("sc1a1coup2.h", "PDF / PNG height:", width = "138px", 
                          min = 4, max = 20, value = 6, step = 0.5)), 
           div(style="display:inline-block", 
             numericInput("sc1a1coup2.w", "PDF / PNG width:", width = "138px", 
                          min = 4, max = 20, value = 8, step = 0.5)) 
         )  # End of column (6 space) 
       )    # End of fluidRow (4 space) 
 
     ),  # End of tab (2 space) 
 
 
 
 
 ### b3 
 tabPanel( 
  HTML("Split-By: Gene Coexpression"), 
  h4("Coexpression of two genes on reduced dimensions"), 
  "In this tab, users can visualise the coexpression of two genes ", 
  "on low-dimensional representions, while splitting between meta.data", 
  "columns which only have two possibilities.", 
  br(),br(), 
  fluidRow( 
    column( 
      3, h4("Dimension Reduction"), 
        selectInput("sc1b3drX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                    selected = sc1def$dimred[1]), 
        selectInput("sc1b3drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                    selected = sc1def$dimred[2]), 
        h4("Gene Expression"), 
        selectInput("sc1b3inp1", "Gene 1:", choices=NULL) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "Gene expression to colour cells by", 
                 content = c("Select gene to colour cells by gene expression", 
                             paste0("- Gene expression are coloured in a ", 
                                    "White-Red colour scheme which can be ", 
                                    "changed in the plot controls"))), 
        selectInput("sc1b3inp2", "Gene 2:", choices=NULL) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "Gene expression to colour cells by", 
                 content = c("Select gene to colour cells by gene expression", 
                             paste0("- Gene expression are coloured in a ", 
                                    "White-Blue colour scheme which can be ", 
                                    "changed in the plot controls"))), 
        #selectInput("sc1b3spl", "Split by:", choices = NULL), 
        selectInput("sc1b3split", "Split by:", choices = sc1conf[sc1conf$split==TRUE]$UI), 
        uiOutput("sc1b3split_select1_out1"), 
        uiOutput("sc1b3split_select2_out1"), 
        actionButton("sc1b3tog1", "Toggle plot controls"), 
        conditionalPanel( 
          condition = "input.sc1b3tog1 % 2 == 1", 
          radioButtons("sc1b3col1", "Colour:", 
                       choices = c("Red (Gene1); Blue (Gene2)", 
                                   "Orange (Gene1); Blue (Gene2)", 
                                   "Red (Gene1); Green (Gene2)", 
                                   "Green (Gene1); Blue (Gene2)"), 
                       selected = "Red (Gene1); Blue (Gene2)"), 
          radioButtons("sc1b3ord1", "Plot order:", 
                       choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                       selected = "Max-1st", inline = TRUE) 
        ), 
        br(),br(), 
        actionButton("sc1b3togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc1b3togL % 2 == 1", 
          selectInput("sc1b3sub1", "Cell information to subset:", 
                      choices = sc1conf[grp == TRUE]$UI, 
                      selected = sc1def$grp1), 
            uiOutput("sc1b3sub1.ui"), 
            actionButton("sc1b3sub1all", "Select all groups", class = "btn btn-primary"), 
            actionButton("sc1b3sub1non", "Deselect all groups", class = "btn btn-primary") 
        ), 
        br(),br(), 
        actionButton("sc1b3tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc1b3tog0 % 2 == 1", 
          sliderInput("sc1b3siz", "Point size:", 
                      min = 0, max = 4, value = 1.25, step = 0.25), 
          radioButtons("sc1b3psz", "Plot size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE), 
          radioButtons("sc1b3fsz", "Font size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Small", inline = TRUE), 
          radioButtons("sc1b3asp", "Aspect ratio:", 
                      choices = c("Square", "Fixed", "Free"), 
                      selected = "Square", inline = TRUE), 
          checkboxInput("sc1b3txt", "Show axis text", value = FALSE), 
          br(), br() 
        ), 
        column(12, 
          br(), br(), 
          uiOutput("sc1b3oup3.ui"), 
          downloadButton("sc1b3oup3.pdf", "Download PDF"), 
          downloadButton("sc1b3oup3.png", "Download PNG") 
        ) 
      ), # End of column (6 space) 
 
      column( 
        4, 
        #style="border-right: 2px solid black", 
        uiOutput("sc1b3oup1.ui"), 
        downloadButton("sc1b3oup1.pdf", "Download PDF"), 
        downloadButton("sc1b3oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc1b3oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc1b3oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 10, step = 0.5)), 
        br(), h4("Cell numbers"), # paste0( plot title, " cell numbers") ? 
        dataTableOutput("sc1b3oup1.dt") 
      ), 
      column ( 
        4, 
        uiOutput("sc1b3oup2.ui"), 
        downloadButton("sc1b3oup2.pdf", "Download PDF"), 
        downloadButton("sc1b3oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc1b3oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc1b3oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 10, step = 0.5)), 
        br(), h4("Cell numbers"), 
        dataTableOutput("sc1b3oup2.dt") 
      ) 
    ) 
  ),     # End of tab (2 space) 
 
 ### Tab1.c1: violinplot / boxplot 
 tabPanel( 
    HTML("Violinplot / Boxplot"),  
   h4("Cell information / gene expression violin plot / box plot"), 
   "In this tab, users can visualise the gene expression or continuous cell information ",  
   "(e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).", 
   br(),br(), 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", 
       selectInput("sc1c1inp1", "Cell information (X-axis):", 
                   choices = sc1conf[grp == TRUE]$UI, 
                   selected = sc1def$grp1) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell information to group cells by",  
                content = c("Select categorical cell information to group cells by",  
                            "- Single cells are grouped by this categorical covariate",  
                            "- Plotted as the X-axis of the violin plot / box plot")),  
       selectInput("sc1c1inp2", "Cell Info / Gene name (Y-axis):", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell Info / Gene to plot", 
                content = c("Select cell info / gene to plot on Y-axis", 
                            "- Can be continuous cell information (e.g. nUMIs / scores)", 
                            "- Can also be gene expression")), 
       radioButtons("sc1c1typ", "Plot type:", 
                    choices = c("violin", "boxplot"), 
                    selected = "violin", inline = TRUE), 
       checkboxInput("sc1c1pts", "Show data points", value = FALSE), 
       actionButton("sc1c1togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.sc1c1togL % 2 == 1", 
         selectInput("sc1c1sub1", "Cell information to subset:", 
                     choices = sc1conf[grp == TRUE]$UI, 
                     selected = sc1def$grp1), 
         uiOutput("sc1c1sub1.ui"), 
         actionButton("sc1c1sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("sc1c1sub1non", "Deselect all groups", class = "btn btn-primary") 
       ), br(), br(), 
       actionButton("sc1c1tog", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.sc1c1tog % 2 == 1", 
         sliderInput("sc1c1siz", "Data point size:",  
                     min = 0, max = 4, value = 1.25, step = 0.25),  
         radioButtons("sc1c1psz", "Plot size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE), 
         radioButtons("sc1c1fsz", "Font size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE)) 
     ), # End of column (6 space) 
     column(9, uiOutput("sc1c1oup.ui"),  
            downloadButton("sc1c1oup.pdf", "Download PDF"),  
            downloadButton("sc1c1oup.png", "Download PNG"), br(), 
            div(style="display:inline-block", 
                numericInput("sc1c1oup.h", "PDF / PNG height:", width = "138px", 
                             min = 4, max = 20, value = 8, step = 0.5)), 
            div(style="display:inline-block", 
                numericInput("sc1c1oup.w", "PDF / PNG width:", width = "138px", 
                             min = 4, max = 20, value = 10, step = 0.5)) 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
 
 
### Tab1.c2: Proportion plot 
tabPanel( 
  HTML("Proportion plot"), 
  h4("Proportion / cell numbers across different cell information"), 
  "In this tab, users can visualise the composition of single cells based on one discrete ", 
  "cell information across another discrete cell information. ",  
  "Usage examples include the library or cellcycle composition across clusters.", 
  br(),br(), 
  fluidRow( 
    column( 
      3, style="border-right: 2px solid black", 
      selectInput("sc1c2inp1", "Cell information to plot (X-axis):", 
                  choices = sc1conf[grp == TRUE]$UI, 
                  selected = sc1def$grp2) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to plot cells by",  
               content = c("Select categorical cell information to plot cells by", 
                           "- Plotted as the X-axis of the proportion plot")), 
      selectInput("sc1c2inp2", "Cell information to group / colour by:", 
                  choices = sc1conf[grp == TRUE]$UI, 
                  selected = sc1def$grp1) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to group / colour cells by", 
               content = c("Select categorical cell information to group / colour cells by", 
                           "- Proportion / cell numbers are shown in different colours")), 
      radioButtons("sc1c2typ", "Plot value:", 
                   choices = c("Proportion", "CellNumbers"), 
                   selected = "Proportion", inline = TRUE), 
      checkboxInput("sc1c2flp", "Flip X/Y", value = FALSE), 
      actionButton("sc1c2togL", "Toggle to subset cells"), 
      conditionalPanel( 
        condition = "input.sc1c2togL % 2 == 1", 
        selectInput("sc1c2sub1", "Cell information to subset:", 
                    choices = sc1conf[grp == TRUE]$UI, 
                    selected = sc1def$grp1), 
        uiOutput("sc1c2sub1.ui"), 
        actionButton("sc1c2sub1all", "Select all groups", class = "btn btn-primary"), 
        actionButton("sc1c2sub1non", "Deselect all groups", class = "btn btn-primary") 
      ), br(), br(), 
      actionButton("sc1c2tog", "Toggle graphics controls"), 
      conditionalPanel( 
        condition = "input.sc1c2tog % 2 == 1", 
        radioButtons("sc1c2psz", "Plot size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE), 
        radioButtons("sc1c2fsz", "Font size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE)) 
    ), # End of column (6 space) 
    column(9, uiOutput("sc1c2oup.ui"),  
           downloadButton("sc1c2oup.pdf", "Download PDF"),  
           downloadButton("sc1c2oup.png", "Download PNG"), br(), 
           div(style="display:inline-block", 
               numericInput("sc1c2oup.h", "PDF / PNG height:", width = "138px", 
                            min = 4, max = 20, value = 8, step = 0.5)), 
           div(style="display:inline-block", 
               numericInput("sc1c2oup.w", "PDF / PNG width:", width = "138px", 
                            min = 4, max = 20, value = 10, step = 0.5)) 
    )  # End of column (6 space) 
  )    # End of fluidRow (4 space) 
),     # End of tab (2 space) 
 
 
 
  ### Tab1.d1: Multiple gene expr 
  tabPanel( 
    HTML("Bubbleplot / Heatmap"), 
    h4("Gene expression bubbleplot / heatmap"), 
    "In this tab, users can visualise the gene expression patterns of ", 
    "multiple genes grouped by categorical cell information (e.g. library / cluster).", br(), 
    "The normalised expression are averaged, log-transformed and then plotted.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, style="border-right: 2px solid black", 
        textAreaInput("sc1d1inp", HTML("List of gene names <br /> 
                                          (Max 50 genes, separated <br /> 
                                           by , or ; or newline):"), 
                      height = "200px", 
                      value = paste0(sc1def$genes, collapse = ", ")) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "List of genes to plot on bubbleplot / heatmap", 
                 content = c("Input genes to plot", 
                             "- Maximum 50 genes (due to ploting space limitations)", 
                             "- Genes should be separated by comma, semicolon or newline")), 
        selectInput("sc1d1grp", "Group by:", 
                    choices = sc1conf[grp == TRUE]$UI, 
                    selected = sc1conf[grp == TRUE]$UI[1]) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "Cell information to group cells by", 
                 content = c("Select categorical cell information to group cells by", 
                             "- Single cells are grouped by this categorical covariate", 
                             "- Plotted as the X-axis of the bubbleplot / heatmap")), 
        radioButtons("sc1d1plt", "Plot type:", 
                     choices = c("Bubbleplot", "Heatmap"), 
                     selected = "Bubbleplot", inline = TRUE), 
        checkboxInput("sc1d1scl", "Scale gene expression", value = TRUE), 
        checkboxInput("sc1d1row", "Cluster rows (genes)", value = TRUE), 
        checkboxInput("sc1d1col", "Cluster columns (samples)", value = FALSE), 
        br(), 
        actionButton("sc1d1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc1d1togL % 2 == 1", 
          selectInput("sc1d1sub1", "Cell information to subset:", 
                      choices = sc1conf[grp == TRUE]$UI, 
                      selected = sc1def$grp1), 
          uiOutput("sc1d1sub1.ui"), 
          actionButton("sc1d1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc1d1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ), br(), br(), 
        actionButton("sc1d1tog", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc1d1tog % 2 == 1", 
          radioButtons("sc1d1cols", "Colour scheme:", 
                       choices = c("White-Red", "Blue-Yellow-Red", 
                                   "Yellow-Green-Purple"), 
                       selected = "Blue-Yellow-Red"), 
          radioButtons("sc1d1psz", "Plot size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE), 
          radioButtons("sc1d1fsz", "Font size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE)) 
      ), # End of column (6 space) 
      column(9, h4(htmlOutput("sc1d1oupTxt")), 
             uiOutput("sc1d1oup.ui"), 
             downloadButton("sc1d1oup.pdf", "Download PDF"), 
             downloadButton("sc1d1oup.png", "Download PNG"), br(), 
             div(style="display:inline-block", 
                 numericInput("sc1d1oup.h", "PDF / PNG height:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)), 
             div(style="display:inline-block", 
                 numericInput("sc1d1oup.w", "PDF / PNG width:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  )      # End of tab (2 space) 
  , 
tabPanel( 
  HTML("Markers Seurat Cluster, All"), 
  h4("This page is meant to display associations between genetic markers and Seurat data clusters through an embedded spreadsheet."), 
  br(),br(), 
  fluidRow( 
    selectInput("sc1m_all_select", "Select cluster:", choices=strsplit(sc1conf[grp==TRUE][UI=="seurat_clusters"]$fID, "\\|")[[1]]), # would you ever want to look at all? 
    column(12, 
      DT::dataTableOutput("sc1m_all.ui") 
    ) 
  ) 
) 
, 
tabPanel( 
  HTML("Markers Seurat Cluster, Top 20"), 
  h4("This page is meant to display associations between the top 20 genetic markers and Seurat data clusters through an embedded spreadsheet."), 
  br(),br(), 
  fluidRow( 
    column(12, 
      DT::dataTableOutput("sc1m_t20") 
    ) 
  ) 
) 
, 
tabPanel( 
  HTML("Diff. Exp. Genes"), 
  h4("This page is meant to display associations between genetic markers and cell types tagged to their Seurat data cluster."), 
  br(),br(), 
  fluidRow( 
    selectInput("sc1de_genes_select", "Select differential expression:", choices=strsplit(sc1conf$DEs[1], "\\|")[[1]]), # would you ever want to look at all? 
    column(12, 
      #DT::dataTableOutput("sc1de_genes.ui") 
      dataTableOutput("sc1de_genes.ui") 
    ) 
  ) 
) 
,
 
tabPanel( 
  HTML("Gene signature"), 
  h4("View gene signature data through custom list of genes; uses AUCell for calculation."), 
  br(),br(), 
  fluidRow( 
    column(3, style="border-right: 2px solid black", 
      textAreaInput("gsig_list", "List of genes", value=paste0(sc1def$genes, collapse = "\n"), height="200px"), 
      selectInput("gsig_group", "Group by:", choices=sc1conf[grp==TRUE]$UI, selected=sc1def$grp1), 
      actionButton("gsig_generate_plot", "Generate / Update plot", class = "btn-success"), 
      br(),br(), 
      actionButton("gsig_subset_toggle", "Toggle to subset cells"), 
      conditionalPanel("input.gsig_subset_toggle % 2 == 1", 
                       selectInput("gsig_subset", "Cell information to subset:", choices=sc1conf[grp==TRUE]$UI, selected=sc1def$grp1), 
                       uiOutput("gsig_subset.ui"), 
                       actionButton("gsig_subset_sel_all", "Select all groups"), 
                       actionButton("gsig_subset_desel_all", "Deselect all groups") 
      ), 
      br(),br(), 
      actionButton("gsig_graphics_toggle", "Toggle graphics controls"), 
      conditionalPanel("input.gsig_graphics_toggle % 2 == 1", 
              #radioButtons("gsig_graphics_plot_size", "Plot size:", choices = c("Small", "Medium", "Large"), inline=TRUE), 
              radioButtons("gsig_graphics_font_size", "Font size:", choices = c("Small", "Medium", "Large"), inline=TRUE), 
              sliderInput("gsig_graphics_point_size", "Point size:", min=0, max=4, value=0.60, step=0.10), 
              checkboxInput("gsig_graphics_cell_info", "Show cell info labels", value=TRUE), 
              radioButtons("gsig_graphics_color", "Color:", choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"), selected="White-Red") 
      ) 
    ), 
    column(9, 
      htmlOutput("gsig_list.html"), 
      plotOutput("gsig_plot", width="100%", height="600px"), 
      #htmlOutput("gsig_graphs.html"), 
      downloadButton("gsig_download.pdf", "Download PDF"), 
      downloadButton("gsig_download.png", "Download PNG"), 
      br(), 
      fillRow( 
        numericInput("gsig_download_height", "PDF / PNG height:", value=6, min=4, max=20, step=0.5, width="138px"), 
        numericInput("gsig_download_width", "PDF / PNG width:", value=14, min=4, max=20, step=0.5, width="140px"), 
        width="45%" 
      ) 
    ) 
  ) 
) 
, 

tabPanel( 
  HTML("Diff. Gene Exp., Volcano"), 
  h4("Generate volcano plots from differential gene expression analysis."), 
  br(),br(), 
  fluidRow( 
    column(3, style="border-right: 2px solid black", 
      #textAreaInput("volc_tags", "Genes to tag") 
      selectInput("volc_de_select", "Select differential expression:", choices=strsplit(sc1conf$DEs[2], "\\|")[[1]]), # would you ever want to look at all? 
      actionButton("volc_subset_toggle", "Toggle to subset cells"), 
      conditionalPanel("input.volc_subset_toggle % 2 == 1", 
                       #selectInput("volc_subset", "Cell information to subset:", choices=sc1conf[grp==TRUE]$UI, selected=sc1def$grp1), 
                       selectInput("volc_subset", "Cell information to subset:", choices=sc1conf$DEs[3]), 
                       uiOutput("volc_subset.ui"), 
                       actionButton("volc_subset_sel_all", "Select all groups"), 
                       actionButton("volc_subset_desel_all", "Deselect all groups") 
      ), 
      br(), 
      checkboxInput("volc_top10_toggle", "Toggle Top 10 Genes", value=FALSE), 
      br(), 
      sliderInput("volc_fold_change", "Fold Change", 0, 2, 0.3, round=-1, step=0.1), # dynamic max? 
      textInput("volc_p_value", "Adjusted p-Value", placeholder="0.05"), 
      sliderInput("volc_xlim_up", "Upperbound range of fold change", 0, 10, 5, round=-1, step=0.1), 
      sliderInput("volc_xlim_down", "Lowerbound range of fold change", -10, 0, -5, round=-1, step=0.1), 
      sliderInput("volc_ylim", "Ceiling of p-value", 30, 600, 300, round=1, step=1) 
    ), 
    column(9, 
      plotOutput("volc_plot", width="100%", height="600px"), 
      downloadButton("volc_download.pdf", "Download PDF"), 
      downloadButton("volc_download.png", "Download PNG"), 
      br(), 
      fillRow( 
        numericInput("volc_download_height", "PDF / PNG height:", value=6, min=4, max=20, step=0.5, width="138px"), 
        numericInput("volc_download_width", "PDF / PNG width:", value=14, min=4, max=20, step=0.5, width="140px"), 
        width="45%" 
      ) 
    ) 
  ) 
) ,    
br(), 
p("", style = "font-size: 125%;"), 
p(em("This webpage was made using "), a("ShinyCell", 
  href = "https://github.com/SGDDNB/ShinyCell",target="_blank")), 
br(),br(),br(),br(),br() 
))) 
 
 
 
 