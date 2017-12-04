#in R studio, click Run App (i suggest Run External in browser)

#required libraries
require(shiny)
require(compositions)
require(zCompositions)
require(markdown)
require(ALDEx2)
require(propr)
require(igraph)
#require(vegan)

#userinterface -------------------------------------------------------
ui <- fluidPage(theme= "bootstrap.css",
  titlePanel("omicplotR: Visual exploration of omic datasets as compositions"),

  navbarPage("omicplotR",
             tabPanel("Getting started",
                      fluidRow(
                        column(6, offset = 3, includeMarkdown("help.Rmd")
                               )
                      )
                      ),
             tabPanel("Input data",
                      sidebarLayout(
                        sidebarPanel(
                          tabsetPanel(
                            type = "tabs",
                            selected = "Data",
                            tabPanel(
                              "Data",
                              fileInput(
                              'file1',
                              label = h3('Choose Data'),
                              accept = c('text/csv',
                                         'text/comma-separated-values,text/plain',
                                         '.csv')
                            ),
                            actionButton("showdata", "Show data"),

                            fileInput(
                              'file2',
                              h3('Choose Metadata'),
                              accept = c('text/csv',
                                         'text/comma-separated-values,text/plain', '.csv')
                            ),
                            actionButton("showmetadata", "Show metadata")),
                            tabPanel(
                              "Example Data",
                              checkboxInput("exampledata", "Vaginal dataset (data and metadata)"),
                              checkboxInput("exampledata2", "Selex dataset (data only)")
                            )
                          )
                        ),#sidebarpanel
                      mainPanel(
                        tabsetPanel(
                          type = "tabs",
                          selected="Data",
                          tabPanel("Data",
                                   dataTableOutput('datatable')
                                   ),
                          tabPanel("Metadata",
                                   dataTableOutput("metadatatable")
                                   )
                        ))
                      )
             ),#tabpanel
             tabPanel("PCA Biplots",
                      sidebarLayout(
                        sidebarPanel(
                          tabsetPanel(
                            type = "tabs",
                            selected = "Filtering",
                            tabPanel("Filtering",
                            h3(textOutput("textTitle")),
                            fluidRow(column(6, numericInput("mincounts", "Minimum count per OTU", value = 0)),
                                     column(6, numericInput("minreads", "Minimum count sum per sample", value = 0))),
                            fluidRow(column(6, numericInput("minprop", "Minimum proportional abundance", value = 0.0)),
                                     column(6, numericInput("maxprop", "Maximum proportional abundance", value = 1))),
                            sliderInput("minsum", "Minimum count sum per OTU", min = 0, max = 10000, value = 1),
                            uiOutput("varianceslider"),
                            sliderInput(
                              inputId = "scale",
                              label = "Adjust scale (0 = between samples, 1 = between OTUs)",
                              value = 0,
                              min = 0,
                              max = 1,
                              step = 1
                            ),

                            textOutput("removedsamples"),
                            textOutput("removedotus")
                            ),
                            tabPanel("Colouring Options",
                                     uiOutput("choose_column"),
                                     radioButtons(
                                       label = h4("Choose data type"),
                                       inputId = "colouringtype",
                                       choices = list(
                                         "Continuous colouring" = 1,
                                         "Quartile" = 2,
                                         "Nonnumeric metadata (up to 10 categories)" = 3
                                       )
                                     ),
                                     actionButton("choosemeta", "Filter by metadata values"),
                                     textOutput("test"),
                                     checkboxInput("removesamplenames", label = "Remove sample names"),
                                     checkboxInput("arrowcheckbox", label = "Show arrows", value = FALSE),
                                     checkboxInput("taxoncheckbox", label = "Show taxonomy", value = FALSE),
                                     uiOutput("taxchoice")#,
                                     #textInput("biplot_title", label = "Biplot title", placeholder = "Title for biplot")
                                     )
                          )
                        ),
                        mainPanel(
                            tabsetPanel(
                          selected = "Biplot",
                          tabPanel(
                            "Biplot",
                            h3(textOutput("nodata")),
                            splitLayout(
                              cellWidths = c("60%", "40%"),
                              plotOutput("biplot", width = 600, height = 600),
                              plotOutput("screeplot")
                            )
                          ),
                          tabPanel(
                            "Colored Biplot",
                            h3(textOutput("nometadata")),
                            splitLayout(
                              cellWidths = c("65%", "35%"),
                              plotOutput("coloredBiplot", width = 600, height = 600),
                              plotOutput("metahist")
                            )
                          ),
                          tabPanel("Removed data",
                                   fluidRow(
                                     column(6, plotOutput("colsums")),
                                     column(6, plotOutput("rowsums")),
                                     tabsetPanel(tabPanel("Samples removed", fluidRow(column(
                                       12, dataTableOutput("removedDT")
                                     )),
                                     actionButton("showremoved", "Show removed samples/OTUs")),
                                     tabPanel("Features removed", fluidRow(column(
                                       12, dataTableOutput("removedDTotu")
                                     ))))
                                   )
                          )
                        )))#sidebarlayout
             ), #tabpanel
            #  tabPanel("Association plots",
            #           sidebarLayout(
            #             sidebarPanel(width = 3,
            #               numericInput("rhocutoff",
            #                            label = "Input rho cutoff",
            #                            value = -0.25),
            #               textOutput("associationtext")),
            #             mainPanel(plotOutput("associationplot"))
            #           )),
             tabPanel("Relative abundance plots",
                      sidebarLayout(
                        sidebarPanel(width = 3,
                          sliderInput("abundcutoffbarplot", "Choose abundance cutoff", min = 0, max = 0.25, value= 0, step = 0.01),
                          selectInput("clustermethod", "Select clustering method",
                                      choices = list("complete" = 1, "single" = 2, "ward.D2" = 3), selected = 3),
                          selectInput("dismethod", "Select distance matrix method",
                                      choices = list("euclidean" = 1, "maximum" = 2, "manhattan" = 3), selected = 1),
                                      textOutput("filter_warning_dendro")
                        ),
                        mainPanel(
                          textOutput("dendrotext"),
                          plotOutput("dendrogram",
                                     dblclick = "dendro_dblclick",
                                     brush = brushOpts(
                                       id = "dendro_brush",
                                       resetOnNew = TRUE
                                     )),
                          plotOutput("barplot",
                                     dblclick = "bp_dblclick",
                                     brush = brushOpts(
                                       id = "bp_brush",
                                       resetOnNew = TRUE
                                     ))
                        )
                      )
                      ),
             tabPanel("Effect plots",
                      sidebarLayout(
                        sidebarPanel(
                          tabsetPanel(
                            tabPanel("Calculate",
                              radioButtons("ep_chooseconds", "How will you select your conditions?", choices = list("Manually" = 1, "From metadata" = 2)),
                              selectInput("denomchoice", "Choose ALDEx2 method", choices = list("all" = 1, "iqlr" = 2, "zero"= 3)),
                              actionButton("effectplot_ab", "Generate effect plot"),
                              # textOutput("anosim"),
                              conditionalPanel("input.ep_chooseconds == '2'",
                                               uiOutput("colselectcond")),
                              uiOutput("conditions"),
                              textOutput("effectwarning"),
                              textOutput("filter_warning_effect")
                            ),
                            tabPanel("Input",

                            fileInput(
                            'effect_file',
                            label = h3('Input ALDEx2 file'),
                            accept = c('text/csv',
                                       'text/comma-separated-values,text/plain',
                                       '.csv')
                            ),
                            actionButton("effectplot_ab2", "Generate effect plot")
                          )
                          )
                        ),
                        mainPanel(
                          tabsetPanel(selected = "Calculate",
                            tabPanel("Calculate",
                          h3(textOutput("nostripchart"),
                          splitLayout(
                            cellWidths = c("50%", "50%"),
                            plotOutput("stripchart"),

                            plotOutput("effectMW",
                                      hover = hoverOpts(id = "mw_hover"))
                          ),
                          fluidRow(column(6, offset = 6,
                                          uiOutput("mw_hovertext"),
                                          plotOutput("effectMA",
                                                     hover = hoverOpts(
                                                       id = "ma_hover"
                                                     )),
                                          uiOutput("ma_hovertext"))
                          ))),
                          tabPanel("ALDEx2 input",
                          splitLayout(
                            cellWidths = c("50%", "50%"),
                            plotOutput("table_effect",
                                       hover = hoverOpts(id = "mw_hover2")),
                            plotOutput("table_bland",
                                      hover = hoverOpts(
                                        id = "ma_hover"
                                      ))
                          ),
                          uiOutput("featurename"),
                          textInput("point.colour", label = "Colour points by name",
                          placeholder = "Input string to search in row names..."
                        ),
                        actionButton("update_points", label = "Update")
                        )
                        )
                      ))),
              navbarMenu("More",
                         tabPanel("Instructions",
                                  fluidRow(
                                    (column(6, offset = 3, includeMarkdown("extrainformation.rmd")))
                                  )
                                  ),
                         tabPanel("Quit omicplotR",
                                  actionButton("stopApp", label="Quit omicplotR")
                                  )
                         )
  )
 )
