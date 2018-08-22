server <- function(input, output) {

################################################################################
    #downloads

    filtering_options <- reactive({
        #get the reactive inputs for downloadable file
        reactive_values <- paste("min.count <- ", input$mincounts, "\n",
        "min.prop <- ", input$minprop, "\n",
        "max.prop <- ", input$maxprop, "\n",
        "min.sum <- ", input$minsum, "\n",
        "min.reads <- ", input$minreads, "\n",
        "taxselect <- ", as.numeric(input$taxlevel), "\n",
        "taxoncheck <- ", as.numeric(input$taxlevel), "\n",
        "arrowcheck <- ", input$arrowcheckbox, "\n",
        "scale.slider <- ", input$scale, "\n",
        "removenames <- ", input$removesamplenames, "\n",
        "var.filt <- ", input$varslider, "\n",
        "abund <- ", input$abundcutoffbarplot, "\n",
        "dist <- ", input$dismethod, "\n",
        "clust <- ", input$clustermethod, "\n",
        "vals <- ", vals$data, "\n",
        "metaval <- ", metaval$data,
        sep = "")
    })

    #combine reactive values with script for PCA biplot
    PCA_script <- reactive({
        y <- filtering_options()
        PCA_script <- c(y, readLines("./PCA_script.R"))
    })

    rab_script <- reactive({
        y <- filtering_options()
        rab_script <- c(y, readLines("./rab_script.R"))
    })

    effect_script <- reactive({
        y <- filtering_options()
        effect_script <- c(y, readLines("./effect_script.R"))
    })

    file_name <- reactive({
        inFile <- input$file1
        file_name <- inFile$name
    })

    #download script for PCA biplot
    output$PCA_download <- downloadHandler(
        filename=function(){paste(file_name(), "_PCA.r", sep = "")},
        content= function(file) {
            writeLines(PCA_script(), file)
        }
    )

    #download script for rab plots
    output$rab_download <- downloadHandler(
        filename=function(){paste(file_name(), "_rab.r", sep = "")},
        content= function(file) {
            writeLines(rab_script(), file)
        }
    )

    #download script for effect plots
    output$effect_download <- downloadHandler(
        filename=function(){paste(file_name(), "_effect.r", sep = "")},
        content= function(file) {
            writeLines(effect_script(), file)
        }
    )

################################################################################
  #observe events
  #quit button
  options(shiny.maxRequestSize=70*1024^2)
  observeEvent(input$stopApp, {
    stopApp(returnValue = invisible())
  })

  #reset filter
  observeEvent(input$reset, {
    vals$data <- NULL
    removeModal()
  })

  #pop up the modal
  observeEvent(input$choosemeta, {
    showModal(metaModal())
  })

  #convert input to reactive object
  observeEvent(input$update, {
    vals$data <-
    c(input$cnname, input$cnname2, input$cnname3, input$cnname4)
    metaval$data <- input$select_biplot_filter
    removeModal()
  })

  #dendrogram brush ranges
  ranges <- reactiveValues(x = NULL, y = NULL)
  bp_ranges <- reactiveValues(x = NULL, y = NULL)

  ###########################################################################
  #custom UIs
  #make the pop up
  metaModal <- function (x) {
    modalDialog(
      renderUI({
        metadata <- metadata()
        #message if there is nothing uploaded
        if (is.null(metadata)) {
          selectInput(
            "select_biplot_filter",
            label = h3("Filter biplot by metadata column"),
            choices = list(
              "Input data" = 1
            ),
            selected = 1
          )
        } else {
          options <- colnames(metadata)
          selectInput(
            "select_biplot_filter",
            label = h3("Filter data by metadata column"),
            choices = options,
            selected = 1
          )
        }
      }),
      textInput("cnname", "Value 1", placeholder = "Enter value from metadata or leave blank"),
      textInput("cnname2", "Value 2", placeholder = "Enter value from metadata or leave blank"),
      textInput("cnname3", "Value 3", placeholder = "Enter value from metadata or leave blank"),
      textInput("cnname4", "Value 4", placeholder = "Enter value from metadata or leave blank"),
      title = "Choose metadata values to filter PCA biplot",
      "Enter metadata values of the metadata variable you have chosen (must be exact match).
      The biplot will be replotted with only the samples identified by
      the metadata values. Leave inputs blank if you don't need them.
      Samples that do not have metadata will be coloured black",
      footer = tagList(actionButton("update", "Update Filter"),
      actionButton("reset", "Reset Filter"),
      modalButton("Cancel")),
      easyClose = TRUE
    )
  }

  output$taxchoice <- renderUI({
    if (input$taxoncheckbox) {
      radioButtons("taxlevel",
      "Choose taxonomic level to display",
      c("Kingdom" = 1,
      "Phylum" = 2,
      "Class" = 3,
      "Order" = 4,
      "Family" = 5,
      "Genus" = 6,
      "Species" = 7),
      selected = 4
    )
  }
})

#dynamic variance slider
output$varianceslider <- renderUI({
  data <- data()

  sliderInput(
    "varslider",
    "Variance cutoff",
    min = 0,
    max = 5,
    value = 0,
    step = 0.05
  )
})

#effect plot conditions
output$conditions<- renderUI({
  if (input$ep_chooseconds == 1) {
    tagList(
      numericInput("group1s", "Number of columns for group 1", value = 0),
      numericInput("group2s", "Number of columns for group 2", value = 0)
    )
  } else {
    meta <- metadata()
    options <- colnames(meta)
    cn <- input$colselect
    choice <- unique(meta[[cn]])
    tagList(
      selectInput("group1", "Choose group 1", choices = choice),
      selectInput("group2", "Choose group 2", choices = choice))
    }
  })

  #choose column from metadata
  output$colselectcond <- renderUI({
    meta <- metadata()
    options <- colnames(meta)
    selectInput("colselect", "Select column from metadata", choices = options)
  })

  ###########################################################################
  #data input
  #get data from uploaded file
  data <- reactive({
    if (input$exampledata) {
      read.table(
        "example_data.txt",
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        quote = "",
        row.names = 1,
        check.names = FALSE,
        comment.char = "",
        na.strings = ""
      )
    } else {
      if (input$exampledata2) {
        read.table(
          "selex.txt",
          header = TRUE,
          sep = "\t",
          stringsAsFactors = FALSE,
          quote = "",
          row.names = 1,
          check.names = FALSE,
          comment.char = "",
          na.strings = ""
        )
      } else{
        #reactive input file
        inFile <- input$file1

        #return NULL when no file is uploaded
        if (is.null(inFile))
        return(NULL)

        #reads the file
        #if there's an error, it doesn't quit omicplotR
        data <- try(read.table(
          inFile$datapath,
          header = TRUE,
          sep = "\t",
          stringsAsFactors = FALSE,
          quote = "",
          row.names = 1,
          check.names = FALSE,
          comment.char = "",
          na.strings = ""), silent = TRUE)

        #if rownames error
        if(grepl("duplicate 'row.names'", data[1], fixed = TRUE)) {
          showModal(rownamesModal())
        } else {
          #if no error, then continue
          data <- data
        }
      }
    }
  })

  #get metadata from uploaded file
  metadata <- reactive({

    if (input$exampledata) {
      read.table(
        "example_metadata.txt",
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        quote = "",
        row.names = 1,
        check.names = FALSE,
        comment.char = ""
      )
    } else {
      #reactive input file
      inFile2 <- input$file2

      #return NULL when no file is uploaded
      if (is.null(inFile2))
      return(NULL)

      #reads the file
      read.table(
        inFile2$datapath,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        quote = "",
        row.names = 1,
        check.names = FALSE,
        comment.char = ""
      )
    }
  })

################################################################################
#tests to make sure input is correct format
formatModal <- function(failed = FALSE) {
    modalDialog(
        title="Incorrect format",
        "There are duplicate column names. Ensure all column names unique, then re-import your file.",
        easyClose = TRUE,
        footer = modalButton("Dismiss")
    )
}

  rownamesModal <- function(failed = FALSE) {
    modalDialog(
      title="Incorrect format",
      "There are duplicate row names. Ensure all row names unique, then re-import your file.",
      easyClose = TRUE,
      footer = modalButton("Dismiss")
    )
  }

#check data format
   observeEvent(input$showdata, {
     inFile <- data()

     #return NULL when no file is uploaded
     if (is.null(inFile)) {
         return(NULL)
     }

     #make frequency table of occurence of colnames
     c_occur <- data.frame(table(colnames(inFile)))

     #if frequency is more than one, show the "your format is wrong" modalDialog
     if (max(c_occur$Freq) > 1) {
         showModal(formatModal())
     } else {
         return(NULL)
     }
 })

#metadata check format
 observeEvent(input$showmetadata, {
   inFile2 <- metadata()

   #return NULL when no file is uploaded
   if (is.null(inFile2)) {
     return(NULL)
   }

   #make frequency table of occurence of colnames
   c_occur <- data.frame(table(colnames(inFile2)))

   #if frequency is more than one, show the "your format is wrong" modalDialog
   if (max(r_occur$Freq) > 1) {
     showModal(formatModal())
   } else {
     return(NULL)
   }
 })

################################################################################


  #input ALDEx2 table
  effect_input <- reactive({

    inFile3 <- input$effect_file

    #return NULL when no file is uploaded
    if (is.null(inFile3))
    return(NULL)

    read.table(
      inFile3$datapath,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      quote = "",
      row.names = 1,
      check.names = FALSE,
      comment.char = "",
      na.strings = ""
    )
  })

  #create reactive object for biplot filtering inputs
  vals <- reactiveValues(data = NULL)
  metaval <- reactiveValues(data = NULL)

  #get column selection for metadata
  column <- reactive({
    cn <- input$selectcolumn
  })

  ###########################################################################
  #data processing
  #generate taxonomy vector for biplots
  taxonomy <- reactive({
    data <- data.t()
    tax.result <- data.check()

    #if no taxonomy column, return NULL for later
    if (isTRUE(tax.result)) {
      tax <- NULL
    } else {
      #remove taxonomy column
      tax <- data$taxonomy
      data$taxonomy <- NULL

      #get genus names
      genus <- vapply(strsplit(as.character(tax), "[[:punct:]]"), "[", taxselected, FUN.VALUE=character(1))
    }
  })

  #if taxonomic column present, return boolean. no tax column returns TRUE
  data.check <- reactive ({
    data.r <- data()
    if (is.null(data.r$taxonomy)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  })

  #do transformations of data for biplot
  data.t <- reactive ({
    #get the data
    x <- data()
    meta <- metadata()
    cn <- column()

    #set to zero, but inputs from sliders/numeric inputs
    min.count <- input$mincounts
    min.prop <- input$minprop
    max.prop <- input$maxprop
    min.sum <- input$minsum
    min.reads <- input$minreads

    #catch the errors
    validate((need(input$minprop, "Processing..")))

    #check for tax
    if (is.null(x$taxonomy)) {
      taxCheck <- TRUE
    } else {
      taxCheck <- FALSE
    }

    #get filtered data if filtered
    if (is.null(vals$data)) {
      x <- x
    } else {

      x <- omicplotr.metadataFilter(x, meta, column = metaval$data, values = vals$data)
    }

    x.filt <- omicplotr.filter(x, min.reads = min.reads, min.count = min.count, min.prop = min.prop, max.prop = max.prop, min.sum = min.sum)

  })

  #prcomp object
  data.prcomp <- reactive({
    #get data
    data.t <- data.t()
    var.filt <- input$varslider
    data <- data()

    validate(need(input$varslider, "Calculating..."))

    #check for tax
    if (is.null(data$taxonomy)) {
      taxCheck <- TRUE
    } else {
      taxCheck <- FALSE
    }

    #catch error if there is no data
    if (is.null(data)) {
      return(NULL)
    } else {
      if (isTRUE(taxCheck)) {
        data.t <- data.t
      } else {
        data.t$taxonomy <- NULL
      }

      data.pr <- omicplotr.clr(data.t, var.filt)

    }
    return(data.pr)
  })

#these are for effect plot choices by metadata
cn <- ""
group1 <- ""
group2 <- ""
denom <- ""

#forces effect plot to update ONLY when generate action button is clicked
#these values are needed for d.clr and ald.obj
observeEvent(input$effectplot_ab, {
          cn <<- input$colselect
          group1 <<- input$group1
          group2 <<- input$group2
          denom <<- input$denomchoice
  })

  #computing aldex object
  d.clr <- reactive({
    x <- data.t()
    g1s <- input$group1s
    g2s <- input$group2s
    meta <- metadata()

    #require user to click action button
    validate(need(input$effectplot_ab, ""))

    #get filtered metadata if it is filtered
    if (is.null(vals$data)) {
      x <- x
    } else {

      x <- omicplotr.metadataFilter(x, meta, column = metaval$data, values = vals$data)

      meta <- meta[which(rownames(meta) %in% colnames(x)),]
  }

    # Separate the taxonomy column from the counts
    if (is.null(x$taxonomy)) {
      d <- x
    } else {
      d <- x[, 0:(dim(x)[2] - 1)]
      taxon <- x[(dim(x)[2])]
    }

    if (input$ep_chooseconds == 1) {
      conds <- c(rep("group1", g1s),
      rep("group2", g2s))
    }

    if (input$ep_chooseconds ==2) {
      #filter the metadata and keep only the data which have been chosen
      group1 <- rownames(meta[which(meta[[cn]] == group1), ])
      group1.filt <- group1[group1 %in% (colnames(x))]
      data1.filt <- x[,which(colnames(x) %in% group1.filt)]
      one <- length(data1.filt)

      group2 <- rownames(meta[which(meta[[cn]] == group2), ])
      group2.filt <- group2[group2 %in% (colnames(x))]
      data2.filt <- x[,which(colnames(x) %in% group2.filt)]
      two <- length(data2.filt)

      conds <- c(rep("Group 1", one),
      rep("Group 2", two))

      #combine
      d <- cbind(data1.filt, data2.filt)

    }

    # Generate ALDEx2 object and get a distribution of values that are centered
    # log ratio transformed at each feature*sample

    choices = c("all", "iqlr", "zero")

    method <- choices[as.numeric(denom)]
    withProgress(
      message = "Calculating",
      detail = "Calculating expected clr values", value = 1/2, {

        d.clr <- aldex.clr(d, mc.samples=128, conds = conds, denom = method, verbose=TRUE)
        incProgress(1/2, message = "clr values calculated")
      }
    )
    #output d.clr
    d.clr <- d.clr

  })

  aldex.obj <- reactive({
    x <- data.t()
    d.clr <- d.clr()
    meta <- metadata()
    g1s <- input$group1s
    g2s <- input$group2s

    #get filtered data if filtered
    if (is.null(vals$data)) {
      x <- x
    } else {

      x <- omicplotr.metadataFilter(x, meta, column = metaval$data, values = vals$data)

      meta <- meta[which(rownames(meta) %in% colnames(x)),]
  }

    if (is.null(x$taxonomy)) {
      d <- x
    } else {
      d <- x[, 0:(dim(x)[2] - 1)]
      taxon <- x[(dim(x)[2])]
    }

    if (input$ep_chooseconds == 1) {
      conds <- c(rep("group1", g1s),
      rep("group2", g2s))
    }

    if (input$ep_chooseconds ==2) {
      #filter the meta and keep only the data which have been chosen

      #make group of samples from metadata which are chosen by user
      group1 <- rownames(meta[which(meta[[cn]] == group1), ])

      #make sure these are actually in the otu table
      group1.filt <- group1[group1 %in% (colnames(x))]

      #subset dataframe to take samples that are only in the filtered group
      data1.filt <- x[,which(colnames(x) %in% group1.filt)]

      #get length for later
      one <- length(data1.filt)

      group2 <- rownames(meta[which(meta[[cn]] == group2), ])
      group2.filt <- group2[group2 %in% (colnames(x))]
      data2.filt <- x[,which(colnames(x) %in% group2.filt)]
      two <- length(data2.filt)

      conds <- c(rep("Group 1", one),
      rep("Group 2", two))
    }

    withProgress(
      message = "Calculating t-test...",
      detail = "Calculating t test and effect size", value = 1/2, {
        # Run a t-test between the two groups specified by conds
        #TODO option to generate effect plots ttest to save time.
        # this will require changing aldex.plot in effect plot output
        x.tt <- aldex.ttest(d.clr, conds, paired.test=FALSE)
        incProgress(1/4, message = "Calculating effect size...")
        # Generate effect sizes
        x.effect <- aldex.effect(d.clr, conds, include.sample.summary=TRUE,
          verbose=TRUE)
          incProgress(1/4, message = "Efect sizes calculated")
        }

      )
      # Combine into one table
      x.all <- data.frame(x.tt, x.effect, stringsAsFactors=FALSE)
    })

    #stripchart data
    clr.strip <- reactive({
      x <- d.clr()
      codaSeq.repl <- function(x){
        # Initialize medians vector to be reshaped later to matrix
        MC.means <- numeric()


        # Iterating all samples
        for(i in seq_len(numConditions(x))) {
          sample.means <- numeric()

          # Get sample's MC Instances for features
          MC.matrix <- getMonteCarloReplicate(x, i)

          # Iterating features
          for(j in seq_len(numFeatures(x))) {
            feature.means <- mean(MC.matrix[j, ])
            sample.means <- append(sample.means, feature.means)
          }

          # Add subject's MC instance medians to master vector
          MC.means <- append(MC.means, sample.means)
        }

        # Turn the vector of medians into a matrix and then data frame with
        result <- matrix(MC.means, nrow = numFeatures(x), ncol = numConditions(x))
        result <- data.frame(result, row.names = getFeatureNames(x))
        colnames(result) <- getSampleIDs(x)

        return(result)
      }

      clr.strip <- codaSeq.repl(x)
    })

    ###########################################################################
    # outputs
    #show removed OTUs/samples
    output$removedDT <- renderDataTable({
      data.pr <- data.prcomp()
      data.in <- data()

      validate(need(input$showremoved, ""))

      omicplotr.getRemovedSamples(data.in, data.pr)
    },
    #force datatable to size of window
    options = list(scrollX = TRUE)
  )

  output$removedDTotu <- renderDataTable({
    data.pr <- data.prcomp()
    data.in <- data()

    validate(need(input$showremoved, "Click 'Show removed samples/OTUs' to view removed OTUs"))

    omicplotr.getRemovedFeatures(data.in, data.pr)
  },
  #force datatable to size of window
  options = list(scrollX = TRUE))

  #choose the metadata column you want to colour by
  output$choose_column <- renderUI({
    metadata <- metadata()

    #message if there is nothing uploaded
    if (is.null(metadata)) {
      selectInput(
        "selectcolumn",
        label = h3("Choose column to colour by"),
        choices = list(
          "No metadata file detected" = 1
        ),
        selected = 1
      )
    } else {
      options <- colnames(metadata)
      selectInput(
        "selectcolumn",
        label = h3("Choose column to colour by"),
        choices = options,
        selected = 1
      )
    }
  })

  output$textTitle <- renderText({
    "Choose filtering options"
  })

  output$nometadata <- renderText({
    metadata <- metadata()
    if (is.null(metadata)) {
      "Select metadata file"
    }
  })

  output$nodata <- renderText({
    metadata <- metadata()
    if (is.null(data)) {
      "Select data file"
    }
  })

  output$nostripchart <- renderText({
    validate(need(input$effectplot_ab, "Select conditions and click 'Generate effect plot'"))
  })

  output$test <- renderText({
    if (!is.null(vals$data)) {
      c("Data filtered by metadata column:", metaval$data,  "value:", vals$data)
    } else {
      return(NULL)
    }
  })

  output$filter_warning_effect <- renderText({
    if (!is.null(vals$data)) {
      c("Data filtered by metadata column:", metaval$data,  "value:", vals$data)
    } else {
      return(NULL)
    }
  })

    output$filter_warning_dendro <- renderText({
      if (!is.null(vals$data)) {
        c("Data filtered by metadata column:", metaval$data,  "value:", vals$data)
      } else {
        return(NULL)
      }
    })

  output$datatable <- renderDataTable({
    validate((need(input$showdata, "Click 'Check data' to check format and display table")))

    if (input$showdata) {
      data <- data()

      output <- cbind(OTUs = rownames(data), data)

    }
  }, options = list(scrollX = TRUE)) #force data table to window size

  output$metadatatable <- renderDataTable({

    validate((need(input$showmetadata, "Click 'Check metadata' to check format and display table")))

    if (input$showmetadata) {
      meta <- metadata()

      output <- cbind(Samples = rownames(meta), meta)
    }
  }, options = list(scrollX = TRUE)) #force data table to window size

  output$contents <- renderTable(include.rownames = TRUE, {
    data()
  })

  ##############################################################################

  # PCA biplots

  #biplot
  output$biplot <- renderPlot({

    #get reactive objects
    data <- data.prcomp()
    tax <- data.t()
    removenames <- input$removesamplenames
    scale.slider <- input$scale
    arrowcheck <- input$arrowcheckbox
    taxoncheck <- input$taxoncheckbox
    taxselect <- as.numeric(input$taxlevel)
    #title <- input$biplot_title
    opacity <- input$opacity_samples_pca
    sample_size <- input$size_samples_pca

    x.var <- sum(data$sdev ^ 2)
    PC1 <- paste("PC 1 Variance: %", round(sum(data$sdev[1] ^ 2) / x.var * 100, 1))
    PC2 <- paste("PC 2 Variance: %", round(sum(data$sdev[2] ^ 2) / x.var*100, 1))

    #if no data don't display anything
    if (is.null(data))
    return(NULL)

    if (input$taxoncheckbox) {
      validate(need(input$taxlevel, "Processing..."))
    }

    taxonomy <- tax$taxonomy

    #if Show taxonomy checkbox is clicked, take chosen level
    if (isTRUE(taxoncheck)) {
        taxselect <- as.numeric(input$taxlevel)
    } else {
        taxselect <- 6
    }

    #get genus (or other level)
    genus <- vapply(strsplit(as.character(taxonomy), "[[:punct:]]"), "[", taxselect, FUN.VALUE=character(1))

    #biplot colouring options
    col = c(rgb(0,0,0,opacity), rgb(0, 0, 0, 0.2))

    #do checks for arrows
    if (isTRUE(arrowcheck)) {
      arrows = FALSE
    } else {
      arrows = TRUE
    }

    #if taxonomy is is null, use periods as points for features..
    #if taxoncheckbox is checked, show genus instead of points, otherwise
      if (isTRUE(taxoncheck)) {
        points <- genus
      } else {
        points <- c(rep(".", length(dimnames(data$rotation)[[1]])))
      }

    #remove sample names
    if (isTRUE(removenames)) {
      xlabs <- c(rep(".", length(dimnames(data$x)[[1]])))
      size = c(5, 0.8)
    } else {
      xlabs = unlist(dimnames(data$x)[1])
      size = c(sample_size, 0.8)
    }

    biplot(
      data,
      main = "Principal Component Analysis",
      cex.main = 1.5,
      cex = size,
      col = col,
      scale = scale.slider,
      var.axes = arrows,
      xlab = PC1,
      ylab = PC2,
      xlabs = xlabs,
      ylabs = points
    )
  })

  #coloured biplot
  output$coloredBiplot <- renderPlot({
    #get reactive objects
    d <- data()
    data <- data.prcomp()
    tax <- data.t()
    tax.check <- data.check()
    meta <- metadata()
    cn <- column()
    taxoncheck <- input$taxoncheckbox
    taxselect <- as.numeric(input$taxlevel)
    cols <- colours()
    opacity <- input$opacity_samples_pca
    sample_size <- input$size_samples_pca

    validate(need(input$selectcolumn, "Click 'Colouring Options' to choose metadata for colouring biplot"))

    if (is.null(vals$data)) {
      #reactive metadata instead
      df <- as.data.frame(unlist(dimnames(data$x)[1]))

    } else {
      #order the metadata
      meta <- meta[order(rownames(meta)), ]

      #get the input from the pop up
      select <- vals$data

      #get rid of blank values
      s <- select[select != ""]

      #filter metatable if it was filtered
      meta <- meta[which(meta[[metaval$data]] %in% s),]

      df <- unlist(dimnames(data$x)[1])

      df <- as.data.frame(df[which(df %in% rownames(meta))])

    }

    if (input$taxoncheckbox) {
      validate(need(input$taxlevel, "Processing..."))
    }

    #if no data don't display anything
    if (is.null(meta))
    return(NULL)

    #get rid of taxonomy column if its there
    if (!isTRUE(tax.check)) {
      d$taxonomy <- NULL
    }

    colourvector <- omicplotr.colvec(data, meta, opacity, cn, type = input$colouringtype)

    omicplotr.colouredPCA(data, colourvector, scale = input$scale, arrows = input$arrowcheckbox, taxonomy = tax$taxonomy, show.taxonomy = taxoncheck, tax.level = taxselect, removenames = input$removesamplenames, names.cex=sample_size)
  })

  #histograms
  output$metahist <- renderPlot({
    #import reactive objects
    cn <- column()
    meta <- metadata()
    data <- data.prcomp()

    validate(need(input$selectcolumn, ""))

    if (is.null(meta))
    return(NULL)

    #get sample names
    df <- unlist(dimnames(data$x)[1])

    #filter rownames meta based on sample names in data
    meta <- meta[rownames(meta) %in% df,]

    #filter meta based on colnames data so it updates

    #colouring options
    if (input$colouringtype == 1) {

      validate(need(is.numeric(meta[[cn]]), "Nonnumeric data detected. Please click 'Nonnumeric metadata'"))

      #get length of unique values
      #dont laugh at my variables
      un <- unique(meta[[cn]])
      uniq <- un[!is.na(un)]
      uni <- length(uniq)

      #ramp it up
      colfunc <- colorRampPalette(c("red", "blue"))

      #make the colour
      col <- colfunc(uni)

      barplot(
        table(meta[[cn]]),
        space = 0,
        col = col,
        main = "Variable frequency",
        xlab = cn
      )
    } else if (input$colouringtype == 2) {

      validate(need(is.numeric(meta[[cn]]), "Nonnumeric data detected. Please click 'Nonnumeric metadata'"))

      #calculate quartile
      q <- quantile(meta[[cn]])

      #get freq/value
      t <- as.data.frame(table(meta[[cn]]))

      #convert factors to numeric
      t[, 1] <- as.numeric(as.character(t[, 1]))

      #add quartile rank column
      tn <-
      within(t, quartile <-
        as.integer(cut(
          t[, 1], quantile(meta[[cn]]), include.lowest = TRUE
        )))

        c <- colorRampPalette(c("red", "blue"))(4)

        #replace with colours
        tn$quartile[tn$quartile == 1] <- c[1]
        tn$quartile[tn$quartile == 2] <- c[2]
        tn$quartile[tn$quartile == 3] <- c[3]
        tn$quartile[tn$quartile == 4] <- c[4]

        barplot(
          table(meta[[cn]]),
          col = tn$quartile,
          space = 0,
          main = "Variable frequency",
          xlab = paste(cn, "by Quartiles")
        )

      } else if (input$colouringtype == 3) {
        #get frequency table
        t <- as.data.frame(table(meta[[cn]], exclude = NULL))

        #make colour column
        t$col <- NA

        #get unique features of columns in metadata
        u <- as.character(unique(meta[[cn]]))

        #colour palette
        colours <- c("indianred1", "steelblue3",  "skyblue1", "mediumorchid",
        "royalblue4", "olivedrab3", "pink", "#FFED6F", "mediumorchid3",
        "ivory2", "tan1", "aquamarine3", "#C0C0C0", "mediumvioletred",
        "#999933", "#666699", "#CC9933", "#006666", "#3399FF",
        "#993300", "#CCCC99", "#666666", "#FFCC66", "#6699CC",
        "#663366", "#9999CC", "#CCCCCC", "#669999", "#CCCC66",
        "#CC6600", "#9999FF", "#0066CC", "#99CCCC", "#999999",
        "#FFCC00", "#009999", "#FF9900", "#999966", "#66CCCC",
        "#339966", "#CCCC33", "#EDEDED"
      )

      #the tricky part here is to order the colouring based on what the coloredBiplot is doing,
      #so you need to order t based on the order of the unique values found in meta
      t <- t[match(u, t$Var1), ]

      #colour NAs black or follow color scheme
      if (NA %in% t$Var1) {
        y <- which(is.na(t$Var1))
        t$col[y] <- "black"
        for (i in (seq_len(length(t$Var1)))) {
          if (is.na(t$Var1[i])) {
            next
          } else {
            t$col[i] <- colours[i - 1]
          }
        }

        #colour NAs black
      } else {
        for (i in (seq_len(length(t$Var1)))) {
          t$col[i] <- colours[i]
        }
      }

      barplot(
        t$Freq,
        space = 0,
        col = t$col,
        main = "Number of occurences",
        names.arg = as.character(t$Var1),
        xlab = cn,
        ylab = "Frequency"
      )
    }
  })

  #generate screeplot
  output$screeplot <- renderPlot({
    data <- data.prcomp()

    #if no data don't display anything
    if (is.null(data.t()))
    return(NULL)

    #generate screeplot
    screeplot(data, main = "Screeplot of Variances", xlab = "Component")
  })

  #removed samples

  #colsums for samples
  output$colsums <- renderPlot({
    data <- data()

    if (is.null(data)) {
      return(NULL)
    }

    ab <- input$minreads

    if (is.null(data$taxonomy)){
      x <- colSums(data)

      plot(x, ylab = "Counts", xlab = "Sample Number", pch = 19, col = ifelse({x > ab}, "gray0", "red"), main = "Samples removed by filtering (count sum)")

      abline(h = ab, col = "red", lty = 2, lwd = 2)
      legend("topright", legend = c("Remaining", "Removed"), col = c("black", "red"), pch = 19)

    } else {
      x <- colSums(data[,(seq_along(data) - 1)])

      plot(x, ylab = "Counts", xlab = "Sample Number", pch = 19, col = ifelse({x > ab}, "gray0", "red"), main = "Samples removed by filtering")

      abline(h = ab, col = "red", lty = 2, lwd = 2)
      legend("topright", legend = c("Remaining", "Removed"), col = c("black", "red"), pch = 19)
    }
  })

  output$rowsums <- renderPlot({
    data <- data()

    if (is.null(data)) {
      return(NULL)
    }

    ab <- input$minsum

    if (is.null(data$taxonomy)){
      x <- rowSums(data)

      plot(x, ylab = "Counts", xlab = "Row Number", pch = 19, col = ifelse({x > ab}, "gray0", "grey"), main = "Rows removed by filtering (count sum)")

      abline(h = ab, col = "grey", lty = 2, lwd = 2)
      legend("topright", legend = c("Remaining", "Removed"), col = c("black", "grey"), pch = 19)

    } else {
      x <- rowSums(data[,(seq_along(data) - 1)])

      plot(x, ylab = "Counts", xlab = "Row Number", pch = 19, col = ifelse({x > ab}, "gray0", "grey"), main = "Rows removed by filtering")

      abline(h = ab, col = "grey", lty = 2, lwd = 2)
      legend("topright", legend = c("Remaining", "Removed"), col = c("black", "grey"), pch = 19)
    }
  })

  ##############################################################################
  #this will be implemented in a later version

  # #association plot
  #
  # output$associationplot <- renderPlot({
  #   d <- data()
  #   rhocut <- input$rhocutoff
  #
  #   if (is.null(d)) {
  #     return(NULL)
  #   }
  #   associationPlot <- function(x, cutoff = rhocut) {
  #     # Separate out the taxonomy column from the counts
  #     d <- x[, 0:(dim(x)[2] - 1)]
  #     taxon <- x[(dim(x)[2])]
  #
  #     # Generate ALDEx2 object so values are centered log ratio transformed,
  #     # then convert to propr object with a symmetric rho statistic matrix
  #     d.clr <- aldex.clr(d, mc.samples = 128, verbose = TRUE)
  #     d.sma.df <- aldex2propr(d.clr, how = "perb")
  #
  #     # Make a reference list with pairs related by a rho statistic
  #     # less than the cutoff
  #     d.sma.lo.rho <- d.sma.df["<", cutoff]
  #
  #     # **** propr and ALDEx stuff done ****
  #
  #     # igraph: Convert the connections into a graphical object using propr's
  #     # cytescape function to first generate a table of indexed pairs and
  #     # proportions
  #     g <- graph.data.frame(cytescape(d.sma.lo.rho), directed = FALSE)
  #
  #     # igraph: Find the clusters
  #     g.clust <- clusters(g)
  #
  #     # Make a table to examine the cluster membership by hand
  #     g.df <-
  #     data.frame(
  #       Systematic.name = V(g)$name,
  #       cluster = g.clust$membership,
  #       cluster.size = g.clust$csize[g.clust$membership]
  #     )
  #
  #     # Generate a set of clusters larger than some size
  #     # Minimum cluster size is 2 (obviously)
  #     big <- g.df[which(g.df$cluster.size >= 2), ]
  #     colnames(big) <- colnames(g.df)
  #
  #     # Get genera
  #     genera <- c()
  #     for (i in 1:dim(taxon)[1]) {
  #       genera <- c(genera, sapply(strsplit(as.character(taxon[i,]), "[[:punct:]]"),
  #       "[", 6))
  #     }
  #
  #     # igraph: Rename the cluster members by their genus name
  #     # gsub(pattern, replacement, strings, perl-syntax)
  #     V(g)$name <- gsub("(^[A-Za-z]{3}).+", "\\1",
  #     as.vector(genera[V(g)]), perl = TRUE)
  #
  #     plot.new()
  #
  #     # igraph:
  #     # vertex.size controls point and text color
  #     # vertex.color controls point color
  #     # vertex.frame controls point outline color
  #     plot(
  #       g,
  #       vertex.size = 5,
  #       vertex.color = rgb(0, 0, 0, 0.2),
  #       vertex.frame.color = "white"
  #     )
  #   }
  #   associationPlot(d)
  # })
  #
  # output$associationtext <- renderText({
  #   "This association plot uses Spearman's rho to measure the strength and direction of association between two variables.
  #   The cutoff values ranges from -1 to 1. If greater than 1000 features are detected, input will be required on the R console."
  # })

  ################################################################################

  # Relative abundance plots

  #dendrogram
  output$dendrogram <- renderPlot({
    x <- data.t()
    meta <- metadata()
    abund <- input$abundcutoffbarplot

    observeEvent(input$dendro_dblclick, {
      brush <- input$dendro_brush
      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)
      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }
    })

    #get filtered data if filtered
    if (is.null(vals$data)) {
      x <- x
    } else {

      x <- omicplotr.metadataFilter(x, meta, column = metaval$data, values = vals$data)
    }

    if (is.null(x$taxonomy)) {
      d <- x
    } else {
      d <- x[, seq_along(x) - 1]
      taxon <- x[(dim(x)[2])]
    }

    validate(need(x$taxonomy, "No taxonomy column detected"))

    # Get genera
    genera <- c()
    for(i in (seq_len(nrow(taxon)))) {
      genera <- c(genera, vapply(strsplit(as.character(taxon[i, ]), "[[:punct:]]"),
      "[", 6, FUN.VALUE=character(1)))
    }

    # sum counts by name
    d.agg <- aggregate(d, by=list(genera), FUN=sum)
    tax.agg <- d.agg$Group.1
    d.agg$Group.1 <- NULL

    # convert to abundances
    d.prop <- apply(d.agg, 2, function(x){x/sum(x)})

    #filters by abundance (slider bar)
    d.abund <- d.agg[apply(d.prop, 1, max) > abund,]
    tax.abund.u <- tax.agg[apply(d.prop, 1, max) > abund]

    if (any(d.abund == 0)) {
    d.abund <- t(cmultRepl(t(d.abund), label = 0, method = "CZM"))
    } else {
        d.abund <- d.abund
    }
    # get proportions of the filtered data for plotting below
    # in log-ratio speak, you are re-closing your dataset
    d.P.u <- apply(d.abund, 2, function(x){x/sum(x)})

    # order by OTU abundances
    new.order <- rownames(d.P.u)[order(apply(d.P.u, 1, sum), decreasing=T)]
    tax.abund <- tax.abund.u[order(apply(d.P.u, 1, sum), decreasing=T)]
    d.P <- d.P.u[new.order, ]
    d.clr <- apply(d.P, 2, function(x){log2(x) - mean(log2(x))})

    #distance matrix
    inp <- as.numeric(input$dismethod)
    method.m <- c("euclidean", "maximum", "manhattan")
    method.mchoice <- method.m[inp]

    dist.d.clr <- dist(t(d.clr), method=method.mchoice)

    #get clustering method for dendrogram
    num <- as.numeric(input$clustermethod)
    method <- c("complete", "single", "ward.D2")
    method.choice <- method[num]

    clust.d <- hclust(dist.d.clr, method=method.choice)

    plot.new()
    par(fig=c(0, 1, 0, 1), new=TRUE)

    plot(as.dendrogram(clust.d), main=NULL, cex=0.8, xlab="", xlim = ranges$x, ylim= ranges$y, xpd = T)
  })

  output$dendrotext <- renderText({
    "Drag mouse and double click to zoom in. \nDouble click again to zoom out"
  })

  #taoxnomic distribution
  output$barplot <- renderPlot({
    x <- data.t()
    meta <- metadata()
    abund <- input$abundcutoffbarplot

    observeEvent(input$bp_dblclick, {
      bp_brush <- input$bp_brush
      if (!is.null(bp_brush)) {
        bp_ranges$x <- c(bp_brush$xmin, bp_brush$xmax)
        bp_ranges$y <- c(bp_brush$ymin, bp_brush$ymax)
      } else {
        bp_ranges$x <- NULL
        bp_ranges$y <- NULL
      }
    })

    #get filtered data if filtered
    if (is.null(vals$data)) {
      x <- x
    } else {

      x <- omicplotr.metadataFilter(x, meta, column = metaval$data, values = vals$data)
    }

    if (is.null(x$taxonomy)) {
      d <- x
    } else {
      d <- x[, 0:(dim(x)[2] - 1)]
      taxon <- x[(dim(x)[2])]
    }

    validate(need(x$taxonomy, ""))

    # Get genera
    genera <- c()
    for(i in (seq_len(nrow(taxon)))) {
      genera <- c(genera, vapply(strsplit(as.character(taxon[i, ]), "[[:punct:]]"),
      "[", 6, FUN.VALUE=character(1)))
    }

    # sum counts by name
    d.agg <- aggregate(d, by=list(genera), FUN=sum)
    tax.agg <- d.agg$Group.1
    d.agg$Group.1 <- NULL

    # convert to abundances
    d.prop <- apply(d.agg, 2, function(x){x/sum(x)})

    d.abund <- d.agg[apply(d.prop, 1, max) > abund,]
    tax.abund.u <- tax.agg[apply(d.prop, 1, max) > abund]

    if (any(d.abund == 0)) {
    d.abund <- t(cmultRepl(t(d.abund), label = 0, method = "CZM"))
    } else {
        d.abund <- d.abund
    }

    # get proportions of the filtered data for plotting below
    # in log-ratio speak, you are re-closing your dataset
    d.P.u <- apply(d.abund, 2, function(x){x/sum(x)})

    # order by OTU abundances
    new.order <- rownames(d.P.u)[order(apply(d.P.u, 1, sum), decreasing=T)]
    tax.abund <- tax.abund.u[order(apply(d.P.u, 1, sum), decreasing=T)]
    d.P <- d.P.u[new.order, ]
    d.clr <- apply(d.P, 2, function(x){log2(x) - mean(log2(x))})

    #distance matrix
    inp <- as.numeric(input$dismethod)
    method.m <- c("euclidean", "maximum", "manhattan")
    method.mchoice <- method.m[inp]

    dist.d.clr <- dist(t(d.clr), method=method.mchoice)

    #clustering
    num <- as.numeric(input$clustermethod)
    method <- c("complete", "single", "ward.D2")
    method.choice <- method[num]
    clust.d <- hclust(dist.d.clr, method=method.choice)

    # standard colour scheme (Jean Macklaim)
    colours <- c("steelblue3", "skyblue1", "indianred", "mediumpurple1",
    "olivedrab3", "pink", "#FFED6F", "mediumorchid3",
    "ivory2", "tan1", "aquamarine3", "#C0C0C0", "royalblue4",
    "mediumvioletred", "#999933", "#666699", "#CC9933", "#006666",
    "#3399FF", "#993300", "#CCCC99", "#666666", "#FFCC66",
    "#6699CC", "#663366", "#9999CC", "#CCCCCC", "#669999",
    "#CCCC66", "#CC6600", "#9999FF", "#0066CC", "#99CCCC",
    "#999999", "#FFCC00", "#009999", "#FF9900", "#999966",
    "#66CCCC", "#339966", "#CCCC33", "#EDEDED")

    plot.new()
    par(fig=c(0,1,0,1), new = TRUE)

    par(fig=c(0,0.80,0, 1), new = TRUE)

    barplot(d.P[,clust.d$order], names.arg = clust.d$labels, space=0, xlim = bp_ranges$x, ylim= bp_ranges$y, col=colours, las=2, axisnames=T, border = NA, xpd = T)

    par(fig=c(0.8,1, 0, 1), new=TRUE)
    par(xpd = TRUE)

    leg.col <- rev(colours[seq_len(nrow(d.P))])

    legend(x="center", legend=rev(tax.abund), col=leg.col, lwd=5, cex=0.8,
    border=NULL)
  })



################################################################################

# effect plots

#find points by string and colour them
point.colour <- eventReactive(input$update_points, {
  effect <- effect_input()
  point.colour <- input$point.colour

  row.num <- grep(point.colour, rownames(effect))

  points(effect$diff.win[row.num], effect$diff.btw[row.num], pch=19, col=rgb(1,0,0,0.8), cex=1)
})

#find points by string and colour them
BA.point.colour <- eventReactive(input$update_points, {
  effect <- effect_input()
  point.colour <- input$point.colour

  row.num <- grep(point.colour, rownames(effect))

  points(effect$rab.all[row.num], effect$diff.btw[row.num], pch=19, col=rgb(1,0,0,0.8), cex=1)
})

#effect plots after calculating
output$effectMW <- renderPlot({
  x.all <- aldex.obj()

 if (input$effectplot_ab) {

        if (is.null(x.all)){
          return(NULL)
        } else {
          aldex.plot(x.all, type="MW", test="welch", all.cex = 1.5, rare.cex = 1.5, called.cex = 1.5, xlab = "Dispersion", ylab = "Difference")
          title(main = "Effect Plot")
        }
}
  })

output$effectMA <- renderPlot({
  x.all <- aldex.obj()
 if (input$effectplot_ab) {
    if (is.null(x.all)){
      return(NULL)
    } else {
      aldex.plot(x.all, type="MA", test="welch", all.cex = 1.5, rare.cex = 1.5, called.cex = 1.5, xlab = "CLR abundance", ylab = "Difference between")
      title(main = "Bland-Altman Plot")
    }
}
})

#effect plots for inputted aldex table
output$table_effect <- renderPlot({
  effect <- effect_input()
  point.colour <- input$point.colour

  if (input$effectplot_ab2) {
    if (is.null(effect)){
      return(NULL)
    } else {
      plot(effect$diff.win, effect$diff.btw, pch = 19, col=rgb(0,0,0,0.1), cex=0.4, xlab = "Difference within", ylab = "Difference between")
      title(main = "Effect Plot")
    }
    if (input$update_points == 0){
      return()
    } else {
      point.colour()
    }
  }
})

#effect plots for inputted aldex table
output$table_bland <- renderPlot({
  effect <- effect_input()
  point.colour <- input$point.colour

  if (input$effectplot_ab2) {
    if (is.null(effect)){
      return(NULL)
    } else {
      plot(effect$rab.all, effect$diff.btw, pch = 19, col=rgb(0,0,0,0.1), cex=0.4, xlab = "CLR abundance", ylab = "Difference within")
      title(main = "Bland-Altman Plot")
    }
    if (input$update_points == 0){
      return()
    } else {
      BA.point.colour()
    }
  }
})

#display information of hovered point
output$ma_hovertext <- renderUI({
  x.all <- aldex.obj()

  row <- nearPoints(x.all, input$ma_hover, xvar = "rab.all", yvar = "diff.btw")

  x <- paste("Feature Id: ", rownames(row))
  y <- paste("Median Log2 relative abundance: ", round(row$rab.all, digits =3))
  z <- paste("Between condition difference size: ", round(row$diff.btw, digits = 3))
  e <- paste("Effect size: ", round(row$effect, digits = 3))
  a <- ""
  HTML(paste(x, y, z, e, a, sep = "<br/>"))
})

#display information of hovered point
output$mw_hovertext <- renderUI({
  x.all <- aldex.obj()
  a <- d.clr()

  row <- nearPoints(x.all, input$mw_hover, xvar = "diff.win", yvar = "diff.btw", maxpoints = 1)

  if (is.null(input$mw_hover)) {
    HTML("Hover over MW plot with mouse")

  } else {

    p.v <- round(row$we.eBH, digits = 4)
    ps <- format(p.v, scientific=TRUE)
    x <- paste("Feature Id: ", rownames(row))
    y <- paste("Within condition difference size: ", round(row$diff.win, digits =3))
    z <- paste("Between condition difference size: ", round(row$diff.btw, digits = 3))
    e <- paste("Effect size: ", round(row$effect, digits = 3))
    p <- paste("Benjami Hochberg corrected p-value:", ps)
    a <- ""
    HTML(paste(x, y, z, e, p, a, sep = "<br/>"))
  }
})

#strip chart of expected CLR values for each sample per condition of hovered
# point
output$stripchart <- renderPlot({
  x <- data.t()
  cond1 <- input$group1s
  cond2 <- input$group2s
  g1s <- input$group1s
  g2s <- input$group2s
  obj <- clr.strip()
  x.all <- aldex.obj()
  meta <- metadata()
  cn <- input$colselect
  group1 <- input$group1
  group2 <- input$group2

  if (is.null(vals$data)) {
    x <- x
  } else {

    x <- omicplotr.metadataFilter(x, meta, column = metaval$data, values = vals$data)

    #meta <- meta[which(rownames(meta) %in% colnames(x)),]
}

  if (is.null(x$taxonomy)) {
    d <- x
  } else {
    d <- x[, 0:(dim(x)[2] - 1)]
    taxon <- x[(dim(x)[2])]
  }

  row <- nearPoints(x.all, input$mw_hover, xvar = "diff.win", yvar = "diff.btw", maxpoints = 1)

  feature <- rownames(row)

  if (length(feature) == 0) {
    return(NULL)
  }

  if (input$ep_chooseconds == 1) {
    cond1.clr <- obj[feature,1:cond1]
    cond2.clr <- obj[feature,(cond1+1):(cond1+cond2)]

    #make conditions vector
    conds <- c(rep("group1", g1s),
    rep("group2", g2s))
  }

  if (input$ep_chooseconds ==2) {
    #filter the meta and keep on the data which have been chosen
    group1 <- rownames(meta[which(meta[[cn]] == group1), ])
    group1.filt <- group1[group1 %in% (colnames(x))]
    data1.filt <- x[,which(colnames(x) %in% group1.filt)]
    one <- length(data1.filt)

    group2 <- rownames(meta[which(meta[[cn]] == group2), ])
    group2.filt <- group2[group2 %in% (colnames(x))]
    data2.filt <- x[,which(colnames(x) %in% group2.filt)]
    two <- length(data2.filt)

    cond1.clr <- obj[feature,1:one]
    cond2.clr <- obj[feature,(one+1):(one+two)]

    #make condition vector
    conds <- c(rep("Group 1", one),
    rep("Group 2", two))

  }

  s <- cbind(cond1.clr, cond2.clr)

  #clr values for each condition for the otu
  s.clr <- list("Condition 1" = cond1.clr, "Condition 2" = cond2.clr)


  stripchart(s.clr, main = c("Difference between conditions for feature ", feature), xlab = "Conditions", ylab = "Expected CLR values",
  col = c("black", "red"), pch =16, vertical = TRUE, method = "jitter", jitter = 0.10)
  #text(0, labels = txt)
})

#displays row name for plots when inputting ALDEx2 table
output$featurename<- renderUI({
  x.all <- effect_input()

  row <- nearPoints(x.all, input$mw_hover2, xvar = "diff.win", yvar = "diff.btw")

  ba.row <- nearPoints(x.all, input$ma_hover, xvar = "rab.all", yvar = "diff.btw")

  feature <- rownames(row)

  feature2 <- rownames(ba.row)

  number <- c(rownames(x.all[feature,]), rownames(x.all[feature2,]))
  return(number)
})

output$effectwarning <- renderText({
  "See the GitHub Wiki for more information on how to select conditions"
})
}
