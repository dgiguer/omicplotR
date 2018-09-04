## ----global_options, include=FALSE---------------------------------------
#makes the figure position be held to its position in the text.
knitr::opts_chunk$set(fig.pos = "h")

## ---- eval=FALSE---------------------------------------------------------
#  install.packages("BiocManager")
#  BiocManager::install("omicplotR")
#  library(omicplotR)
#  omicplotr.run()

## ---- echo = FALSE, out.width='80%', fig.align='center', fig.cap="Figure 1: Screenshot of input data page. The 'Example data' tab on the sidebar panel provides access to the provided datasets within the Shiny app."----
knitr::include_graphics("inputdata.png")

## ---- echo = FALSE, out.width='80%', fig.align='center', fig.cap="Figure 2: Example data. If taxonomy column is present, it must use the column name 'taxonomy'. Image taken from modified version of Vaginal dataset."----
knitr::include_graphics("./example_data.png")

## ---- echo = FALSE, out.width='70%', fig.align='center', fig.cap="Figure 3: Example metadata file. Metadata maybe be numerical or categorical. Any blank spaces will be replaced as NA when importing the file. Any values of T or F will be read as TRUE or FALSE. Image taken from modified version of Vaginal dataset."----
knitr::include_graphics("./example_metadata.png")

## ---- echo = FALSE, fig.align='center', fig.cap="Figure 4: Screenshot of coloured PCA biplot."----
knitr::include_graphics("col_PCA.png")

## ---- echo = FALSE, fig.align='center', out.width='50%', fig.cap="Figure 5: Pop-up generated from clicking 'Filter by metadata values'. Select a column and enter values from the column to re-plot only those values. They must match exactly. The filter can be reset or updated."----
knitr::include_graphics("metafilter.png")

## ---- echo = FALSE, fig.align='center', fig.cap="Figure 6: Screenshot of effect plots. Hovering over a feature's point in the effect plot generates a stripchart to compare the relative abundances calculated by ALDEx2 for each sample."----
knitr::include_graphics("effect.png")

