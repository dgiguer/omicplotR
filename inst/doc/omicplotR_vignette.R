## ---- echo =FALSE--------------------------------------------------------
#set the height and width globally
knitr::opts_chunk$set(fig.width=6, fig.height=6, fig.align = 'center') 

## ------------------------------------------------------------------------
library(omicplotR)
data(otu_table)
data(metadata)

## ------------------------------------------------------------------------
d.mf <- omicplotr.metadataFilter(data = otu_table, meta = metadata, 
                                 column = "probio", values = c("y", "n"))

## ------------------------------------------------------------------------
d.f <- omicplotr.filter(data = d.mf, min.reads = 10)

## ------------------------------------------------------------------------
d.clr <- omicplotr.clr(data = d.f)

## ------------------------------------------------------------------------
biplot(d.clr, cex.main = 1.5, cex = c(0.8, 0.4))


## ------------------------------------------------------------------------
colvec <- omicplotr.colvec(data = d.clr, meta = metadata, column = "probio", type = 3)

## ------------------------------------------------------------------------
omicplotr.colouredPCA(d.clr, colourvector = colvec, main = "")

