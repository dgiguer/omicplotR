#script for reproducing the effect plots. values chosen for filtering in omicplotR are listed above

library(ALDEx2)
library(omicplotR)

##################################################################
#                                                                #
# You must change the file variable to your own file name        #
# and set conditions to calculate ALDEx2                         #
#                                                                #
##################################################################
#change filename.txt to your filename, make sure it follow input requirements
file <- "filename.txt"

#set which columns correspond to what conditions conditions.
#your columns must be ordered according to these conditions (view wiki).
#i.e., first three columns are time 0, last three columns are time 1.
#these conditions are for the selex dataset
conds <- c(rep("0", 7), rep("1", 7))

#read file
data <- read.table(
  file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  row.names = 1,
  check.names = FALSE,
  comment.char = "",
  na.strings = "")

#filter data
x.filt <- omicplotr.filter(data, min.reads = min.reads, min.count = min.count, min.prop = min.prop, max.prop = max.prop, min.sum = min.sum)

#check for taxonomy column
if (is.null(data$taxonomy)) {
   taxCheck <- TRUE
} else {
   taxCheck <- FALSE
}

#get rid of taxonomy column if present
if (isTRUE(taxCheck)) {
   x.filt <- x.filt
} else {
   x.filt$taxonomy <- NULL
}

#use ALDEx2 to for differential abundance analysis
#change MC or denominator if desired. More (1000-2000) MC samples is
# recommended, but takes longer to calculate.
d.clr <- aldex.clr(x.filt, mc.samples=128, conds = conds, denom = "all", verbose=TRUE)

#effect plot
aldex.plot(d.clr, type="MW", test="welch", all.cex = 1.5, rare.cex = 1.5, called.cex = 1.5, xlab = "Dispersion", ylab = "Difference")
title(main = "Effect Plot")

#Bland-Altman plot
aldex.plot(d.clr, type="MA", test="welch", all.cex = 1.5, rare.cex = 1.5, called.cex = 1.5, xlab = "CLR abundance", ylab = "Difference")
title(main = "Bland-Altman Plot")

######### save both in one pdf
# pdf("omicplotR_effect.pdf", width = 10)
# par(mfrow=c(1,2))
# aldex.plot(d.clr, type="MW", test="welch", all.cex = 1.5, rare.cex = 1.5, called.cex = 1.5, xlab = "Dispersion", ylab = "Difference")
# title(main = "Effect Plot")
# #Bland-Altman plot
# aldex.plot(d.clr, type="MA", test="welch", all.cex = 1.5, rare.cex = 1.5, called.cex = 1.5, xlab = "CLR abundance", ylab = "Difference")
# title(main = "Bland-Altman Plot")
# dev.off()
#########





















a
