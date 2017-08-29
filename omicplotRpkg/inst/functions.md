### Functions/workflow for omicplotR

###### Filter the data table by counts
`omicplotr.filter()`
* Use functions from CoDaSeq to filter the data table. This will pretty much be the codaSeq.filter function.
* Input: data table, min.reads, min.count, min.prop, min.occurence
* output: filtered data frame with features or samples removed

###### Convert to CLR and generate prcomp object for biplot
`omicplotr.clr()`
* gets CLR object, filters by variance, generate prcomp object for biplot (later)
* input: data frame from `omicplotr.filter`, variance
* output: prcomp class object for biplot.

##### Create vector of colours for colouredBiplot
`omicplotr.colvec()`
* takes data + metadata file and creates vector of colour for colouredbiplot
* input: data, metadata, column, type of data (1-3)
* output: list of colours for use in coloured biplot that will colour the sample names according to metadata.

##### Create PCA plot
`omicplotr.colouredPCA()`
* creates colouredBiplot from data, colourvector
* calculates variance of the components
* input: data, vector from `omicplotr.colvec`
* output: PCA biplot

##### Show removed features and samples
`omicplotr.getRemovedFeatures()` and `omicplotr.getRemovedSamples()`
* generates table of removed features and/or samples
* input: unfiltered data, filtered data by `omicplotr.filter`
* output: table of removed data (either features or samples)

##### Filter data based on metadata values
`omicplotr.metadataFilter()`
* filters data based on selected metadata
* input: unfiltered data, column, values
* output: data frame of filtered data

##### Generate downloadable report for ALDEx2 results
`omicplotr.report()`
* generates report of anosim p-value, PERMANOVA

##### Calculate anosim p-value
`omciplotr.anosim()`
* calculate anosim p-value using `vegan` between conditions
* input: (filtered) data from `omicplotr.filter`, vector of conditions
* output: anosim p-value for different groups
