## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(DeadOrAlive)
library(ggplot2)

## ----out.width="100%", echo=FALSE---------------------------------------------
knitr::include_graphics("Robots.png")

## ----out.width="100%", echo=FALSE---------------------------------------------
knitr::include_graphics("Plates.png")

## ----eval=FALSE---------------------------------------------------------------
#  #Get the directory of files to be analysed
#  dir <- system.file("extdata", "images", package="DeadOrAlive")
#  
#  #View the files
#  list.files(dir)
#  
#  #Get the reference image
#  reference <- system.file("extdata", "reference.jpg", package="DeadOrAlive")
#  
#  #Analyse the files to identify whether there is a colony or not in each position
#  #Note: This will create a new directory called 'Image_Analysis'
#  colonyThreshold(dir=dir, reference=reference)

## ----eval=FALSE---------------------------------------------------------------
#  #Get a csv file showing the identity of each sample on each plate
#  sampleReferenceFile <- system.file("extdata", "sampleReferenceFile.csv", package="DeadOrAlive")
#  View(read.csv(sampleReferenceFile))
#  
#  #Get a csv file showing the identity of each plate to be analysed
#  plateReferenceFile <- system.file("extdata", "plateReferenceFile.csv", package="DeadOrAlive")
#  View(read.csv(plateReferenceFile))
#  
#  #Get the patterns of colonies from the files (from left to right across the plate)
#  myColonyVectors <- extractColonyVectors("Image_Analysis", plateReferenceFile, sampleReferenceFile)
#  myColonyVectors

## ----eval=FALSE---------------------------------------------------------------
#  #Perform a maximum likelihood estimation of the number of viable cells
#  #Note: This will save a csv and markdown file in the current working directory
#  CFUsMLE <- analyseColonyVectors(myColonyVectors)

## ----echo=FALSE---------------------------------------------------------------
#Get a data.frame showing maximum likelihood estimation as produced by analyseColonyVectors
CFUsMLE <- system.file("extdata", "CFUsMLE.csv", package="DeadOrAlive")
CFUsMLE <- read.csv(CFUsMLE)

## -----------------------------------------------------------------------------
#Remove low quality data points
CFUsMLE <- CFUsMLE[CFUsMLE$TotalExclusions<=1,]

## -----------------------------------------------------------------------------
#Plot proxy calculation for the wt (972 h-)
g1 <- plotProxy(CFUsMLE, "972 h-")
print(g1)

## -----------------------------------------------------------------------------
#Calculate a proxy for all samples
proxy <- calculateProxy(CFUsMLE)

#Add proxies to CFUsMLE
CFUsMLE$Proxy <- proxy$Proxy[match(CFUsMLE$Sample, proxy$Sample)]

#Plot all lifespan curves and color by proxy
g2 <- g2 <- ggplot(CFUsMLE, aes(Time, ColonyFormingUnitsPerDroplet + 1, group=Sample, color=Proxy))
g2 <- g2 + geom_point() + geom_line()
g2 <- g2 +  scale_y_log10() + xlab("Time / days") + ylab("CFUs / droplet + 1") + scale_color_viridis_c()
print(g2)

