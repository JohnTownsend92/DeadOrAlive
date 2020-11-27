#A function to concatenate a matrix of present/absent into nested lists of vectors of present/absent for each sample at each dilution factor
matrixConcatenation <- function(matrix, well.plate.format) {
  well.plate.format <- rev(well.plate.format)
  colPerWell <- ncol(matrix)/well.plate.format[1]
  rowPerWell <- nrow(matrix)/well.plate.format[2]
  plate <- lapply(1:well.plate.format[2], function(j) {
    
    pattern <- matrix[c(((j-1)*rowPerWell+1):(j*rowPerWell)),]
    if(is.vector(pattern)) {
      pattern <- t(as.matrix(pattern))
    }
    pattern <- lapply(1:ncol(pattern), function(x) {pattern[,x]})
    pattern <- t(matrix(lapply(1:well.plate.format[1], function(i) {
      do.call("c", pattern[c(((i-1)*colPerWell+1):(i*colPerWell))])
    })))
    
  })
  plate <- do.call("rbind", plate)
  return(plate)
}

#A function to create a list of concatenated matricies for a group of .dat files
deadAliveMatrix <- function(datFiles, well.plate.format) {
  lapply(1:length(datFiles), function(i) {
    x <- utils::read.csv(datFiles[i], stringsAsFactors=F)
    matrixConcatenation(gitter2matrix(x), well.plate.format)
  })
}

#A function to concatenate multiple repeat plates into an object the same dimensions as a single plate
concatenatePlates <- function(plates, well.plate.format) {
  matrix(lapply(1:prod(well.plate.format), function(i) {do.call("c", lapply(plates, `[[`, i))}), nrow=well.plate.format[1], ncol=well.plate.format[2])
}

#Constructor function for class colonyVectors
colonyVectors <- function(sample, time, colonies, grid) {
  structure(list(sample=sample, time=time, colonies=colonies, grid=grid), class="colonyVectors")
}

#A function to take a directory of .dat files and rearrange the data into vectors of present/absent colonies for each unique combination of sample and time based on information in supporting files provided by the user.
#A list is returned containing 4 elements of the same length (these are vectors or lists themselves), with the information about each sample
#'@export
#'@title Extraction of colony patterns from analysed images
#'@description \code{extractColonyVectors} analyses a batch of DAT files produced by \code{\link{colonyThreshold}} in order to extract colony patterns for each sample.
#'@param dir Directory of DAT files. Defaults to current working directory.
#'@param plateReferenceFile Path to a csv file containing information on the identity of each plate. This file must have 3 columns as follows:
#'\itemize{
#'\item{\code{time}: Numerical values indicating the time in days.}
#'\item{\code{plate}: Characters indicating the name of each plate.}
#'\item{\code{fileName}: Full or relative file paths of each DAT file to analyse. Only DAT files listed here will be analysed.}
#'}
#'@param sampleReferenceFile Path to a csv file containing information on the samples present on each plate. Each column must correspond to a different plate, with the column names corresponding to the plates listed in the \code{plate} column of the \code{plateReferenceFile}.
#'Each sample name should be unique. The number of samples in each column should correspond to the number of samples on each plate.
#'In the case that \code{byRow} is \code{TRUE}, this should mean that there is one sample for each column on the plate. For example, if colonies were pinned from a 96-well plate, this would mean 8 samples.
#'Conversely, in the case that \code{byRow} is \code{FALSE}, this should mean that there is one sample for each row on the plate. For example, if colonies were pinned from a 96-well plate, this would mean 12 samples.
#'@param byRow Do dilution series run across the rows of the plate (i.e left to right)? Defaults to \code{TRUE}. Set to \code{FALSE} if dilution series run across columns (i.e. top to bottom). 
#'@inheritParams colonyThreshold
#'@return An object of class \code{colonyVectors}, containing the patterns extracted from each plate.
#'@template extractColonyVectorsExample
extractColonyVectors <- function(dir=getwd(), plateReferenceFile, sampleReferenceFile, plate.format=c(16, 24), well.plate.format=c(8, 12), byRow=T) {
  
  #Read in the reference and sample files
  plateReferenceFile <- utils::read.csv(plateReferenceFile, stringsAsFactors=F, header=T)
  plateReferenceFile$plate <- make.names(plateReferenceFile$plate)
  sampleReferenceFile <- utils::read.csv(sampleReferenceFile, stringsAsFactors=F, header=T)
  
  #Set the working directory
  oldDir <- getwd()
  on.exit(setwd(oldDir))
  setwd(dir)
  
  #Subset the plate reference file to contain only files present in the directory
  plateReferenceFile <- plateReferenceFile[plateReferenceFile$fileName %in% list.files(),]
  
  #Get unique combinations of sample and time
  sampletimes <- unique(plateReferenceFile[,c("plate", "time")])
  
  #Read in the .dat files and rearrange into nested lists, grouped by plates containing identical samples and timpoints
  deadAliveMatricies <- lapply(1:nrow(sampletimes), function(a) {
    deadAliveMatrix(plateReferenceFile$fileName[plateReferenceFile$plate==sampletimes$plate[a] & plateReferenceFile$time==sampletimes$time[a]], well.plate.format)
  })
  
  #Concatinate the multiple repeats into a single plate for each sample at each time
  deadAliveMatricies <- lapply(deadAliveMatricies, concatenatePlates, well.plate.format)
  
  #Create empty vectors and lists to store information about the sample, time, colony pattern and number of grid positions
  sample <- vector()
  time <- vector()
  colonies <- list()
  grid <- vector()
  
  #Create a counter to track indexing
  counter <- 0
  
  #Set the number of loops which will run for each plate based on whether the samples run by col or by row
  samplesPerPlate <- ifelse(byRow==F, well.plate.format[2], well.plate.format[1])
  
  #Split each plate based on sample, to give vectors/lists the same length as the number of unqiue samples
  for(i in 1:length(deadAliveMatricies)) {
    
    for(j in 1:samplesPerPlate) {
      
      counter <- counter + 1
      
      sample[counter] <- sampleReferenceFile[[sampletimes$plate[i]]][j]
      time[counter] <- sampletimes$time[i]
      if(byRow) {
        colonies[[counter]] <- unlist(deadAliveMatricies[[i]][j,])
      } else {
        colonies[[counter]] <- unlist(deadAliveMatricies[[i]][,j])
      }
      grid[counter] <- prod(plate.format/well.plate.format)*nrow(plateReferenceFile[plateReferenceFile$plate==sampletimes$plate[i] & plateReferenceFile$time==sampletimes$time[i],])
      
    }
    
  }
  
  #Return a list with the extracted data
  return(colonyVectors(sample, time, colonies, grid))
  
}

#'@export
length.colonyVectors <- function(x) {
  length(x$colonies)
}

#'@export
`[.colonyVectors` <- function(x, i) {
  colonyVectors(x$sample[i], x$time[i], x$colonies[i], x$grid[i])
}

colonyVectorsPrint <- function(x, n) {
  for(i in n) {
    cat("[")
    cat(i)
    cat("]\tSample: ")
    cat(x$sample[i])
    cat(", Time: ")
    cat(x$time[i])
    cat("\n\tColony Pattern: ")
    cat(compactColonies(x$colonies[[i]], x$grid[i]))
    cat("\n")
  }
}

#'@export
print.colonyVectors <- function(x, ...) {
  cat("colonyVectors with ")
  cat(length(x))
  cat(" colony pattern(s).\n")
  if(length(x) <= 11) {
    colonyVectorsPrint(x, 1:length(x))
  } else {
    colonyVectorsPrint(x, 1:5)
    cat("---\t---\n")
    colonyVectorsPrint(x, (length(x)-4):length(x))
  }
}

#'@export
c.colonyVectors <- function(...) {
  args <- list(...)
  sample <- do.call("c", lapply(1:length(args), function(i) {args[[i]]$sample}))
  time <- do.call("c", lapply(1:length(args), function(i) {args[[i]]$time}))
  colonies <- do.call("c", lapply(1:length(args), function(i) {args[[i]]$colonies}))
  grid <- do.call("c", lapply(1:length(args), function(i) {args[[i]]$grid}))
  return(colonyVectors(sample, time, colonies, grid))
}

