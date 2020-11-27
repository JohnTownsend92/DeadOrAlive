#A function to turn a vector of the binary colonies (present/absent) into the number present at each dilution factor
compactColonies <- function(colonies, grid) {
  colonies <- matrix(colonies, nrow=grid)
  apply(colonies, 2, function(x) {ifelse(all(is.na(x)==T), NA, sum(x, na.rm=T))})
}

#A function to calculate the expected number of cells per droplet for each grid position based on the number of cells per droplet of undiluted culture, the dilution factor and number of dilutions
nExpectedCells <- function(cellsPerDroplet, dilutionFactor, length, initialDilution) {
  sapply(1:length, function(x) {cellsPerDroplet/(initialDilution*dilutionFactor^(x-1))})
}

#A function to calculate the probability of having a colony at a particular grid position
colonyProbability <- function(cellsPerDroplet, dilutionFactor, length, initialDilution) {
  1-stats::dpois(0,nExpectedCells(cellsPerDroplet, dilutionFactor, length, initialDilution))
}

#A function to calculate the probability of observing the real data based on the estimated data for a particular cells per droplet value
binomialProbability <- function(realData, gridSize, cellsPerDroplet, dilutionFactor, length, initialDilution) {
  stats::dbinom(realData, gridSize, colonyProbability(cellsPerDroplet, dilutionFactor, length, initialDilution))
}

#A function to calculate the total likelihood of observing a pattern
totalLikelihood <- function(realData, gridSize, cellsPerDroplet, dilutionFactor, length, initialDilution, positive=T, logCellsPerDroplet=F, logLikelihood=F) {
  sign <- ifelse(positive, 1, -1)
  if(logCellsPerDroplet==T) {
    cellsPerDroplet <- exp(cellsPerDroplet)
  }
  likelihood <- prod(binomialProbability(realData, gridSize, cellsPerDroplet, dilutionFactor, length, initialDilution), na.rm=T)
  if(logLikelihood==T) {
    likelihood <- log(likelihood)
  }
  sign <- ifelse(positive, 1, -1)
  return(sign*likelihood)
}

#A function to return the total likelihood plus a value, to allow the likelihood to be solved to a particular value
solveTotalLikelihood <- function(cellsPerDroplet, realData, gridSize, dilutionFactor, length, initialDilution, positive=T, logCellsPerDroplet=F, logLikelihood=F, c) {
  totalLikelihood(realData, gridSize, cellsPerDroplet, dilutionFactor, length, initialDilution, positive, logCellsPerDroplet, logLikelihood) + c
}

#A function to produce axes labels for the dilution factors
axesLabel <- function(colonyLength, dilution) {
  graphics::axis(1, at=1:colonyLength, labels=sapply(1:colonyLength, function(x) {parse(text=paste0(dilution, "^", x-1, "*x"))}))
}

#Create an environment to store data used to plot the markdown
snapshotEnvir <- new.env(parent=emptyenv())

#A funciton to create a list of parameters which will be used to plot in the markdown. These will be stored in an internal package environment
snapshot <- function(sample, time, attempt, envir) {
  newList <- list(
    Likelihoods = get("Likelihoods", envir=envir),
    bestEstimate = get("bestEstimate", envir=envir),
    maxLikelihood = get("maxLikelihood", envir=envir),
    loglikelihoodCutOff = get("loglikelihoodCutOff", envir=envir),
    CI = get("CI", envir=envir),
    expectedCells = get("expectedCells", envir=envir),
    expectedProportion = get("expectedProportion", envir=envir),
    informativeRegion = get("informativeRegion", envir=envir),
    colonies2plot = get("colonies2plot", envir=envir),
    grid2plot = get("grid2plot", envir=envir),
    probReal = get("probReal", envir=envir),
    excludedPositions = get("excludedPositions", envir=envir),
    excludedGridPositions = get("excludedGridPositions", envir=envir),
    exclusionOrder = get("exclusionOrder", envir=envir)
  )
  snapshotEnvir[[paste0(sample, ".", time, ".", attempt)]] <- newList
  return()
}

#A function to return an error message if a vector does not match the required characteristics
checkVector <- function(Vector, Class, Length) {
  if(is.vector(Vector) == F) {
    return(T)
  }
  if(Class == "numeric") {
    if(is.numeric(Vector) == F) {
      return(T)
    }
  } else {
    if(class(Vector) != Class) {
      return(T)
    }
  }
  if(missing(Length)==F) {
    if(is.na(Length)==F) {
      if(length(Vector) != Length) {
        return(T)
      }
    }
  }
  return(F)
}

#A function to pass an error message is there is a problem with any of the parameters
errorCheck <- function(sample, time, colonies, grid, tolerance, dilution, initialDilution, minCellsPerDroplet, maxCellsPerDroplet, CIprob) {
  
  #Check the sample name
  if(missing(sample)) {
    return(paste0("Please specify 'sample' as a character vector."))
  }
  if(checkVector(sample, "character", NA)) {
    return(paste0("Please specify 'sample' as a character vector."))
  }
  
  #Check the time
  if(missing(time)) {
    return(paste0("Please specify 'time' as a vector."))
  }
  if(is.vector(time)==F) {
    return(paste0("Please specify 'time' as a vector."))
  }
  
  #Check the colonies
  if(missing(colonies)) {
    return("Please specify 'colonies' as a list of numeric vectors.")
  }
  if(is.list(colonies)==F) {
    return("Please specify 'colonies' as a list of numeric vectors.")
  }
  if(all(sapply(colonies, checkVector, "numeric")==F)==F) {
    return("Please specify colonies as a list of numeric vectors.")
  }
  
  #Check the colonies are a binary output
  if(any(sapply(colonies, function(i) {any(!(i==0 | i == 1))}))) {
    return("Please specify 'colonies' as a list of numeric vectors, with each vector being a binary input consisting of 0s and 1s.")
  }
  
  #Check the grid
  if(missing(grid)) {
    return(paste0("Please specify 'grid' as a numeric vector."))
  }
  if(checkVector(grid, "numeric", NA)) {
    return(paste0("Please specify 'grid' as a numeric vector."))
  }
  
  #Check the grid is a positive integer
  if(any(grid<=0 | grid%%1!=0)) {
    return(paste0("Please specify 'grid' as a vector of positive integers."))
  }
  
  #Check that the colonies is a multiple of grid
  if(any(sapply(colonies, length) %% grid != 0)) {
    problemElement <- min(which(sapply(colonies, length) %% grid != 0))
    return(paste0("Please make sure that the length of each vector of 'colonies' is a multiple of the respective 'grid' value.\nSee Element ", problemEmelent, " (Sample: ", sample[problemEmelent], ", Day: ", day[problemEmelent], ")."))
  }
  
  #Check that the lenghts of the uploaded data are all equal
  if(identical(length(sample), length(time)) & identical(length(sample), length(colonies)) & identical(length(sample), length(grid))==F) {
    return(paste0(
      "Please make sure that the lengths of the data are all equal ('sample', 'time', 'colonies' and 'grid'). ",
      "'sample' is a vector of length ", length(sample), ". ",
      "'time' is a vector of length ", length(time), ". ",
      "'colonies' is a list of length ", length(colonies), ". ",
      "'grid' is a vector of length ", length(grid), "."
    ))
  }
  
  #Check that each sample/day combination is unique
  if(nrow(unique(data.frame(sample, time))) != length(sample)) {
    return("Please make sure that each data point has a unqiue combination of sample name and time.")
  }
  
  #Check the tolerance
  if(checkVector(tolerance, "numeric", 1)) {
    return("Please specify 'tolerance' as a length 1 numeric vector.")
  }
  
  #Check the tolerance is a probability between 0 and 1
  if(tolerance<=0 | tolerance>=1) {
    return("Please specify 'tolerance' as a probability between 0 and 1.")
  }
  
  #Check the dilution
  if(checkVector(dilution, "numeric", 1)) {
    return("Please specify 'dilution' as a length 1 numeric vector.")
  }
  
  #Check the dilution is greater than 1
  if(dilution<=1) {
    return("Please specify 'dilution' as a number greater than 1.")
  }
  
  #Check the initial dilution
  if(checkVector(initialDilution, "numeric", 1)) {
    return("Please specify 'initialDilution' as a length 1 numeric vector.")
  }
  
  #Check the minCellsPerDroplet
  if(checkVector(minCellsPerDroplet, "numeric", 1)) {
    return("Please specify 'minCellsPerDroplet' as a length 1 numeric vector.")
  }
  
  #Check the minCellsPerDroplet is greater than 0
  if(minCellsPerDroplet<=0) {
    return("Please specify 'minCellsPerDroplet' as a number greater than 0.")
  }
  
  if(!is.null(maxCellsPerDroplet)) {
    
    #Check the maxCellsPerDroplet
    if(checkVector(maxCellsPerDroplet, "numeric", 1)) {
      return("Please specify 'maxCellsPerDroplet' as a length 1 numeric vector.")
    }
    
    #Check the minCellsPerDroplet is greater than the minimum
    if(maxCellsPerDroplet<=minCellsPerDroplet) {
      return("Please make sure 'maxCellsPerDroplet' is greater than 'minCellsPerDroplet'.")
    }
    
  }
  
  #Check the CIprob
  if(checkVector(CIprob, "numeric", 1)) {
    return("Please specify 'CIprob' as a length 1 numeric vector.")
  }
  
  #Check the CIprob is a probability between 0 and 1
  if(CIprob<=0 | CIprob>=1) {
    return("Please specify 'CIprob' as a probability between 0 and 1.")
  }
  
  return()
  
}

#A function to generate a default maximum value of cells per droplet
genMaxCells <- function(dilution, colonyLength, initialDilution) {
  max <- initialDilution*(20*(dilution^colonyLength))
  return(max)
}

#A function to generate the default series of cells per droplet values
genCellsPerDrop <- function(dilution, minCellsPerDroplet, maxCellsPerDroplet, nSteps, initialDilution) {
  return(dilution^(seq(from=log(minCellsPerDroplet, dilution), to=log(maxCellsPerDroplet, dilution), length.out=nSteps)))
}

#A function to exclude a particular well
wellExclude <- function(colonies, colonyLength, excludedPositions, grid) {
  colonies[which(rep(1:colonyLength, each=grid) %in% excludedPositions)] <- NA
  colonies
}

#A function to exclude a particular grid position
gridExclude <- function(colonies, colonyLength, excludedGridPositions, grid) {
  colonies[which(rep(1:grid, colonyLength) %in% excludedGridPositions)] <- NA
  colonies
}

#This function will provide the MLE, and create an object containing graphical parameters for the markdown file
MLE <- function(sample, time, colonies, grid, tolerance=0.001, dilution=3, initialDilution=1, minCellsPerDroplet=0.001, maxCellsPerDroplet, nSteps=1000, cellsPerDropletValues=NULL, CIprob=0.95) {
  
  #Get the length of the pattern
  colonyLength <- length(colonies)/grid
  
  #Assign a value to 'maxCellsPerDroplet' if the user has not specified one
  if(missing(maxCellsPerDroplet)) {
    maxCellsPerDroplet <- genMaxCells(dilution, colonyLength, initialDilution)
  }
  
  #Assign a value to 'cellsPerDropletValues' if the user has not specified one
  if(missing(cellsPerDropletValues)) {
    cellsPerDropletValues <- genCellsPerDrop(dilution, minCellsPerDroplet, maxCellsPerDroplet, nSteps, initialDilution)
  }
  
  #Create a vector to store the excluded positions
  excludedPositions <- numeric()
  excludedGridPositions <- numeric()
  exclusionOrder <- numeric()
  
  #Set the counter
  counter <- 0
  
  repeat{
    
    counter <- counter+1
    
    #Update the colonies and grid based on the excluded data
    colonies2plot <- colonies
    colonies2plot <- wellExclude(colonies2plot, colonyLength, excludedPositions, grid)
    colonies2plot <- gridExclude(colonies2plot, colonyLength, excludedGridPositions, grid)
    grid2plot <- grid
    grid2plot <- grid - length(excludedGridPositions)
    
    #Compact the colony vector to show the number at each position
    coloniesCompacted <- compactColonies(colonies2plot, grid)
    
    #Perform the optimisation
    if(all(coloniesCompacted==0, na.rm=T)) {
      bestEstimate <- 0
      maxLikelihood <- 1
    }
    if(all(coloniesCompacted==grid2plot, na.rm=T)) {
      bestEstimate <- Inf
      maxLikelihood <- 1
    }
    if(all(coloniesCompacted==0, na.rm=T)==F & all(coloniesCompacted==grid2plot, na.rm=T)==F) {
      #Set a sensible starting value for the optimisation
      sensibleStart <- initialDilution*dilution^(min(which(coloniesCompacted!=grid2plot))-1)
      optimisedLikelihood <- suppressWarnings(stats::optim(par=log(sensibleStart), totalLikelihood, realData=coloniesCompacted, gridSize=grid2plot, dilutionFactor=dilution, length=colonyLength, initialDilution=initialDilution, positive=F, logCellsPerDroplet=T, logLikelihood=T, method="Brent", lower=log(minCellsPerDroplet), upper=log(maxCellsPerDroplet)))
      bestEstimate <- exp(optimisedLikelihood$par)
      maxLikelihood <- exp(-optimisedLikelihood$value)
    }
    
    #Get the appropriate quantile of the chisq distribution for the desired confidence interval
    X2 <- stats::qchisq(CIprob, 1)
    
    #Find the cutoff loglikelihood
    if(all(coloniesCompacted==grid2plot, na.rm=T) | all(coloniesCompacted==0, na.rm=T)) {
      loglikelihoodCutOff <- -X2/2
    } else {
      loglikelihoodCutOff <- -optimisedLikelihood$value -X2/2
    }
    
    #Get the confidence interval
    CI <- vector()
    if(all(coloniesCompacted==grid2plot, na.rm=T) | all(coloniesCompacted==0, na.rm=T)) {
      if(all(coloniesCompacted==grid2plot, na.rm=T)) {
        CI[1] <- exp(suppressWarnings(stats::uniroot(solveTotalLikelihood, lower=log(minCellsPerDroplet), upper=log(maxCellsPerDroplet), realData=coloniesCompacted, gridSize=grid2plot, dilutionFactor=dilution, length=colonyLength, initialDilution=initialDilution, positive=T, logCellsPerDroplet=T, logLikelihood=T, c=-loglikelihoodCutOff))$root)
        CI[2] <- Inf
      }
      if(all(coloniesCompacted==0, na.rm=T)) {
        CI[1] <- 0
        CI[2] <- exp(suppressWarnings(stats::uniroot(solveTotalLikelihood, lower=log(minCellsPerDroplet), upper=log(maxCellsPerDroplet), realData=coloniesCompacted, gridSize=grid2plot, dilutionFactor=dilution, length=colonyLength, initialDilution=initialDilution, positive=T, logCellsPerDroplet=T, logLikelihood=T, c=-loglikelihoodCutOff))$root)
      }
    } else {
      CI[1] <- exp(suppressWarnings(stats::uniroot(solveTotalLikelihood, lower=log(minCellsPerDroplet), upper=log(bestEstimate), realData=coloniesCompacted, gridSize=grid2plot, dilutionFactor=dilution, length=colonyLength, initialDilution=initialDilution, positive=T, logCellsPerDroplet=T, logLikelihood=T, c=-loglikelihoodCutOff))$root)
      CI[2] <- exp(suppressWarnings(stats::uniroot(solveTotalLikelihood, lower=log(bestEstimate), upper=log(maxCellsPerDroplet), realData=coloniesCompacted, gridSize=grid2plot, dilutionFactor=dilution, length=colonyLength, initialDilution=initialDilution, positive=T, logCellsPerDroplet=T, logLikelihood=T, c=-loglikelihoodCutOff))$root)
    }
    
    #Create a vector of likelihoods as a function of 'cellsPerDroplet'
    Likelihoods <- sapply(cellsPerDropletValues, totalLikelihood, realData=coloniesCompacted, gridSize=grid2plot, dilutionFactor=dilution, length=colonyLength, initialDilution=initialDilution)
    
    #Create a vector of the expected number of cells per droplet based on the MLE
    expectedCells <- nExpectedCells(bestEstimate, dilution, colonyLength, initialDilution)
    
    #Create a vector of expected proportion of grid positions to be filled for each row based on the best estimate of cells per droplet
    expectedProportion <- colonyProbability(bestEstimate, dilution, colonyLength, initialDilution)
    
    #Create a vector containing the grid positions which are likely to contain informative informaiton
    informativeRegion <- which(expectedProportion>=0.05 & expectedProportion<=0.95)
    if(length(informativeRegion>0)) {
      informativeRegion <- c(min(informativeRegion)-1, informativeRegion, max(informativeRegion)+1)
    }
    
    #Create a vector of the probability of observing the real data against the simulated data for the best estimate for each grid position
    probReal <- binomialProbability(coloniesCompacted, grid2plot, bestEstimate, dilution, colonyLength, initialDilution)
    
    #If there are grid positions for which the likelihoods are below the tolerance, then exclude the most unlikely datapoint and retry. Otherwise, exit the loop
    if(any(probReal<tolerance, na.rm=T)) {
      
      #Create a vector to store new MLEs
      newMLE <- vector()
      
      #Calculate a new MLE by dropping each of the wells
      for(i in 1:colonyLength) {
        
        #Drop a well and compact
        coloniesCompactedNew <- compactColonies(wellExclude(colonies2plot, colonyLength, i, grid), grid)
        
        #Perform the optimisation
        if(all(coloniesCompactedNew==0, na.rm=T)) {
          newMLE[i] <- 1
        }
        if(all(coloniesCompactedNew==grid2plot, na.rm=T)) {
          newMLE[i] <- 1
        }
        if(all(coloniesCompactedNew==0, na.rm=T)==F & all(coloniesCompactedNew==grid2plot, na.rm=T)==F) {
          #Set a sensible starting value for the optimisation
          sensibleStart <- initialDilution*dilution^(min(which(coloniesCompacted!=grid2plot))-1)
          newMLE[i] <- -suppressWarnings(stats::optim(par=log(sensibleStart), totalLikelihood, realData=coloniesCompactedNew, gridSize=grid2plot, dilutionFactor=dilution, length=colonyLength, initialDilution=initialDilution, positive=F, logCellsPerDroplet=T, logLikelihood=F, method="Brent", lower=log(minCellsPerDroplet), upper=log(maxCellsPerDroplet)))$value
        }
        
      }
      
      #Calculate a new MLE by dropping each of the grid positions
      for(j in 1:grid) {
        
        #Drop a grid position and compact
        coloniesCompactedNew <- compactColonies(gridExclude(colonies2plot, colonyLength, j, grid), grid)
        grid2plotNew <- grid2plot-1
        
        #Perform the optimisation
        if(all(coloniesCompactedNew==0, na.rm=T)) {
          newMLE[i+j] <- 1
        }
        if(all(coloniesCompactedNew==grid2plotNew, na.rm=T)) {
          newMLE[i+j] <- 1
        }
        if(all(coloniesCompactedNew==0, na.rm=T)==F & all(coloniesCompactedNew==grid2plotNew, na.rm=T)==F) {
          #Set a sensible starting value for the optimisation
          sensibleStart <- initialDilution*dilution^(min(which(coloniesCompacted!=grid2plotNew))-1)
          newMLE[i+j] <- -suppressWarnings(stats::optim(par=log(sensibleStart), totalLikelihood, realData=coloniesCompactedNew, gridSize=grid2plotNew, dilutionFactor=dilution, length=colonyLength, initialDilution=initialDilution, positive=F, logCellsPerDroplet=T, logLikelihood=F, method="Brent", lower=log(minCellsPerDroplet), upper=log(maxCellsPerDroplet)))$value
        }
        
      }
      
      #Find the best improvement
      bestImprovement <- which.max(newMLE)
      if(bestImprovement<=colonyLength) {
        excludedPositions <- c(excludedPositions, bestImprovement)
        exclusionOrder <- c(exclusionOrder, 1)
      } else {
        excludedGridPositions <- c(excludedGridPositions, (bestImprovement-colonyLength))
        exclusionOrder <- c(exclusionOrder, 0)
      }
      
      #Create a list of parameters which will be passed to the markdown via an internal package environment
      snapshot(sample, time, counter, environment())
      
    } else {
      
      #Create a list of parameters which will be passed to the markdown via an internal package environment
      snapshot(sample, time, counter, environment())
      
      break
      
    }
    
  }
  
  #Return a row of a data frame
  return(data.frame(Sample = sample, Time=time, ColonyFormingUnitsPerDroplet=bestEstimate, Likelihood=maxLikelihood, LowerCI=CI[1], UpperCI=CI[2], GridSize=grid2plot, ExcludedWellPositions=paste(excludedPositions, collapse=", "), ExcludedGridPositions=paste(excludedGridPositions, collapse=", "), TotalExclusions=length(exclusionOrder), stringsAsFactors=F))
  
}

#A function to peform the MLE on multiple samples and produce a html report
#'@export
#'@title Maximum likelihood estimation of colony forming units per droplet from colony patterns
#'@description \code{analyseColonyVectors} analyses the colony patterns extracted using the \code{\link{extractColonyVectors}} function in order to estimate the number of colony forming units per droplet from the original liquid culture.
#'This is done by performing a maximum likelihood estimation on the pattern using Brent optimisation (see \code{\link[stats]{optim}}).
#'@param colonyVectors An object of class \code{colonyVectors} created by \code{\link{extractColonyVectors}} containing the extracted colony patterns.
#'@param tolerance The minimum probability of observing data at each dilution factor to be considered acceptable.
#'If the probability of observing the data at any of the dilution factors falls below this value, the pattern will be considered anomalous.
#'In this case, an algorithm will re-attempt the maximum likelihood estimation with excluded data points, until a suitable solution is found.
#'Any data points which are removed are shown in the output table and markdown files.
#'@param dilution The dilution factor used in the the experiment. Defaults to \code{3}.
#'@param initialDilution The initial dilution factor. Defaults to \code{1} (i.e. undiluted).
#'@param minCellsPerDroplet The minimum number of viable cells per droplet to test in the optimisation. Defaults to \code{0.001}.
#'@param maxCellsPerDroplet The maximum number of viable cells per droplet to test in the optimisation. In the default case (\code{NULL}), a number is selected which should cover the range expected in the experiment based on the \code{minCellsPerDroplet}, \code{dilution} and length of the colony pattern.
#'@param CIprob The probabilty of the confidence intervals to be returned. Defaults to \code{0.95}.
#'@param save.table Should a csv file of the results be saved. Defaults to \code{TRUE}.
#'@param table.name Name of the csv file.
#'@param save.markdown Should a markdown (html) document showing the results of the maximum likelihood estimation be saved. Defaults to \code{TRUE}.
#'@param markdown.name Name of the markdown file.
#'@param save.directory Directory to save the csv and markdown file. Defaults to the current working directory.
#'@return \code{analyseColonyVectors} returns a \code{data.frame} containing the results of the maximum likelihood estimation for each sample at each time.
#'The table contains the following columns:
#'\item{Sample}{The sample name.}
#'\item{Time}{The time in days.}
#'\item{ColonyFormingUnitsPerDroplet}{The maximum likelihood estimation of the number of viable cells per droplet.}
#'\item{Likelihood}{The value of the maximum likelihood (i.e. the probability of observing the pattern provided in the case that the maximum likelihood estimation is correct).}
#'\item{LowerCI}{The lower confidence interval for the estimation of viable cells per droplet.}
#'\item{UpperCI}{The upper confidence interval for the estimation of viable cells per droplet.}
#'\item{GridSize}{The size of the grid (i.e. the number of grid positions for each dilution factor).}
#'\item{ExcludedWellPositions}{Any well positions which were excluded. This is as a result of an attempted optimisation not passing the \code{tolerance} threshold.}
#'\item{ExcludedGridPositions}{Any grid positions which were excluded. This is as a result of an attempted optimisation not passing the \code{tolerance} threshold.}
#'\item{TotalExclusions}{The total number of excluded well and grid positions.}
#'In addition, a csv file containing this data and a markdown (html) report showing the maximum likelihood optimisation for each colony pattern are saved by default.
#'@template extractColonyVectorsExample
#'@examples
#'
#'#Perform a maximum likelihood estimation of the number of viable cells
#'#Note: This will save a csv and markdown file in the current working directory
#'\dontrun{
#'CFUsMLE <- analyseColonyVectors(myColonyVectors)
#'}
analyseColonyVectors <- function(colonyVectors, tolerance=0.001, dilution=3, initialDilution=1, minCellsPerDroplet=0.001, maxCellsPerDroplet=NULL, CIprob=0.95, save.table=T, table.name="CFUsMLE.csv", save.markdown=T, markdown.name="CFUsMLE.html", save.directory=getwd()) {
  
  #Clear snapshot environment upon exit
  on.exit({
    rm(list=ls(envir = snapshotEnvir), envir = snapshotEnvir)
  })
  
  #Get information from the colonyVectors
  sample=colonyVectors$sample
  time=colonyVectors$time
  colonies=colonyVectors$colonies
  grid=colonyVectors$grid
  
  #Get the file path to the markdown template
  markdownLocation <- system.file("rmd", "MLE.rmd", package="DeadOrAlive")
  
  #Check the arguments
  errorMessage <- errorCheck(sample, time, colonies, grid, tolerance, dilution, initialDilution, minCellsPerDroplet, maxCellsPerDroplet, CIprob)
  if(is.null(errorMessage)==F) {
    stop(errorMessage)
  }
  
  #Get the length of the pattern
  colonyLength <- max(sapply(colonies, length)/grid)
  
  #Assign a value to 'maxCellsPerDroplet' if the user has not specified one
  if(missing(maxCellsPerDroplet)) {
    maxCellsPerDroplet <- genMaxCells(dilution, colonyLength, initialDilution)
  }
  
  #Assign a value to 'cellsPerDropletValues' (for plotting purposes)
  cellsPerDropletValues <- genCellsPerDrop(dilution, minCellsPerDroplet, maxCellsPerDroplet, 1000, initialDilution)
  
  #Perform the MLE
  Viability <- data.frame()
  for(i in 1:length(sample)) {
    Viability <- rbind(Viability, MLE(sample[i], time[i], colonies[[i]], grid[i], tolerance, dilution, initialDilution, minCellsPerDroplet, maxCellsPerDroplet, nSteps, cellsPerDropletValues, CIprob))
    cat("Analysed ", sample[i], " day ", time[i], ". Analysis ", 100*(round(i/length(sample), 2)), "% complete.\n", sep="")
  }
  
  #Remember to return to old working directory
  oldDir <- getwd()
  on.exit(setwd(oldDir), add = T)
  
  if(save.table) {
    
    #Save the table
    setwd(save.directory)
    utils::write.csv(Viability, table.name, row.names=F)
    
  }
  
  if(save.markdown) {
    
    setwd(save.directory)
    
    cat("\nGenerating html report - this may take a while...")
    
    #Create a html report
    rmarkdown::render(markdownLocation, envir=new.env(parent=environment()), output_dir=save.directory, output_file=markdown.name)
    
  }
  
  return(Viability)
  
}
