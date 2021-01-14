#Function to create fit
createFit <- function(Time, CFUs, constraint, ...) {
  
  cobs::cobs(Time, log10(CFUs), constraint, print.mesg = F, print.warn = F, ...)
  
}

#Function to find roots
solveProxy <- function(fit, maxTime, viability) {
  
  initialCFUs <- 10^stats::predict(fit, 0)[2] - 1
  proxyCFUs <- initialCFUs*viability
  proxyCFUslog <- log10(proxyCFUs + 1)
  
  predictCFUs <- Vectorize(function(Time) {
    unname(stats::predict(fit, Time)[2] - proxyCFUslog)
  })
  
  rootSolve::uniroot.all(predictCFUs, c(0, maxTime))
  
}

#'@export
#'@title Calculation of Lifespan Proxies
#'@description \code{calculateProxy} will calculate a proxy for each sample, providing a single number which encapsulates lifespan.
#'@param CFUsMLE A \code{data.frame} produced by \code{\link{analyseColonyVectors}} containing CFU data from which proxies are to be calculated.
#'Must contain the column names \code{"Sample"}, \code{"Time"} and \code{"ColonyFormingUnitsPerDroplet"}.
#'@param maxTime The maximum time. Proxy solutions will be searched for in the interval between \code{0} and \code{"maxTime"}. Defaults to the maximum value in the \code{"Time"} column of \code{CFUsMLE}.
#'@param viability The proportion of viable cells. See 'Details'. Defaults to \code{0.05}.
#'@param constraint The constraint used when fitting a spline to CFU data. Must be one of \code{"none"}, \code{"increase"}, \code{"decrease"}, \code{"convex"}, \code{"concave"} or \code{"periodic"}. See \code{\link[cobs]{cobs}}. Defaults to \code{"decrease"}.
#'@param ... Additional arguments passed to \code{\link[cobs]{cobs}}.
#'@details Proxies are calculated as the square root of the amount of time take for a sample to reach a specified viability.
#'This is achieved by fitting a median spline (see \code{\link[cobs]{cobs}}) to the CFU data.
#'Specifically, \code{CFUsMLE$Time} is fitted against \code{log10(CFUsMLE$ColonyFormingUnitsPerDroplet + 1)}.
#'The proxy is the square root of the amount of time taken for the fitted values to reach the indicated viability.
#'In the case that no solutions are found within the specified time interval, the proxy is set as \code{NA}.
#'In the case that multiple solutions are found, the proxy is set as the first solution.
#'@return \code{calculateProxy} will return a \code{data.frame} containing the proxy for each sample.
#'@examples
#'#Get a data.frame showing maximum likelihood estimation as produced by analyseColonyVectors
#'CFUsMLE <- system.file("extdata", "CFUsMLE.csv", package="DeadOrAlive")
#'CFUsMLE <- read.csv(CFUsMLE)
#'
#'#Remove low quality data points
#'CFUsMLE <- CFUsMLE[CFUsMLE$TotalExclusions<=1,]
#'
#'#Plot proxy calculation for the wt (972 h-)
#'g1 <- plotProxy(CFUsMLE, "972 h-")
#'print(g1)
#'
#'#Calculate a proxy for all samples
#'proxy <- calculateProxy(CFUsMLE)
#'
#'#Add proxies to CFUsMLE
#'CFUsMLE$Proxy <- proxy$Proxy[match(CFUsMLE$Sample, proxy$Sample)]
#'
#'#Plot all lifespan curves and color by proxy
#'require(ggplot2)
#'g2 <- g2 <- ggplot(CFUsMLE, aes(Time, ColonyFormingUnitsPerDroplet + 1, group=Sample, color=Proxy))
#'g2 <- g2 + geom_point() + geom_line()
#'g2 <- g2 +  scale_y_log10() + xlab("Time (days)") + ylab("CFUs / droplet + 1")
#'print(g2)
calculateProxy <- function(CFUsMLE, maxTime=max(CFUsMLE$Time), viability=0.05, constraint="decrease", ...) {
  
  #Add 1 to CFUs - for log transformation
  CFUsMLE$ColonyFormingUnitsPerDroplet <- CFUsMLE$ColonyFormingUnitsPerDroplet + 1
  
  #Get sample names
  samples <- unique(CFUsMLE$Sample)
  
  #For each sample, subset, fit cobs, solve roots and return sqrt of solution
  proxies <- lapply(samples, function(i) {
    CFUsMLEsub <- CFUsMLE[CFUsMLE$Sample==i,]
    fit <- createFit(CFUsMLEsub$Time, CFUsMLEsub$ColonyFormingUnitsPerDroplet, constraint, ...)
    solution <- solveProxy(fit, maxTime, viability)
    return(sqrt(solution))
  })
  
  #Check how many solutions were found for each sample
  n <- sapply(proxies, length)
  
  #Take the first solution of each proxy
  proxies <- sapply(proxies, `[`, 1)
  
  #Create data.frame with results
  res <- data.frame(Sample=samples, Proxy=proxies)
  rownames(res) <- NULL
  
  #Return warning messages for cases where no solution or more than 1 solution was found
  if(any(n==0)) {
    missing <- res$Sample[n==0]
    warning(paste0(
      "No proxy solution has been found in the specified time interval for the following samples:\n",
      paste0(missing, collapse = "\n")
    ))
  }
  if(any(n>1)) {
    missing <- res$Sample[n>1]
    warning(paste0(
      "Multiple proxy solutions have been found within the specified time interval for some samples. The first solution will be used for these samples:\n",
      paste0(missing, collapse = "\n")
    ))
  }
  
  #Return
  return(res)
  
}

#'@export
#'@rdname calculateProxy
#'@description \code{plotProxy} will create a plot showing how the proxy is calculated for a specified sample.
#'@param Sample Length 1 character vector specifying which sample to plot.
#'@return \code{plotProxy} will return a \code{\link[ggplot2]{ggplot}} object.
plotProxy <- function(CFUsMLE, Sample, maxTime=max(CFUsMLE$Time), viability=0.05, constraint="decrease", ...) {
  
  #Check arguments
  if(length(Sample)!=1) {
    stop("'Sample' must be a vector of length 1.")
  }
  if(Sample %in% CFUsMLE$Sample==F) {
    stop("'Sample' must match one of the samples in the 'Sample' column of 'CFUsMLE'.")
  }
  
  #Add 1 to CFUs - for log transformation
  CFUsMLE$ColonyFormingUnitsPerDroplet <- CFUsMLE$ColonyFormingUnitsPerDroplet + 1
  
  #Subset for sample
  CFUsMLEsub <- CFUsMLE[CFUsMLE$Sample==Sample,]
  
  #Create fit
  fit <- createFit(CFUsMLEsub$Time, CFUsMLEsub$ColonyFormingUnitsPerDroplet, constraint, ...)
  
  #Find the solutions
  solution <- solveProxy(fit, maxTime, viability)
  
  #If there are multiple solutions, use only the first
  if(length(solution)>1) {
    solution <- solution[1]
    warning("Multiple proxy solutions have been found within the specified time interval. The first solution will be used.")
  }
  
  #Plot CFUs
  g <- ggplot2::ggplot(CFUsMLEsub, ggplot2::aes(Time, ColonyFormingUnitsPerDroplet)) + ggplot2::geom_point()
  g <- g + ggplot2::xlab("Time (days)") + ggplot2::ylab("CFUs / droplet + 1") + ggplot2::scale_y_log10()
  
  #Add line for  cobs fit
  fittedValues <- as.data.frame(stats::predict(fit, seq(from=0, to=maxTime, by=0.01)))
  fittedValues$fit <- 10^fittedValues$fit
  g <- g + ggplot2::geom_line(ggplot2::aes(z, fit), fittedValues, color="dodgerblue1")
  
  #Plot horizontal lines showing initial CFUs
  initialCFUs <- 10^stats::predict(fit, 0)[2]
  g <- g + ggplot2::geom_hline(yintercept = initialCFUs, lty="dashed", color="firebrick1")
  
  #Subtract 1 from initial CFUs to return to original scale
  #This is necessary to properly calculate the desired viability
  #Then add 1 to desired viability for plotting purposes and plot
  initialCFUs <- initialCFUs - 1
  proxyCFUs <- initialCFUs*viability + 1
  g <- g + ggplot2::geom_hline(yintercept = proxyCFUs, lty="dashed", color="firebrick1")
  
  #If a solution has been found, add it
  if(length(solution)==0) {
    warning("No proxy solution has been found within the specified time interval.")
  } else {
    g <- g + ggplot2::geom_vline(xintercept = solution, lty="dashed", color="firebrick1")
  }
  
  return(g)
  
}
