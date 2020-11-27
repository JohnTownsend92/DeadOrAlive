#A function to create a matrix from a gitter file
gitter2matrix <- function(gitterFile) {
  matrix(gitterFile$present, ncol=max(gitterFile$col), byrow=T)
}

#A function to plot a heatmap of colonies
colonyHeatmap <- function(matrix) {
  gplots::heatmap.2(matrix, Rowv = NULL, Colv = NULL, breaks=c(0, 0.5, 1), col=c("dodgerblue3", "firebrick2"), dendrogram = "none", trace="none", notecol = "White", sepcolor="lightblue", colsep=1:ncol(matrix), rowsep=1:nrow(matrix), key=F, margins=c(0,0), lwid=c(1,1000), lhei=c(1,1000))
}

#A function to take a matrix of present/absent and collapse into how many present at each well
matrix2countMatrix <- function(matrix, outputDimensions) {
  outputDimensions <- rev(outputDimensions)
  colPerWell <- ncol(matrix)/outputDimensions[1]
  rowPerWell <- nrow(matrix)/outputDimensions[2]
  t(sapply(1:outputDimensions[2], function(j) {
    
    pattern <- matrix[c(((j-1)*rowPerWell+1):(j*rowPerWell)),]
    if(is.vector(pattern)) {
      pattern <- t(as.matrix(pattern))
    }
    pattern <- apply(pattern, 2, sum)
    pattern <- sapply(1:outputDimensions[1], function(i) {
      sum(pattern[c(((i-1)*colPerWell+1):(i*colPerWell))])
    })
    
  }))
}

#A function to plot a heatmap of colony counts for each dilution factor
count_heatmap <- function(countMatrix, nPins, boxWidth, colonyMatrix) {
  gplots::heatmap.2(countMatrix, Rowv = NULL, Colv = NULL, breaks=seq(0, nPins, length.out=nPins+2), scale="none", dendrogram = "none", trace="none", col=grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral")[-c(1, 6, 10:11)])(nPins+1), notecol = "white", sepcolor="white", colsep=1:ncol(countMatrix), rowsep=1:nrow(countMatrix), sepwidth=0.05*rev(dim(countMatrix)/dim(colonyMatrix)), key=F, margins=c(0,0), lwid=c(1,1000), lhei=c(1,1000))
}
count_heatmap_grid <- function(countMatrix, nPins, boxWidth, colonyMatrix) {
  gplots::heatmap.2(countMatrix, Rowv = NULL, Colv = NULL, breaks=seq(0, nPins, length.out=nPins+2), scale="none", dendrogram = "none", trace="none", col=rep("white", nPins+1), notecol = "darkgray", sepcolor="darkgray", colsep=1:ncol(countMatrix), rowsep=1:nrow(countMatrix), sepwidth=0.05*rev(dim(countMatrix)/dim(colonyMatrix)), cellnote = countMatrix, notecex=0.04*boxWidth, key=F, margins=c(0,0), lwid=c(1,1000), lhei=c(1,1000))
}

#A function to calulate the total linear distance of how well the equally spaced gird with a particular offset aligns with the real grid
gridFit <- function(realGrid, newGrid, offset) {sum(abs(realGrid - (newGrid+offset)))}

#A function to prefix a plus on any number greater than or equal to 0
with_plus <- function(x, ...) {
  if (x >= 0)
  {
    sprintf(
      fmt = "+%s", 
      format(x, ...)
    )
  }
  else
  {
    as.character(x)
  }
}

#Cross correlation 2d registration - change from gitter default
.register2d <- function(im, r2, c2, lag.max) {
  return(im)
}

#A function to process images through gitter, threshold colony sizes and produce plots of plates overlayed with heatmaps
#'@export
#'@title Binary thresholding of gridded colony images
#'@description \code{colonyThreshold} analyses a batch of images of gridded colonies and quantifies whether or not there is a colony in each position. Grid identification and colony size quantification is done using the \code{\link[gitter]{gitter.batch}} function, which is then followed by additional thresholding to classify whether a colony is present or absent at each position.
#'@param dir Directory of images. Defaults to current working directory.
#'@param fileFormat File format of the images. Defaults to \code{"jpg"}.
#'@param files Character vector of file names to be analysed. Defaults to all files of the specified \code{fileFormat} in the specified \code{dir}. Relative file paths to images in sub directories or full paths to images in other directories can also be specified.
#'@param reference File path of the image to be used as a grid reference.
#'@param outputFolder Logical. Should output files from the analysis be saved in a new directory. Defaults to \code{TRUE}.
#'@param plate.format Integer vector of length 2 specifying the number of rows and columns on the agar target plate, respectively. Defaults to \code{c(16, 24)}.
#'@param inverse Logical. Is the image colour inverted (i.e. colonies are darker compared to background). Defaults to \code{TRUE}.
#'@param grid.save Logical. Should a gridded image be saved. Defaults to \code{TRUE}. This image is produced by \code{\link[gitter]{gitter}}.
#'@param threshold.save Logical. Should an image be saved showing which grid positions have been thesholded as colony present. Defaults to \code{TRUE}.
#'@param count.save Logical. Should an image be saved showing how many colonies are present at each dilution factor for each sample. Defaults to \code{TRUE}.
#'@param well.plate.format Integer vector of length 2 specifying the number of wells and columns of the liquid source plate, respectively. Dimensions must be factors of \code{plate.format} Defaults to \code{c(8, 12)}.
#'@inheritParams gitter::gitter
#'@return
#'\code{colonyThreshold} does not recturn any objects, but instead creates up to 4 files for each image analysed, as described below.
#'By default (\code{outputFolder=T}), files will be saved in a sub-directory called \code{Image_Analysis} within the current working directory.
#'\item{DAT file}{
#'  A tab delimited file containing quantified colony sizes, as described in \code{\link[gitter]{gitter}}.
#'  An sixth column marking whether a colony has been classified as present (1) or absent (0) is added.
#'}
#'\item{Gridded image}{
#'  As described in \code{\link[gitter]{gitter}}.
#'}
#'\item{Threshold image}{
#'  Thresholded image showing whether a colony has been classified as present or absent for each position on the plate.
#'}
#'\item{Count image}{
#'  Image showing how many colonies have been classified as present for each sample at each dilution factor.
#'}
#'@examples
#'#Get the directory of files to be analysed
#'dir <- system.file("extdata", "images", package="DeadOrAlive")
#'
#'#Get the reference image
#'reference <- system.file("extdata", "reference.jpg", package="DeadOrAlive")
#'
#'#Analyse the files to identify whether there is a colony or not in each position
#'#Note: This will create a new directory called 'Image_Analysis'
#'\dontrun{
#'colonyThreshold(dir=dir, reference=reference)
#'}
colonyThreshold <- function(dir=getwd(), fileFormat="jpg", files=list.files(path=dir, pattern=paste0("\\.", fileFormat), full.names = T), reference, outputFolder=T, plate.format=c(16, 24), verbose="p", remove.noise=F, autorotate=F, inverse=T, contrast=NULL, grid.save=T, threshold.save=T, count.save=T, well.plate.format=c(8, 12)) {
  
  #Check the file format
  if(fileFormat %in% c("jpg", "tiff") == F) {stop("Error: Please specify 'fileFormat' as 'jpg' or 'tiff'.")}
  
  #Check that the dimensions of the plate format are multiples of the dimensions of the well plate format
  if(!all(plate.format %% well.plate.format == 0)) {stop("Error: The dimensiosn of 'plate.format' are not multiples of 'well.plate.format'.")}
  
  #If files are to be saved in a new directory, create the directory
  if(outputFolder) {
    if("./Image_Analysis" %in% list.dirs()==F) {
      dir.create("Image_Analysis")
    }
    dir.output <- paste0(getwd(), "/Image_Analysis")
  } else {
    dir.output <- getwd()
  }
  
  #Set if the grid will be saved
  if(grid.save) {
    grid.save <- dir.output
  } else {
    grid.save <- NULL
  }
  
  #Set new version of .regester2d in gitter. Put old version back on exit
  #This will alter the way gitter::gitter.back processes the images and provide more stability by removing the ability to re-align the grid
  .register2dOld <- gitter:::.register2d
  utils::assignInNamespace(".register2d", .register2d, "gitter")
  on.exit({
    .register2d <- .register2dOld
    utils::assignInNamespace(".register2d", .register2d, "gitter")
  }, add = T)
  
  #Run the gitter analysis
  gitter::gitter.batch(files, reference, verbose=verbose, plate.format=plate.format, grid.save=grid.save, dat.save=dir.output, remove.noise=remove.noise, autorotate=autorotate, inverse=inverse, contrast=contrast)
  
  #Get colony coordinates for the reference plate
  coords <- gitter::gitter(reference, plate.format=plate.format, verbose="n", grid.save=NULL, dat.save=NULL, remove.noise=remove.noise, autorotate=autorotate, inverse=inverse, contrast=contrast, .is.ref=T)
  coords <- attr(coords, "params")
  coords <- list(x=unique(coords$coords$x), y=unique(coords$coords$y), dims=c(length(coords$col.sums), length(coords$row.sums)))
  
  #Get the average distance between each coordinate for the x and y axes
  xDistAvg <- mean(sapply(1:(length(coords$x)-1), function(i) {coords$x[i+1]-coords$x[i]}))
  yDistAvg <- mean(sapply(1:(length(coords$y)-1), function(i) {coords$y[i+1]-coords$y[i]}))
  
  #Get the pixel dimensions of the heatmap to be plotted
  heatmapDims <- round(c(xDistAvg*(length(coords$x)), yDistAvg*(length(coords$y))))
  
  #Create coordinates of evenly spaced grid starting at (1,1)
  xStarting <- seq(from=1, by=xDistAvg, length.out=length(coords$x))
  yStarting <- seq(from=1, by=yDistAvg, length.out=length(coords$y))
  
  #Optimise the offset to minimise the linear distance between the equally spaced grid and the real grid for the x and y directions
  xOffset <- stats::optimize(gridFit, c(0, coords$dims[1]-xStarting[length(xStarting)]), realGrid=coords$x, newGrid=xStarting)$minimum
  yOffset <- stats::optimize(gridFit, c(0, coords$dims[2]-yStarting[length(yStarting)]), realGrid=coords$y, newGrid=yStarting)$minimum
  
  #Calculate the offsets for the heatmaps
  heatmapOffsets <- round(c(xOffset, yOffset) - 0.5*c(xDistAvg, yDistAvg))
  heatmapOffsets <- sapply(heatmapOffsets, with_plus)
  heatmapOffsets <- paste(heatmapOffsets, collapse="")
  
  #Split file names by folder
  filesSplit <- strsplit(files, "/")
  
  #Get file names without directories
  filesWithoutDirectories <- sapply(filesSplit, function(x) {x[length(x)]})
  
  #Create a list of the files create in the analysis
  datFiles <- paste0(filesWithoutDirectories, ".dat")
  thresholdFiles <- paste0("threshold_", filesWithoutDirectories)
  countFiles <- paste0("count_", filesWithoutDirectories)
  
  #Edit the file names to include the sub folder if necessary
  if(outputFolder) {
    thresholdFiles <- paste0("Image_Analysis/", thresholdFiles)
    countFiles <- paste0("Image_Analysis/", countFiles)
  }
  
  #Read all .dat files
  datFilesPath <- datFiles
  if(outputFolder) {
    datFilesPath <- paste0("Image_Analysis/", datFilesPath)
  }
  datFilesRead <- lapply(datFilesPath, gitter::gitter.read)
  
  #Create vector of all colony sizes
  colonySizes <- unlist(lapply(datFilesRead, `[`, "size"))
  
  #Log transform data
  logColonySizes <- log(colonySizes+1)
  
  #Create a vector of breaks (thresholds) to test
  uniqueSizes <- unique(logColonySizes)
  uniqueSizes <- uniqueSizes[order(uniqueSizes)]
  breaks <- zoo::rollmean(uniqueSizes, 2)
  
  #Get the total number of data points
  counts <- length(logColonySizes)
  
  #Calculate the between class variance for each possible threshold
  betweenClassVarianceVector <- sapply(breaks, function(i) {
    absent <- logColonySizes<i
    nAbsent <- sum(absent)
    weightAbsent <- nAbsent/counts
    weightPresent <- 1-weightAbsent
    meanAbsent <- mean(logColonySizes[absent])
    meanPresent <- mean(logColonySizes[!absent])
    betweenClassVariance <- weightAbsent*weightPresent*((meanAbsent-meanPresent)^2)
  })
  
  #Set the threshold as the break value which maximises the between class variance
  otsuThreshold <- breaks[which.max(betweenClassVarianceVector)]
  otsuThreshold <- exp(otsuThreshold) + 1
  
  #For each file
  for(i in 1:length(datFiles)) {
    
    #Read the gitter files and skip to the next loop if not there
    if(datFiles[i] %in% list.files(dir.output)) {
      gitterData <- datFilesRead[[i]]
    } else {
      next
    }
    
    #Add a column to indicate whether the colony is classified as preset or absent
    gitterData$present <- ifelse(gitterData$size>otsuThreshold, 1, 0)
    
    #Quality check - for those with "S" flags (spillage), remove colonies with low circularity
    gitterData$present <- ifelse(grepl("S", gitterData$flags), ifelse(gitterData$size*gitterData$circularity<otsuThreshold, 0, gitterData$present), gitterData$present)
    
    #Save the updated gitter file
    datFile <- datFiles[i]
    datFile <- ifelse(outputFolder, paste0("Image_Analysis/", datFile), datFile)
    utils::write.csv(gitterData, datFile, row.names=F)
    
    #If any images are being saved
    if(threshold.save | count.save) {
      
      #Convert the gitter file to a matrix showing present or absent
      colonyMatrix <- gitter2matrix(gitterData)
      
      #Read in the raw image of the plate
      plateImage <- magick::image_read(files[i])
      
      #If image is not already inverted, invert it
      if(inverse==F) {
        plateImage <- magick::image_negate(plateImage)
      }
      
    }
    
    #If the thresholded image is being saved
    if(threshold.save) {
      
      #Plot a heatmap of the thresholded values
      grDevices::jpeg(thresholdFiles[i], width=heatmapDims[1], height=heatmapDims[2])
      colonyHeatmap(colonyMatrix)
      grDevices::dev.off()
      
      #Read in the image of the threshold heatmap
      thresholdHeatmap <- magick::image_read(thresholdFiles[i])
      
      #Overlay the heatmap and the original image
      thresholdOutput <- magick::image_composite(plateImage, composite_image=thresholdHeatmap, operator = "color-burn", offset=heatmapOffsets)
      magick::image_write(thresholdOutput, thresholdFiles[i], fileFormat)
      
      #Remove
      rm(list=c("thresholdOutput", "thresholdHeatmap"))

    }
    
    #If count image is being saved
    if(count.save) {
      
      #Invert the image and increase contrast
      plateImage <- magick::image_negate(plateImage)
      plateImage <- magick::image_contrast(plateImage, 1)
      
      #Convert to a matrix of summed values
      countMatrix <- matrix2countMatrix(colonyMatrix, well.plate.format)
      
      #Calculate the number of times each well plate was pinned
      nPins <- prod(plate.format/well.plate.format)
      
      #Plot a heatmap of the thresholded values
      grDevices::jpeg(countFiles[i], width=heatmapDims[1], height=heatmapDims[2])
      count_heatmap(countMatrix, nPins, min(yDistAvg*(plate.format/well.plate.format)[1], xDistAvg*(plate.format/well.plate.format)[2]), colonyMatrix)
      grDevices::dev.off()
      
      #Read in the image of the count heatmap
      countHeatmap <- magick::image_read(countFiles[i])
      
      #Plot a heatmap of the grid overlay
      grDevices::jpeg(countFiles[i], width=heatmapDims[1], height=heatmapDims[2])
      count_heatmap_grid(countMatrix, nPins, min(yDistAvg*(plate.format/well.plate.format)[1], xDistAvg*(plate.format/well.plate.format)[2]), colonyMatrix)
      grDevices::dev.off()
      
      #Read in the image of the grid overlay
      countHeatmapGrid <- magick::image_read(countFiles[i])
      
      #Add a border
      borderDist <- round(0.025*min(yDistAvg*(plate.format/well.plate.format)[1], xDistAvg*(plate.format/well.plate.format)[2]))
      countHeatmapGrid <- magick::image_crop(countHeatmapGrid, paste0(heatmapDims[1]-2*borderDist, "x", heatmapDims[2]-2*borderDist, "+", borderDist, "+", borderDist))
      countHeatmapGrid <- magick::image_border(countHeatmapGrid, "darkgray", paste0(borderDist, "x", borderDist))
      
      #Make white channel transparent and invert (and disable output to console from 'magick::image_transparent')
      invisible(utils::capture.output(countHeatmapGrid <- magick::image_transparent(countHeatmapGrid, "white", 10)))
      countHeatmapGrid <- magick::image_negate(countHeatmapGrid)
      
      #Overlay the heatmap and the original image
      countOutput <- magick::image_composite(plateImage, composite_image=countHeatmap, operator = "multiply", offset=heatmapOffsets)
      countOutput <- magick::image_composite(countOutput, composite_image=countHeatmapGrid, operator = "Atop", offset=heatmapOffsets)
      magick::image_write(countOutput, countFiles[i], fileFormat)
      
      #Remove
      rm(list=c("countOutput", "countHeatmap", "countHeatmapGrid"))
      
    }
    
    #If any images are being saved
    if(threshold.save | count.save) {
      
      #Remove
      rm("plateImage")
      
      #Empty garbage
      gc(full=T)
      
    }
    
  }
  
  return()
  
}
