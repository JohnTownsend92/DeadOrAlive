#'@examples
#'#Get a csv file showing the identity of each plate to be analysed
#'plateReferenceFile <- system.file("extdata", "plateReferenceFile.csv", package="DeadOrAlive")
#'
#'#Get a csv file showing the identity of each sample on each plate
#'sampleReferenceFile <- system.file("extdata", "sampleReferenceFile.csv", package="DeadOrAlive")
#'
#'#Get the directory of files to be analysed
#'dir <- system.file("extdata", "Image_Analysis", package="DeadOrAlive")
#'
#'#Get the patterns of colonies from the files (from left to right across the plate)
#'myColonyVectors <- extractColonyVectors(dir=dir, plateReferenceFile, sampleReferenceFile)
