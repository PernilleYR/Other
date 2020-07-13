
myPath <- "/Users/pernillerainer/Projects/ImageRadi/RawImages"
myOutputPath <- "/Users/pernillerainer/Projects/ImageRadi/OutputImage/"

library(magick); library(grid)

#' myAmazingFunctionToArrangeTiff
#' @param pathOfImages full path to the initial image
#' @param saveTo where to put result file
#' @param nameImage name of the output image, MergeImage.pdf by default
#' @param nrow number of rows of output image, defined such as image is close to square by default
#' @param ncol number of col of output image, defined such as image is close to square by default
myAmazingFunctionToArrangeTiff <- function(pathOfImages, saveTo, 
                                           nameImage = "MergeImage.pdf", 
                                           nrow = NULL, ncol = NULL){
  #Create list of image names in pathOfImages
  filenames <- list.files(path = pathOfImages , pattern=".tif", full.names = T)
  #Import Images in a list
  myImages <- lapply(filenames, function(x) readImage(mypath = x))
  
  #If the number of row or number of col is not set then there are defined
  #such as the output is the closest to a square
  if(is.null(nrow)){
    nrow = round(sqrt(length(filenames)))
  }
  if(is.null(ncol)){
    ncol = round(sqrt(length(filenames)))
  }
  
  #Create a pdf
  pdf(paste0(saveTo,nameImage))
    #print() is only necessary since we are plotting inside a function 
    print(marrangeGrob(myImages, nrow = nrow, ncol = ncol))
  dev.off()
}
readImage <- function(mypath){
  return(rasterGrob(image_read(path = mypath)))
}

myAmazingFunctionToArrangeTiff(pathOfImages = myPath, saveTo = myOutputPath)

