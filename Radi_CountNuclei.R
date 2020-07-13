################################################################
#                                                              #
#                         Count Nuclei                         #
#                                                              #
################################################################

##---------------------------------------------##
##---------------Load & set data---------------##
##---------------------------------------------##

## -- Table
Nuclei <- data.table::fread("/Users/pernillerainer/Projects/CodeRadi/CountNuclei/Nuclei.txt", 
                            stringsAsFactors = F, sep = "\t", data.table = F)$Label

## -- Values to grep
to_grep <- c()
for(i in 2:7){
  to_grep <- c(to_grep, paste0("R0",i, "-C0", 2:9), paste0("R0",i,"-C", 10:11))
};rm(i)

##---------------------------------------------##
##------------------Functions------------------##
##---------------------------------------------##
#' count
#' @description count number of time togrep appears in vector
#' @param vector the vector to count from
#' @param togrep the value to grep
#' @return a numeric, # of time togrep appears in vector
count <- function(vector, togrep){
  return(length(grep(togrep, vector)))
}

##---------------------------------------------##
##--------------------Main---------------------##
##---------------------------------------------##

## -- count
output <- sapply(to_grep, function(x) count(Nuclei, x))

## -- output as a data.frame
output <- reshape2::melt(output)

write.table(output, file = "/Users/pernillerainer/Projects/CodeRadi/CountNuclei/CountNuclei_output.txt", sep = "\t" )
