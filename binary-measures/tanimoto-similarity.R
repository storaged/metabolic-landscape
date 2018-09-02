######################################################
## TANIMOTO SIMILARITY MEASURE                      ##
## x, y - logical vectors                           ##
######################################################
tanimoto.similarity <- function(x,y){
  l <- length(x)
  z <- sum(x!=y); 
  (l - z)/(l + z)
}
