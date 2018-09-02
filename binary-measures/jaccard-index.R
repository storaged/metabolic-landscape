######################################################
## JACCARD INDEX FUNCTION                           ##
## reac - logical vector of activity states         ##
## clustering.logical - two clusters representation ##
## clust.sizes - sizes of the given clusters        ##
######################################################
jaccard.index <- function(reac, clustering.logical, 
                          clust.sizes = c(sum(clustering.logical), sum(!clustering.logical))){
    reac <- as.logical(reac)
    #message(length(clustering.logical))
    #message(length(reac))
    cl1.1 <- sum(reac[clustering.logical])
    cl1.0 <- clust.sizes[1] - cl1.1
          
    cl2.1 <- sum(reac[!clustering.logical])
    cl2.0 <- clust.sizes[2] - cl2.1
          
    #message(cl1.1, cl1.0, cl2.1, cl2.0)
    un.num <- min(cl1.1, cl2.1) + min(cl1.0, cl2.0)
    su.num <- max(cl1.1, cl2.1) + max(cl1.0, cl2.0)
    
    un.num/su.num
}
