###
## Landscape generation from activity matrices
###

generate.landscapes <- function(gene.expresion.matrix){
    message("Number of landscapes calculated:")
    no.samples <- ncol(gene.expresion.matrix)
    i <<- 0
    landscapes <- apply(gene.expresion.matrix, 2, function(act){
                i <<- i + 1
                act        <- as.numeric(as.character(act)); 
                names(act) <- rownames(gene.expresion.matrix);
                if(i %% 10 == 1)
                    message(paste(i, "of", no.samples, sep = "\t"))
                as.numeric(determine.landscape(act))
              })
    colnames(landscapes) <- paste("X", colnames(landscapes), sep="")
    rownames(landscapes) <- paste("i", which(recon$rules != ""), sep="")
    landscapes
}
