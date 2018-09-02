#####################################################################################
## L function, that selects the features introducing the highest variance to data. ##
#####################################################################################

select.attributes.by.rotation <- function(rotation.mat, top.from.pc = 1, number.of.pcs = 2){
    lapply(1:number.of.pcs, function(pc.num){
        rot.vec <- abs(rotation.mat[, pc.num])
        threshold <- sort(unique(rot.vec), decreasing = T)[top.from.pc]
        which(rot.vec >= threshold)
    })
}

