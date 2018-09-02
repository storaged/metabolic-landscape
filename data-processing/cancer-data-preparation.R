##
## For a given dataset prefix prepare a data structure for the furhter analysis
##
prepare_data <- function(cancer_name, .unique.genes = unique.genes){
    enzymes.expression <- rse_data_to_counts(
        paste("rse_data/rse_gene_", cancer_name, ".Rdata", sep=""), 
        gene_list = .unique.genes)
    colnames(enzymes.expression) <- sapply(
        colnames(enzymes.expression), 
        function(name) 
            recover.name(name, isXpref = F)
        )
    landscape.dat <- all_landscapes[[cancer_name]]
    print(names(landscape.dat))
    colnames(landscape.dat$act_mat) <- sapply(
        colnames(landscape.dat$act_mat), 
        function(name) 
            recover.name(name, isXpref = F)
        )
    colnames(landscape.dat$landscapes) <- sapply(
        colnames(landscape.dat$landscapes), 
        function(name) 
            recover.name(name, isXpref = T)
        )
    common.colnames <- intersect(intersect(colnames(enzymes.expression), 
                                 colnames(landscape.dat$landscapes)),
                                 colnames(landscape.dat$act_mat))
    
    print(colnames(enzymes.expression)[1:5])
    print(colnames(landscape.dat$landscapes)[1:5])
    print(colnames(landscape.dat$act_mat)[1:5])

    enzymes.expression       <- enzymes.expression[, common.colnames]
        print(ncol(enzymes.expression))
    
    landscape.dat$act_mat    <- landscape.dat$act_mat[, common.colnames]
        print(ncol(landscape.dat$act_mat))

    landscape.dat$landscapes <- landscape.dat$landscapes[, common.colnames]
        print(ncol(landscape.dat$landscapes))

    landscape.dat$aggregate.comp <- aggregate.by.factor(landscape.dat$landscapes, 
                                                        compartments.factor$group)
    landscape.dat$aggregate.rule <- aggregate.by.factor(landscape.dat$landscapes, 
                                                        compartments.factor$genetic.rule.rules)
    
    list(expression=enzymes.expression, landscape=landscape.dat)
}
