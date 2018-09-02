################################
### CONVERT RSE_GENE.RDATA   ###
### FROM RECOUNT PROJECT     ###         
### INTO THE ENZYMATIC GENES ###
### ACTIVITY MATRIX          ###
################################
ensembl.to.hgnc.table <- read.table("ensembl_to_hgnc_table.tab", sep="\t", stringsAsFactors = F)

get.counts.matrix <- function(path.to.rse.gene){
    require(recount)
    load(path.to.rse.gene)
    rse.s <- scale_counts(rse_gene)
    rm(rse_gene); gc(verbose = F)
    assay(rse.s)
}

select_enzymatic_genes <- function(count_matrix, gene_list, mapping_tab = ensembl.to.hgnc.table){
    unique.genes.ensemblID <- mapping_tab[as.numeric(mapping_tab$hgnc) %in% gene_list,]
    trimmed.rownames       <- do.call(rbind, strsplit(rownames(count_matrix), "\\."))[,1]
    enzymatic.genes        <- count_matrix[trimmed.rownames %in% unique.genes.ensemblID$ensembl,]
    enzymatic.genes
}

aggregated_counts <- function(count_matrix, mapping_tab = ensembl.to.hgnc.table){
    list.to.aggr <- do.call(rbind, strsplit(rownames(count_matrix), "\\."))[,1]
    head(list.to.aggr)
    counts.aggr  <- aggregate(count_matrix, by=list(list.to.aggr), FUN=mean)
    head(counts.aggr)
    target            <- counts.aggr$Group.1
    names.by.hgnc     <- as.character(mapping_tab[match(target, mapping_tab$ensembl),]$hgnc)
    head(names.by.hgnc)
    rownames(counts.aggr) <- names.by.hgnc
    head(names.by.hgnc)
    counts.aggr[,-1] # remove the Group.1 column of the ensembl name
}

create_activity_matrix <- function(count_matrix, threshold = 10){
    activity_mat_logical <- count_matrix > threshold
    activity_mat <- apply(activity_mat_logical, 2, as.numeric)
    rownames(activity_mat) <- rownames(activity_mat_logical)
    activity_mat
}

rse_data_to_activity_matrix <- function(path_to_rse, gene_list, mapping_tab = ensembl.to.hgnc.table, threshold = 10){
    rse_counts            <- get.counts.matrix(path_to_rse)
    enzymatic_counts      <- select_enzymatic_genes(rse_counts, gene_list, mapping_tab = mapping_tab)
    enzymatic_counts_hgnc <- aggregated_counts(enzymatic_counts, mapping_tab= mapping_tab)
    final_act_gene        <- create_activity_matrix(enzymatic_counts_hgnc, threshold = threshold)
    final_act_gene
}

rse_data_to_counts <- function(path_to_rse, gene_list, mapping_tab = ensembl.to.hgnc.table, threshold = 10){
    rse_counts            <- get.counts.matrix(path_to_rse)
    enzymatic_counts      <- select_enzymatic_genes(rse_counts, gene_list, mapping_tab = mapping_tab)
    enzymatic_counts_hgnc <- aggregated_counts(enzymatic_counts, mapping_tab= mapping_tab)
    enzymatic_counts_hgnc
}
