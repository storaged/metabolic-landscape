# compartments.factor is a factor that groups genetic rules with respect to compartments they are applied
# Rules of grouping: 
# a : reaction takes place in compartment a
# a>b : reaction goes from a to b
# a.b>c : reaction goed from compartments a and b to compartment c
compartments.factor <- apply(enzymatic.reactions, 1, function(reaction){
        reactants_left <- metabolites_from_reaction(as.character(reaction[2]), side='left') 
        reactants_right <- metabolites_from_reaction(as.character(reaction[2]), side='right') 
        from.m <- sort(unique(sapply(reactants_left, function(r){ q <- unlist(strsplit(x = r, "_")); q[length(q)] })))
        to.m   <- sort(unique(sapply(reactants_right, function(r){ q <- unlist(strsplit(x = r, "_")); q[length(q)] })))
        compartments <- paste(unique(c(paste(from.m, collapse='.'), paste(to.m, collapse='.'))), collapse = ">")
        c(genetic.rule = reaction[5], group = paste(paste(str_replace_all(reaction[5],pattern = " ",replacement = "-" )), compartments, sep="-"))
})

# Aggregation of data by compartments structure
aggregate.by.factor <- function(lands, group.vector=compartments.factor$group){
    df <- cbind(lands,group=group.vector)
    df.by.group <- df %>% 
        group_by(group) %>%
        summarise_all("mean")
    res <- as.matrix(df.by.group[,-1])
    rownames(res) <- as.character(df.by.group$group)
    list(activities = res, grouping = df.by.group$group)
}
