#####################################
## EXTRACT THE METABOLITES USED    ##
## IN THE REACTION AS:             ##
##   - SUBSTRATES (side='left'),   ##
##   - PRODUCTS (side='right')     ##
##   - ALL REACTANTS (side='both') ##
##################################### 
#' Extract the metabolites from the reaction given as character.
#' 
#' @param reaction A character describing reaction (e.g 'M_1 + M_2 = M_3').
#' @param side Which metabolites should be extracted: all, substrates or products. 
#'   Corresponding values are: c("both",left","right").
#' @return A vector of characters from the given \code{side}.
#' @examples
#' metabolites_from_reaction("M_1 + M_2 = "M_3")
#' metabolites_from_reaction("M_1 + M_2 = "M_3", "left")
metabolites_from_reaction <- function(reaction, side='both'){
    if(side=='left')  which.metabolites <- 1
    if(side=='right') which.metabolites <- 2
    if(side=='both')  which.metabolites <- c(1,2)
    res <- sapply(which.metabolites, function(s){
        as.character(reaction) %>% strsplit(., "=") %>% 
        lapply(., "[[", s) %>% unlist() %>% trimws() %>%
        strsplit(., "\\+") %>% unlist() %>% trimws() %>%
        strsplit(., " ") %>% unlist() %>% trimws() %>% 
        .[!grepl(pattern = "^[0-9]+", x = .)] 
    })
    unlist(res)
}

############################################
## PLOT THE FULL METABOLIC NETWORK        ##
## AS A BIPARTITE GRAPH                   ##
## REQUIRES:                              ##
## - recon.matrix                         ##
##   (matrix representation,              ##
##    rows: reactions, cols: metabolites) ##
############################################
run_me = F
if(run_me){
    no.reac  <- nrow(recon.matrix)

    out.m    <- colSums(recon.matrix > 0)
    in.m     <- colSums(recon.matrix < 0 )
    sel.cols <- (out.m > 1 | in.m > 1)

    sel.cols.num <- which((out.m > 1 | in.m > 1))

    y.r      <- nrow(recon.matrix[1:no.reac, ])
    y.m      <- ncol(recon.matrix[ , ])
    my.scale <- y.r/y.m

    pdf(file = "recon_viz_PC1.pdf", width = 10, height=300)
    psz(w = 10, h = 300)

    plot(x = rep(0, y.r), y=1:y.r, xlim=c(0,1), pch = 20)
    points(x = rep(1, y.m), y=(1:y.m)*my.scale, pch = 20)

    message("Drawing *in* reactinos.")

    idxs <- 1:nrow(recon.matrix[1:no.reac, ])
    #idxs <- which(rownames(recon.matrix) %in% rownames(recon[qq,]))
    r.out <- sapply(idxs, function(i){ 
        sapply(sel.cols.num, function(j){
            if(recon.matrix[i,j] == 1){
                c(i,j)
            }
        })
    })
    left <- (unlist(r.out))[c(1:length(unlist(r.out))) %% 2 == 1]
    right <- my.scale * (unlist(r.out))[c(1:length(unlist(r.out))) %% 2 == 0]
    arrows(x0 = rep(0, length(left)), y0 = left, x1 = rep(1, length(right)), y1=right, lty=2, angle = 20, lwd = 0.5)

    message("Drawing *out* reactions.")

    r.in <- sapply(idxs, function(i){ 
        sapply(sel.cols.num, function(j){
            if(recon.matrix[i,j] == -1){
                c(i,j)
            }
        })
    })
    left <- (unlist(r.in))[c(1:length(unlist(r.in))) %% 2 == 1]
    right <- my.scale * (unlist(r.in))[c(1:length(unlist(r.in))) %% 2 == 0]
    arrows(x1 = rep(0, length(left)), y1 = left, x0 = rep(1, length(right)), y0=right, lty=2, angle = 20, lwd = 0.5)
    dev.off()
}

#############
## The code below can help in an investigation 
## of metabolic composition of RECON, e.g. what is the composition of reactions
## how often a metabolite appears in specific compartment and so on
#############
out.m    <- colSums(recon.matrix > 0)
in.m     <- colSums(recon.matrix < 0)

sapply(
    as.character(enzymatic.reactions$reaction), 
    function(reaction) metabolites_from_reaction(reaction, "left")
    ) %>% unlist() %>% table() -> table_unique_enzymatic_substrates
sapply(
    as.character(enzymatic.reactions$reaction), 
    function(reaction) metabolites_from_reaction(reaction, "right")
    ) %>% unlist() %>% table() -> table_unique_enzymatic_products
sapply(
    as.character(enzymatic.reactions$reaction), 
    function(reaction) metabolites_from_reaction(reaction, "both")
    ) %>% unlist() %>% table() -> table_unique_enzymatic_metabolites


q  <- cbind(in.m, out.m)[order(in.m + out.m, decreasing = T), ]
qq <- cbind(q, 
          no.enz.reac = sapply(rownames(q), function(reac) table_unique_enzymatic_metabolites[reac]), 
          consumed    = sapply(rownames(q), function(reac) table_unique_enzymatic_substrates[reac]), 
          produced    = sapply(rownames(q), function(reac) table_unique_enzymatic_products[reac]))

qq[grepl(pattern = "^M_h2o_.*", x = rownames(qq)),]
#cbind(out.m, in.m)[order(in.m, decreasing = T), ][1:30,]
                        
