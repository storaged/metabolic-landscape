######################################################
## CALCULATION OF GRAPH ENTROPY                     ##
## x, y - logical vectors                           ##
######################################################

## Probability of a vetex given: 
## stoch.mat: a verted for a given stochastic matrix
## my.subset: a subset to which the vertex belongs
## inner_p: boolean if we calculate inner or outer probability
## row_n/col_n: the number of vertex row-/col-wise
p.vert <- function(stoch.mat, my.subset, inner_p = T, row.n = NULL, col.n = NULL){
  if(is.null(my.subset$cols)){print("NULL is coming")}
  if(!is.null(row.n)){
    if(is.null(my.subset$cols)){return(ifelse(inner_p, 0, 1))}
    q <- sum(abs(stoch.mat[row.n, my.subset$cols]))/sum(abs(stoch.mat[row.n, ]))
    return(ifelse(inner_p, q, 1-q))
  } else if(!is.null(col.n)) {
    if(is.null(my.subset$rows)){return(ifelse(inner_p, 0, 1))}
    q <- sum(abs(stoch.mat[my.subset$rows, col.n]))/sum(abs(stoch.mat[ ,col.n]))
    return(ifelse(inner_p, q, 1-q))
  } else {
    return(NULL)
  }
}

## Entropy of a vertex given:
## stoch.mat: a verted for a given stochastic matrix
## my.subset: a subset to which the vertex belongs
## row_n/col_n: the number of vertex row-/col-wise
entr.vert <- function(stoch.mat, my.subset, row.n = NULL, col.n = NULL){
  pi <- p.vert(stoch.mat, my.subset, inner_p = T, row.n, col.n)
  po <- p.vert(stoch.mat, my.subset, inner_p = F, row.n, col.n)
  qi <- ifelse(pi == 0, 0, -pi * log2(pi))
  qo <- ifelse(po == 0, 0, -po * log2(po))
  q <- qi + qo
  return(q)
}

## Entropy of a graph given:
## stoch.mat: a verted for a given stochastic matrix
## my.subset: a subset to which the vertex belongs
entr.graph <- function(stoch.mat, my.subset){
  e1 <- sum(sapply(my.subset$rows, function(x) entr.vert(stoch.mat, my.subset, row.n = x)))
  e2 <- ifelse(is.null(my.subset$cols), 0,
               sum(sapply(my.subset$cols, function(x) entr.vert(stoch.mat, my.subset, col.n = x))))
  return(e1+e2)
}