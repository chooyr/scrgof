#a to b distance matrix
eucDist <- function(xy1, xy2){
    i <- sort(rep(1:nrow(xy2), nrow(xy1)))
    dvec <- sqrt((xy1[, 1] - xy2[i, 1])^2 + (xy1[, 2] - xy2[i, 2])^2)
    matrix(dvec, nrow=nrow(xy1), ncol=nrow(xy2), byrow=F)
}