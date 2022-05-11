#' @useDynLib gyro, .registration=TRUE
#' @importFrom Rcpp evalCpp
NULL

hdelaunay <- function(points, isolations = FALSE, exact = FALSE){
  stopifnot(is.matrix(points))
  stopifnot(ncol(points) == 2L)
  stopifnot(nrow(points) >= 3L)
  stopifnot(isBoolean(isolations))
  stopifnot(isBoolean(exact))
  if(exact){
    hdel <- hdelaunay_EK(t(points), isolations)
  }else{
    hdel <- hdelaunay_K(t(points), isolations)
  }
  hdel[["vertices"]] <- t(hdel[["vertices"]])
  hdel[["edges"]] <- t(hdel[["edges"]])
  hdel[["faces"]] <- t(hdel[["faces"]])
  if(isolations){
    mv <- apply(hdel[["edges"]], 1L, function(edge){
      (edge[1L] %in% hdel[["mvertices"]]) && (edge[2L] %in% hdel[["mvertices"]])
    })
    hdel[["medges"]] <- hdel[["edges"]][mv, ]
  }
  hdel
}
