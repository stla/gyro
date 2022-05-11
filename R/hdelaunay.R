#' @useDynLib gyro, .registration=TRUE
#' @importFrom Rcpp evalCpp
NULL

#' @title Hyperbolic Delaunay triangulation
#' @description Computes the hyperbolic Delaunay triangulation of a set of
#'   points in the Poincaré disk.
#'
#' @encoding UTF-8
#'
#' @param points points in the unit disk given as a numeric matrix with
#'   two columns
#' @param isolations Boolean, whether to identify isolated vertices and edges
#' @param exact Boolean, whether to perform exact calculations; slower but
#'   more accurate
#'
#' @return A list with three fields \code{vertices}, \code{edges} and
#'   \code{faces}, and two additional fields \code{mvertices} and
#'   \code{medges} if \code{isolations=TRUE}, giving the non-isolated
#'   vertices and edges. The input \code{points} matrix and the output
#'   \code{vertices} matrix are the same up to the order of the rows.
#' @export
#'
#' @examples
#' library(gyro)
#' library(uniformly)
#' set.seed(666)
#' points <- runif_in_sphere(10L, d = 2)
#' hdelaunay(points)
hdelaunay <- function(points, isolations = FALSE, exact = FALSE){
  stopifnot(is.matrix(points))
  stopifnot(ncol(points) == 2L)
  stopifnot(nrow(points) >= 3L)
  stopifnot(is.numeric(points))
  if(anyNA(points)){
    stop("Missing values in `points`.")
  }
  sqnorms <- apply(points, 1L, dotprod)
  if(any(sqnorms >= 1)){
    stop("The points must be in the unit disk.")
  }
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
  class(hdel) <- "hdelaunay"
  hdel
}
