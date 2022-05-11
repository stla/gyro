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
#'   \code{triangles}, and two additional fields \code{mvertices} and
#'   \code{medges} if \code{isolations=TRUE}, giving the non-isolated
#'   vertices and edges ("m" for "multivalent"). The input \code{points}
#'   matrix and the output \code{vertices} matrix are the same up to the
#'   order of the rows.
#' @export
#'
#' @seealso \code{\link{plotHdelaunay}}
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
  hdel[["triangles"]] <- t(hdel[["faces"]])
  if(isolations){
    mv <- apply(hdel[["edges"]], 1L, function(edge){
      (edge[1L] %in% hdel[["mvertices"]]) && (edge[2L] %in% hdel[["mvertices"]])
    })
    hdel[["medges"]] <- hdel[["edges"]][mv, ]
  }
  class(hdel) <- "hdelaunay"
  hdel
}

#' @title Plot hyperbolic Delaunay triangulation
#' @description Plot a hyperbolic Delaunay triangulation obtained
#'   with \code{\link{hdelaunay}}.
#'
#' @param hdel an output of \code{\link{hdelaunay}}
#' @param remove \code{NULL}, \code{"ivertices"}, \code{"iedges"} or \code{c("ivertices", "iedges")}
#' @param vertices Boolean
#' @param edges Boolean
#' @param circle Boolean,
#' @param color either
#' @param hue,luminosity passed to \code{\link[randomcoloR]{randomColor}}
#'
#' @return No returned value.
#' @export
#'
#' @importFrom plotrix draw.circle
#' @importFrom graphics par polypath lines points
#' @importFrom randomcoloR randomColor distinctColorPalette
#'
#' @examples
#' library(gyro)
#' library(uniformly)
#' set.seed(666)
#' points <- runif_in_sphere(25L, d = 2)
#' hdel <- hdelaunay(points)
#' plotHdelaunay(hdel)
plotHdelaunay <- function(
  hdel, remove = NULL, vertices = TRUE, edges = TRUE, circle = TRUE,
  color = "random", hue = "random", luminosity = "random"
){
  opar <- par(mar = c(0, 0, 0, 0))
  plot(
    NULL, type = "n", asp = 1, xlim = c(-1, 1), ylim = c(-1, 1),
    xlab = NA, ylab = NA, axes = FALSE
  )
  draw.circle(0, 0, radius = 1, border = "black")
  pts <- hdel[["vertices"]]
  hedges <- hdel[["edges"]]
  triangles <- hdel[["triangles"]]
  colors <- randomColor(nrow(triangles), hue = hue, luminosity = luminosity)
  for(i in 1L:nrow(triangles)){
    trgl <- triangles[i, ]
    hpolypath <- rbind(
      Mgyrosegment(pts[trgl[1L], ], pts[trgl[3L], ], s = 1, n = 60)[-1L, ],
      Mgyrosegment(pts[trgl[3L], ], pts[trgl[2L], ], s = 1, n = 60)[-1L, ],
      Mgyrosegment(pts[trgl[2L], ], pts[trgl[1L], ], s = 1, n = 60)[-1L, ]
    )
    polypath(hpolypath, border = NA, col = colors[i])
  }
  for(i in 1L:nrow(hedges)){
    hedge <- hedges[i, ]
    hsegment <- Mgyrosegment(pts[hedge[1L], ], pts[hedge[2L], ], s = 1, n = 60)
    lines(hsegment, lty = "solid", col = "black", lwd = 1.5)
  }
  points(pts, pch = 19, cex = 0.9)
  par(opar)
  invisible(NULL)
}
