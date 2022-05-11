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
#' @return A list with four fields \code{vertices}, \code{edges},
#'   \code{triangles} and \code{ntriangles}, and two additional fields
#'   \code{mvertices} and \code{medges} if \code{isolations=TRUE}, giving
#'   the non-isolated vertices and edges ("m" for "multivalent"). The input
#'   \code{points} matrix and the output \code{vertices} matrix are the same
#'   up to the order of the rows.
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
  hdel[["ntriangles"]] <- ncol(hdel[["faces"]])
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
#' @param remove what you want to remove, \code{NULL} for nothing,
#'   \code{"ivertices"} for the isolated vertices, \code{"iedges"}
#'   for the isolated edges or \code{c("ivertices", "iedges")} for
#'   the isolated vertices and the isolated edges; if not \code{NULL},
#'   this assumes you ran \code{\link{hdelaunay}} with the
#'   option \code{isolations=TRUE}
#' @param vertices Boolean, whether to plot the vertices
#' @param edges Boolean, whether to plot the edges
#' @param circle Boolean, whether to plot the unit circle
#' @param color this argument controls the colors of the triangles; it can be
#'   \code{NA} for no color, \code{"random"} for random colors generated
#'   with \code{\link[randomcoloR]{randomColor}}, \code{"distinct"} for
#'   distinct colors generated with
#'   \code{\link[randomcoloR]{distinctColorPalette}}, a single color or
#'   a vector of colors
#' @param hue,luminosity passed to \code{\link[randomcoloR]{randomColor}}
#'   if \code{color="random"}
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
#' points <- runif_in_sphere(35L, d = 2)
#' hdel <- hdelaunay(points)
#' plotHdelaunay(hdel)
plotHdelaunay <- function(
  hdel, remove = NULL, vertices = TRUE, edges = TRUE, circle = TRUE,
  color = "distinct", hue = "random", luminosity = "random"
){
  if(!inherits(hdel, "hdelaunay")){
    stop("The `hdel` argument must be an output of the `hdelaunay` function.")
  }
  if(!is.null(remove)){
    isolations <- "medges" %in% names(hdel)
    if(!isolations){
      stop(
        "In order to use the `remove` argument you have to run the ",
        "`hdelaunay` function with `isolations=TRUE`."
      )
    }
  }
  opar <- par(mar = c(0, 0, 0, 0))
  plot(
    NULL, type = "n", asp = 1, xlim = c(-1, 1), ylim = c(-1, 1),
    xlab = NA, ylab = NA, axes = FALSE
  )
  if(circle){
    draw.circle(0, 0, radius = 1, border = "black")
  }
  pts <- hdel[["vertices"]]
  if(length(color) > 1L || !is.na(color)){
    triangles <- hdel[["triangles"]]
    ntriangles <- nrow(triangles)
    if(length(color) > 1L){
      colors <- color
    }else{
      if(color == "random"){
        colors <- randomColor(ntriangles, hue = hue, luminosity = luminosity)
      }else if(color == "distinct"){
        colors <- distinctColorPalette(ntriangles)
      }else{
        colors <- rep(color, ntriangles)
      }
    }
    for(i in 1L:ntriangles){
      trgl <- triangles[i, ]
      hpolypath <- rbind(
        Mgyrosegment(pts[trgl[1L], ], pts[trgl[3L], ], s = 1, n = 50)[-1L, ],
        Mgyrosegment(pts[trgl[3L], ], pts[trgl[2L], ], s = 1, n = 50)[-1L, ],
        Mgyrosegment(pts[trgl[2L], ], pts[trgl[1L], ], s = 1, n = 50)[-1L, ]
      )
      polypath(hpolypath, border = NA, col = colors[i])
    }
  }
  if(edges){
    if("iedges" %in% remove){
      hedges <- hdel[["medges"]]
    }else{
      hedges <- hdel[["edges"]]
    }
    for(i in 1L:nrow(hedges)){
      hedge <- hedges[i, ]
      hseg <- Mgyrosegment(pts[hedge[1L], ], pts[hedge[2L], ], s = 1, n = 50)
      lines(hseg, lty = "solid", col = "black", lwd = 1.5)
    }
  }
  if(vertices){
    if("ivertices" %in% remove){
      pts <- hdel[["mvertices"]]
    }
    points(pts, pch = 19, cex = 0.9)
  }
  par(opar)
  invisible(NULL)
}
