gammaF <- function(A, s){
  1 / sqrt(1 - dotprod(A)/(s*s))
}

# gyromidpointE <- function(A, B, s=1){
#   gA <- gamm(A, s=s); gB <- gamm(B, s=s)
#   (gA*A + gB*B) / (gA+gB)
# }

gyrocentroidE <- function(A, B, C, s){
  gA <- gammaF(A, s); gB <- gammaF(B, s); gC <- gammaF(C, s)
  (gA*A + gB*B + gC*C) / (gA + gB + gC)
}

PhiEU <- function(A, s){
  gammaF(A, s) * A
}

betaF <- function(A, s) 1 / sqrt(1 + dotprod(A)/(s*s))

PhiUE <- function(A, s){
  betaF(A, s) * A
}

# gyromidpointU <- function(A, B, s=1){
#   PhiEU(gyromidpointE(PhiUE(A,s=s),PhiUE(B,s=s),s=s),s=s)
# }

gyrocentroid <- function(A, B, C, s){
  PhiEU(gyrocentroidE(PhiUE(A, s), PhiUE(B, s), PhiUE(C, s), s), s)
}

gyroadd <- function(A, B, s){
  betaA <- betaF(A, s)
  betaB <- betaF(B, s)
  (1 + betaA/(1+betaA) * dotprod(A, B)/(s*s) + (1-betaB)/betaB) * A + B
}

gyroscalar <- function(r, A, s){
  h <- sqrt(dotprod(A)) / s
  sinh(r*asinh(h)) * A / h
}

.gyroABt <- function(A, B, t, s){
  gyroadd(A, gyroscalar(t, gyroadd(-A, B, s), s), s)
}

#' @title Point on a gyroline
#' @description Point of coordinate \code{t} on the gyroline passing through
#'   two given points \code{A} and \code{B}. This is \code{A} for \code{t=0}
#'   and this is \code{B} for \code{t=1}. For \code{t=1/2} this is the
#'   gyromidpoint of the gyrosegment joining \code{A} and \code{B}.
#'
#' @param A,B two distinct points
#' @param t a number
#' @param s positive number, the parameter defining the hyperbolic curvature
#'
#' @return A point.
#' @export
gyroABt <- function(A, B, t, s){
  stopifnot(isPositiveNumber(s))
  stopifnot(isPoint(A))
  stopifnot(isPoint(B))
  stopifnot(isNumber(t))
  stopifnot(areDistinct(A, B))
  .gyroABt(A, B, t, s)
}

gyromidpoint <- function(A, B, s){
  .gyroABt(A, B, 0.5, s)
}

.gyrosegment <- function(A, B, s, n){
  stopifnot(isPositiveNumber(s))
  stopifnot(isPositiveInteger(n))
  stopifnot(areDistinct(A, B))
  t(vapply(seq(0, 1, length.out = n), function(t){
    gyroABt(A, B, t, s)
  }, numeric(length(A))))
}


#' @title Gyrosegment
#' @description Gyrosegment joining two given points.
#'
#' @param A,B two distinct points (of the same dimension)
#' @param s positive number, the curvature
#' @param n number of points forming the gyrosegment from \code{A} to \code{B}
#'
#' @return A numeric matrix with \code{n} rows. Each row is a point on the
#'   gyrosegment from \code{A} (the first row) to \code{B} (the last row).
#' @export
#'
#' @examples library(gyro)
#' # a 2D example ####
#' A <- c(1, 2); B <- c(1, 1)
#' plot(rbind(A, B), type = "p", pch = 19, xlab = NA, ylab = NA,
#'      xlim = c(0, 2), ylim = c(0, 2), asp = 1)
#' AB <- gyrosegment(A, B, s = 0.2)
#' lines(AB) # this is a piece of an hyperboloid
#' text(t(A), expression(italic(A)), pos = 1)
#' text(t(B), expression(italic(B)), pos = 3)
#'
#' # a 3D hyperbolic triangle
#' library(rgl)
#' A <- c(1, 0, 0); B <- c(0, 1, 0); C <- c(0, 0, 1)
#' s <- 0.3
#' AB <- gyrosegment(A, B, s)
#' AC <- gyrosegment(A, C, s)
#' BC <- gyrosegment(B, C, s)
#' view3d(30, 30, zoom = 0.75)
#' lines3d(AB, lwd = 3); lines3d(AC, lwd = 3); lines3d(BC, lwd = 3)
gyrosegment <- function(A, B, s = 1, n = 100){
  stopifnot(isPoint(A))
  stopifnot(isPoint(B))
  stopifnot(length(A) == length(B))
  .gyrosegment(A, B, s, n)
}

#' @title Gyrotube (tubular gyrosegment)
#' @description Tubular gyrosegment joining two given 3D points.
#'
#' @param A,B distinct 3D points
#' @param s positive number, the curvature (higher value, less curved)
#' @param n number of points forming the gyrosegment
#' @param radius radius of the tube around the gyrosegment
#' @param sides number of sides in the polygon cross section
#' @param caps Boolean, whether to put caps on the ends of the tube
#'
#' @return A \code{\link[rgl]{mesh3d}} object.
#' @export
#'
#' @importFrom rgl cylinder3d
#'
#' @examples library(gyro)
#' library(rgl)
#' A <- c(1, 2, 0); B <- c(1, 1, 0)
#' tube <- gyrotube(A, B, s = 0.2, radius = 0.02)
#' shade3d(tube, color = "orangered")
#'
#' # a 3D hyperbolic triangle ####
#' library(rgl)
#' A <- c(1, 0, 0); B <- c(0, 1, 0); C <- c(0, 0, 1)
#' s <- 0.3
#' r <- 0.03
#' AB <- gyrotube(A, B, s, radius = r)
#' AC <- gyrotube(A, C, s, radius = r)
#' BC <- gyrotube(B, C, s, radius = r)
#' view3d(30, 30, zoom = 0.75)
#' shade3d(AB, color = "gold")
#' shade3d(AC, color = "gold")
#' shade3d(BC, color = "gold")
#' spheres3d(rbind(A, B, C), radius = 0.04, color = "gold")
gyrotube <- function(A, B, s = 1, n = 100, radius, sides = 90, caps = FALSE){
  stopifnot(isPositiveNumber(s))
  stopifnot(is3dPoint(A))
  stopifnot(is3dPoint(B))
  stopifnot(isPositiveInteger(n))
  stopifnot(isPositiveInteger(sides))
  stopifnot(isBoolean(caps))
  points <- .gyrosegment(A, B, s, n)
  closed <- ifelse(caps, -2, 0)
  cylinder3d(points, radius = radius, sides = sides, closed = closed)
}

gyrosubdiv <- function(A1, A2, A3, s){
  M12 <- gyromidpoint(A1, A2, s)
  M13 <- gyromidpoint(A1, A3, s)
  M23 <- gyromidpoint(A2, A3, s)
  list(
    list(A1, M12, M13),
    list(A2, M23, M12),
    list(A3, M13, M23),
    list(M12, M13, M23)
  )
}

#' @title Gyrotriangle in 3D space
#' @description 3D gyrotriangle as a mesh.
#'
#' @param A,B,C three distinct 3D points
#' @param s positive number, the curvature (the smaller, the more curved)
#' @param iterations the gyrotriangle is constructed by iterated subdivisions,
#'   this argument is the number of iterations
#' @param palette a vector of colors to decorate the triangle, or \code{NULL}
#'   if you don't want to use a color palette
#' @param bias,interpolate if \code{palette} is not \code{NULL}, these
#'   arguments are passed to \code{\link[grDevices]{colorRamp}}
#' @param g a function from [0,1] to [0,1]; if \code{palette} is not
#'   \code{NULL}, this function is applied to the scalars defining the colors
#'   (the normalized gyrodistances to the gyrocentroid of the gyrotriangle)
#'
#' @return A \code{\link[rgl]{mesh3d}} object.
#' @export
#'
#' @importFrom purrr flatten
#' @importFrom rgl tmesh3d
#' @importFrom Rvcg vcgClean
#' @importFrom grDevices colorRamp rgb
#'
#' @examples library(gyro)
#' library(rgl)
#' A <- c(1, 0, 0); B <- c(0, 1, 0); C <- c(0, 0, 1)
#' ABC <- gyrotriangle(A, B, C, s = 0.3)
#' open3d(windowRect = c(50, 50, 562, 562))
#' view3d(30, 30, zoom = 0.75)
#' shade3d(ABC, color = "navy", specular = "cyan")
#'
#' # using a color palette ####
#' library(trekcolors)
#' ABC <- gyrotriangle(
#'   A, B, C, s = 0.5,
#'   palette = trek_pal("klingon"), bias = 1.5, interpolate = "spline"
#' )
#' open3d(windowRect = c(50, 50, 562, 562))
#' view3d(zoom = 0.75)
#' shade3d(ABC)
#'
#' # hyperbolic icosahedron ####
#' library(rgl)
#' library(Rvcg) # to get the edges with the `vcgGetEdge` function
#' icosahedron <- icosahedron3d() # mesh with 12 vertices, 20 triangles
#' vertices <- t(icosahedron$vb[-4, ])
#' triangles <- t(icosahedron$it)
#' edges <- as.matrix(vcgGetEdge(icosahedron)[, c("vert1", "vert2")])
#' s <- 0.3
#' open3d(windowRect = c(50, 50, 562, 562))
#' view3d(zoom = 0.75)
#' for(i in 1:nrow(triangles)){
#'   triangle <- triangles[i, ]
#'   A <- vertices[triangle[1], ]
#'   B <- vertices[triangle[2], ]
#'   C <- vertices[triangle[3], ]
#'   gtriangle <- gyrotriangle(A, B, C, s)
#'   shade3d(gtriangle, color = "midnightblue")
#' }
#' for(i in 1:nrow(edges)){
#'   edge <- edges[i, ]
#'   A <- vertices[edge[1], ]
#'   B <- vertices[edge[2], ]
#'   gtube <- gyrotube(A, B, s, radius = 0.03)
#'   shade3d(gtube, color = "lemonchiffon")
#' }
#' spheres3d(vertices, radius = 0.05, color = "lemonchiffon")
gyrotriangle <- function(
  A, B, C, s = 1, iterations = 5,
  palette = NULL, bias = 1, interpolate = "linear", g = identity
){
  subd <- gyrosubdiv(A, B, C, s)
  for(i in seq_len(iterations-1)){
    subd <- flatten(lapply(subd, function(triplet){
      gyrosubdiv(triplet[[1L]], triplet[[2L]], triplet[[3L]], s)
    }))
  }
  vertices <-
    do.call(cbind, lapply(subd, function(triplet) do.call(cbind, triplet)))
  indices <- matrix(1L:ncol(vertices), nrow = 3L)
  mesh0 <- tmesh3d(
    vertices = vertices,
    indices = indices
  )
  mesh <- vcgClean(mesh0, sel = c(0, 7), silent = TRUE)
  mesh[["remvert"]] <- NULL
  mesh[["remface"]] <- NULL
  if(!is.null(palette)){
    fpalette <- colorRamp(palette, bias = bias, interpolate = interpolate)
    gyroG <- gyrocentroid(A, B, C, s)
    dists <- sqrt(apply(mesh$vb[-4L, ], 2L, function(v){
      dotprod(gyroadd(-gyroG, v, s))
    }))
    dists <- (dists - min(dists))/diff(range(dists))
    RGB <- fpalette(g(dists))
    colors <- rgb(RGB[, 1L], RGB[, 2L], RGB[, 3L], maxColorValue = 255)
    mesh[["material"]] = list(color = colors)
  }
  mesh
}

#' @importFrom cxhull cxhull VerticesXYZ EdgesXYZ TrianglesXYZ
#' @noRd
.cxhull <- function(points){
  hull <- cxhull(points, triangulate = TRUE)
  Vertices <- VerticesXYZ(hull)
  Edges <- EdgesXYZ(hull)
  Edges <- Edges[Edges[, "border"] == 1, c("x", "y", "z")]
  nedges <- nrow(Edges) %/% 2L
  Edges <- lapply(split(as.data.frame(Edges), gl(nedges, 2L)), as.matrix)
  Triangles <- TrianglesXYZ(hull)
  ntriangles <- length(hull[["facets"]])
  Triangles <-
    lapply(split(as.data.frame(Triangles), gl(ntriangles, 3L)), as.matrix)
  return(list(vertices = Vertices, edges = Edges, triangles = Triangles))
}

#' @title Hyperbolic convex hull
#' @description Plot the hyperbolic convex hull of a set of 3D points.
#'
#' @param points matrix of 3D points, one point per row
#' @param s curvature parameter
#' @param iterations argument passed to \code{\link{gyrotriangle}}
#' @param n argument passed to \code{\link{gyrotube}} or
#'   \code{\link{gyrosegment}}, the number of points for each edge
#' @param edgesAsTubes Boolean, whether to represent tubular edges
#' @param verticesAsSpheres Boolean, whether to represent the vertices as
#'   spheres
#' @param edgesColor a color for the edges
#' @param spheresColor a color for the spheres, if
#'   \code{verticesAsSpheres = TRUE}
#' @param tubesRadius radius of the tubes, if \code{edgesAsTubes = TRUE}
#' @param spheresRadius radius of the spheres,
#'   if \code{verticesAsSpheres = TRUE}
#' @param facesColor this argument sets the color of the faces; it can be
#'   either a single color or a color palette, i.e. a vector of colors; if it
#'   is a color palette, it will be passed to the argument \code{palette} of
#'   \code{\link{gyrotriangle}}
#' @param bias,interpolate,g these arguments are passed to
#'   \code{\link{gyrotriangle}} in the case when \code{facesColor} is a color
#'   palette
#'
#' @return No value, called for plotting.
#' @export
#'
#' @importFrom rgl shade3d spheres3d
#' @importFrom Morpho mergeMeshes
#' @importFrom Rvcg vcgClean
#'
#' @examples library(gyro)
#' library(rgl)
#' # Triangular orthobicopula ####
#' points <- rbind(
#'   c(1, -1/sqrt(3), sqrt(8/3)),
#'   c(1, -1/sqrt(3), -sqrt(8/3)),
#'   c(-1, -1/sqrt(3), sqrt(8/3)),
#'   c(-1, -1/sqrt(3), -sqrt(8/3)),
#'   c(0, 2/sqrt(3), sqrt(8/3)),
#'   c(0, 2/sqrt(3), -sqrt(8/3)),
#'   c(1, sqrt(3), 0),
#'   c(1, -sqrt(3), 0),
#'   c(-1, sqrt(3), 0),
#'   c(-1, -sqrt(3), 0),
#'   c(2, 0, 0),
#'   c(-2, 0, 0)
#' )
#' open3d(windowRect = c(50, 50, 562, 562))
#' view3d(zoom = 0.7)
#' plotGyrohull3d(points, s = 0.4)
#'
#' # a non-convex polyhedron with triangular faces ####
#' vertices <- rbind(
#'   c(-2.1806973249, -2.1806973249, -2.1806973249),
#'   c(-3.5617820682, 0.00000000000, 0.00000000000),
#'   c(0.00000000000, -3.5617820682, 0.00000000000),
#'   c(0.00000000000, 0.00000000000, -3.5617820682),
#'   c(-2.1806973249, -2.1806973249, 2.18069732490),
#'   c(0.00000000000, 0.00000000000, 3.56178206820),
#'   c(-2.1806973249, 2.18069732490, -2.1806973249),
#'   c(0.00000000000, 3.56178206820, 0.00000000000),
#'   c(-2.1806973249, 2.18069732490, 2.18069732490),
#'   c(2.18069732490, -2.1806973249, -2.1806973249),
#'   c(3.56178206820, 0.00000000000, 0.00000000000),
#'   c(2.18069732490, -2.1806973249, 2.18069732490),
#'   c(2.18069732490, 2.18069732490, -2.1806973249),
#'   c(2.18069732490, 2.18069732490, 2.18069732490))
#' triangles <- 1 + rbind(
#'   c(3, 2, 0),
#'   c(0, 1, 3),
#'   c(2, 1, 0),
#'   c(4, 2, 5),
#'   c(5, 1, 4),
#'   c(4, 1, 2),
#'   c(6, 7, 3),
#'   c(3, 1, 6),
#'   c(6, 1, 7),
#'   c(5, 7, 8),
#'   c(8, 1, 5),
#'   c(7, 1, 8),
#'   c(9, 2, 3),
#'   c(3, 10, 9),
#'   c(9, 10, 2),
#'   c(5, 2, 11),
#'   c(11, 10, 5),
#'   c(2, 10, 11),
#'   c(3, 7, 12),
#'   c(12, 10, 3),
#'   c(7, 10, 12),
#'   c(13, 7, 5),
#'   c(5, 10, 13),
#'   c(13, 10, 7))
#' edges0 <- do.call(c, lapply(1:nrow(triangles), function(i){
#'   face <- triangles[i, ]
#'   list(
#'     sort(c(face[1], face[2])),
#'     sort(c(face[1], face[3])),
#'     sort(c(face[2], face[3]))
#'   )
#' }))
#' edges <- do.call(rbind, edges0)
#' edges <- edges[!duplicated(edges), ]
#' s <- 2
#' library(rgl)
#' open3d(windowRect = c(50, 50, 1074, 562))
#' mfrow3d(1, 2)
#' view3d(zoom = 0.65)
#' for(i in 1:nrow(triangles)){
#'   triangle <- triangles[i, ]
#'   A <- vertices[triangle[1], ]
#'   B <- vertices[triangle[2], ]
#'   C <- vertices[triangle[3], ]
#'   gtriangle <- gyrotriangle(A, B, C, s)
#'   shade3d(gtriangle, color = "violetred")
#' }
#' for(i in 1:nrow(edges)){
#'   edge <- edges[i, ]
#'   A <- vertices[edge[1], ]
#'   B <- vertices[edge[2], ]
#'   gtube <- gyrotube(A, B, s, radius = 0.06)
#'   shade3d(gtube, color = "darkviolet")
#' }
#' spheres3d(vertices, radius = 0.09, color = "deeppink")
#' # now plot the hyperbolic convex hull
#' next3d()
#' view3d(zoom = 0.65)
#' plotGyrohull3d(vertices, s)
#'
#' # an example of color palette ####
#' library(trekcolors)
#' library(uniformly)
#' set.seed(666)
#' points <- runif_on_sphere(50, d = 3)
#' open3d(windowRect = c(50, 50, 562, 562))
#' plotGyrohull3d(
#'   points, edgesColor = "brown",
#'   facesColor = trek_pal("lcars_series"), g = function(u) 1-u^2
#' )
plotGyrohull3d <- function(
  points, s = 1, iterations = 5, n = 100, edgesAsTubes = TRUE,
  verticesAsSpheres = edgesAsTubes, edgesColor = "yellow",
  spheresColor = edgesColor, tubesRadius = 0.03, spheresRadius = 0.05,
  facesColor = "navy", bias = 1, interpolate = "linear", g = identity
){
  stopifnot(isBoolean(edgesAsTubes))
  stopifnot(isBoolean(verticesAsSpheres))
  hull <- .cxhull(points)
  Triangles <- hull[["triangles"]]
  Edges <- hull[["edges"]]
  Vertices <- hull[["vertices"]]
  ntriangles <- length(Triangles)
  Gtriangles <- vector("list", ntriangles)
  palette <- if(length(facesColor) > 1L) facesColor
  for(i in 1L:ntriangles){
    triangle <- Triangles[[i]]
    Gtriangles[[i]] <- gyrotriangle(
      triangle[1L, ], triangle[2L, ], triangle[3L, ],
      s = s, iterations = iterations, palette = palette,
      bias = bias, interpolate = interpolate, g = g
    )
  }
  mesh <- vcgClean(mergeMeshes(Gtriangles), sel = 0, silent = TRUE)
  if(is.null(palette)){
    shade3d(mesh, color = facesColor)
  }else{
    shade3d(mesh)
  }
  if(edgesAsTubes){
    for(edge in Edges){
      gtube <- gyrotube(
        edge[1L, ], edge[2L, ], s = s, n = n, radius = tubesRadius
      )
      shade3d(gtube, color = edgesColor)
    }
  }else{
    for(edge in Edges){
      gsegment <- gyrosegment(
        edge[1L, ], edge[2L, ], s = s, n = n
      )
      lines3d(gsegment, color = edgesColor, lwd = 2)
    }
  }
  if(verticesAsSpheres){
    spheres3d(Vertices, radius = spheresRadius, color = spheresColor)
  }
  invisible(NULL)
}


