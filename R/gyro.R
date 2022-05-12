gammaF <- function(A, s){
  1 / sqrt(1 - dotprod(A)/(s*s))
}

betaF <- function(A, s) 1 / sqrt(1 + dotprod(A)/(s*s))

Ugyroadd <- function(A, B, s){
  betaA <- betaF(A, s)
  betaB <- betaF(B, s)
  (1 + betaA/(1+betaA) * dotprod(A, B)/(s*s) + (1-betaB)/betaB) * A + B
}

Mgyroadd <- function(X, Y, s){
  s2 <- s * s
  x <- dotprod(X) / s2
  y <- dotprod(Y) / s2
  xy <- 2 * dotprod(X, Y) / s2
  ((1 + xy + y) * X + (1 - x) * Y) / (1 + xy + x*y)
}

Ugyroscalar <- function(r, A, s){
  h <- sqrt(dotprod(A)) / s
  sinh(r*asinh(h)) * A / h
}

Mgyroscalar <- function(r, X, s){
  Xnorm <- sqrt(dotprod(X))
  s / Xnorm * tanh(r * atanh(Xnorm / s)) * X
}

UgyroABt <- function(A, B, t, s){
  Ugyroadd(A, Ugyroscalar(t, Ugyroadd(-A, B, s), s), s)
}

MgyroABt <- function(A, B, t, s){
  Mgyroadd(A, Mgyroscalar(t, Mgyroadd(-A, B, s), s), s)
}

#' @title Point on a gyroline
#' @description Point of coordinate \code{t} on the gyroline passing through
#'   two given points \code{A} and \code{B}. This is \code{A} for \code{t=0}
#'   and this is \code{B} for \code{t=1}. For \code{t=1/2} this is the
#'   gyromidpoint of the gyrosegment joining \code{A} and \code{B}.
#'
#' @encoding UTF-8
#'
#' @param A,B two distinct points
#' @param t a number
#' @param s positive number, the radius of the Poincaré ball if
#'   \code{model="M"}, otherwise, if \code{model="U"}, this number
#'   defines the hyperbolic curvature
#' @param model the hyperbolic model, either \code{"M"} (Möbius model, i.e.
#'   Poincaré model) or \code{"U"} (Ungar model, i.e. hyperboloid model)
#'
#' @return A point.
#' @export
gyroABt <- function(A, B, t, s = 1, model = "U"){
  model <- match.arg(model, c("M", "U"))
  stopifnot(isPositiveNumber(s))
  stopifnot(isPoint(A))
  stopifnot(isPoint(B))
  stopifnot(isNumber(t))
  stopifnot(areDistinct(A, B))
  if(model == "M"){
    if(dotprod(A) >= s || dotprod(B) >= s){
      stop(
        "In the M\u00f6bius gyrovector space, points must be ",
        "strictly inside the centered ball of radius `s`.",
        call. = TRUE
      )
    }
    MgyroABt(A, B, t, s)
  }else{
    UgyroABt(A, B, t, s)
  }
}

PhiUE <- function(A, s){
  gammaF(A, s) * A
}

PhiEU <- function(A, s){
  betaF(A, s) * A
}

PhiEM <- function(A, s){
  x <- dotprod(A) / s / s
  2 * A / (1 + x)
}

#' @title Isomorphism from Möbius gyrovector space to Ungar gyrovector space
#' @description Isomorphism from the Möbius gyrovector space to
#'   the Ungar gyrovector space.
#'
#' @encoding UTF-8
#'
#' @param A a point whose norm is lower than \code{s}
#' @param s positive number, the radius of the Poincaré ball
#'
#' @return The point of the Ungar gyrovector space corresponding to \code{A}
#'   by isomorphism.
#' @export
PhiUM <- function(A, s = 1){
  stopifnot(isPoint(A))
  stopifnot(s > 0)
  if(dotprod(A) >= s){
    stop(
      "In the M\u00f6bius gyrovector space, points must be ",
      "strictly inside the centered ball of radius `s`.",
      call. = TRUE
    )
  }
  PhiUE(PhiEM(A, s), s)
}

PhiME <- function(A, s){
  gamm <- gammaF(A, s)
  gamm * A / (1 + gamm)
}

#' @title Isomorphism from Ungar gyrovector space to Möbius gyrovector space
#' @description Isomorphism from the Ungar gyrovector space to
#'   the Möbius gyrovector space.
#'
#' @encoding UTF-8
#'
#' @param A a point in the Ungar vector space with curvature \code{s}
#' @param s a positive number, the hyperbolic curvature of the Ungar
#'   vector space
#'
#' @return The point of the Poincaré ball of radius \code{s} corresponding
#'   to \code{A} by isomorphism.
#' @export
PhiMU <- function(A, s = 1){
  stopifnot(isPoint(A))
  stopifnot(s > 0)
  PhiME(PhiEU(A, s), s)
}

Egyromidpoint <- function(A, B, s){
  gA <- gammaF(A, s); gB <- gammaF(B, s)
  (gA*A + gB*B) / (gA + gB)
}

Ugyromidpoint <- function(A, B, s){
  PhiUE(Egyromidpoint(PhiEU(A, s), PhiEU(B, s), s), s)
}

Mgyromidpoint <- function(A, B, s){
  PhiME(Egyromidpoint(PhiEM(A, s), PhiEM(B, s), s), s)
}

# Ugyromidpoint <- function(A, B, s){
#   UgyroABt(A, B, 0.5, s)
# }
#
# Mgyromidpoint <- function(A, B, s){
#   MgyroABt(A, B, 0.5, s)
# }

Ugyrosegment <- function(A, B, s, n){
  stopifnot(isPositiveNumber(s))
  stopifnot(isPositiveInteger(n))
  stopifnot(areDistinct(A, B))
  t(vapply(seq(0, 1, length.out = n), function(t){
    UgyroABt(A, B, t, s)
  }, numeric(length(A))))
}

Mgyrosegment <- function(A, B, s, n){
  stopifnot(isPositiveNumber(s))
  stopifnot(isPositiveInteger(n))
  stopifnot(areDistinct(A, B))
  t(vapply(seq(0, 1, length.out = n), function(t){
    MgyroABt(A, B, t, s)
  }, numeric(length(A))))
}

#' @title Gyrosegment
#' @description Gyrosegment joining two given points.
#'
#' @encoding UTF-8
#'
#' @param A,B two distinct points (of the same dimension)
#' @param s positive number, the radius of the Poincaré ball if
#'   \code{model="M"}, otherwise, if \code{model="U"}, this number
#'   defines the hyperbolic curvature
#' @param model the hyperbolic model, either \code{"M"} (Möbius model, i.e.
#'   Poincaré model) or \code{"U"} (Ungar model, i.e. hyperboloid model)
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
gyrosegment <- function(A, B, s = 1, model = "U", n = 100){
  model <- match.arg(model, c("M", "U"))
  stopifnot(isPoint(A))
  stopifnot(isPoint(B))
  stopifnot(length(A) == length(B))
  if(model == "M"){
    if(dotprod(A) >= s || dotprod(B) >= s){
      stop(
        "In the M\u00f6bius gyrovector space, points must be ",
        "strictly inside the centered ball of radius `s`.",
        call. = TRUE
      )
    }
    Mgyrosegment(A, B, s, n)
  }else{
    Ugyrosegment(A, B, s, n)
  }
}

#' @title Gyrotube (tubular gyrosegment)
#' @description Tubular gyrosegment joining two given 3D points.
#'
#' @encoding UTF-8
#'
#' @param A,B distinct 3D points
#' @param s positive number, the radius of the Poincaré ball if
#'   \code{model="M"}, otherwise, if \code{model="U"}, this number
#'   defines the hyperbolic curvature (higher value, less curved)
#' @param model the hyperbolic model, either \code{"M"} (Möbius model, i.e.
#'   Poincaré model) or \code{"U"} (Ungar model, i.e. hyperboloid model)
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
gyrotube <- function(
    A, B, s = 1, model = "U", n = 100, radius, sides = 90, caps = FALSE
){
  model <- match.arg(model, c("M", "U"))
  stopifnot(isPositiveNumber(s))
  stopifnot(is3dPoint(A))
  stopifnot(is3dPoint(B))
  stopifnot(isPositiveInteger(n))
  stopifnot(isPositiveInteger(sides))
  stopifnot(isBoolean(caps))
  if(model == "M"){
    if(dotprod(A) >= s || dotprod(B) >= s){
      stop(
        "In the M\u00f6bius gyrovector space, points must be ",
        "strictly inside the centered ball of radius `s`.",
        call. = TRUE
      )
    }
    points <- Mgyrosegment(A, B, s, n)
  }else{
    points <- Ugyrosegment(A, B, s, n)
  }
  closed <- ifelse(caps, -2, 0)
  cylinder3d(points, radius = radius, sides = sides, closed = closed)
}

Ugyrosubdiv <- function(A1, A2, A3, s){
  M12 <- Ugyromidpoint(A1, A2, s)
  M13 <- Ugyromidpoint(A1, A3, s)
  M23 <- Ugyromidpoint(A2, A3, s)
  list(
    list(A1, M12, M13),
    list(A2, M23, M12),
    list(A3, M13, M23),
    list(M12, M13, M23)
  )
}

Mgyrosubdiv <- function(A1, A2, A3, s){
  M12 <- Mgyromidpoint(A1, A2, s)
  M13 <- Mgyromidpoint(A1, A3, s)
  M23 <- Mgyromidpoint(A2, A3, s)
  list(
    list(A1, M12, M13),
    list(A2, M23, M12),
    list(A3, M13, M23),
    list(M12, M13, M23)
  )
}

Egyrocentroid <- function(A, B, C, s){
  gA <- gammaF(A, s); gB <- gammaF(B, s); gC <- gammaF(C, s)
  (gA*A + gB*B + gC*C) / (gA + gB + gC)
}

Ugyrocentroid <- function(A, B, C, s){
  PhiUE(Egyrocentroid(PhiEU(A, s), PhiEU(B, s), PhiEU(C, s), s), s)
}

Mgyrocentroid <- function(A, B, C, s){
  s2 <- s * s
  gA2 <- 1 / (1 - dotprod(A)/s2)
  gB2 <- 1 / (1 - dotprod(B)/s2)
  gC2 <- 1 / (1 - dotprod(C)/s2)
  if(
    gA2 < 0 || gB2 < 0 || gC2 < 0 ||
    is.infinite(gA2) || is.infinite(gB2) || is.infinite(gC2)
  ){
    stop(
      "In the M\u00f6bius gyrovector space, points must be ",
      "strictly inside the centered ball of radius `s`.",
      call. = FALSE
    )
  }
  Mgyroscalar(0.5, (gA2*A + gB2*B + gC2*C) / (gA2 + gB2 + gC2 - 1.5), s)
}

#' @title Gyrotriangle in 3D space
#' @description 3D gyrotriangle as a mesh.
#'
#' @encoding UTF-8
#'
#' @param A,B,C three distinct 3D points
#' @param s positive number, the radius of the Poincaré ball if
#'   \code{model="M"}, otherwise, if \code{model="U"}, this number
#'   defines the hyperbolic curvature (the smaller, the more curved)
#' @param model the hyperbolic model, either \code{"M"} (Möbius model, i.e.
#'   Poincaré model) or \code{"U"} (Ungar model, i.e. hyperboloid model)
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
#' \donttest{open3d(windowRect = c(50, 50, 562, 562))
#' view3d(30, 30, zoom = 0.75)
#' shade3d(ABC, color = "navy", specular = "cyan")}
#'
#' # using a color palette ####
#' library(trekcolors)
#' ABC <- gyrotriangle(
#'   A, B, C, s = 0.5,
#'   palette = trek_pal("klingon"), bias = 1.5, interpolate = "spline"
#' )
#' \donttest{open3d(windowRect = c(50, 50, 562, 562))
#' view3d(zoom = 0.75)
#' shade3d(ABC)}
#'
#' # hyperbolic icosahedron ####
#' library(rgl)
#' library(Rvcg) # to get the edges with the `vcgGetEdge` function
#' icosahedron <- icosahedron3d() # mesh with 12 vertices, 20 triangles
#' vertices <- t(icosahedron$vb[-4, ])
#' triangles <- t(icosahedron$it)
#' edges <- as.matrix(vcgGetEdge(icosahedron)[, c("vert1", "vert2")])
#' s <- 0.3
#' \donttest{open3d(windowRect = c(50, 50, 562, 562))
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
#' spheres3d(vertices, radius = 0.05, color = "lemonchiffon")}
gyrotriangle <- function(
    A, B, C, s = 1, model = "U", iterations = 5,
    palette = NULL, bias = 1, interpolate = "linear", g = identity
){
  model <- match.arg(model, c("M", "U"))
  if(model == "M"){
    if(dotprod(A) >= s || dotprod(B) >= s || dotprod(C) >= s){
      stop(
        "In the M\u00f6bius gyrovector space, points must be ",
        "strictly inside the centered ball of radius `s`.",
        call. = TRUE
      )
    }
    subd <- Mgyrosubdiv(A, B, C, s)
    for(i in seq_len(iterations-1)){
      subd <- flatten(lapply(subd, function(triplet){
        Mgyrosubdiv(triplet[[1L]], triplet[[2L]], triplet[[3L]], s)
      }))
    }
  }else{
    subd <- Ugyrosubdiv(A, B, C, s)
    for(i in seq_len(iterations-1)){
      subd <- flatten(lapply(subd, function(triplet){
        Ugyrosubdiv(triplet[[1L]], triplet[[2L]], triplet[[3L]], s)
      }))
    }
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
    if(model == "M"){
      gyroG <- Mgyrocentroid(A, B, C, s)
      dists <- sqrt(apply(mesh$vb[-4L, ], 2L, function(v){
        dotprod(Mgyroadd(-gyroG, v, s))
      }))
    }else{
      gyroG <- Ugyrocentroid(A, B, C, s)
      dists <- sqrt(apply(mesh$vb[-4L, ], 2L, function(v){
        dotprod(Ugyroadd(-gyroG, v, s))
      }))
    }
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
#' @encoding UTF-8
#'
#' @param points matrix of 3D points, one point per row
#' @param s positive number, the radius of the Poincaré ball if
#'   \code{model="M"}, otherwise, if \code{model="U"}, this number
#'   defines the hyperbolic curvature (the smaller, the more curved)
#' @param model the hyperbolic model, either \code{"M"} (Möbius model, i.e.
#'   Poincaré model) or \code{"U"} (Ungar model, i.e. hyperboloid model)
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
#' @importFrom rgl shade3d spheres3d lines3d
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
#' \donttest{open3d(windowRect = c(50, 50, 562, 562))
#' view3d(zoom = 0.7)
#' plotGyrohull3d(points, s = 0.4)}
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
#' \donttest{library(rgl)
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
#' plotGyrohull3d(vertices, s)}
#'
#' # an example of color palette ####
#' library(trekcolors)
#' library(uniformly)
#' set.seed(666)
#' points <- runif_on_sphere(50, d = 3)
#' \donttest{open3d(windowRect = c(50, 50, 562, 562))
#' plotGyrohull3d(
#'   points, edgesColor = "brown",
#'   facesColor = trek_pal("lcars_series"), g = function(u) 1-u^2
#' )}
plotGyrohull3d <- function(
    points, s = 1, model = "U", iterations = 5, n = 100, edgesAsTubes = TRUE,
    verticesAsSpheres = edgesAsTubes, edgesColor = "yellow",
    spheresColor = edgesColor, tubesRadius = 0.03, spheresRadius = 0.05,
    facesColor = "navy", bias = 1, interpolate = "linear", g = identity
){
  model <- match.arg(model, c("M", "U"))
  if(model == "M"){
    snorms <- apply(points, 1L, dotprod)
    if(any(snorms >= s)){
      stop(
        "In the M\u00f6bius gyrovector space, points must be ",
        "strictly inside the centered ball of radius `s`.",
        call. = TRUE
      )
    }
  }
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
      s = s, model = model, iterations = iterations, palette = palette,
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
        edge[1L, ], edge[2L, ], s = s, model = model,
        n = n, radius = tubesRadius
      )
      shade3d(gtube, color = edgesColor)
    }
  }else{
    if(model == "M"){
      for(edge in Edges){
        gsegment <- Mgyrosegment(
          edge[1L, ], edge[2L, ], s = s, n = n
        )
        lines3d(gsegment, color = edgesColor, lwd = 2)
      }
    }else{
      for(edge in Edges){
        gsegment <- Ugyrosegment(
          edge[1L, ], edge[2L, ], s = s, n = n
        )
        lines3d(gsegment, color = edgesColor, lwd = 2)
      }
    }
  }
  if(verticesAsSpheres){
    spheres3d(Vertices, radius = spheresRadius, color = spheresColor)
  }
  invisible(NULL)
}

#' @title Examples of the 'gyro' package
#' @description Some examples of hyperbolic polyhedra realized with the 'gyro'
#'   package.
#'
#' @return No value. The function firstly copies the demo files in a
#'   temporary directory. If you use RStudio, the function opens these files.
#'   Otherwise it prints a message giving the instructions to access to these
#'   files.
#'
#' @note The \emph{BarthLike} file has this name because the figure it
#'   generates looks like the Barth sextic (drawing by Patrice Jeener):
#'
#' \if{html}{
#'
#'   \out{<div style="text-align: center">}\figure{SextiqueDeBarth.png}{options: style="max-width:75\%;"}\out{</div>}
#'
#' }
#' \if{latex}{
#'
#'   \out{\begin{center}}\figure{SextiqueDeBarth.png}\out{\end{center}}
#'
#' }
#' @export
#'
#' @importFrom rstudioapi isAvailable navigateToFile
#' @importFrom clipr clipr_available write_clip
gyrodemos <- function(){
  folder <- system.file("gyrodemos", package = "gyro")
  tmpdir <- file.path(tempdir(), "gyrodemos")
  dir.create(tmpdir)
  files <- list.files(folder)
  tmpfiles <- file.path(tmpdir, files)
  invisible(file.copy(file.path(folder, files), tmpfiles))
  if(isAvailable(version_needed = "0.99.719")){
    invisible(sapply(tmpfiles, navigateToFile))
    message(sprintf("Opened files: %s.", toString(files)))
  }else{
    line <- sprintf(
      'wd <- setwd("%s")\n',
      normalizePath(tmpdir, winslash = "/")
    )
    message(
      "Copy the following line to go to the demos folder:\n",
      line,
      "Type `setwd(wd)` to come back."
    )
    if(clipr_available()){
      write_clip(line)
      message("The line has been copied to the clipboard.")
    }
  }
}

#' @title Changes of sign
#' @description Sometimes, the coordinates of the vertices of a polyhedron are
#'   given with changes of sign (with a symbol \strong{+/-}). This function
#'   performs the changes of sign.
#'
#' @param M a numeric matrix of coordinates of some points (one point per row)
#' @param changes either the indices of the columns of \code{M} where the
#'   changes of sign must be done, or \code{"all"} to select all the indices
#'
#' @return A numeric matrix, \code{M} transformed by the changes of sign.
#' @export
#'
#' @importFrom purrr imap
#'
#' @examples library(gyro)
#' library(rgl)
#' ## ~~ rhombicosidodecahedron ~~##
#' phi <- (1 + sqrt(5)) / 2
#' vs1 <- rbind(
#'   c(1, 1, phi^3),
#'   c(phi^2, phi, 2 * phi),
#'   c(2 + phi, 0, phi^2)
#' )
#' vs2 <- rbind(vs1, vs1[, c(2, 3, 1)], vs1[, c(3, 1, 2)]) # even permutations
#' vs <- changesOfSign(vs2)
#' \donttest{open3d(windowRect = c(50, 50, 562, 562), zoom = 0.65)
#' plotGyrohull3d(vs)}
changesOfSign <- function(M, changes = "all"){
  if(!is.matrix(M)) M <- rbind(M)
  if(identical(changes, "all")) changes <- 1L:ncol(M)
  `colnames<-`(as.matrix(do.call(rbind, apply(M, 1L, function(row){
    expand.grid(
      imap(row, ~ (if(.x == 0 || !.y %in% changes) .x else c(-.x, .x)))
    )
  }))), NULL)
}
