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
  stopifnot(isPositiveNumber(s))
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
  stopifnot(isPositiveNumber(s))
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

#' @title Gyromidpoint
#' @description The gyromidpoint of a \code{\link{gyrosegment}}.
#'
#' @encoding UTF-8
#'
#' @param A,B two distinct points (of the same dimension)
#' @param s positive number, the radius of the Poincaré ball if
#'   \code{model="M"}, otherwise, if \code{model="U"}, this number
#'   defines the hyperbolic curvature
#' @param model the hyperbolic model, either \code{"M"} (Möbius model, i.e.
#'   Poincaré model) or \code{"U"} (Ungar model, i.e. hyperboloid model)
#'
#' @return A point, the gyromidpoint of a the \code{\link{gyrosegment}}
#'   joining \code{A} and \code{B}.
#' @export
#' @note This is the same as \code{gyroABt(A, B, 1/2, s)} but the
#'   calculation is more efficient.
gyromidpoint <- function(A, B, s = 1, model = "U"){
  model <- match.arg(model, c("M", "U"))
  stopifnot(isPoint(A))
  stopifnot(isPoint(B))
  stopifnot(length(A) == length(B))
  stopifnot(isPositiveNumber(s))
  stopifnot(areDistinct(A, B))
  if(model == "M"){
    if(dotprod(A) >= s || dotprod(B) >= s){
      stop(
        "In the M\u00f6bius gyrovector space, points must be ",
        "strictly inside the centered ball of radius `s`.",
        call. = TRUE
      )
    }
    Mgyromidpoint(A, B, s)
  }else{
    Ugyromidpoint(A, B, s)
  }
}

Ugyrosegment <- function(A, B, s, n){
  t(vapply(seq(0, 1, length.out = n), function(t){
    UgyroABt(A, B, t, s)
  }, numeric(length(A))))
}

Mgyrosegment <- function(A, B, s, n){
  t(Mgyrosegment_cpp(A, B, s, n))
  # t(vapply(seq(0, 1, length.out = n), function(t){
  #   MgyroABt(A, B, t, s)
  # }, numeric(length(A))))
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
#' @return A numeric matrix with \code{n} rows. Each row is a point of the
#'   gyrosegment from \code{A} (the first row) to \code{B} (the last row).
#' @export
#'
#' @note The gyrosegment is obtained from \code{\link{gyroABt}} by varying
#' \code{t} from \code{0} to \code{1}.
#'
#' @examples
#' library(gyro)
#' # a 2D example ####
#' A <- c(1, 2); B <- c(1, 1)
#' opar <- par(mfrow = c(1, 2), mar = c(2, 2, 2, 0.5))
#' plot(rbind(A, B), type = "p", pch = 19, xlab = NA, ylab = NA,
#'      xlim = c(0, 2), ylim = c(0, 2), main = "s = 0.2")
#' s <- 0.2
#' AB <- gyrosegment(A, B, s)
#' lines(AB, col = "blue", lwd = 2)
#' text(t(A), expression(italic(A)), pos = 2)
#' text(t(B), expression(italic(B)), pos = 3)
#' # this is an hyperbola whose asymptotes meet at the origin
#' # approximate asymptotes
#' lines(rbind(c(0, 0), gyroABt(A, B, t = -20, s)), lty = "dashed")
#' lines(rbind(c(0, 0), gyroABt(A, B, t = 20, s)), lty = "dashed")
#' # plot the gyromidoint
#' points(
#'  rbind(gyromidpoint(A, B, s)),
#'  type = "p", pch = 19, col = "red"
#' )
#' # another one, with a different `s`
#' plot(rbind(A, B), type = "p", pch = 19, xlab = NA, ylab = NA,
#'      xlim = c(0, 2), ylim = c(0, 2), main = "s = 0.1")
#' s <- 0.1
#' AB <- gyrosegment(A, B, s)
#' lines(AB, col = "blue", lwd = 2)
#' text(t(A), expression(italic(A)), pos = 2)
#' text(t(B), expression(italic(B)), pos = 3)
#' # approximate asymptotes
#' lines(rbind(c(0, 0), gyroABt(A, B, t = -20, s)), lty = "dashed")
#' lines(rbind(c(0, 0), gyroABt(A, B, t = 20, s)), lty = "dashed")
#' # plot the gyromidoint
#' points(
#'  rbind(gyromidpoint(A, B, s)),
#'  type = "p", pch = 19, col = "red"
#' )
#'
#' # a 3D hyperbolic triangle ####
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
  stopifnot(isPositiveNumber(s))
  stopifnot(isPositiveInteger(n))
  stopifnot(n >= 2)
  stopifnot(areDistinct(A, B))
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
