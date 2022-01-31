betaF <- function(A, s) 1 / sqrt(1 + dotprod(A)/(s*s))

gyroadd <- function(A, B, s){
  betaA <- betaF(A, s)
  betaB <- betaF(B, s)
  (1 + betaA/(1+betaA) * dotprod(A, B)/(s*s) + (1-betaB)/betaB) * A + B
}

gyroscalar <- function(r, A, s){
  h <- sqrt(dotprod(A)) / s
  sinh(r*asinh(h)) * A / h
}

gyroABt <- function(A, B, t, s){
  gyroadd(A, gyroscalar(t, gyroadd(-A, B, s), s), s)
}

gyromidpoint <- function(A, B, s){
  gyroABt(A, B, 0.5, s)
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

