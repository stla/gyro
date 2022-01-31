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

#' @title Gyrosegment
#' @description Gyrosegment joining two given points.
#'
#' @param A,B two points (of the same dimension)
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
  stopifnot(isPositiveNumber(s))
  stopifnot(isPoint(A))
  stopifnot(isPoint(B))
  stopifnot(length(A) == length(B))
  stopifnot(isPositiveInteger(n))
  d <- length(A)
  t(vapply(seq(0, 1, length.out = n), function(t){
    gyroABt(A, B, t, s)
  }, numeric(d)))
}


