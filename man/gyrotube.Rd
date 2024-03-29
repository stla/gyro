% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gyro.R
\encoding{UTF-8}
\name{gyrotube}
\alias{gyrotube}
\title{Gyrotube (tubular gyrosegment)}
\usage{
gyrotube(A, B, s = 1, model = "U", n = 100, radius, sides = 90, caps = FALSE)
}
\arguments{
\item{A, B}{distinct 3D points}

\item{s}{positive number, the radius of the Poincaré ball if
\code{model="M"}, otherwise, if \code{model="U"}, this number
defines the hyperbolic curvature (higher value, less curved)}

\item{model}{the hyperbolic model, either \code{"M"} (Möbius model, i.e.
Poincaré model) or \code{"U"} (Ungar model, i.e. hyperboloid model)}

\item{n}{number of points forming the gyrosegment}

\item{radius}{radius of the tube around the gyrosegment}

\item{sides}{number of sides in the polygon cross section}

\item{caps}{Boolean, whether to put caps on the ends of the tube}
}
\value{
A \code{\link[rgl]{mesh3d}} object.
}
\description{
Tubular gyrosegment joining two given 3D points.
}
\examples{
library(gyro)
library(rgl)
A <- c(1, 2, 0); B <- c(1, 1, 0)
tube <- gyrotube(A, B, s = 0.2, radius = 0.02)
shade3d(tube, color = "orangered")

# a 3D hyperbolic triangle ####
library(rgl)
A <- c(1, 0, 0); B <- c(0, 1, 0); C <- c(0, 0, 1)
s <- 0.3
r <- 0.03
AB <- gyrotube(A, B, s, radius = r)
AC <- gyrotube(A, C, s, radius = r)
BC <- gyrotube(B, C, s, radius = r)
view3d(30, 30, zoom = 0.75)
shade3d(AB, color = "gold")
shade3d(AC, color = "gold")
shade3d(BC, color = "gold")
spheres3d(rbind(A, B, C), radius = 0.04, color = "gold")
}
