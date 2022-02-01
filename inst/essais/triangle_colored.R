library(gyro)
library(rgl)

s <- 0.5

A <- c(1, 0, 0); B <- c(0, 1, 0); C <- c(0, 0, 1)
ABC <- gyrotriangle(A, B, C, s)

gyroG <- gyro:::gyrocentroid(A, B, C, s)

dists <- apply(ABC$vb[-4L, ], 2L, function(v){
  c(crossprod(gyro:::gyroadd(-gyroG, v, s)))
})

ndists <- (dists - min(dists))/diff(range(dists))

library(trekcolors)

fpalette <- colorRamp(trek_pal("klingon"), bias = 1.5, interpolate = "spline")
RGB <- fpalette(ndists)
colors <- rgb(RGB[, 1L], RGB[, 2L], RGB[, 3L], maxColorValue = 255)


ABC$material = list(color = colors)

shade3d(ABC)
