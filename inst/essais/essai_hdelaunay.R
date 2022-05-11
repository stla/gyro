library(gyro)
library(uniformly)
library(rsvg)

path <- "./inst/trash"

set.seed(314)

pts <- rbind(
  runif_in_annulus(45L, c(0, 0), 2/3, 0.98),
  runif_in_annulus(10L, c(0, 0), 1/3, 0.45),
  runif_in_sphere(5L, d = 2, r = 0.25)
)

hdel <- hdelaunay(pts)
svg(file.path(path, "hdelaunay.svg"))
plotHdelaunay(hdel, circle = FALSE)
dev.off()

rsvg_png(file.path(path, "hdelaunay.svg"), file.path(path, "hdelaunay.png"))
