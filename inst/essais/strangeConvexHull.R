library(gyro)

pt <- function(x){
  c(
    sin(x) * cos(2 * x),
    sin(x) * sin(2 * x),
    cos(x)
  )
}
pts <- t(vapply(seq(0, pi, length.out = 50), pt, numeric(3L)))

open3d(windowRect = c(50, 50, 562, 562))
view3d(zoom = 0.9)
plotGyrohull3d(pts, s = 0.8, tubesRadius = 0.02, spheresRadius = 0.04)

# animation ####
M <- par3d("userMatrix")
movie3d(
  par3dinterp(
    time = seq(0, 1, len = 9),
    userMatrix = list(
      M,
      rotate3d(M, pi, 1, 0, 0),
      rotate3d(M, pi, 1, 1, 0),
      rotate3d(M, pi, 1, 1, 1),
      rotate3d(M, pi, 0, 1, 1),
      rotate3d(M, pi, 0, 1, 0),
      rotate3d(M, pi, 1, 0, 1),
      rotate3d(M, pi, 0, 0, 1),
      M
    )
  ),
  fps = 120,
  duration = 1,
  dir = ".",
  movie = "zzpic",
  convert = FALSE,
  clean = FALSE,
  webshot = FALSE
)

command <- "gifski --fps=10 --frames=zzpic*.png -o strangeConvexHull.gif"
system(command)

pngfiles <- list.files(pattern = "^zzpic?.*png$")
file.remove(pngfiles)
