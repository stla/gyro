# Triakis octahedron
library(gyro)
library(rgl)
library(trekcolors)

C0 <- 1 + sqrt(2)
vertices <- rbind(
  V0  = c( 0.0,  0.0,   C0),
  V1  = c( 0.0,  0.0,  -C0),
  V2  = c(  C0,  0.0,  0.0),
  V3  = c( -C0,  0.0,  0.0),
  V4  = c( 0.0,   C0,  0.0),
  V5  = c( 0.0,  -C0,  0.0),
  V6  = c( 1.0,  1.0,  1.0),
  V7  = c( 1.0,  1.0, -1.0),
  V8  = c( 1.0, -1.0,  1.0),
  V9  = c( 1.0, -1.0, -1.0),
  V10 = c(-1.0,  1.0,  1.0),
  V11 = c(-1.0,  1.0, -1.0),
  V12 = c(-1.0, -1.0,  1.0),
  V13 = c(-1.0, -1.0, -1.0)
)

faces <- 1 + rbind(
  c(  6,  0,  2 ),
  c(  6,  2,  4 ),
  c(  6,  4,  0 ),
  c(  7,  1,  4 ),
  c(  7,  4,  2 ),
  c(  7,  2,  1 ),
  c(  8,  0,  5 ),
  c(  8,  5,  2 ),
  c(  8,  2,  0 ),
  c(  9,  1,  2 ),
  c(  9,  2,  5 ),
  c(  9,  5,  1 ),
  c( 10,  0,  4 ),
  c( 10,  4,  3 ),
  c( 10,  3,  0 ),
  c( 11,  1,  3 ),
  c( 11,  3,  4 ),
  c( 11,  4,  1 ),
  c( 12,  0,  3 ),
  c( 12,  3,  5 ),
  c( 12,  5,  0 ),
  c( 13,  1,  5 ),
  c( 13,  5,  3 ),
  c( 13,  3,  1 )
)

edges <- do.call(rbind, lapply(1:nrow(faces), function(i){
  f <- sort(faces[i, ])
  rbind(
    c(f[1], f[2]),
    c(f[1], f[3]),
    c(f[2], f[3])
  )
}))
edges <- edges[!duplicated(edges), ]

s <- 0.7

open3d(windowRect = c(50, 50, 562, 562), zoom = 0.65)
bg3d(rgb(54, 57, 64, maxColorValue = 255))
light3d(diffuse = "orangered")
for(i in 1L:nrow(faces)){
  v1 <- vertices[faces[i, 1], ]
  v2 <- vertices[faces[i, 2], ]
  v3 <- vertices[faces[i, 3], ]
  mesh <- gyrotriangle(
    v1, v2, v3, s = s,
    palette = hcl.colors(256, palette = "Spectral"),#palette.colors(palette = "Classic Tableau", alpha = 1),
    bias = 1.5, interpolate = "spline", g = function(u) u^2
  )
  shade3d(mesh, specular = "violetred")
}
for(i in 1L:nrow(edges)){
  idx <- edges[i, ]
  A <- vertices[idx[1], ]
  B <- vertices[idx[2], ]
  edge <- gyrotube(A, B, s = s, radius = 0.03)
  shade3d(edge, color = "deeppink")
}
spheres3d(vertices, radius = 0.04, color = "deeppink")


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

command <- "gifski --fps=10 --frames=zzpic*.png -o triakisOctahedron.gif"
system(command)

pngfiles <- list.files(pattern = "^zzpic?.*png$")
file.remove(pngfiles)
