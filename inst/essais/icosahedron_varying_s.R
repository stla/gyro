# hyperbolic icosahedron vayring s ####
library(gyro)
library(rgl)
library(Rvcg) # to get the edges with the `vcgGetEdge` function
library(trekcolors)
icosahedron <- icosahedron3d() # mesh with 12 vertices, 20 triangles
vertices <- t(icosahedron$vb[-4, ])
triangles <- t(icosahedron$it)
edges <- as.matrix(vcgGetEdge(icosahedron)[, c("vert1", "vert2")])
s_ <- seq(0.2, 1, length.out = 40)
for(k in 1:length(s_)){
  s <- s_[k]
  open3d(windowRect = c(50, 50, 562, 562))
  view3d(30, 30, zoom = 0.7)
  for(i in 1:nrow(triangles)){
    triangle <- triangles[i, ]
    A <- vertices[triangle[1], ]
    B <- vertices[triangle[2], ]
    C <- vertices[triangle[3], ]
    gtriangle <- gyrotriangle(
      A, B, C, s,
      palette = trek_pal("klingon"), bias = 1.5, interpolate = "spline"
    )
    shade3d(gtriangle)
  }
  for(i in 1:nrow(edges)){
    edge <- edges[i, ]
    A <- vertices[edge[1], ]
    B <- vertices[edge[2], ]
    gtube <- gyrotube(A, B, s, radius = 0.03)
    shade3d(gtube, color = "orangered")
  }
  spheres3d(vertices, radius = 0.05, color = "orangered")
  rgl.snapshot(sprintf("zzpic%03d.png", k))
  close3d()
}


command <-
  "gifski --fps=10 --frames=zzpic*.png -b -o icosahedron_varying_s.gif"
system(command)

pngfiles <- list.files(pattern = "^zzpic?.*png$")
file.remove(pngfiles)
