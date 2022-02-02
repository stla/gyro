library(gyro)
library(rgl)

## ~~ rhombicosidodecahedron ~~##

changesOfSign <- function(M, changes = "all") {
  if (!is.matrix(M)) M <- rbind(M)
  if (identical(changes, "all")) changes <- 1L:ncol(M)
  `colnames<-`(as.matrix(do.call(rbind, apply(M, 1L, function(row) {
    expand.grid(
      purrr::imap(row, ~ (if (.x == 0 || !.y %in% changes) .x else c(-.x, .x)))
    )
  }))), NULL)
}

phi <- (1 + sqrt(5)) / 2
vs1 <- rbind(
  c(1, 1, phi^3),
  c(phi^2, phi, 2 * phi),
  c(2 + phi, 0, phi^2)
)
vs2 <- rbind(vs1, vs1[, c(2, 3, 1)], vs1[, c(3, 1, 2)]) # even permutations
vs <- changesOfSign(vs2)

plotGyrohull3d(vs)

# square faces and their centers
library(cxhull)
h <- cxhull(vs)
nvertices <- sapply(h[["facets"]], function(f) length(f[["vertices"]]))
indices <- which(nvertices == 4)
vertices <- vapply(indices, function(i) h[["facets"]][[i]][["vertices"]], integer(4))
centers <- apply(vertices, 2, function(ijkl) {
  (vs[ijkl[1], ] + vs[ijkl[2], ] + vs[ijkl[3], ] + vs[ijkl[4], ]) / 4
})

polygonize <- function(edges) { # function which orders the vertices
  nedges <- nrow(edges)
  es <- edges[1, ]
  i <- es[2]
  edges <- edges[-1, ]
  for (. in 1:(nedges - 2)) {
    j <- which(apply(edges, 1, function(e) i %in% e))
    i <- edges[j, ][which(edges[j, ] != i)]
    es <- c(es, i)
    edges <- edges[-j, ]
  }
  es
}
squares <-
  vapply(indices, function(r) polygonize(h[["facets"]][[r]][["edges"]]), integer(4))


# mesh
faces <- matrix(integer(3 * 6 * length(indices)), ncol = 3, nrow = 6 * length(indices))
for (j in 1L:length(indices)) {
  v1 <- squares[1, j]
  v2 <- squares[2, j]
  v3 <- squares[3, j]
  v4 <- squares[4, j]
  v5 <- j + nrow(vs)
  faces[6 * (j - 1) + 1, ] <- c(v1, v2, v3)
  faces[6 * (j - 1) + 2, ] <- c(v1, v3, v4)
  faces[6 * (j - 1) + 3, ] <- c(v1, v2, v5)
  faces[6 * (j - 1) + 4, ] <- c(v2, v3, v5)
  faces[6 * (j - 1) + 5, ] <- c(v3, v4, v5)
  faces[6 * (j - 1) + 6, ] <- c(v4, v1, v5)
}

# vertices and edges including pyramids
Vertices <- rbind(vs, t(centers) * 1.5)
edges <- do.call(rbind, lapply(1:nrow(faces), function(i) {
  f <- faces[i, ]
  rbind(
    c(f[1], f[2]),
    c(f[1], f[3]),
    c(f[2], f[3])
  )
}))
edges <- edges[!duplicated(edges), ]
edgeLengths <-
  apply(edges, 1, function(ij) crossprod(Vertices[ij[1], ] - Vertices[ij[2], ]))
edges <- edges[-which(edgeLengths == 8), ]


# plot
s <- 1.8
theta <- 1
h <- 0
v <- c(h * cos(theta), h * sin(theta), 0)
library(rgl)
open3d(windowRect = c(50, 50, 562, 562), zoom = 0.65)
bg3d(rgb(54, 57, 64, maxColorValue = 255))
spheres3d(t(t(Vertices) + v), radius = 0.06, color = "orange")
for (i in 1L:nrow(edges)) {
  idx <- edges[i, ]
  A <- Vertices[idx[1], ]
  B <- Vertices[idx[2], ]
  edge <- gyrotube(A + v, B + v, s = s, radius = 0.05)
  shade3d(edge, color = "yellow")
}
for (i in 1L:nrow(faces)) {
  v1 <- Vertices[faces[i, 1], ] + v
  v2 <- Vertices[faces[i, 2], ] + v
  v3 <- Vertices[faces[i, 3], ] + v
  mesh <- gyrotriangle(v1, v2, v3, s = s)
  shade3d(mesh, color = "midnightblue")
}
