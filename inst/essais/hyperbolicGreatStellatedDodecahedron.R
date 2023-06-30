library(rgl)
library(Rvcg)
library(gyro)

mesh <- vcgPlyRead(
  "greatStellatedDodecahedron.ply", updateNormals = FALSE, clean = FALSE
)
Vertices <- mesh[["vb"]][-4L, ]
Triangles <- mesh[["it"]]
Edges <- as.matrix(vcgGetEdge(mesh)[, c(1L, 2L)])

s <- 1 
iterations <- 3L
facesColor <- "firebrick1"
n <- 50
edgesColor <- "black"

ntriangles <- ncol(Triangles)
Gtriangles <- vector("list", ntriangles)
palette <- NULL
for(i in 1L:ntriangles) {
  triangle <- Vertices[, Triangles[, i]]
  Gtriangles[[i]] <- gyrotriangle(
    triangle[, 1L], triangle[, 2L], triangle[, 3L], s = s, 
    model = "U", iterations = iterations
  )
  #   palette = palette, bias = bias, interpolate = interpolate, 
  # g = g)
}
mesh <- vcgClean(Morpho::mergeMeshes(Gtriangles), sel = 0, silent = TRUE)
if(is.null(palette)) {
  shade3d(mesh, color = facesColor, polygon_offset = 1)
} else {
  shade3d(mesh)
}
if(edgesAsTubes) {
  for(i in 1L:nrow(Edges)) {
    edge <- Edges[i, ]
    gtube <- gyrotube(
      Vertices[, edge[1L]], Vertices[, edge[2L]], 
      s = s, model = model, n = n, radius = tubesRadius
    )
    shade3d(gtube, color = edgesColor)
  }
} else {
  if(model == "M") {
    for(i in 1L:nrow(Edges)) {
      edge <- Edges[i, ]
      gsegment <- Mgyrosegment(
        Vertices[, edge[1L]], Vertices[, edge[2L]], s = s, n = n
      )
      lines3d(gsegment, color = edgesColor, lwd = 2, line_antialias = TRUE)
    }
  } else {
    for(i in 1L:nrow(Edges)) {
      edge <- Edges[i, ]
      gsegment <- gyro:::Ugyrosegment(
        Vertices[, edge[1L]], Vertices[, edge[2L]], s = s, n = n
      )
      lines3d(gsegment, color = edgesColor, lwd = 2, line_antialias = TRUE)
    }
  }
}
if(verticesAsSpheres) {
  spheres3d(t(Vertices), radius = spheresRadius, color = spheresColor)
}
