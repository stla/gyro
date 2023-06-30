library(AlphaHull3D)
library(rgl)

rho <- sqrt((5 - sqrt(5))/10)
vs1 <- t(vapply(0:4, function(i){
  c(rho*cos(2*i*pi/5), rho*sin(2*i*pi/5), 0)
}, numeric(3L)))
R <- sqrt((25 - 11*sqrt(5))/10)
vs2 <- t(vapply(0:4, function(i){
  c(R*cos(2*i*pi/5 + pi/5), R*sin(2*i*pi/5 + pi/5), 0)
}, numeric(3L)))
pts <- rbind(vs1, vs2, c(0, 0, 0.1), c(0, 0, -0.1))



ahull <- fullAhull3d(pts)
mesh <- setAlpha(ahull, alpha = 0.2)
mesh$normals <- NULL

open3d(windowRect = c(50, 50, 512, 512))
view3d(15, -15, zoom = 0.7)
shade3d(mesh, color = "darkorange")
wire3d(mesh)

library(gyro)
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(15, 15, zoom = 0.85)
plotGyroMesh(mesh, s = 0.2, tubesRadius = 0.005, spheresRadius = 0.01)

