library(gyro)
library(trekcolors)

phi <- (1 + sqrt(5)) / 2
theta <- head(seq(0, pi/2, length.out = 11), -1L)
a <- phi^(sqrt(2*theta/pi) - 1)
a <- phi^((2*theta/pi)^0.8 - 1)
u <- a * cos(theta)
v <- a * sin(theta)
x <- c(0, u, -v, -u, v)
y <- c(0, v, u, -v, -u)
pts <- cbind(x, y) / 1.03

hdel <- hdelaunay(pts, exact = TRUE)

fcolor <- function(t){
  RGB <- colorRamp(trek_pal("klingon"))(t)
  rgb(RGB[, 1L], RGB[, 2L], RGB[, 3L], maxColorValue = 255)
}

plotHdelaunay(
  hdel, vertices = FALSE, circle = FALSE, color = fcolor
)


theta = seq(0, 8*pi, length.out = 16)
x = cos(theta) * phi^(theta/pi)
y = sin(theta) * phi^(theta/pi)

theta = seq(0, pi/2, length.out = 11)[-11]
u = cos(theta) * phi^(sqrt(2*theta/pi) - 1)
v = sin(theta) * phi^(sqrt(2*theta/pi) - 1)
x <- c(u, -v, -u, v)
y <- c(v, u, -v, -u)

pts = cbind(c(0, x), c(0, y))/1.03
plot(pts, asp = 1, # xlim = c(-15, 15), ylim = c(-15, 15),
     xlab = NA, ylab = NA, axes = FALSE)

hdel = hdelaunay(pts, exact = T)
plotHdelaunay(hdel, color = "random", luminosity = "bright")

plot(pts, asp = 1, # xlim = c(-15, 15), ylim = c(-15, 15),
     xlab = NA, ylab = NA, axes = FALSE)

opar <- par(mar = c(0, 0, 0, 0), bg = "black")
plot(pts, asp = 1, # xlim = c(-15, 15), ylim = c(-15, 15),
     xlab = NA, ylab = NA, axes = FALSE)

del <- delaunay(pts)
v <- voronoi(del)
plotVoronoiDiagram(v , colors = viridisLite::turbo(281))
plotVoronoiDiagram(v , colors = randomcoloR::distinctColorPalette(281))
plotVoronoiDiagram(v, luminosity = "dark")


fplot <- function(){
  opar <- par(mar = c(0, 0, 0, 0), bg = "black")
  plot(NULL, asp = 1, xlim = c(-15, 15), ylim = c(-15, 15),
       xlab = NA, ylab = NA, axes = FALSE)
  plotVoronoiDiagram(v , colors = viridisLite::turbo(281))
}
