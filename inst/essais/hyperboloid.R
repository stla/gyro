library(gyro)
library(rgl)

s <- 0.5
z <- 1
r <- 1
alpha_ <- seq(0, pi, length.out=100)[-1]
for(alpha in alpha_){
  A <- c(r * c(cos(alpha), -sin(alpha)), z)
  B <- c(r * c(cos(alpha+pi), -sin(alpha+pi)), z)
  shade3d(
    gyrotube(A, B, s=s, radius= 0.06, caps = TRUE),
    color = "yellow"
  )
}


s <- 0.5
z <- 5
x <- 4
y_ <- seq(-4, 4, length.out=100)
for(y in y_){
  A <- c(x, y, z)
  B <- c(-x, -y, z)
  shade3d(
    gyrotube(A, B, s=s, radius= 0.06),
    color = "yellow"
  )
}
