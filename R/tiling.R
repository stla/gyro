circumcircle <- function(p1, p2, p3){
  x1 <- p1[1L]
  y1 <- p1[2L]
  x2 <- p2[1L]
  y2 <- p2[2L]
  x3 <- p3[1L]
  y3 <- p3[2L]
  a <- det(cbind(rbind(p1, p2, p3), 1))
  q1 <- dotprod(p1)
  q2 <- dotprod(p2)
  q3 <- dotprod(p3)
  q <- c(q1, q2, q3)
  x <- c(x1, x2, x3)
  y <- c(y1, y2, y3)
  Dx <- det(cbind(q, y, 1))
  Dy <- -det(cbind(q, x, 1))
  c <- det(cbind(q, x, y))
  center <- 0.5 * c(Dx, Dy) / a
  r <- sqrt(dotprod(center - p1))
  list("center" = center, "radius" = r)
}

inversion <- function(circle, M) {
  v <- M - circle[["center"]]
  circle[["center"]] + circle[["radius"]]^2 * v / dotprod(v)
}

