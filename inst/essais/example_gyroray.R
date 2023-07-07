library(gyro)
# a 2D example ####
O <- c(1, 2); A <- c(1, 1)
opar <- par(mar = c(2, 2, 2, 0.5))
plot(rbind(O, A), type = "p", pch = 19, xlab = NA, ylab = NA,
     xlim = c(0, 2), ylim = c(0, 3), main = "s = 0.3")
s <- 0.3
ray <- gyroray(O, A, s)
lines(ray, col = "blue", lwd = 2)
text(t(O), expression(italic(O)), pos = 2)
text(t(A), expression(italic(A)), pos = 3)
# opposite gyroray
yar <- gyroray(O, A, s, OtoA = FALSE)
lines(yar, col = "red", lwd = 2)
par(opar)
