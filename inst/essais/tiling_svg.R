library(gyro)
library(rsvg)

path <- "./inst/trash"

svg(file.path(path, "htiling.svg"))
tiling(3, 7, 9, circle = FALSE)
dev.off()

rsvg_png(file.path(path, "htiling.svg"), file.path(path, "htiling.png"))
