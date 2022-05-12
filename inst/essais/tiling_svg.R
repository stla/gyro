library(gyro)
library(rsvg)

path <- "./inst/trash"

svg(file.path(path, "htiling_4-6.svg"))
tiling(4, 6, 4, circle = FALSE, colors = c("#FFFA0C", "#590000"))
dev.off()

rsvg_png(file.path(path, "htiling_4-6.svg"), file.path(path, "htiling_4-6.png"))
