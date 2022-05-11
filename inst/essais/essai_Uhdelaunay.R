library(gyro)
library(randomcoloR)
library(plotrix)
library(uniformly)

s <- 1
phiEM <- function(A){
  gyro:::Mgyroscalar(2, A, s)
  # gamm <- gyro:::gammaF(A, s)
  # gamm*A / (1+gamm)
}
phiMU <- function(A){
  gyro:::PhiEU(phiEM(A), s)
}

plotHdelaunayU <- function(
    hdel, remove = NULL, vertices = TRUE, edges = TRUE, circle = TRUE,
    color = "distinct", hue = "random", luminosity = "random"
){
  if(!inherits(hdel, "hdelaunay")){
    stop("The `hdel` argument must be an output of the `hdelaunay` function.")
  }
  if(!is.null(remove)){
    isolations <- "medges" %in% names(hdel)
    if(!isolations){
      stop(
        "In order to use the `remove` argument you have to run the ",
        "`hdelaunay` function with `isolations=TRUE`."
      )
    }
  }
  opar <- par(mar = c(0, 0, 0, 0))
  pts <- t(apply(s*hdel[["vertices"]], 1L, phiMU))
  plot(
    pts, type = "p", asp = 1,
    xlab = NA, ylab = NA, axes = FALSE
  )
  if(circle){
    draw.circle(0, 0, radius = 1, border = "black")
  }
  if(length(color) > 1L || !is.na(color)){
    triangles <- hdel[["triangles"]]
    ntriangles <- nrow(triangles)
    if(length(color) > 1L){
      colors <- color
    }else{
      if(color == "random"){
        colors <- randomColor(ntriangles, hue = hue, luminosity = luminosity)
      }else if(color == "distinct"){
        colors <- distinctColorPalette(ntriangles)
      }else{
        colors <- rep(color, ntriangles)
      }
    }
    for(i in 1L:ntriangles){
      trgl <- triangles[i, ]
      hpolypath <- rbind(
        gyro:::Ugyrosegment(pts[trgl[1L], ], pts[trgl[3L], ], s = s, n = 50)[-1L, ],
        gyro:::Ugyrosegment(pts[trgl[3L], ], pts[trgl[2L], ], s = s, n = 50)[-1L, ],
        gyro:::Ugyrosegment(pts[trgl[2L], ], pts[trgl[1L], ], s = s, n = 50)[-1L, ]
      )
      polypath(hpolypath, border = NA, col = colors[i])
    }
  }
  if(edges){
    if("iedges" %in% remove){
      hedges <- hdel[["medges"]]
    }else{
      hedges <- hdel[["edges"]]
    }
    for(i in 1L:nrow(hedges)){
      hedge <- hedges[i, ]
      hseg <- gyro:::Ugyrosegment(pts[hedge[1L], ], pts[hedge[2L], ], s = s, n = 50)
      lines(hseg, lty = "solid", col = "black", lwd = 1.5)
    }
  }
  if(vertices){
    if("ivertices" %in% remove){
      pts <- hdel[["mvertices"]]
    }
    points(pts, pch = 19, cex = 0.9)
  }
  par(opar)
  invisible(NULL)
}

set.seed(314)

pts <- rbind(
  runif_in_annulus(20L, c(0, 0), 0.88, 0.9),
  runif_in_annulus(10L, c(0, 0), 0.85, 0.87),
  runif_in_sphere(5L, d = 2, r = 0.83)
)

hdel <- hdelaunay(pts, exact = TRUE)
plotHdelaunayU(hdel, circle = FALSE)

path <- "./inst/trash"
svg(file.path(path, "Uhdelaunay.svg"))
plotHdelaunayU(hdel, circle = FALSE)
dev.off()

rsvg::rsvg_png(file.path(path, "Uhdelaunay.svg"), file.path(path, "Uhdelaunay.png"))
