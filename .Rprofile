dllunload <- function(){
  dyn.unload(
    system.file("libs", "x64", "gyro.dll", package = "gyro")
  )
}

makedoc <- function(){
  roxygen2::roxygenise(load_code = roxygen2::load_installed)
}

myinstall <- function(){
  withr::with_libpaths("C:/SL/Rloclib", devtools::install(quick = TRUE))
}
