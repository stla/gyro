dllunload <- function(){
  dyn.unload(
    system.file("libs", "x64", "gyro.dll", package = "gyro")
  )
}

myinstall <- function() {
  try(dllunload())
  if(rstudioapi::isAvailable()) {
    rstudioapi::restartSession(
      "devtools::install(quick = TRUE, keep_source = TRUE)"
    )
  } else {
    devtools::install(quick = TRUE, keep_source = TRUE)
  }
}

mydocument <- function() {
  if(rstudioapi::isAvailable()) {
    rstudioapi::restartSession(
      "roxygen2::roxygenise(load_code = roxygen2::load_installed)"
    )
  } else {
    roxygen2::roxygenise(load_code = roxygen2::load_installed)
  }
}
