isPositiveNumber <- function(x){
  is.numeric(x) && length(x) == 1L && !is.na(x) && x > 0
}

isPoint <- function(A){
  is.atomic(A) && is.numeric(A) && length(A) >= 2L && !anyNA(A)
}

isPositiveInteger <- function(m){
  isPositiveNumber(m) && floor(m) == m
}

dotprod <- function(x, y = NULL){
  c(crossprod(x, y))
}

