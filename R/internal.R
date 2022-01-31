isBoolean <- function(x){
  is.atomic(x) && is.logical(x) && length(x) == 1L && !is.na(x)
}

isPositiveNumber <- function(x){
  is.numeric(x) && length(x) == 1L && !is.na(x) && x > 0
}

isPoint <- function(A){
  is.atomic(A) && is.numeric(A) && length(A) >= 2L && !anyNA(A)
}

is3dPoint <- function(A){
  is.atomic(A) && is.numeric(A) && length(A) == 3L && !anyNA(A)
}

isPositiveInteger <- function(m){
  isPositiveNumber(m) && floor(m) == m
}

areDistinct <- function(A, B){
  !isTRUE(all.equal(A, B))
}

dotprod <- function(x, y = NULL){
  c(crossprod(x, y))
}

