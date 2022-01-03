##' @title Statistical description before using function ispls
##'
##' @description The function describes the basic statistical information of the data, including sample mean, sample variance of X and Y, the first direction of partial least squares method, etc.
##'
##' @param x list of data matrices, L datasets of explanatory variables.
##' @param y list of data matrices, L datasets of dependent variables.
##' @param L numeric, number of datasets.
##' @param scale.x character, "TRUE" or "FALSE", whether or not to scale the variables x. The default is TRUE.
##' @param scale.y character, "TRUE" or "FALSE", whether or not to scale the variables y. The default is TRUE.
##'
##' @return A 'preview.pls' object that contains the list of the following items.
##' \itemize{
##' \item{x:}{ list of data matrices, L datasets of explanatory variables with centered columns. If scale.x is TRUE, the columns of L datasets are standardized to have mean 0 and standard deviation 1.}
##' \item{y:}{ list of data matrices, L datasets of dependent variables with centered columns. If scale.y is TRUE, the columns of L datasets are standardized to have mean 0 and standard deviation 1.}
##' \item{loading:}{ the estimated first direction vector.}
##' \item{meanx:}{ list of numeric vectors, column mean of the original datasets x.}
##' \item{normx:}{ list of numeric vectors, column standard deviation of the original datasets x.}
##' \item{meany:}{ list of numeric vectors, column mean of the original datasets y.}
##' \item{normy:}{ list of numeric vectors, column standard deviation of the original datasets y.}
##' }
##'
##' @seealso See Also as \code{\link{ispls}}.
##'
##' @import caret
##' @import irlba
##' @import graphics
##' @import stats
##' @importFrom grDevices rainbow
##' @export
##' @examples
##' library(iSFun)
##' data("simData.pls")
##' x <- simData.pls$x
##' y <- simData.pls$y
##' L <- length(x)
##'
##' prev_pls <- preview.pls(x = x, y = y, L = L, scale.x = TRUE, scale.y = TRUE)

preview.pls <- function(x, y, L, scale.x = TRUE, scale.y = TRUE){

  if (class(x) != "list") { stop("x should be of list type.") }
  if (class(y) != "list") { stop("y should be of list type.") }

  # initialization
  x  <- lapply(x, as.matrix)
  y  <- lapply(y, as.matrix)
  nl <- as.numeric(lapply(x, nrow))
  pl <- as.numeric(lapply(x, ncol))
  ql <- as.numeric(lapply(y, ncol))
  p  <- unique(pl)
  q  <- unique(ql)
  if(length(p) > 1){ stop("The dimension of data x should be consistent among different datasets.")}
  if(length(q) > 1){ stop("The dimension of data y should be consistent among different datasets.")}
  ip <- c(1:p)

  # center & scale x & y
  meanx <- lapply(1:L, function(l) drop( matrix(1, 1, nl[l]) %*% x[[l]] / nl[l] ) )
  meany <- lapply(1:L, function(l) drop( matrix(1, 1, nl[l]) %*% y[[l]] / nl[l] ) )
  x <- lapply(1:L, function(l) scale(x[[l]], meanx[[l]], FALSE) )
  y <- lapply(1:L, function(l) scale(y[[l]], meany[[l]], FALSE) )

  x.scale <- function(l){
    one <- matrix(1, 1, nl[l])
    normx <- sqrt(drop(one %*% (x[[l]]^2)) / (nl[l] - 1))
    if (any(normx < .Machine$double.eps)) {
      stop("Some of the columns of the predictor matrix have zero variance.")
    }
    return(normx)
  }
  y.scale <- function(l){
    one <- matrix(1, 1, nl[l])
    normy <- sqrt(drop(one %*% (y[[l]]^2)) / (nl[l] - 1))
    if (any(normy < .Machine$double.eps)) {
      stop("Some of the columns of the response matrix have zero variance.")
    }
    return(normy)
  }

  if (scale.x) { normx <- lapply(1:L, x.scale ) } else { normx <- rep(list(rep(1,p)), L) }
  if (scale.y) { normy <- lapply(1:L, y.scale ) } else { normy <- rep(list(rep(1,q)), L) }
  if (scale.x) { x <- lapply(1:L, function(l) scale(x[[l]], FALSE, normx[[l]]) ) }
  if (scale.y) { y <- lapply(1:L, function(l) scale(y[[l]], FALSE, normy[[l]]) ) }

  fun.1 <- function(l) {
    Z_l <- t(x[[l]]) %*% y[[l]]
    Z_l <- Z_l / nl[l]
  }
  Z <- matrix(mapply(fun.1, c(1:L)), nrow = p)

  c <- mapply(function(l) svd(Z[, ((l - 1) * q + 1):(l * q)] %*% t(Z[, ((l - 1) * q + 1):(l * q)] ), nu = 1)$u, 1:L)#/ (nl[l]^2)

  # normalization
  what <- c

  listname <- mapply(function(l) paste("Dataset ", l), 1:L)
  names(meanx) <- listname
  names(meany) <- listname
  names(normx) <- listname
  names(normy) <- listname
  names(x) <- listname
  names(y) <- listname

  for (l in 1:L) {
    plot(x = 1:p, y = as.numeric(what[,l]),
         main = paste("Dataset ", l, "\n", "The first direction vector"),
         xlab = "Dimension", ylab = "Value", pch = 15)
  }

  # return objects
  object <- list(
    x = x, y = y, loading = what,
    meanx = meanx, normx = normx, meany = meany, normy = normy
  )
  class(object) <- "preview.spls"

  return(object)
}
