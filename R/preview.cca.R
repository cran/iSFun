##' @title Statistical description before using function iscca
##'
##' @description The function describes the basic statistical information of the data, including sample mean, sample variance of X and Y, and the first pair of canonical vectors.
##'
##' @param x list of data matrices, L datasets of explanatory variables.
##' @param y list of data matrices, L datasets of dependent variables.
##' @param L numeric, number of datasets.
##' @param scale.x character, "TRUE" or "FALSE", whether or not to scale the variables x. The default is TRUE.
##' @param scale.y character, "TRUE" or "FALSE", whether or not to scale the variables y. The default is TRUE.
##'
##' @return An 'preview.cca' object that contains the list of the following items.
##' \itemize{
##' \item{x:}{ list of data matrices, L datasets of explanatory variables with centered columns. If scale.x is TRUE, the columns of L datasets are standardized to have mean 0 and standard deviation 1.}
##' \item{y:}{ list of data matrices, L datasets of dependent variables with centered columns. If scale.y is TRUE, the columns of L datasets are standardized to have mean 0 and standard deviation 1.}
##' \item{loading.x:}{ the estimated canonical vector of variables x.}
##' \item{loading.y:}{ the estimated canonical vector of variables y.}
##' \item{meanx:}{ list of numeric vectors, column mean of the original datasets x.}
##' \item{normx:}{ list of numeric vectors, column standard deviation of the original datasets x.}
##' \item{meany:}{ list of numeric vectors, column mean of the original datasets y.}
##' \item{normy:}{ list of numeric vectors, column standard deviation of the original datasets y.}
##' }
##' @seealso See Also as \code{\link{iscca}}.
##'
##' @import caret
##' @import irlba
##' @import graphics
##' @import stats
##' @importFrom grDevices rainbow
##' @export
##' @examples
##' # Load a list with 3 data sets
##' library(iSFun)
##' data("simData.cca")
##' x <- simData.cca$x
##' y <- simData.cca$y
##' L <- length(x)
##'
##' prev_cca <- preview.cca(x = x, y = y, L = L, scale.x = TRUE, scale.y = TRUE)
##'


preview.cca <- function(x, y, L, scale.x = TRUE, scale.y = TRUE) {

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
  iq <- c(1:q)

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

  # define Z
  fun.1 <- function(l) {
    Z_l <- irlba( t(x[[l]]) %*% y[[l]] , nu =1, nv = 1)
    u_l <- Z_l$u
    return(u_l)
  }
  U <- matrix(mapply(fun.1, 1:L), nrow = p)

  fun.2 <- function(l) {
    Z_l <- irlba( t(x[[l]]) %*% y[[l]] , nu =1, nv = 1)
    v_l <- Z_l$v
    return(v_l)
  }
  V <- mapply(fun.2, 1:L)

  what.u <- U
  what.v <- V

  listname <- mapply(function(l) paste("Dataset ", l), 1:L)
  names(meanx) <- listname
  names(meany) <- listname
  names(normx) <- listname
  names(normy) <- listname
  names(x) <- listname
  names(y) <- listname
  colnames(what.u) <- listname
  rownames(what.u) <- c(1 : p)
  colnames(what.v) <- listname
  rownames(what.v) <- c(1 : q)

  plot_loading <- function(order){
    opar <- par(mfrow = c(1,2))
    on.exit(par(opar))
    for (l in order) {
      plot(x = 1:p, y = U[, l],
           main = paste("Dataset ", l, "\n", "The first canonical vector u"),
           xlab = "Dimension", ylab = "Value", pch = 15)
      plot(x = 1:q, y = V[, l],
           main = paste("Dataset ", l, "\n", "The first canonical vector v"),
           xlab = "Dimension", ylab = "Value", pch = 15)
    }
  }

  plot_loading(order = 1:L)

  # return objects
  object <- list(
    x = x, y = y, loading.x = what.u, loading.y = what.v,
    meanx = meanx, normx = normx, meany = meany, normy = normy)
  class(object) <- "preview.cca"
  return(object)
}


