##' @title Statistical description before using function ispca
##'
##' @description The function describes the basic statistical information of the data, including sample mean, sample co-variance of X and Y, the first eigenvector, eigenvalue and principal component, etc.
##'
##' @param x list of data matrices, L datasets of explanatory variables.
##' @param L numeric, number of data sets.
##' @param scale.x character, "TRUE" or "FALSE", whether or not to scale the variables x. The default is TRUE.
##'
##' @return An 'preview.pca' object that contains the list of the following items.
##' \itemize{
##' \item{x:}{ list of data matrices, L datasets of explanatory variables with centered columns. If scale.x is TRUE, the columns of L datasets are standardized to have mean 0 and standard deviation 1.}
##' \item{eigenvalue:}{ the estimated first eigenvalue.}
##' \item{eigenvector:}{ the estimated first eigenvector.}
##' \item{component:}{ the estimated first component.}
##' \item{meanx:}{ list of numeric vectors, column mean of the original datasets x.}
##' \item{normx:}{ list of numeric vectors, column standard deviation of the original datasets x.}
##' }
##' @seealso See Also as \code{\link{ispca}}.
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
##' data("simData.pca")
##' x <- simData.pca$x
##' L <- length(x)
##'
##' prev.pca <- preview.pca(x = x, L = L, scale.x = TRUE)
##'

preview.pca <- function(x, L, scale.x = TRUE) {

  if (class(x) != "list") { stop("x should be of list type.") }

  # initialization
  x  <- lapply(x, as.matrix)
  nl <- as.numeric(lapply(x, nrow))
  pl <- as.numeric(lapply(x, ncol))
  p  <- unique(pl)
  if(length(p) > 1){ stop("The dimension of data x should be consistent among different datasets.")}
  ip <- c(1:p)

  # center & scale x
  meanx <- lapply(1:L, function(l) drop( matrix(1, 1, nl[l]) %*% x[[l]] / nl[l] ) )
  x <- lapply(1:L, function(l) scale(x[[l]], meanx[[l]], FALSE) )

  x.scale <- function(l){
    one <- matrix(1, 1, nl[l])
    normx <- sqrt(drop(one %*% (x[[l]]^2)) / (nl[l] - 1))
    if (any(normx < .Machine$double.eps)) {
      stop("Some of the columns of the predictor matrix have zero variance.")
    }
    return(normx)
  }

  if (scale.x) { normx <- lapply(1:L, x.scale ) } else { normx <- rep(list(rep(1,p)), L) }
  if (scale.x) { x <- lapply(1:L, function(l) scale(x[[l]], FALSE, normx[[l]]) ) }

  # initilize objects
  what <- matrix(0, p, L)

  fun.1 <- function(l) {
    Z_l <- irlba(x[[l]], nu = 1, nv = 1)
    u_l <- Z_l$v
    return(u_l)
  }
  U <- matrix(mapply(fun.1, 1:L), nrow = p)

  what <- U

  eigenvalue <- mapply(function(l) (t(what[, l])%*%cov(x[[l]]) %*% what[, l]), 1:L)
  comp <- lapply(1:L, function(l) x[[l]] %*% what[, l])

  listname <- mapply(function(l) paste("Dataset ", l), 1:L)
  names(meanx) <- listname
  names(normx) <- listname
  names(x) <- listname
  colnames(what) <- listname
  rownames(what) <- c(1 : p)
  names(comp) <- listname
  names(eigenvalue) <- listname

  for (l in 1:L) {
    plot(x = 1:p, y = what[,l],
         main = paste("The first principal component of dataset", l),
         xlab = "Dimension", ylab = "Value", pch = 15)
  }

  heat <- lapply(1:L, function(l) cov(x[[l]]) )

  rc <- rainbow(nrow(heat[[1]]), start = 0, end = .3)
  cc <- rainbow(ncol(heat[[1]]), start = 0, end = .3)
  for (l in 1:L) {
    heatmap(x = heat[[l]], RowSideColors = cc, ColSideColors = rc,
            Rowv = NA, Colv = NA, scale = "row", xlab = paste("Dataset ", l, ": X"), ylab = "Y",
            main = paste("Heatmap of covariance"))
  }

  # return objects
  object <- list(
    x = x, eigenvalue = eigenvalue, eigenvector = what, component = comp,
    meanx = meanx, normx = normx)
  class(object) <- "preview.pca"
  return(object)
}


