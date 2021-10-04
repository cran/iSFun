##' @title Sparse principal component analysis
##'
##' @description This function provides penalty-based integrative sparse principal component analysis to obtain the direction of first principal component of a given dataset with high dimensions.
##'
##' @param x data matrix of explanatory variables.
##' @param mu1 numeric, sparsity penalty parameter.
##' @param eps numeric, the threshold at which the algorithm terminates.
##' @param scale.x character, "TRUE" or "FALSE", whether or not to scale the variables x. The default is TRUE.
##' @param maxstep numeric, maximum iteration steps. The default value is 50.
##' @param trace character, "TRUE" or "FALSE". If TRUE, prints out its screening results of variables.
##'
##' @return An 'spca' object that contains the list of the following items.
##' \itemize{
##' \item{x:}{ data matrix of explanatory variables with centered columns. If scale.x is TRUE, the columns of data matrix are standardized to have mean 0 and standard deviation 1.}
##' \item{eigenvalue:}{ the estimated first eigenvalue.}
##' \item{eigenvector:}{ the estimated first eigenvector.}
##' \item{component:}{ the estimated first principal component.}
##' \item{variable:}{ the screening results of variables.}
##' \item{meanx:}{ column mean of the original dataset x.}
##' \item{normx:}{ column standard deviation of the original dataset x.}
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
##' library(iSFun)
##' data("simData.pca")
##' x <- simData.pca$x[[1]]
##'
##' res <- spca(x = x, mu1 = 0.08, trace = TRUE)

spca <- function(x, mu1, eps =1e-4, scale.x = TRUE, maxstep = 50, trace = FALSE) {

  # initialization
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  ip <- c(1:p)

  # center & scale x
  one <- matrix(1, 1, n)
  meanx <- drop(one %*% x) / n
  x <- scale(x, meanx, FALSE)

  if (scale.x) {
    normx <- sqrt(drop(one %*% (x^2)) / (n - 1))
    if (any(normx < .Machine$double.eps)) {
      stop("Some of the columns of the predictor matrix have zero variance.")
    }
    x <- scale(x, FALSE, normx)
  } else {
    normx <- rep(1, p)
  }

  ro <- function(x, mu, alpha) {
    f <- function(x) mu * (1 > x / (mu * alpha)) * (1 - x / (mu * alpha))
    r <- integrate(f, 0, x)
    return(r)
  }

  ro_d1st <- function(x, mu, alpha) {
    r <- mu * (1 > x / (mu * alpha)) * (1 - x / (mu * alpha))
    return(r)
  }

  # initilize objects
  what <- matrix(0, p, 1)
  Z <- irlba(x, nu = 1, nv = 1)
  U <- Z$v * Z$d[1]
  V <- Z$u

  u <- U
  v <- V
  iter <-1
  dis.u <- 10
  loading_trace <- matrix(0, nrow = p, ncol = maxstep)

  # main iteration
  if (trace) {
    cat("The variables that join the set of selected variables at final step:\n")
  }
  while (dis.u > eps & iter <= maxstep) {
    u.old <- u
    v.old <- v
    s <- colSums( x * matrix(rep(v, p), ncol=p, nrow = n) ) / n
    ro_d <- mapply(function(j) ro_d1st(s[j], mu1, 6), 1:p)
    fun.c <- function(j) {
      c <- sign(s[j]) * (abs(s[j]) > ro_d[j]) * (abs(s[j]) - ro_d[j])
      return(c)
    }
    u <- mapply(fun.c, 1:p)
    if ( sum(abs(u) <= 1e-4 ) == p ) {
      cat("The value of mu1 is too large");
      what_cut <- u
      break}
    v <- x %*% u / sqrt( sum( (x %*% u)^2 ) )

    u_norm <- sqrt(sum(u^2))
    u_norm <- ifelse(u_norm == 0, 0.0001, u_norm)
    u.scale <- u / u_norm
    dis.u <- sqrt(sum((u - u.old)^2)) / sqrt(sum(u.old^2))

    what <- u.scale
    what_cut <- ifelse(abs(what) > 1e-4, what, 0)
    what_cut_norm <- sqrt(sum(what_cut^2))
    what_cut_norm <- ifelse(what_cut_norm == 0, 0.0001, what_cut_norm)
    what_dir <- what_cut / what_cut_norm
    what_cut <- ifelse(abs(what_dir) > 1e-4, what_dir, 0)
    loading_trace[, iter] <- as.numeric(what_cut)

    iter <- iter + 1
    if ( sum(abs(u.scale) <= 1e-4 ) == p ) {
      cat("The value of mu1 is too large");
      break}
  }
  loading_trace <- loading_trace[,1:(iter-1)]

  # normalization
  what <- what_cut

  # selected variables
  new2A <- which(what != 0)
  if (trace) {
    if (length(new2A) <= 10) {
      cat(paste("DataSet: \n", sep = ""))
      cat(paste("X", new2A, ", ", sep = " "))
      cat("\n")
    } else {
      cat(paste("DataSet: \n", sep = ""))
      nlines <- ceiling(length(new2A) / 10)
      for (i in 0:(nlines - 2))
      {
        cat(paste("X", new2A[(10 * i + 1):(10 * (i + 1))], ", ", sep = " "))
        cat("\n")
      }
      cat(paste("X", new2A[(10 * (nlines - 1) + 1):length(new2A)], ", ", sep = " "))
      cat("\n")
    }
  }

  eigenvalue <- t(what) %*% cov(x) %*% what
  comp <- x %*% what

  # return objects
  object <- list(
    x = x, eigenvalue = eigenvalue, eigenvector = what, component = comp, variable = new2A,
    meanx = meanx, normx = normx, mu1 = mu1
  )
  class(object) <- "spca"
  return(object)
}


