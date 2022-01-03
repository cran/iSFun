##' @title Sparse canonical correlation analysis
##'
##' @description This function provides penalty-based sparse canonical correlation analysis to get the first pair of canonical vectors.
##'
##' @param x data matrix of explanatory variables
##' @param y data matrix of dependent variables.
##' @param mu1 numeric, sparsity penalty parameter for vector u.
##' @param mu2 numeric, sparsity penalty parameter for vector v.
##' @param eps numeric, the threshold at which the algorithm terminates.
##' @param scale.x character, "TRUE" or "FALSE", whether or not to scale the variables x. The default is TRUE.
##' @param scale.y character, "TRUE" or "FALSE", whether or not to scale the variables y. The default is TRUE.
##' @param maxstep numeric, maximum iteration steps. The default value is 50.
##' @param trace character, "TRUE" or "FALSE". If TRUE, prints out its screening results of variables.
##'
##' @return An 'scca' object that contains the list of the following items.
##' \itemize{
##' \item{x:}{ data matrix of explanatory variables with centered columns. If scale.x is TRUE, the columns of data matrix are standardized to have mean 0 and standard deviation 1.}
##' \item{y:}{ data matrix of dependent variables with centered columns. If scale.y is TRUE, the columns of data matrix are standardized to have mean 0 and standard deviation 1.}
##' \item{loading.x:}{ the estimated canonical vector of variables x.}
##' \item{loading.y:}{ the estimated canonical vector of variables y.}
##' \item{variable.x:}{ the screening results of variables x.}
##' \item{variable.y:}{ the screening results of variables y.}
##' \item{meanx:}{ column mean of the original dataset x.}
##' \item{normx:}{ column standard deviation of the original dataset x.}
##' \item{meany:}{ column mean of the original dataset y.}
##' \item{normy:}{ column standard deviation of the original dataset y.}
##' }
##' @seealso See Also as \code{\link{iscca}}, \code{\link{meta.scca}}.
##'
##' @import caret
##' @import irlba
##' @import graphics
##' @import stats
##' @importFrom grDevices rainbow
##' @export
##' @examples
##' library(iSFun)
##' data("simData.cca")
##' x.scca <- do.call(rbind, simData.cca$x)
##' y.scca <- do.call(rbind, simData.cca$y)
##' res_scca <- scca(x = x.scca, y = y.scca, mu1 = 0.1, mu2 = 0.1, eps = 1e-3,
##'                  scale.x = TRUE, scale.y = TRUE, maxstep = 50, trace = FALSE)

scca <- function(x, y, mu1, mu2, eps =1e-4, scale.x = TRUE, scale.y = TRUE, maxstep = 50, trace = FALSE) {

  # initialization
  x <- as.matrix(x)
  y <- as.matrix(y)
  nx <- nrow(x)
  ny <- nrow(y)
  p <- ncol(x)
  q <- ncol(y)
  if(nx != ny){ stop("The rows of data x and data y should be consistent.")}
  n <- nx
  ip <- c(1:p)
  iq <- c(1:q)

  # center & scale x
  one <- matrix(1, 1, n)
  meany <- drop(one %*% y / n)
  y <- scale(y, meany, FALSE)
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

  if (scale.y) {
    normy <- sqrt(drop(one %*% (y^2)) / (n - 1))
    if (any(normy < .Machine$double.eps)) {
      stop("Some of the columns of the response matrix have zero variance.")
    }
    y <- scale(y, FALSE, normy)
  } else {
    normy <- rep(1, q)
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

  Z <- irlba(t(x) %*% y, nu = 1, nv = 1)
  U <- Z$u
  V <- Z$v

  u <- U
  v <- V
  iter <- 1
  dis.u <- 10
  dis.v <- 10

  if (trace) {
    cat("The variables that join the set of selected variables at final step:\n")
  }

  while (dis.u > eps & dis.v > eps & iter <= maxstep){
    u.old <- u
    v.old <- v

    s <- t(v) %*% t(y) %*% x / (n)
    ro_d <- mapply(function(j) ro_d1st(s[j], mu1, 6), 1:p)
    fun.c <- function(j) {
      c <- sign(s[j]) * (abs(s[j]) > ro_d[j]) * (abs(s[j]) - ro_d[j])
      return(c)
    }
    u <- mapply(fun.c, 1:p)
    u_norm <- sqrt(sum(u^2))
    u_norm <- ifelse(u_norm == 0, 0.0001, u_norm)
    u <- u / u_norm
    dis.u <- sqrt(sum((u - u.old)^2)) / sqrt(sum(u.old^2))

    what <- u
    what_cut <- ifelse(abs(what) > 1e-4, what, 0)
    what_cut_norm <- sqrt(sum(what_cut^2))
    what_cut_norm <- ifelse(what_cut_norm == 0, 0.0001, what_cut_norm)
    what_dir <- what_cut / what_cut_norm
    what_cut <- ifelse(abs(what_dir) > 1e-4, what_dir, 0)

    s <- t(u) %*% t(x) %*% y / (n)
    ro_d <- mapply(function(j) ro_d1st(s[j], mu2, 6), 1:q)
    fun.c <- function(j) {
      c <- sign(s[j]) * (abs(s[j]) > ro_d[j]) * (abs(s[j]) - ro_d[j])
      return(c)
    }
    v <- mapply(fun.c, 1:q)
    v_norm <- sqrt(sum(v^2))
    v_norm <- ifelse(v_norm == 0, 0.0001, v_norm)
    v <- v / v_norm
    dis.v <- sqrt(sum((v - v.old)^2)) / sqrt(sum(v.old^2))

    what_v <- v
    what_cut_v <- ifelse(abs(what_v) > 1e-4, what_v, 0)
    what_cut_norm <- sqrt(sum(what_cut_v^2))
    what_cut_norm <- ifelse(what_cut_norm == 0, 0.0001, what_cut_norm)
    what_dir <- what_cut_v / what_cut_norm
    what_cut_v <- ifelse(abs(what_dir) > 1e-4, what_dir, 0)

    dis.u <- sqrt(sum((u - u.old)^2)) / sqrt(sum(u.old^2))
    dis.v <- sqrt(sum((v - v.old)^2)) / sqrt(sum(v.old^2))
    if ( sum(abs(v) <= 1e-4 ) == q & sum(abs(u) <= 1e-4 ) == p) {
      cat("Stop! The values of mu1 and mu2 are too large");
      break }
    if ( sum(abs(v) <= 1e-4 ) != q & sum(abs(u) <= 1e-4 ) == p) {
      cat("Stop! The value of mu1 is too large");
      break }
    if ( sum(abs(v) <= 1e-4 ) == q & sum(abs(u) <= 1e-4 ) != p) {
      cat("Stop! The value of mu2 is too large");
      break }

    if (trace) {
      new2A_u <- ip[what_cut != 0]
      new2A_v <- iq[what_cut_v != 0]
      cat("\n")
      cat(paste("--------------------", "\n"))
      cat(paste("----- Step", iter, " -----\n", sep = " "))
      cat(paste("--------------------", "\n"))
      new2A_l <- new2A_u
      if (length(new2A_l) <= 10) {
        cat(paste("X", new2A_l, ", ", sep = " "))
        cat("\n")
      } else {
        nlines <- ceiling(length(new2A_l) / 10)
        for (i in 0:(nlines - 2))
        {
          cat(paste("X", new2A_l[(10 * i + 1):(10 * (i + 1))], ", ", sep = " "))
          cat("\n")
        }
        cat(paste("X", new2A_l[(10 * (nlines - 1) + 1):length(new2A_l)], ", ", sep = " "))
        cat("\n")
      }

      new2A_l <- new2A_v
      if (length(new2A_l) <= 10) {
        cat(paste("Y", new2A_l, ", ", sep = " "))
        cat("\n")
      } else {
        nlines <- ceiling(length(new2A_l) / 10)
        for (i in 0:(nlines - 2))
        {
          cat(paste("Y", new2A_l[(10 * i + 1):(10 * (i + 1))], ", ", sep = " "))
          cat("\n")
        }
        cat(paste("Y", new2A_l[(10 * (nlines - 1) + 1):length(new2A_l)], ", ", sep = " "))
        cat("\n")
      }
    }
    iter <- iter + 1
  }

  # normalization

  what <- what_cut
  what_v <- what_cut_v

  # selected variables

  new2A <- which(what != 0)
  new2A_v <- which(what_v != 0)

  # return objects
  object <- list(
    x = x, y = y, loading.x = what, loading.y = what_v,
    variable.x = new2A, variable.y = new2A_v, meanx = meanx,
    normx = normx, meany = meany, normy = normy, mu1 = mu1, mu2 = mu2
  )
  class(object) <- "scca"
  return(object)
}


