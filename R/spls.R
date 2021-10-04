##' @title Sparse partial least squares
##'
##' @description This function provides penalty-based sparse partial least squares analysis for single dataset with high dimensions., which aims to have the direction of the first loading.
##'
##' @param x matrix of explanatory variables.
##' @param y matrix of dependent variables.
##' @param mu1 numeric, sparsity penalty parameter.
##' @param eps numeric, the threshold at which the algorithm terminates.
##' @param kappa numeric, 0 < kappa < 0.5 and the parameter reduces the effect of the concave part of objective function.
##' @param scale.x character, "TRUE" or "FALSE", whether or not to scale the variables x. The default is TRUE.
##' @param scale.y character, "TRUE" or "FALSE", whether or not to scale the variables y. The default is TRUE.
##' @param maxstep numeric, maximum iteration steps. The default value is 50.
##' @param trace character, "TRUE" or "FALSE". If TRUE, prints out its screening results of variables.
##'
##' @return An 'spls' object that contains the list of the following items.
##' \itemize{
##' \item{x:}{ data matrix of explanatory variables with centered columns. If scale.x is TRUE, the columns of data matrix are standardized to have mean 0 and standard deviation 1.}
##' \item{y:}{ data matrix of dependent variables with centered columns. If scale.y is TRUE, the columns of data matrix are standardized to have mean 0 and standard deviation 1.}
##' \item{betahat:}{ the estimated regression coefficients.}
##' \item{loading:}{ the estimated first direction vector.}
##' \item{variable:}{ the screening results of variables.}
##' \item{meanx:}{ column mean of the original dataset x.}
##' \item{normx:}{ column standard deviation of the original dataset x.}
##' \item{meany:}{ column mean of the original dataset y.}
##' \item{normy:}{ column standard deviation of the original dataset y.}
##' }
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
##' x <- simData.pls$x[[1]]
##' y <- simData.pls$y[[1]]
##'
##' res <- spls(x = x, y = y, mu1 = 0.25, trace = TRUE)

spls <- function(x, y, mu1, eps = 1e-4, kappa = 0.05, scale.x = TRUE, scale.y = TRUE, maxstep = 50, trace = FALSE) {

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

  # initilize objects

  what <- matrix(0, p, 1)

  # main iteration

  if (trace) {
    cat("The variables that join the set of selected variables at final step:\n")
  }

  # define Z
  Z <- t(x) %*% y
  Z <- Z / (n)

  # fit direction vector
  ro <- function(x, mu, alpha) {
    f <- function(x) mu * (1 > x / (mu * alpha)) * (1 - x / (mu * alpha))
    r <- integrate(f, 0, x)
    return(r)
  }

  ro_d1st <- function(x, mu, alpha) {
    r <- mu * (1 > x / (mu * alpha)) * (1 - x / (mu * alpha))
    return(r)
  }

  # main iteration: optimize u and v iteratively

  # initial value for u(l) (outside the unit circle)

  c <- svd(Z %*% t(Z) , nu = 1)$u
  a <- c
  #u <- U
  #v <- V
  iter <-1
  dis <- 10
  kappa2 <- (1 - kappa) / (1 - 2 * kappa)
  loading_trace <- matrix(0, nrow = p, ncol = maxstep)

  while (dis > eps & iter <= maxstep) {
    # optimize u(l) for fixed v(l)
    c.old <- c
    a.old <- a

    M <- Z %*% t(Z) / q
    h <- function(lambda){
      alpha <- solve(M + lambda * diag(p)) %*% M %*% c
      obj <- t(alpha) %*% alpha - 1 / kappa2^2
      return(obj)
    }
    while (h(1e-4) * h(1e+12) > 0) { # while( h(1e-4) <= 1e+5 )
      {
        M <- 2 * M
        c <- 2 * c
      }
    }
    lambdas <- uniroot(h, c(1e-4, 1e+12))$root
    a <- kappa2 * solve(M + lambdas * diag(p)) %*% M %*% c

    s <- t(a) %*% Z %*% t(Z) / (n)
    ro_d <- mapply(function(j) ro_d1st(s[j], mu1, 6), 1:p)
    fun.c <- function(j) {
      c <- sign(s[j]) * (abs(s[j]) > ro_d[j]) * (abs(s[j]) - ro_d[j])
      return(c)
    }
    c <- mapply(fun.c, 1:p)

    # calculate discrepancy between a & c
    c_norm <- sqrt(sum(c^2))
    c_norm <- ifelse(c_norm == 0, 0.0001, c_norm)
    c <- c / c_norm
    #dis <- max(abs(c - c.old))
    dis <- sqrt(sum((c - c.old)^2)) / sqrt(sum(c.old^2))

    what <- c
    what_cut <- ifelse(abs(what) > 1e-4, what, 0)
    what_cut_norm <- sqrt(sum(what_cut^2))
    what_cut_norm <- ifelse(what_cut_norm == 0, 0.0001, what_cut_norm)
    what_dir <- what_cut / what_cut_norm
    what_cut <- ifelse(abs(what_dir) > 1e-4, what_dir, 0)
    loading_trace[, iter] <- as.numeric(what_cut)

    iter <- iter + 1
    if ( sum(abs(c) <= 1e-4 ) == p ) {
      cat("The value of mu1 is too large");
      break} # exists an l such that c(l)=0
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

  # print out variables that join the active set

  betahat <- matrix(0, nrow = p, ncol = q)
  w <- what
  t_hat <- x %*% w

  if (sum(w == 0) != p) {
    fit <- lm(y ~ t_hat - 1)
    betahat <- matrix(w %*% coef(fit), nrow = p, ncol = q)
  }else {
    betahat <- matrix(0, nrow = p, ncol = q)
  }
  if (!is.null(colnames(x))) {
    rownames(betahat) <- c(1 : p)
  }
  if (q > 1 & !is.null(colnames(y))) {
    colnames(betahat) <- c(1 : q)
  }


  # return objects
  object <- list(
    x = x, y = y, betahat = betahat, loading = what, variable = new2A,
    meanx = meanx, normx = normx, meany = meany, normy = normy, mu1 = mu1
  )
  class(object) <- "spls"
  return(object)
}


