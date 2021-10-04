##' @title Meta-analytic sparse partial least squares method in integrative study
##'
##' @description This function provides penalty-based sparse canonical correlation meta-analytic method to handle the multiple datasets with high dimensions generated under similar protocols, which is based on the principle of maximizing the summary statistics.
##'
##' @param x list of data matrices, L datasets of explanatory variables.
##' @param y list of data matrices, L datasets of dependent variables.
##' @param L numeric, number of datasets.
##' @param eps numeric, the threshold at which the algorithm terminates.
##' @param mu1 numeric, sparsity penalty parameter.
##' @param kappa numeric, 0 < kappa < 0.5 and the parameter reduces the effect of the concave part of objective function.
##' @param scale.x character, "TRUE" or "FALSE", whether or not to scale the variables x. The default is TRUE.
##' @param scale.y character, "TRUE" or "FALSE", whether or not to scale the variables y. The default is TRUE.
##' @param maxstep numeric, maximum iteration steps. The default value is 50.
##' @param trace character, "TRUE" or "FALSE". If TRUE, prints out its screening results of variables.
##'
##' @return A 'meta.spls' object that contains the list of the following items.
##' \itemize{
##' \item{x:}{ list of data matrices, L datasets of explanatory variables with centered columns. If scale.x is TRUE, the columns of L datasets are standardized to have mean 0 and standard deviation 1.}
##' \item{y:}{ list of data matrices, L datasets of dependent variables with centered columns. If scale.y is TRUE, the columns of L datasets are standardized to have mean 0 and standard deviation 1.}
##' \item{betahat:}{ the estimated regression coefficients.}
##' \item{loading:}{ the estimated first direction vector.}
##' \item{variable:}{ the screening results of variables x.}
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
##' res <- meta.spls(x = x, y = y, L = L, mu1 = 0.03, trace = TRUE)

meta.spls <- function(x, y, L, mu1, eps = 1e-4, kappa = 0.05, scale.x = TRUE, scale.y = TRUE, maxstep = 50, trace = FALSE){

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

  if (trace) {
    cat("The variables that join the set of selected variables at final step:\n")
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

  fun.1 <- function(l) {
    Z_l <- t(x[[l]]) %*% y[[l]]
  }
  ZZ <- lapply(1:L, fun.1)
  Z  <- matrix(0, nrow = p, ncol = q)
  for (l in 1:L) { Z <- Z + ZZ[[l]] }
  Z  <- Z / ( sum(nl) )

  # main iteration: optimize u and v iteratively

  # initial value for u(l) (outside the unit circle)

  c <- svd(Z %*% t(Z), nu = 1)$u
  a <- c
  #u <- U
  #v <- V
  iter <-1
  dis <- 10
  kappa2 <- (1 - kappa) / (1 - 2 * kappa)

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

    s <- t(a) %*% Z %*% t(Z) / sum(nl)
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

    iter <- iter + 1
    if ( sum(abs(c) <= 1e-4 ) == p ) {
      cat("The value of mu1 is too large");
      break} # exists an l such that c(l)=0
  }

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

  # fit y with component t=xw

  betahat <- matrix(0, nrow = p, ncol = q * L)

  fun.fit <- function(l) {
    x_l <- x[[l]]
    w_l <- what
    t_l <- x_l %*% w_l

    if (sum(w_l == 0) != p) {
      y_l <- y[[l]]
      fit_l <- lm(y_l ~ t_l - 1)
      betahat_l <- matrix(w_l %*% coef(fit_l), nrow = p, ncol = q)
    }
    else {
      betahat_l <- matrix(0, nrow = p, ncol = q)
    }
    if (!is.null(colnames(x[[l]]))) {
      rownames(betahat_l) <- c(1 : p)
    }
    if (q > 1 & !is.null(colnames(y[[l]]))) {
      colnames(betahat_l) <- c(1 : q)
    }
    return(betahat_l)
  }

  betahat <- lapply(1:L, fun.fit)

  listname <- mapply(function(l) paste("Dataset ", l), 1:L)
  names(betahat) <- listname
  names(meanx) <- listname
  names(meany) <- listname
  names(normx) <- listname
  names(normy) <- listname
  names(x) <- listname
  names(y) <- listname

  # return objects
  object <- list(
    x = x, y = y, betahat = betahat, loading = what, variable = new2A,
    meanx = meanx, normx = normx, meany = meany, normy = normy, mu1 = mu1
  )
  class(object) <- "meta.spls"

  return(object)
}
