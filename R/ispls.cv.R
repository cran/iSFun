##' @title Cross-validation for ispls
##'
##' @description Performs K-fold cross validation for the integrative sparse partial least squares over a grid of values for the regularization parameter mu1 and mu2.
##'
##' @param x list of data matrices, L datasets of explanatory variables.
##' @param y list of data matrices, L datasets of dependent variables.
##' @param L numeric, number of datasets.
##' @param K numeric, number of cross-validation folds. Default is 5.
##' @param mu1 numeric, the feasible set of sparsity penalty parameter.
##' @param mu2 numeric, the feasible set of contrasted penalty parameter.
##' @param eps numeric, the threshold at which the algorithm terminates.
##' @param kappa numeric, 0 < kappa < 0.5 and the parameter reduces the effect of the concave part of objective function.
##' @param pen1 character, "homogeneity" or "heterogeneity" type of the sparsity structure. If not specified, the default is homogeneity.
##' @param pen2 character, "magnitude" or "sign" based contrasted penalty. If not specified, the default is magnitude.
##' @param scale.x character, "TRUE" or "FALSE", whether or not to scale the variables x. The default is TRUE.
##' @param scale.y character, "TRUE" or "FALSE", whether or not to scale the variables y. The default is TRUE.
##' @param maxstep numeric, maximum iteration steps. The default value is 50.
##'
##' @return An 'ispls.cv' object that contains the list of the following items.
##' \itemize{
##' \item{x:}{ list of data matrices, L datasets of explanatory variables with centered columns. If scale.x is TRUE, the columns of L datasets are standardized to have mean 0 and standard deviation 1.}
##' \item{y:}{ list of data matrices, L datasets of dependent variables with centered columns. If scale.y is TRUE, the columns of L datasets are standardized to have mean 0 and standard deviation 1.}
##' \item{mu1:}{ the sparsity penalty parameter selected from the feasible set of parameter mu1 provided by users.}
##' \item{mu2:}{ the contrasted penalty parameter selected from the feasible set of parameter mu2 provided by users.}
##' \item{fold:}{ The fold assignments for cross-validation for each observation.}
##' \item{betahat:}{ the estimated regression coefficients with selected tuning parameters mu1 and mu2.}
##' \item{loading:}{ the estimated first direction vector with selected tuning parameters mu1 and mu2.}
##' \item{variable:}{ the screening results of variables x.}
##' \item{meanx:}{ list of numeric vectors, column mean of the original datasets x.}
##' \item{normx:}{ list of numeric vectors, column standard deviation of the original datasets x.}
##' \item{meany:}{ list of numeric vectors, column mean of the original datasets y.}
##' \item{normy:}{ list of numeric vectors, column standard deviation of the original datasets y.}
##' }
##' @references
##' \itemize{
##' \item{Liang W, Ma S, Zhang Q, et al. Integrative sparse partial least squares[J]. Statistics in Medicine, 2021, 40(9): 2239-2256.}
##' }
##' @seealso See Also as \code{\link{ispls}}.
##' @import caret
##' @import irlba
##' @import graphics
##' @import stats
##' @importFrom grDevices rainbow
##' @export
##' @examples
##' \donttest{
##' # Load a list with 3 data sets
##' library(iSFun)
##' data("simData.pls")
##' x <- simData.pls$x
##' y <- simData.pls$y
##' L <- length(x)
##' mu1 <- c(0.04, 0.05)
##' mu2 <- 0.25
##'
##' res_homo_m <- ispls.cv(x = x, y = y, L = L, K = 5, mu1 = mu1, mu2 = mu2, eps = 1e-2,
##'                        kappa = 0.05, pen1 = "homogeneity", pen2 = "magnitude",
##'                        scale.x = TRUE, scale.y = TRUE, maxstep = 50)
##'
##' res_homo_s <- ispls.cv(x = x, y = y, L = L, K = 5, mu1 = mu1, mu2 = mu2, eps = 1e-2,
##'                        kappa = 0.05, pen1 = "homogeneity", pen2 = "sign",
##'                        scale.x = TRUE, scale.y = TRUE, maxstep = 50)
##'
##' res_hete_m <- ispls.cv(x = x, y = y, L = L, K = 5, mu1 = mu1, mu2 = mu2, eps = 1e-2,
##'                        kappa = 0.05, pen1 = "heterogeneity", pen2 = "magnitude",
##'                        scale.x = TRUE, scale.y = TRUE, maxstep = 50)
##'
##' res_hete_s <- ispls.cv(x = x, y = y, L = L, K = 5, mu1 = mu1, mu2 = mu2, eps = 1e-2,
##'                        kappa = 0.05, pen1 = "heterogeneity", pen2 = "sign",
##'                        scale.x = TRUE, scale.y = TRUE, maxstep = 50)
##' }

ispls.cv <- function(x, y, L, K, mu1, mu2, eps = 1e-4, kappa = 0.05, pen1 = "homogeneity", pen2 = "magnitude", scale.x = TRUE, scale.y = TRUE, maxstep = 50) {

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

  # initilize objects
  folds <- lapply(1:L, function(l) createFolds(1:nl[l], K))

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

  c_value_homo <- function(Z, a, c,  p, q, L, mu1, mu2, pen2, nl) {
    # compute s[j,l]
    fun.s <- function(j, l) {
      Z_l <- Z[, ((l - 1) * q + 1):(l * q)]
      a_l <- matrix(a[, l], ncol = 1)
      c_j <- c[j, ]
      s1 <- t(a_l) %*% Z_l %*% t(matrix(Z_l[j, ], nrow = 1)) / (nl[l]^2)
      if (pen2 == "magnitude") s2 <- mu2 * sum(c_j[-l])
      if (pen2 == "sign") {
        s2 <- mu2 * sum(mapply(function(x) {
          r <- x / sqrt(x^2 + 0.5)
          return(r)
        }, c_j[-l])) / sqrt(c_j[l]^2 + 0.5)
      }
      s <- s1 + s2
      return(s)
    }
    result.s <- mapply(fun.s, rep(c(1:p), times = L), rep(c(1:L), each = p))
    s <- matrix(result.s, nrow = p, ncol = L)

    # compute ro'(||c_j||,mu1,a)

    norm_c_j <- apply(c, 1, function(x) {
      return(sqrt(sum(x^2)))
    })
    ro_d <- ro_d1st(norm_c_j, mu1, 6)

    # compute c[j,l]
    fun.c <- function(j, l) {
      s_norm <- sqrt(sum(s[j, ]^2))
      if (pen2 == "magnitude") c <- (s_norm > ro_d[j]) * s[j, l] * (s_norm - ro_d[j]) / ( (1 / q / (nl[l]^2) + mu2 * (L - 1) ) * s_norm )

      if (pen2 == "sign") c <- (s_norm > ro_d[j]) * s[j, l] * (s_norm - ro_d[j]) / ( (1 / q / (nl[l]^2) + mu2 * (L - 1) / (c[j, l]^2 + 0.5) ) * s_norm)

      return(c)
    }
    c <- matrix(mapply(fun.c, rep(c(1:p), times = L), rep(c(1:L), each = p)), nrow = p, ncol = L)
    return(c)
  }

  c_value_hetero <- function(Z, a, c, p, q, L, mu1, mu2, pen2, nl) {
    # compute mu[j,l]
    fun.mu <- function(j, l) {
      c_j <- c[j, ]
      ro_j <- mapply(ro, abs(c_j), mu1, 6)
      s_ro <- sum(as.data.frame(ro_j[1, ]))
      mu_jl <- ro_d1st(s_ro, 1, 1 / 2 * L * 6 * mu1^2) * ro_d1st(abs(c_j[l]), mu1, 6)
      return(mu_jl)
    }
    result.mu <- mapply(fun.mu, rep(c(1:p), times = L), rep(c(1:L), each = p))
    mu <- matrix(result.mu, nrow = p, ncol = L)

    # compute s[j,l]

    fun.s <- function(j, l) {
      Z_l <- Z[, ((l - 1) * q + 1):(l * q)]
      a_l <- matrix(a[, l], ncol = 1)
      c_j <- c[j, ]
      s1 <- t(a_l) %*% Z_l %*% t(matrix(Z_l[j, ], nrow = 1)) / (nl[l]^2)
      if (pen2 == "magnitude") s2 <- mu2 * sum(c_j[-l])
      if (pen2 == "sign") {
        s2 <- mu2 * sum(mapply(function(x) {
          r <- x / sqrt(x^2 + 0.5)
          return(r)
        }, c_j[-l])) / sqrt(c_j[l]^2 + 0.5)
      }
      s <- s1 + s2
      return(s)
    }
    result.s <- mapply(fun.s, rep(c(1:p), times = L), rep(c(1:L), each = p))
    s <- matrix(result.s, nrow = p, ncol = L)

    # compute c[j,l]

    fun.c <- function(j, l) {
      if (pen2 == "magnitude") c <- sign(s[j, l]) * (abs(s[j, l]) > mu[j, l]) * (abs(s[j, l]) - mu[j, l]) / ( 1 / q / (nl[l]^2) + mu2 * (L - 1) )

      if (pen2 == "sign") c <- sign(s[j, l]) * (abs(s[j, l]) > mu[j, l]) * (abs(s[j, l]) - mu[j, l]) / ( 1 / q / (nl[l]^2) + mu2 * (L - 1) / (c[j, l]^2 + 0.5) )

      return(c)
    }
    c <- matrix(mapply(fun.c, rep(c(1:p), times = L), rep(c(1:L), each = p)), nrow = p, ncol = L)
    return(c)
  }

  fun.k <- function(k){
    x.train <- lapply(1:L, function(l) x[[l]][-folds[[l]][[k]],])
    y.train <- lapply(1:L, function(l) y[[l]][-folds[[l]][[k]],])
    x.test  <- lapply(1:L, function(l) x[[l]][folds[[l]][[k]],])
    y.test  <- lapply(1:L, function(l) y[[l]][folds[[l]][[k]],])
    nl.train <- as.numeric(lapply(x.train, nrow))
    nl.test <- as.numeric(lapply(x.test, nrow))
    return(list(x.train = x.train, y.train = y.train, x.test = x.test,
                y.test = y.test, nl.train = nl.train, nl.test = nl.test))
  }
  data.cv = lapply(1:K, fun.k)

  cv.fun <- function(k){
    x.train <- data.cv[[k]]$x.train
    y.train <- data.cv[[k]]$y.train
    x.test  <- data.cv[[k]]$x.test
    y.test  <- data.cv[[k]]$y.test
    nl.train <- data.cv[[k]]$nl.train
    nl.test <- data.cv[[k]]$nl.test

    what <- matrix(0, p, L)

    # define Z

    fun.1 <- function(l) {
      Z_l <- t(x.train[[l]]) %*% y.train[[l]]
      Z_l <- Z_l / nl.train[l]
    }
    #ZZ <- lapply(1:L, fun.1)
    Z <- matrix(mapply(fun.1, c(1:L)), nrow = p)

    fun.2 <- function(l) M_l <- Z[, ((l - 1) * q + 1):(l * q)] %*% t(Z[, ((l - 1) * q + 1):(l * q)]) / q # / (nl.train[l]^2)
    M     <- matrix(mapply(fun.2, c(1:L)), nrow = p)
    dis   <- 10
    iter  <- 1

    # main iteration: optimize c and a iteratively

    kappa2 <- (1 - kappa) / (1 - 2 * kappa)

    # initial value for c(l) (outside the unit circle)

    c <- mapply(function(l) svd(Z[, ((l - 1) * q + 1):(l * q)] %*% t(Z[, ((l - 1) * q + 1):(l * q)]), nu = 1)$u, 1:L)
    a <- c

    while (dis > eps & iter <= maxstep) {
      # optimize a(l) for fixed c(l)
      c.old <- c
      fun.3 <- function(l) {
        h <- function(lambda) {
          alpha <- solve(M[, ((l - 1) * p + 1):(l * p)] + lambda * diag(p)) %*% M[, ((l - 1) * p + 1):(l * p)] %*% c[, l]
          obj <- t(alpha) %*% alpha - 1 / kappa2^2
          return(obj)
        }

        # control size of M_l & c_l if too small

        while (h(1e-4) * h(1e+12) > 0) { # while( h(1e-4) <= 1e+5 )
          {
            M[, ((l - 1) * p + 1):(l * p)] <- 2 * M[, ((l - 1) * p + 1):(l * p)]
            c[, l] <- 2 * c[, l]
          }
        }

        # optimization

        lambdas <- uniroot(h, c(1e-4, 1e+12))$root
        a_l <- kappa2 * solve(M[, ((l - 1) * p + 1):(l * p)] + lambdas * diag(p)) %*% M[, ((l - 1) * p + 1):(l * p)] %*% c[, l]
        return(a_l)
      }

      a <- mapply(fun.3, c(1:L))

      # optimize c(l) for fixed a(l)

      if (pen1 == "homogeneity") c <- c_value_homo(Z, a, c, p, q, L, mu1, mu2, pen2, nl = nl.train)
      if (pen1 == "heterogeneity") c <- c_value_hetero(Z, a, c, p, q, L, mu1, mu2, pen2, nl = nl.train)

      # calculate discrepancy between a & c
      c_norm <- sqrt(colSums(c^2))
      c_norm <- ifelse(c_norm == 0, 0.0001, c_norm)
      c <- t(t(c) / c_norm)
      #dis <- max(abs(c - c.old))
      dis <- max( sqrt(colSums((c - c.old)^2)) / sqrt(colSums(c.old^2)) )

      iter <- iter + 1
      if (sum(apply(c, 2, function(x) sum(abs(x) <= 1e-4) == p)) > 0) {
        cat("The value of mu1 is too large");
        break} # exists an l such that c(l)=0
    }
    what <- c
    what_cut <- ifelse(abs(what) > 1e-4, what, 0)
    what_cut_norm <- sqrt(colSums(what_cut^2))
    what_cut_norm <- ifelse(what_cut_norm == 0, 0.0001, what_cut_norm)
    what_dir <- t(t(what_cut) / what_cut_norm)
    what_cut <- ifelse(abs(what_dir) > 1e-4, what_dir, 0)

    fun.1 <- function(l) {
      Z_l <- ( t(x.test[[l]]) %*% y.test[[l]] ) / nl.test[l]
      rho <- t(what_cut[, l]) %*% Z_l %*% t(Z_l) %*% what_cut[, l]
      return(rho)
    }
    cv.rho <- sum( as.numeric( lapply(1:L, fun.1) ) )
    return(cv.rho)
  }

  mu <- as.matrix(expand.grid(mu1, mu2))
  rho <- c()
  for (loop in 1:nrow(mu)) {
    mu1 = mu[loop, 1]
    mu2 = mu[loop, 2]
    rho[loop] <- mean( mapply(cv.fun, 1:K) )
  }
  index <- which.max(rho)[1]
  mu1.final <- mu[index, 1]
  mu2.final <- mu[index, 2]

  result <- ispls(x, y, L, mu1 = mu1.final, mu2 = mu2.final, eps, kappa, pen1, pen2,
                  scale.x, scale.y, maxstep, trace = FALSE, draw = FALSE)

  # return objects
  object <- list(
    x = x, y = y, mu1 = mu1.final, mu2 = mu2.final, fold = folds,
    betahat = result$betahat, loading = result$loading,
    variable = result$variable,meanx = meanx, normx = normx,
    meany = meany, normy = normy, pen1 = pen1, pen2 = pen2,
    kappa = kappa, loading_trace = result$loading_trace
  )

  class(object) <- "ispls.cv"
  return(object)
}


