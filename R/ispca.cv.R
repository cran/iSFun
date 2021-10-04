##' @title Cross-validation for ispca
##'
##' @description Performs K-fold cross validation for the integrative sparse principal component analysis over a grid of values for the regularization parameter mu1 and mu2.
##'
##' @param x list of data matrices, L datasets of explanatory variables.
##' @param L numeric, number of datasets.
##' @param K numeric, number of cross-validation folds. Default is 5.
##' @param mu1 numeric, the feasible set of sparsity penalty parameter.
##' @param mu2 numeric, the feasible set of contrasted penalty parameter.
##' @param eps numeric, the threshold at which the algorithm terminates.
##' @param pen1 character, "homogeneity" or "heterogeneity" type of the sparsity structure. If not specified, the default is homogeneity.
##' @param pen2 character, "magnitude" or "sign" based contrasted penalty. If not specified, the default is magnitude.
##' @param scale.x character, "TRUE" or "FALSE", whether or not to scale the variables x. The default is TRUE.
##' @param maxstep numeric, maximum iteration steps. The default value is 50.
##'
##' @return An 'isca.cv' object that contains the list of the following items.
##' \itemize{
##' \item{x:}{ list of data matrices, L datasets of explanatory variables with centered columns. If scale.x is TRUE, the columns of L datasets are standardized to have mean 0 and standard deviation 1.}
##' \item{y:}{ list of data matrices, L datasets of dependent variables with centered columns. If scale.y is TRUE, the columns of L datasets are standardized to have mean 0 and standard deviation 1.}
##' \item{mu1:}{ the sparsity penalty parameter selected from the feasible set of parameter mu1 provided by users.}
##' \item{mu2:}{ the contrasted penalty parameter selected from the feasible set of parameter mu2 provided by users.}
##' \item{fold:}{ The fold assignments for cross-validation for each observation.}
##' \item{eigenvalue:}{ the estimated first eigenvalue with selected tuning parameters mu1 and mu2.}
##' \item{eigenvector:}{ the estimated first eigenvector with selected tuning parameters mu1 and mu2.}
##' \item{component:}{ the estimated first component with selected tuning parameters mu1 and mu2.}
##' \item{variable:}{ the screening results of variables.}
##' \item{meanx:}{ list of numeric vectors, column mean of the original datasets x.}
##' \item{normx:}{ list of numeric vectors, column standard deviation of the original datasets x.}
##' }
##' @references
##' \itemize{
##' \item{Fang K, Fan X, Zhang Q, et al. Integrative sparse principal component analysis[J]. Journal of Multivariate Analysis, 2018, 166: 1-16.}
##' }
##' @seealso See Also as \code{\link{ispca}}.
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
##' data("simData.pca")
##' x <- simData.pca$x
##' L <- length(x)
##' mu1 <- c(0.3, 0.5)
##' mu2 <- 0.002
##'
##' res_homo_m <- ispca.cv(x = x, L = L, K = 5, mu1 = mu1, mu2 = mu2,
##'                        pen1 = "homogeneity", pen2 = "magnitude", scale.x = TRUE, maxstep = 50)
##'
##' res_homo_s <- ispca.cv(x = x, L = L, K = 5, mu1 = mu1, mu2 = mu2,
##'                        pen1 = "homogeneity", pen2 = "sign", scale.x = TRUE, maxstep = 50)
##'
##' mu1 <- c(0.1, 0.2)
##' mu2 <- 0.05
##' res_hete_m <- ispca.cv(x = x, L = L, K = 5, mu1 = mu1, mu2 = mu2,
##'                        pen1 = "heterogeneity", pen2 = "magnitude", scale.x = TRUE, maxstep = 50)
##'
##' res_hete_s <- ispca.cv(x = x, L = L, K = 5, mu1 = mu1, mu2 = mu2,
##'                        pen1 = "heterogeneity", pen2 = "sign", scale.x = TRUE, maxstep = 50)
##' }

ispca.cv <- function(x, L, K = 5, mu1, mu2, eps = 1e-4, pen1 = "homogeneity", pen2 = "magnitude", scale.x = TRUE, maxstep = 50) {

  if (class(x) != "list") { stop("x should be of list type.") }

  # initialization

  x  <- lapply(x, as.matrix)
  nl <- as.numeric(lapply(x, nrow))
  pl <- as.numeric(lapply(x, ncol))
  p  <- unique(pl)
  if(length(p) > 1){ stop("The dimension of data x should be consistent among different datasets.")}
  ip <- c(1:p)
  n <- mean(unique(nl))

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

  u_value_homo <- function(u, v, p, L, mu1, mu2, pen2, x, nl) {
    # compute s[j,l]
    fun.s <- function(j, l) {
      s1 <- colSums( x[[l]] * matrix(rep(v[[l]], p), ncol = p, nrow = nl[[l]]) ) / nl[l]
      if (pen2 == "magnitude") s2 <- mu2 * sum(u[j, -l])
      if (pen2 == "sign") {
        s2 <- mu2 * sum(mapply(function(x) {
          r <- x / sqrt(x^2 + 0.5)
          return(r)
        }, u[j, -l])) / sqrt(u[j, l]^2 + 0.5)
      }
      s <- s1[j] + s2
      return(s)
    }
    result.s <- mapply(fun.s, rep(c(1:p), times = L), rep(c(1:L), each = p))
    s <- matrix(result.s, nrow = p, ncol = L)

    # compute ro'(||c_j||,mu1,a)

    norm_u_j <- apply(u, 1, function(x) {
      return(sqrt(sum(x^2)))
    })
    ro_d <- ro_d1st(norm_u_j, mu1, 6)

    # compute c[j,l]
    fun.c <- function(j, l) {
      s_norm <- sqrt(sum(s[j, ]^2))
      if (pen2 == "magnitude") c <- nl[l] * (s_norm > ro_d[j]) * s[j, l] * (s_norm - ro_d[j]) / ( (1 + mu2 * nl[l] * (L - 1) ) * s_norm )

      if (pen2 == "sign") c <- nl[l] * (s_norm > ro_d[j]) * s[j, l] * (s_norm - ro_d[j]) / ( (1 + mu2 * nl[l] * (L - 1) / (u[j, l]^2 + 0.5) ) * s_norm)

      return(c)
    }
    c <- matrix(mapply(fun.c, rep(c(1:p), times = L), rep(c(1:L), each = p)), nrow = p, ncol = L)
    return(c)
  }

  u_value_hetero <- function(u, v, p, L, mu1, mu2, pen2, x, nl) {
    # compute mu[j,l]
    fun.mu <- function(j, l) {
      u_j <- u[j, ]
      ro_j <- mapply(ro, abs(u_j), mu1, 6)
      s_ro <- sum(as.data.frame(ro_j[1, ]))
      mu_jl <- ro_d1st(s_ro, 1, 1 / 2 * L * 6 * mu1^2) * ro_d1st(abs(u_j[l]), mu1, 6)
      return(mu_jl)
    }
    result.mu <- mapply(fun.mu, rep(c(1:p), times = L), rep(c(1:L), each = p))
    mu <- matrix(result.mu, nrow = p, ncol = L)

    # compute s[j,l]

    fun.s <- function(j, l) {
      s1 <- colSums( x[[l]] * matrix(rep(v[[l]], p), ncol = p, nrow = nl[[l]]) ) / nl[l]
      if (pen2 == "magnitude") s2 <- mu2 * sum(u[j, -l])
      if (pen2 == "sign") {
        s2 <- mu2 * sum(mapply(function(x) {
          r <- x / sqrt(x^2 + 0.5)
          return(r)
        }, u[j, -l])) / sqrt(u[j, l]^2 + 0.5)
      }
      s <- s1[j] + s2
      return(s)
    }
    result.s <- mapply(fun.s, rep(c(1:p), times = L), rep(c(1:L), each = p))
    s <- matrix(result.s, nrow = p, ncol = L)

    # compute c[j,l]

    fun.c <- function(j, l) {
      if (pen2 == "magnitude") c <- nl[l] * sign(s[j, l]) * (abs(s[j, l]) > mu[j, l]) * (abs(s[j, l]) - mu[j, l]) / ( 1 + mu2 *nl[l] * (L - 1) )

      if (pen2 == "sign") c <- nl[l] * sign(s[j, l]) * (abs(s[j, l]) > mu[j, l]) * (abs(s[j, l]) - mu[j, l]) / ( 1 + mu2 * nl[l] * (L - 1) / (u[j, l]^2 + 0.5) )

      return(c)
    }
    c <- matrix(mapply(fun.c, rep(c(1:p), times = L), rep(c(1:L), each = p)), nrow = p, ncol = L)
    return(c)
  }

  fun.k <- function(k){
    x.train <- lapply(1:L, function(l) x[[l]][-folds[[l]][[k]],])
    #y.train <- lapply(1:L, function(l) y[[l]][-folds[[l]][[k]],])
    x.test  <- lapply(1:L, function(l) x[[l]][folds[[l]][[k]],])
    #y.test  <- lapply(1:L, function(l) y[[l]][folds[[l]][[k]],])
    nl.train <- as.numeric(lapply(x.train, nrow))
    nl.test <- as.numeric(lapply(x.test, nrow))
    return(list(x.train = x.train, x.test = x.test, nl.train = nl.train, nl.test = nl.test))
  }
  data.cv = lapply(1:K, fun.k)

  cv.fun <- function(k){
    x.train <- data.cv[[k]]$x.train
    x.test  <- data.cv[[k]]$x.test
    nl.train <- data.cv[[k]]$nl.train
    nl.test <- data.cv[[k]]$nl.test

    what <- matrix(0, p, L)
    fun.1 <- function(l) {
      Z_l <- irlba(x.train[[l]], nu = 1, nv = 1)
      u_l <- Z_l$v * Z_l$d[1]
      return(u_l)
    }
    U <- matrix(mapply(fun.1, 1:L), nrow = p)
    fun.2 <- function(l) {
      Z_l <- irlba(x.train[[l]], nu =1, nv = 1)
      v_l <- Z_l$u
      return(v_l)
    }
    V <- lapply(1:L, fun.2)
    sgn <- mapply(function(l) sign((U[,1]%*%U[,l])/(sqrt(sum(U[,1]^2))*sqrt(sum(U[,l]^2)))), 2:L )
    for (l in 2:L) {
      U[, l] <- sgn[l-1]*U[, l]
      V[[l]] <- sgn[l-1]*V[[l]]
    }

    u <- U
    v <- V
    iter <-1
    dis.u <- 10
    while (dis.u > eps & iter <= maxstep) {
      u.old <- u
      v.old <- v
      if (pen1 == "homogeneity") u <- u_value_homo(u, v, p, L, mu1, mu2, pen2, x = x.train, nl = nl.train)
      if (pen1 == "heterogeneity") u <- u_value_hetero(u, v, p, L, mu1, mu2, pen2, x = x.train, nl = nl.train)
      if( sum(apply(u, 2, function(x) sum(abs(x) <= 1e-4) == p)) > 0 ){
        cat("The value of mu1 is too large");
        what_cut <- u
        break }
      v <- lapply(1:L, function(l) x.train[[l]] %*% u[, l] / sqrt( sum( (x.train[[l]] %*% u[, l])^2 ) ))

      u_norm <- sqrt(colSums(u^2))
      u_norm <- ifelse(u_norm == 0, 0.0001, u_norm)
      u.scale <- t(t(u) / u_norm)
      dis.u <- max( sqrt(colSums((u - u.old)^2)) / sqrt(colSums(u.old^2)) )

      what <- u.scale
      what_cut <- ifelse(abs(what) > 1e-4, what, 0)
      what_cut_norm <- sqrt(colSums(what_cut^2))
      what_cut_norm <- ifelse(what_cut_norm == 0, 0.0001, what_cut_norm)
      what_dir <- t(t(what_cut) / what_cut_norm)
      what_cut <- ifelse(abs(what_dir) > 1e-4, what_dir, 0)

      iter <- iter + 1
      if (sum(apply(u.scale, 2, function(x) sum(abs(x) <= 1e-4) == p)) > 0) {
        cat("The value of mu1 is too large");
        break}
    }

    v.test <- lapply(1:L, function(l) x.test[[l]] %*% what_cut[, l])
    fun.1 <- function(l) {
      rho <- norm(t(x.test[[l]]) - what_cut[,l] %*% t(v.test[[l]]), "F")^2 / nl.test[l]
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
  index <- which.min(rho)[1]
  mu1.final <- mu[index, 1]
  mu2.final <- mu[index, 2]

  result <- ispca(x, L, mu1 = mu1.final, mu2 = mu2.final, eps, pen1, pen2,
                  scale.x, maxstep, trace = FALSE, draw = FALSE)

  # return objects
  #object <- list( mu1 = mu1.final, mu2 = mu2.final, fold = folds )
  object <- list(
    x = x, mu1 = mu1.final, mu2 = mu2.final, fold = folds,
    eigenvalue = result$eigenvalue, eigenvector = result$eigenvector,
    component = result$component, variable = result$variable, meanx = meanx,
    normx = normx, pen1 = pen1, pen2 = pen2, loading_trace = result$loading_trace
  )

  class(object) <- "ispca.cv"
  return(object)
}


