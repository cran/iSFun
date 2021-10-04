##' @title Cross-validation for iscca
##'
##' @description Performs K-fold cross validation for the integrative sparse canonical correlation analysis over a grid of values for the regularization parameter mu1, mu2, mu3 and mu4.
##'
##' @param x list of data matrices, L datasets of explanatory variables.
##' @param y list of data matrices, L datasets of dependent variables.
##' @param L numeric, number of datasets.
##' @param K numeric, number of cross-validation folds. Default is 5.
##' @param mu1 numeric, the feasible set of sparsity penalty parameter for vector u.
##' @param mu2 numeric, the feasible set of contrasted penalty parameter for vector u.
##' @param mu3 numeric, the feasible set of sparsity penalty parameter for vector v.
##' @param mu4 numeric, the feasible set of contrasted penalty parameter for vector v.
##' @param eps numeric, the threshold at which the algorithm terminates.
##' @param pen1 character, "homogeneity" or "heterogeneity" type of the sparsity structure. If not specified, the default is homogeneity.
##' @param pen2 character, "magnitude" or "sign" based contrasted penalty. If not specified, the default is magnitude.
##' @param scale.x character, "TRUE" or "FALSE", whether or not to scale the variables x. The default is TRUE.
##' @param scale.y character, "TRUE" or "FALSE", whether or not to scale the variables y. The default is TRUE.
##' @param maxstep numeric, maximum iteration steps. The default value is 50.
##'
##' @return An 'iscca.cv' object that contains the list of the following items.
##' \itemize{
##' \item{x:}{ list of data matrices, L datasets of explanatory variables with centered columns. If scale.x is TRUE, the columns of L datasets are standardized to have mean 0 and standard deviation 1.}
##' \item{y:}{ list of data matrices, L datasets of dependent variables with centered columns. If scale.y is TRUE, the columns of L datasets are standardized to have mean 0 and standard deviation 1.}
##' \item{mu1:}{ the sparsity penalty parameter selected from the feasible set of parameter mu1 provided by users.}
##' \item{mu2:}{ the contrasted penalty parameter selected from the feasible set of parameter mu2 provided by users.}
##' \item{mu3:}{ the sparsity penalty parameter selected from the feasible set of parameter mu3 provided by users.}
##' \item{mu4:}{ the contrasted penalty parameter selected from the feasible set of parameter mu4 provided by users.}
##' \item{fold:}{ The fold assignments for cross-validation for each observation.}
##' \item{loading.x:}{ the estimated canonical vector of variables x with selected tuning parameters.}
##' \item{loading.y:}{ the estimated canonical vector of variables y with selected tuning parameters.}
##' \item{variable.x:}{ the screening results of variables x.}
##' \item{variable.y:}{ the screening results of variables y.}
##' \item{meanx:}{ list of numeric vectors, column mean of the original datasets x.}
##' \item{normx:}{ list of numeric vectors, column standard deviation of the original datasets x.}
##' \item{meany:}{ list of numeric vectors, column mean of the original datasets y.}
##' \item{normy:}{ list of numeric vectors, column standard deviation of the original datasets y.}
##' }
##' @seealso See Also as \code{\link{iscca}}.
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
##' data("simData.cca")
##' x <- simData.cca$x
##' y <- simData.cca$y
##' L <- length(x)
##' mu1 <- c(0.2, 0.4)
##' mu3 <- 0.4
##' mu2 <- mu4 <- 2.5
##'
##' res_homo_m <- iscca.cv(x = x, y = y, L = L, K = 5, mu1 = mu1, mu2 = mu2, mu3 = mu3,
##'                        mu4 = mu4, eps = 1e-2, pen1 = "homogeneity", pen2 = "magnitude",
##'                        scale.x = TRUE, scale.y = TRUE, maxstep = 100)
##'
##' res_homo_s <- iscca.cv(x = x, y = y, L = L, K = 5, mu1 = mu1, mu2 = mu2, mu3 = mu3,
##'                        mu4 = mu4, eps = 1e-2, pen1 = "homogeneity", pen2 = "sign",
##'                        scale.x = TRUE, scale.y = TRUE, maxstep = 100)
##'
##' mu1 <- mu3 <- c(0.1, 0.3)
##' mu2 <- mu4 <- 2
##' res_hete_m <- iscca.cv(x = x, y = y, L = L, K = 5, mu1 = mu1, mu2 = mu2, mu3 = mu3,
##'                        mu4 = mu4, eps = 1e-2, pen1 = "heterogeneity", pen2 = "magnitude",
##'                        scale.x = TRUE, scale.y = TRUE, maxstep = 100)
##'
##' res_hete_s <- iscca.cv(x = x, y = y, L = L, K = 5, mu1 = mu1, mu2 = mu2, mu3 = mu3,
##'                        mu4 = mu4, eps = 1e-2, pen1 = "heterogeneity", pen2 = "sign",
##'                        scale.x = TRUE, scale.y = TRUE, maxstep = 100)
##' }

iscca.cv <- function(x, y, L, K = 5, mu1, mu2, mu3, mu4, eps = 1e-4, pen1 = "homogeneity", pen2 = "magnitude", scale.x = TRUE, scale.y = TRUE, maxstep = 50) {

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

  u_value_homo <- function(u, v, p, q, L, mu1, mu2, pen2, x, y, nl) {
    # compute s[j,l]
    fun.s <- function(j, l) {
      s1 <- t(v[, l]) %*% t(y[[l]]) %*% t(matrix(x[[l]][, j], nrow = 1)) / (nl[l])
      if (pen2 == "magnitude") s2 <- mu2 * sum(u[j, -l])
      if (pen2 == "sign") {
        s2 <- mu2 * sum(mapply(function(x) {
          r <- x / sqrt(x^2 + 0.5)
          return(r)
        }, u[j, -l])) / sqrt(u[j, l]^2 + 0.5)
      }
      s <- s1 + s2
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
      if (pen2 == "magnitude") c <- (s_norm > ro_d[j]) * s[j, l] * (s_norm - ro_d[j]) / ( ( mu2 * (L - 1) ) * s_norm )

      if (pen2 == "sign") c <- (s_norm > ro_d[j]) * s[j, l] * (s_norm - ro_d[j]) / ( ( mu2 * (L - 1) / (u[j, l]^2 + 0.5) ) * s_norm)

      return(c)
    }
    c <- matrix(mapply(fun.c, rep(c(1:p), times = L), rep(c(1:L), each = p)), nrow = p, ncol = L)
    return(c)
  }

  u_value_hetero <- function(u, v, p, q, L, mu1, mu2, pen2, x, y, nl) {
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
      s1 <- t(v[, l]) %*% t(y[[l]]) %*% t(matrix(x[[l]][, j], nrow = 1)) / (nl[l])
      if (pen2 == "magnitude") s2 <- mu2 * sum(u[j, -l])
      if (pen2 == "sign") {
        s2 <- mu2 * sum(mapply(function(x) {
          r <- x / sqrt(x^2 + 0.5)
          return(r)
        }, u[j, -l])) / sqrt(u[j, l]^2 + 0.5)
      }
      s <- s1 + s2
      return(s)
    }
    result.s <- mapply(fun.s, rep(c(1:p), times = L), rep(c(1:L), each = p))
    s <- matrix(result.s, nrow = p, ncol = L)

    # compute c[j,l]

    fun.c <- function(j, l) {
      if (pen2 == "magnitude") c <- sign(s[j, l]) * (abs(s[j, l]) > mu[j, l]) * (abs(s[j, l]) - mu[j, l]) / ( mu2 * (L - 1) )

      if (pen2 == "sign") c <- sign(s[j, l]) * (abs(s[j, l]) > mu[j, l]) * (abs(s[j, l]) - mu[j, l]) / ( mu2 * (L - 1) / (u[j, l]^2 + 0.5) )

      return(c)
    }
    c <- matrix(mapply(fun.c, rep(c(1:p), times = L), rep(c(1:L), each = p)), nrow = p, ncol = L)
    return(c)
  }

  v_value_homo <- function(u, v, p, q, L, mu1, mu2, pen2, x, y, nl) {
    # compute s[j,l]
    fun.s <- function(j, l) {
      s1 <- t(u[, l]) %*% t(x[[l]]) %*% t(matrix(y[[l]][, j], nrow = 1)) / (nl[l])
      if (pen2 == "magnitude") s2 <- mu2 * sum(v[j, -l])
      if (pen2 == "sign") {
        s2 <- mu2 * sum(mapply(function(x) {
          r <- x / sqrt(x^2 + 0.5)
          return(r)
        }, v[j, -l])) / sqrt(v[j, l]^2 + 0.5)
      }
      s <- s1 + s2
      return(s)
    }
    result.s <- mapply(fun.s, rep(c(1:q), times = L), rep(c(1:L), each = q))
    s <- matrix(result.s, nrow = q, ncol = L)

    # compute ro'(||c_j||,mu1,a)

    norm_v_j <- apply(v, 1, function(x) {
      return(sqrt(sum(x^2)))
    })
    ro_d <- ro_d1st(norm_v_j, mu1, 6)

    # compute c[j,l]
    fun.c <- function(j, l) {
      s_norm <- sqrt(sum(s[j, ]^2))
      if (pen2 == "magnitude") c <- (s_norm > ro_d[j]) * s[j, l] * (s_norm - ro_d[j]) / ( ( mu2 * (L - 1) ) * s_norm )

      if (pen2 == "sign") c <- (s_norm > ro_d[j]) * s[j, l] * (s_norm - ro_d[j]) / ( ( mu2 * (L - 1) / (v[j, l]^2 + 0.5) ) * s_norm)

      return(c)
    }
    c <- matrix(mapply(fun.c, rep(c(1:q), times = L), rep(c(1:L), each = q)), nrow = q, ncol = L)
    return(c)
  }

  v_value_hetero <- function(u, v, p, q, L, mu1, mu2, pen2, x, y, nl) {
    # compute mu[j,l]
    fun.mu <- function(j, l) {
      v_j <- v[j, ]
      ro_j <- mapply(ro, abs(v_j), mu1, 6)
      s_ro <- sum(as.data.frame(ro_j[1, ]))
      mu_jl <- ro_d1st(s_ro, 1, 1 / 2 * L * 6 * mu1^2) * ro_d1st(abs(v_j[l]), mu1, 6)
      return(mu_jl)
    }
    result.mu <- mapply(fun.mu, rep(c(1:q), times = L), rep(c(1:L), each = q))
    mu <- matrix(result.mu, nrow = q, ncol = L)

    # compute s[j,l]

    fun.s <- function(j, l) {
      s1 <- t(u[, l]) %*% t(x[[l]]) %*% t(matrix(y[[l]][, j], nrow = 1)) / (nl[l])
      if (pen2 == "magnitude") s2 <- mu2 * sum(v[j, -l])
      if (pen2 == "sign") {
        s2 <- mu2 * sum(mapply(function(x) {
          r <- x / sqrt(x^2 + 0.5)
          return(r)
        }, v[j, -l])) / sqrt(v[j, l]^2 + 0.5)
      }
      s <- s1 + s2
      return(s)
    }
    result.s <- mapply(fun.s, rep(c(1:q), times = L), rep(c(1:L), each = q))
    s <- matrix(result.s, nrow = q, ncol = L)

    # compute c[j,l]

    fun.c <- function(j, l) {
      if (pen2 == "magnitude") c <- sign(s[j, l]) * (abs(s[j, l]) > mu[j, l]) * (abs(s[j, l]) - mu[j, l]) / ( mu2 * (L - 1) )

      if (pen2 == "sign") c <- sign(s[j, l]) * (abs(s[j, l]) > mu[j, l]) * (abs(s[j, l]) - mu[j, l]) / ( mu2 * (L - 1) / (u[j, l]^2 + 0.5) )

      return(c)
    }
    c <- matrix(mapply(fun.c, rep(c(1:q), times = L), rep(c(1:L), each = q)), nrow = q, ncol = L)
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

    fun.1 <- function(l) {
      Z_l <- irlba(t(x.train[[l]]) %*% y.train[[l]], nu = 1, nv = 1)
      u_l <- Z_l$u
      return(u_l)
    }
    U <- matrix(mapply(fun.1, 1:L), nrow = p)

    fun.2 <- function(l) {
      Z_l <- irlba(t(x.train[[l]]) %*% y.train[[l]], nu = 1, nv = 1)
      v_l <- Z_l$v
      return(v_l)
    }
    V <- mapply(fun.2, 1:L)

    sgn <- mapply(function(l) sign((U[,1]%*%U[,l])/(sqrt(sum(U[,1]^2))*sqrt(sum(U[,l]^2)))), 2:L )
    for (l in 2:L) {
      U[, l] <- sgn[l-1]*U[, l]
      V[, l] <- sgn[l-1]*V[, l]
    }

    u <- U
    v <- V
    iter <- 1
    dis.u <- 10

    while (dis.u > eps & iter <= maxstep) {
      u.old <- u
      if (pen1 == "homogeneity") u <- u_value_homo(u, v, p, q, L, mu1, mu2, pen2, x = x.train, y = y.train, nl = nl.train)
      if (pen1 == "heterogeneity") u <- u_value_hetero(u, v, p, q, L, mu1, mu2, pen2, x = x.train, y = y.train, nl = nl.train)
      u_norm <- sqrt(colSums(u^2))
      u_norm <- ifelse(u_norm == 0, 0.0001, u_norm)
      u <- t(t(u) / u_norm)
      dis.u <- max( sqrt(colSums((u - u.old)^2)) / sqrt(colSums(u.old^2)) )
      if (sum(apply(u, 2, function(x) sum(abs(x) <= 1e-4) == p)) > 0) {
        cat("The value of mu1 is too large");
        break}
      iter <- iter + 1
    }
    what <- u
    what_cut <- ifelse(abs(what) > 1e-4, what, 0)
    what_cut_norm <- sqrt(colSums(what_cut^2))
    what_cut_norm <- ifelse(what_cut_norm == 0, 0.0001, what_cut_norm)
    what_dir <- t(t(what_cut) / what_cut_norm)
    what_cut <- ifelse(abs(what_dir) > 1e-4, what_dir, 0)

    u <- U
    v <- V
    iter <- 1
    dis.v <- 10

    while (dis.v > eps & iter <= maxstep) {
      v.old <- v
      if (pen1 == "homogeneity") v <- v_value_homo(u, v, p, q, L, mu1 = mu3, mu2 = mu4, pen2, x = x.train, y = y.train, nl = nl.train)
      if (pen1 == "heterogeneity") v <- v_value_hetero(u, v, p, q, L, mu1 = mu3, mu2 = mu4, pen2, x = x.train, y = y.train, nl = nl.train)
      v_norm <- sqrt(colSums(v^2))
      v_norm <- ifelse(v_norm == 0, 0.0001, v_norm)
      v <- t(t(v) / v_norm)
      dis.v <- max( sqrt(colSums((v - v.old)^2)) / sqrt(colSums(v.old^2)) )
      if (sum(apply(v, 2, function(x) sum(abs(x) <= 1e-4) == q)) > 0) {
        cat("The value of mu1 is too large");
        break}
      iter <- iter + 1
    }
    what_v <- v
    what_cut_v <- ifelse(abs(what_v) > 1e-4, what_v, 0)
    what_cut_norm <- sqrt(colSums(what_cut_v^2))
    what_cut_norm <- ifelse(what_cut_norm == 0, 0.0001, what_cut_norm)
    what_dir <- t(t(what_cut_v) / what_cut_norm)
    what_cut_v <- ifelse(abs(what_dir) > 1e-4, what_dir, 0)

    fun.1 <- function(l) {
      rho <- ( t(what_cut[, l]) %*% t(x.test[[l]]) %*% y.test[[l]] %*% what_cut_v[, l] ) / ( nl.test[l] )
      return(rho)
    }
    cv.rho <- sum( as.numeric( lapply(1:L, fun.1) ) )
    return(cv.rho)
  }

  mu <- as.matrix(expand.grid(mu1, mu2, mu3, mu4))
  rho <- c()
  for (loop in 1:nrow(mu)) {
    mu1 = mu[loop, 1]
    mu2 = mu[loop, 2]
    mu3 = mu[loop, 3]
    mu4 = mu[loop, 4]
    rho[loop] <- mean( mapply(cv.fun, 1:K) )
  }
  index <- which.max(rho)[1]
  mu1.final <- mu[index, 1]
  mu2.final <- mu[index, 2]
  mu3.final <- mu[index, 3]
  mu4.final <- mu[index, 4]

  # return objects
  result <- iscca(x, y, L, mu1 = mu1.final, mu2 = mu2.final, mu3 = mu3.final,
                  mu4 = mu4.final, eps, pen1, pen2, scale.x, scale.y,
                  maxstep, trace = FALSE, draw = FALSE)

  object <- list(
    x = x, y = y, mu1 = mu1.final, mu2 = mu2.final, mu3 = mu3.final,
    mu4 = mu4.final, fold = folds, loading.x = result$loading.x, loading.y = result$loading.y,
    variable.x = result$variable.x, variable.y = result$variable.y,  meanx = meanx,
    normx = normx, meany = meany, normy = normy, pen1 = pen1, pen2 = pen2,
    loading_trace_u = result$loading_trace_u, loading_trace_v = result$loading_trace_v
  )

  class(object) <- "iscca.cv"
  return(object)
}


