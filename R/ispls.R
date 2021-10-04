##' @title Integrative sparse partial least squares
##'
##' @description This function provides a penalty-based integrative sparse partial least squares method to handle the multiple datasets with high dimensions generated under similar protocols, which consists of two built-in penalty items for selecting the important variables for users to choose, and two contrasted penalty functions for eliminating the diffierence (magnitude or sign) between estimators within each group.
##'
##' @param x list of data matrices, L datasets of explanatory variables.
##' @param y list of data matrices, L datasets of dependent variables.
##' @param L numeric, number of datasets.
##' @param mu1 numeric, sparsity penalty parameter.
##' @param mu2 numeric, contrasted penalty parameter.
##' @param eps numeric, the threshold at which the algorithm terminates.
##' @param kappa numeric, 0 < kappa < 0.5 and the parameter reduces the effect of the concave part of objective function.
##' @param pen1 character, "homogeneity" or "heterogeneity" type of the sparsity structure. If not specified, the default is homogeneity.
##' @param pen2 character, "magnitude" or "sign" based contrasted penalty. If not specified, the default is magnitude.
##' @param scale.x character, "TRUE" or "FALSE", whether or not to scale the variables x. The default is TRUE.
##' @param scale.y character, "TRUE" or "FALSE", whether or not to scale the variables y. The default is TRUE.
##' @param maxstep numeric, maximum iteration steps. The default value is 50.
##' @param trace character, "TRUE" or "FALSE". If TRUE, prints out its screening results of variables.
##' @param draw character, "TRUE" or "FALSE". If TRUE, plot the convergence path of loadings.
##'
##' @return An 'ispls' object that contains the list of the following items.
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
##' @references
##' \itemize{
##' \item{Liang W, Ma S, Zhang Q, et al. Integrative sparse partial least squares[J]. Statistics in Medicine, 2021, 40(9): 2239-2256.}
##' }
##' @seealso See Also as \code{\link{preview.pls}}, \code{\link{ispca}}, \code{\link{iscca}}.
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
##' data("simData.pls")
##' x <- simData.pls$x
##' y <- simData.pls$y
##' L <- length(x)
##'
##' prev_pls <- preview.pls(x, y, L, scale.x = TRUE, scale.y = TRUE)
##' res_homo_m <- ispls(x = x, y = y, L = L, mu1 = 0.05, mu2 = 0.25,
##'                     eps = 5e-2, trace = TRUE, draw = TRUE)
##'
##' \donttest{
##' res_homo_s <- ispls(x = x, y = y, L = L, mu1 = 0.05, mu2 = 0.25,
##'                     eps = 5e-2, kappa = 0.05, pen1 = "homogeneity",
##'                     pen2 = "sign", scale.x = TRUE, scale.y = TRUE,
##'                     maxstep = 50, trace = FALSE, draw = FALSE)
##'
##' res_hete_m <- ispls(x = x, y = y, L = L, mu1 = 0.05, mu2 = 0.25,
##'                     eps = 5e-2, kappa = 0.05, pen1 = "heterogeneity",
##'                     pen2 = "magnitude", scale.x = TRUE, scale.y = TRUE,
##'                     maxstep = 50, trace = FALSE, draw = FALSE)
##'
##' res_hete_s <- ispls(x = x, y = y, L = L, mu1 = 0.05, mu2 = 0.25,
##'                     eps = 5e-2, kappa = 0.05, pen1 = "heterogeneity",
##'                     pen2 = "sign", scale.x = TRUE, scale.y = TRUE,
##'                     maxstep = 50, trace = FALSE, draw = FALSE)
##' }

ispls <- function(x, y, L, mu1, mu2, eps = 1e-4, kappa = 0.05, pen1 = "homogeneity", pen2 = "magnitude", scale.x = TRUE, scale.y = TRUE, maxstep = 50, trace = FALSE, draw = FALSE) {

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
  #n <- mean(unique(nl))

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

  what <- matrix(0, p, L)

  # main iteration

  if (trace) {
    cat("The variables that join the set of selected variables at each step:\n")
  }

  # define Z

  fun.1 <- function(l) {
    Z_l <- t(x[[l]]) %*% y[[l]]
    Z_l <- Z_l / nl[l]
  }
  #ZZ <- lapply(1:L, fun.1)
  Z <- matrix(mapply(fun.1, c(1:L)), nrow = p)

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

  c_value_homo <- function(Z, a, c,  p, q, L, mu1, mu2, pen2) {
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

  c_value_hetero <- function(Z, a, c, p, q, L, mu1, mu2, pen2) {
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

  fun.2 <- function(l) M_l <- Z[, ((l - 1) * q + 1):(l * q)] %*% t(Z[, ((l - 1) * q + 1):(l * q)]) / q #/ (nl[l]^2)
  M <- matrix(mapply(fun.2, c(1:L)), nrow = p)
  dis <- 10
  iter <- 1

  # main iteration: optimize c and a iteratively

  kappa2 <- (1 - kappa) / (1 - 2 * kappa)

  # initial value for c(l) (outside the unit circle)

  c <- mapply(function(l) svd(Z[, ((l - 1) * q + 1):(l * q)] %*% t(Z[, ((l - 1) * q + 1):(l * q)] ), nu = 1)$u, 1:L)#/ (nl[l]^2)
  a <- c
  loading_trace <- matrix(0, nrow = p * L, ncol = maxstep)

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

    if (pen1 == "homogeneity") c <- c_value_homo(Z, a, c, p, q, L, mu1, mu2, pen2)
    if (pen1 == "heterogeneity") c <- c_value_hetero(Z, a, c, p, q, L, mu1, mu2, pen2)

    # calculate discrepancy between a & c
    c_norm <- sqrt(colSums(c^2))
    c_norm <- ifelse(c_norm == 0, 0.0001, c_norm)
    c <- t(t(c) / c_norm)
    #dis <- max(abs(c - c.old))
    dis <- max( sqrt(colSums((c - c.old)^2)) / sqrt(colSums(c.old^2)) )


    what <- c
    what_cut <- ifelse(abs(what) > 1e-4, what, 0)
    what_cut_norm <- sqrt(colSums(what_cut^2))
    what_cut_norm <- ifelse(what_cut_norm == 0, 0.0001, what_cut_norm)
    what_dir <- t(t(what_cut) / what_cut_norm)
    what_cut <- ifelse(abs(what_dir) > 1e-4, what_dir, 0)
    loading_trace[, iter] <- as.numeric(what_cut)

    if (trace) {
      # selected variables
      new2A <- list()
      new2A <- mapply(function(l) {
        A <- ip[what_cut[, l] != 0]
        return(A)
      }, c(1:L), SIMPLIFY = FALSE)
      cat("\n")
      cat(paste("--------------------", "\n"))
      cat(paste("----- Step", iter, " -----\n", sep = " "))
      cat(paste("--------------------", "\n"))
      for (l in 1:L)
      {
        new2A_l <- new2A[[l]]
        if (length(new2A_l) <= 10) {
          cat(paste("DataSet ", l, ":\n", sep = ""))
          cat(paste("X", new2A_l, ", ", sep = " "))
          cat("\n")
        } else {
          cat(paste("DataSet ", l, ":\n", sep = ""))
          nlines <- ceiling(length(new2A_l) / 10)
          for (i in 0:(nlines - 2))
          {
            cat(paste("X", new2A_l[(10 * i + 1):(10 * (i + 1))], ", ", sep = " "))
            cat("\n")
          }
          cat(paste("X", new2A_l[(10 * (nlines - 1) + 1):length(new2A_l)], ", ", sep = " "))
          cat("\n")
        }
      }
    }

    iter <- iter + 1
    if (sum(apply(c, 2, function(x) sum(abs(x) <= 1e-4) == p)) > 0) {
      cat("The value of mu1 is too large");
      break} # exists an l such that c(l)=0
  }
  loading_trace <- loading_trace[,1:(iter-1)]

  # normalization

  what <- what_cut

  # selected variables

  new2A <- list()
  new2A <- mapply(function(l) {
    A <- ip[what[, l] != 0]
    return(A)
  }, c(1:L), SIMPLIFY = FALSE)

  # fit y with component t=xw

  betahat <- matrix(0, nrow = p, ncol = q * L)

  fun.fit <- function(l) {
    x_l <- x[[l]]
    w_l <- what[, l]
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
  names(new2A) <- listname
  names(meanx) <- listname
  names(meany) <- listname
  names(normx) <- listname
  names(normy) <- listname
  names(x) <- listname
  names(y) <- listname
  colnames(what) <- listname
  rownames(what) <- c(1 : p)

  if(draw){
    rc <- rainbow(nrow(betahat[[1]]), start = 0, end = .3)
    cc <- rainbow(ncol(betahat[[1]]), start = 0, end = .3)
    for (l in 1:L) {
      matplot(t(loading_trace[((l-1)*p+1):(l*p),]), type = 'l',
              main = paste("Dataset ", l, "\n", "Convergence path of elements in vector w"),
              xlab = "Number of iterations", ylab = "Weight")
    }
  }

  # return objects
  object <- list(
    x = x, y = y, betahat = betahat, loading = what, variable = new2A,
    meanx = meanx, normx = normx, meany = meany, normy = normy,
    pen1 = pen1, pen2 = pen2, mu1 = mu1, mu2 = mu2, kappa = kappa,
    loading_trace = loading_trace
  )
  class(object) <- "ispls"
  return(object)
}


