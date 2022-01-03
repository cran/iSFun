##' @title Integrative sparse principal component analysis
##'
##' @description This function provides a penalty-based integrative sparse principal component analysis method to obtain the direction of first principal component of the multiple datasets with high dimensions generated under similar protocols, which consists of two built-in penalty items for selecting the important variables for users to choose, and two contrasted penalty functions for eliminating the diffierence (magnitude or sign) between estimators within each group.
##'
##' @param x list of data matrices, L datasets of explanatory variables.
##' @param L numeric, number of data sets.
##' @param mu1 numeric, sparsity penalty parameter.
##' @param mu2 numeric, contrasted penalty parameter.
##' @param eps numeric, the threshold at which the algorithm terminates.
##' @param pen1 character, "homogeneity" or "heterogeneity" type of the sparsity structure. If not specified, the default is homogeneity.
##' @param pen2 character, "magnitude" or "sign" based contrasted penalty. If not specified, the default is magnitude.
##' @param scale.x character, "TRUE" or "FALSE", whether or not to scale the variables x. The default is TRUE.
##' @param maxstep numeric, maximum iteration steps. The default value is 50.
##' @param submaxstep numeric, maximum iteration steps in the sub-iterations. The default value is 10.
##' @param trace character, "TRUE" or "FALSE". If TRUE, prints out its screening results of variables.
##' @param draw character, "TRUE" or "FALSE". If TRUE, plot the convergence path of loadings.
##'
##' @return An 'ispca' object that contains the list of the following items.
##' \itemize{
##' \item{x:}{ list of data matrices, L datasets of explanatory variables with centered columns. If scale.x is TRUE, the columns of L datasets are standardized to have mean 0 and standard deviation 1.}
##' \item{eigenvalue:}{ the estimated first eigenvalue.}
##' \item{eigenvector:}{ the estimated first eigenvector.}
##' \item{component:}{ the estimated first component.}
##' \item{variable:}{ the screening results of variables.}
##' \item{meanx:}{ list of numeric vectors, column mean of the original datasets x.}
##' \item{normx:}{ list of numeric vectors, column standard deviation of the original datasets x.}
##' }
##' @references
##' \itemize{
##' \item{Fang K, Fan X, Zhang Q, et al. Integrative sparse principal component analysis[J]. Journal of Multivariate Analysis, 2018, 166: 1-16.}
##' }
##' @seealso See Also as \code{\link{preview.pca}}, \code{\link{ispca.cv}}, \code{\link{meta.spca}}, \code{\link{spca}}.
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
##' prev_pca <- preview.pca(x = x, L = L, scale.x = TRUE)
##' res_homo_m <- ispca(x = x, L = L, mu1 = 0.5, mu2 = 0.002, trace = TRUE, draw = TRUE)
##'
##' \donttest{
##' res_homo_s <- ispca(x = x, L = L, mu1 = 0.5, mu2 = 0.002,
##'                     pen1 = "homogeneity", pen2 = "sign", scale.x = TRUE,
##'                     maxstep = 50, submaxstep = 10, trace = FALSE, draw = FALSE)
##'
##' res_hete_m <- ispca(x = x, L = L, mu1 = 0.1, mu2 = 0.05,
##'                     pen1 = "heterogeneity", pen2 = "magnitude", scale.x = TRUE,
##'                     maxstep = 50, submaxstep = 10, trace = FALSE, draw = FALSE)
##'
##' res_hete_s <- ispca(x = x, L = L, mu1 = 0.1, mu2 = 0.05,
##'                     pen1 = "heterogeneity", pen2 = "sign", scale.x = TRUE,
##'                     maxstep = 50, submaxstep = 10, trace = FALSE, draw = FALSE)
##' }

ispca <- function(x, L, mu1, mu2, eps =1e-4, pen1 = "homogeneity", pen2 = "magnitude", scale.x = TRUE, maxstep = 50, submaxstep = 10, trace = FALSE, draw = FALSE) {

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

  ro <- function(x, mu, alpha) {
    f <- function(x) mu * (1 > x / (mu * alpha)) * (1 - x / (mu * alpha))
    r <- integrate(f, 0, x)
    return(r)
  }

  ro_d1st <- function(x, mu, alpha) {
    r <- mu * (1 > x / (mu * alpha)) * (1 - x / (mu * alpha))
    return(r)
  }

  u_value_homo <- function(u, v, p, L, mu1, mu2, pen2) {
    fun.s <- function(j, l) {
      s1 <- colSums( x[[l]] * matrix(rep(v[[l]], p), ncol=p, nrow = nl[l]) ) / nl[l]
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

    norm_u_j <- apply(u, 1, function(x) {
      return(sqrt(sum(x^2)))
    })
    ro_d <- ro_d1st(norm_u_j, mu1, 6)

    fun.c <- function(j, l) {
      s_norm <- sqrt(sum(s[j, ]^2))
      if (pen2 == "magnitude") c <- nl[l] * (s_norm > ro_d[j]) * s[j, l] * (s_norm - ro_d[j]) / ( (1 + mu2 * nl[l] * (L - 1) ) * s_norm )

      if (pen2 == "sign") c <- nl[l] * (s_norm > ro_d[j]) * s[j, l] * (s_norm - ro_d[j]) / ( (1 + mu2 * nl[l] * (L - 1) / (u[j, l]^2 + 0.5) ) * s_norm)

      return(c)
    }
    c <- matrix(mapply(fun.c, rep(c(1:p), times = L), rep(c(1:L), each = p)), nrow = p, ncol = L)
    return(c)
  }

  u_value_hetero <- function(u, v, p, L, mu1, mu2, pen2) {
    fun.mu <- function(j, l) {
      u_j <- u[j, ]
      ro_j <- mapply(ro, abs(u_j), mu1, 6)
      s_ro <- sum(as.data.frame(ro_j[1, ]))
      mu_jl <- ro_d1st(s_ro, 1, 1 / 2 * L * 6 * mu1^2) * ro_d1st(abs(u_j[l]), mu1, 6)
      return(mu_jl)
    }
    result.mu <- mapply(fun.mu, rep(c(1:p), times = L), rep(c(1:L), each = p))
    mu <- matrix(result.mu, nrow = p, ncol = L)

    fun.s <- function(j, l) {
      s1 <- colSums( x[[l]] * matrix(rep(v[[l]], p), ncol=p, nrow = nl[l]) ) / nl[l]
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

    fun.c <- function(j, l) {
      if (pen2 == "magnitude") c <- nl[l] * sign(s[j, l]) * (abs(s[j, l]) > mu[j, l]) * (abs(s[j, l]) - mu[j, l]) / ( 1 + mu2 *nl[l] * (L - 1) )

      if (pen2 == "sign") c <- nl[l] * sign(s[j, l]) * (abs(s[j, l]) > mu[j, l]) * (abs(s[j, l]) - mu[j, l]) / ( 1 + mu2 * nl[l] * (L - 1) / (u[j, l]^2 + 0.5) )

      return(c)
    }
    c <- matrix(mapply(fun.c, rep(c(1:p), times = L), rep(c(1:L), each = p)), nrow = p, ncol = L)
    return(c)
  }

  # initilize objects
  what <- matrix(0, p, L)
  fun.1 <- function(l) {
    Z_l <- irlba(x[[l]], nu = 1, nv = 1)
    u_l <- Z_l$v * Z_l$d[1]
    return(u_l)
  }
  U <- matrix(mapply(fun.1, 1:L), nrow = p)

  fun.2 <- function(l) {
    Z_l <- irlba(x[[l]], nu =1, nv = 1)
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
  dis.u.iter <- 10
  loading_trace <- matrix(0, nrow = p * L, ncol = maxstep)

  if (trace) {
    cat("The variables that join the set of selected variables at each step:\n")
  }

  while (dis.u.iter > eps & iter <= maxstep) {
    u.iter <- u
    subiter_u <- 1
    dis.u <- 10
    while (dis.u > eps & subiter_u <= submaxstep) {
      u.old <- u
      if (pen1 == "homogeneity") u <- u_value_homo(u, v, p, L, mu1, mu2, pen2)
      if (pen1 == "heterogeneity") u <- u_value_hetero(u, v, p, L, mu1, mu2, pen2)

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
      loading_trace[, iter] <- as.numeric(what_cut)
      if (sum(apply(u.scale, 2, function(x) sum(abs(x) <= 1e-4) == p)) > 0) { break }
      subiter_u <- subiter_u + 1
    }

    if( sum(apply(u, 2, function(x) sum(abs(x) <= 1e-4) == p)) > 0 ){
      cat("Stop! The value of mu1 is too large");
      what_cut <- u
      break }

    v <- lapply(1:L, function(l) x[[l]] %*% u[, l] / sqrt( sum( (x[[l]] %*% u[, l])^2 ) ))

    dis.u.iter <- max( sqrt(colSums((u - u.iter)^2)) / sqrt(colSums(u.iter^2)) )
    if (sum(apply(u, 2, function(x) sum(abs(x) <= 1e-4) == p)) > 0 ) { break }

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

  #print out variables that join the active set
  if(draw){
    for (l in 1:L) {
      matplot(t(loading_trace[((l-1)*p+1):(l*p),]), type = 'l',
              main = paste("Dataset ", l, "\n", "Convergence path of elements in vector u"),
              xlab = "Number of iterations", ylab = "Weight")
    }
  }

  eigenvalue <- mapply(function(l) (t(what[, l])%*%cov(x[[l]]) %*% what[, l]), 1:L)
  comp <- lapply(1:L, function(l) x[[l]] %*% what[, l])

  listname <- mapply(function(l) paste("Dataset ", l), 1:L)
  names(new2A) <- listname
  names(meanx) <- listname
  names(normx) <- listname
  names(x) <- listname
  colnames(what) <- listname
  rownames(what) <- c(1 : p)
  names(comp) <- listname
  names(eigenvalue) <- listname

  # return objects
  object <- list(
    x = x, eigenvalue = eigenvalue, eigenvector = what, component = comp, variable = new2A,
    meanx = meanx, normx = normx, pen1 = pen1, pen2 = pen2,
    mu1 = mu1, mu2 = mu2, loading_trace = loading_trace
  )
  class(object) <- "ispca"
  return(object)
}


