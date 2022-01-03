##' @title Meta-analytic sparse principal component analysis method in integrative study
##'
##' @description This function provides penalty-based sparse principal component meta-analytic method to handle the multiple datasets with high dimensions generated under similar protocols, which is based on the principle of maximizing the summary statistics S.
##'
##' @param x list of data matrices, L datasets of explanatory variables.
##' @param L numeric, number of datasets.
##' @param mu1 numeric, sparsity penalty parameter.
##' @param eps numeric, the threshold at which the algorithm terminates.
##' @param scale.x character, "TRUE" or "FALSE", whether or not to scale the variables x. The default is TRUE.
##' @param maxstep numeric, maximum iteration steps. The default value is 50.
##' @param trace character, "TRUE" or "FALSE". If TRUE, prints out its screening results of variables.
##'
##' @return A 'meta.spca' object that contains the list of the following items.
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
##' \item{Kim S H, Kang D, Huo Z, et al. Meta-analytic principal component analysis in integrative omics application[J]. Bioinformatics, 2018, 34(8): 1321-1328.}
##' }
##' @seealso See Also as \code{\link{ispca}}, \code{\link{spca}}.
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
##' x <- simData.pca$x
##' L <- length(x)
##'
##' res <- meta.spca(x = x, L = L, mu1 = 0.5, trace = TRUE)

meta.spca <- function(x, L, mu1, eps = 1e-4, scale.x = TRUE, maxstep = 50, trace = FALSE){

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

  Z  <- matrix(0, nrow = p, ncol = p)
  for (l in 1:L) { Z <- Z + (nl[l] - 1) * cov(x[[l]]) }
  Z  <- Z / ( sum(nl) - L )

  u <- irlba(Z, nu = 1, nv = 1)$u
  v <- u
  lambda <- rep(0, p)
  dis.u <- 10
  iter <- 1

  if (trace) {
    cat("The variables that join the set of selected variables at each step:\n")
  }

  while (dis.u > eps & iter <= maxstep ) {
    # optimize u(l) for fixed v(l)
    u.old <- u
    v.old <- v

    s <- t(v) %*% Z #/ sum(nl)
    ro_d <- mapply(function(j) ro_d1st(s[j], mu1, 6), 1:p)
    fun.c <- function(j) {
      c <- sign(s[j]) * (abs(s[j]) > ro_d[j]) * (abs(s[j]) - ro_d[j])
      return(c)
    }
    u <- mapply(fun.c, 1:p)

    if ( sum(abs(u) <= 1e-4 ) == p ) {
      cat("The value of mu1 is too large");
      break}
    v <- Z %*% u / sqrt( sum( (Z %*% u)^2 ) )

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

    if ( sum(abs(u) <= 1e-4 ) == p ) {
      cat("The value of mu1 is too large");
      break}

    if (trace) {
      new2A <- ip[what_cut != 0]
      cat("\n")
      cat(paste("--------------------", "\n"))
      cat(paste("----- Step", iter, " -----\n", sep = " "))
      cat(paste("--------------------", "\n"))
      new2A_l <- new2A
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
    }
    iter <- iter + 1
  }

  # normalization
  what <- what_cut

  # selected variables
  new2A <- ip[what != 0]

  eigenvalue <- mapply(function(l) (t(what)%*%cov(x[[l]]) %*% what), 1:L)
  comp <- lapply(1:L, function(l) x[[l]] %*% what)

  listname <- mapply(function(l) paste("Dataset ", l), 1:L)
  names(meanx) <- listname
  names(normx) <- listname
  names(x) <- listname
  names(comp) <- listname
  names(eigenvalue) <- listname

  # return objects
  object <- list(
    x = x, eigenvalue = eigenvalue, eigenvector = what, component = comp, variable = new2A,
    meanx = meanx, normx = normx, mu1 = mu1
  )
  class(object) <- "meta.spca"

  return(object)
}
