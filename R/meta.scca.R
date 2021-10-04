##' @title Meta-analytic sparse canonical correlation analysis method in integrative study
##'
##' @description This function provides penalty-based sparse canonical correlation meta-analytic method to handle the multiple datasets with high dimensions generated under similar protocols, which is based on the principle of maximizing the summary statistics S.
##'
##' @param x list of data matrices, L datasets of explanatory variables.
##' @param y list of data matrices, L datasets of dependent variables.
##' @param L numeric, number of datasets.
##' @param mu1 numeric, sparsity penalty parameter for vector u.
##' @param mu2 numeric, sparsity penalty parameter for vector v.
##' @param eps numeric, the threshold at which the algorithm terminates.
##' @param scale.x character, "TRUE" or "FALSE", whether or not to scale the variables x. The default is TRUE.
##' @param scale.y character, "TRUE" or "FALSE", whether or not to scale the variables y. The default is TRUE.
##' @param maxstep numeric, maximum iteration steps. The default value is 50.
##' @param trace character, "TRUE" or "FALSE". If TRUE, prints out its screening results of variables.
##'
##' @return A 'meta.scca' object that contains the list of the following items.
##' \itemize{
##' \item{x:}{ list of data matrices, L datasets of explanatory variables with centered columns. If scale.x is TRUE, the columns of L datasets are standardized to have mean 0 and standard deviation 1.}
##' \item{y:}{ list of data matrices, L datasets of dependent variables with centered columns. If scale.y is TRUE, the columns of L datasets are standardized to have mean 0 and standard deviation 1.}
##' \item{loading.x:}{ the estimated canonical vector of variables x.}
##' \item{loading.y:}{ the estimated canonical vector of variables y.}
##' \item{variable.x:}{ the screening results of variables x.}
##' \item{variable.y:}{ the screening results of variables y.}
##' \item{meanx:}{ list of numeric vectors, column mean of the original datasets x.}
##' \item{normx:}{ list of numeric vectors, column standard deviation of the original datasets x.}
##' \item{meany:}{ list of numeric vectors, column mean of the original datasets y.}
##' \item{normy:}{ list of numeric vectors, column standard deviation of the original datasets y.}
##' }
##' @references
##' \itemize{
##' \item{Cichonska A, Rousu J, Marttinen P, et al. metaCCA: summary statistics-based multivariate meta-analysis of genome-wide association studies using canonical correlation analysis[J]. Bioinformatics, 2016, 32(13): 1981-1989.}
##' }
##' @seealso See Also as \code{\link{iscca}}.
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
##' data("simData.cca")
##' x <- simData.cca$x
##' y <- simData.cca$y
##' L <- length(x)
##' mu1 <- 0.08
##' mu2 <- 0.08
##'
##' res <- meta.scca(x = x, y = y, L = L, mu1 = mu1, mu2 = mu2, trace = TRUE)

meta.scca <- function(x, y, L, mu1, mu2, eps = 1e-4, scale.x = TRUE, scale.y = TRUE, maxstep = 50, trace = FALSE){

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
  fun.1 <- function(l) {
    Z_l <- t(x[[l]]) %*% y[[l]]
  }
  ZZ <- lapply(1:L, fun.1)
  Z  <- matrix(0, nrow = p, ncol = q)
  for (l in 1:L) { Z <- Z + ZZ[[l]] }
  Z  <- Z / ( sum(nl) - L )
  irlba.Z <- irlba(Z, nu = 1, nv = 1)
  U <- irlba.Z$u
  V <- irlba.Z$v

  u <- U
  v <- V
  iter <-1
  dis.u <- 10
  dis.v <- 10

  # main iteration: optimize u and v iteratively
  if (trace) {
    cat("The variables that join the set of selected variables at final step:\n")
  }
  while (dis.u > eps & dis.v > eps & iter <= maxstep) {
    u.old <- u
    v.old <- v

    s <- t(v) %*% t(Z) #/ (n)
    ro_d <- mapply(function(j) ro_d1st(s[j], mu1, 6), 1:p)
    fun.c <- function(j) {
      c <- sign(s[j]) * (abs(s[j]) > ro_d[j]) * (abs(s[j]) - ro_d[j])
      return(c)
    }
    u <- mapply(fun.c, 1:p)
    u_norm <- sqrt(sum(u^2))
    u_norm <- ifelse(u_norm == 0, 0.0001, u_norm)
    u <- u / u_norm
    if ( sum(abs(u) <= 1e-4 ) == p ) {
      cat("The value of mu1 is too large");
      break}

    s <- t(u) %*% Z #/ (n)
    ro_d <- mapply(function(j) ro_d1st(s[j], mu2, 6), 1:q)
    fun.c <- function(j) {
      c <- sign(s[j]) * (abs(s[j]) > ro_d[j]) * (abs(s[j]) - ro_d[j])
      return(c)
    }
    v <- mapply(fun.c, 1:q)
    v_norm <- sqrt(sum(v^2))
    v_norm <- ifelse(v_norm == 0, 0.0001, v_norm)
    v <- v / v_norm
    if ( sum(abs(v) <= 1e-4 ) == q ) {
      cat("The value of mu2 is too large");
      break}

    dis.u <- sqrt(sum((u - u.old)^2)) / sqrt(sum(u.old^2))
    dis.v <- sqrt(sum((v - v.old)^2)) / sqrt(sum(v.old^2))

    what <- u
    what_cut <- ifelse(abs(what) > 1e-4, what, 0)
    what_cut_norm <- sqrt(sum(what_cut^2))
    what_cut_norm <- ifelse(what_cut_norm == 0, 0.0001, what_cut_norm)
    what_dir <- what_cut / what_cut_norm
    what_cut <- ifelse(abs(what_dir) > 1e-4, what_dir, 0)

    what_v <- v
    what_cut_v <- ifelse(abs(what_v) > 1e-4, what_v, 0)
    what_cut_norm <- sqrt(sum(what_cut_v^2))
    what_cut_norm <- ifelse(what_cut_norm == 0, 0.0001, what_cut_norm)
    what_dir <- what_cut_v / what_cut_norm
    what_cut_v <- ifelse(abs(what_dir) > 1e-4, what_dir, 0)

    iter <- iter + 1
  }

  # normalization

  what <- what_cut
  what_v <- what_cut_v

  # selected variables

  new2A <- which(what != 0)
  new2A_v <- which(what_v != 0)

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

    if (length(new2A_v) <= 10) {
      cat(paste("DataSet: \n", sep = ""))
      cat(paste("X", new2A_v, ", ", sep = " "))
      cat("\n")
    } else {
      cat(paste("DataSet: \n", sep = ""))
      nlines <- ceiling(length(new2A_v) / 10)
      for (i in 0:(nlines - 2))
      {
        cat(paste("X", new2A_v[(10 * i + 1):(10 * (i + 1))], ", ", sep = " "))
        cat("\n")
      }
      cat(paste("X", new2A_v[(10 * (nlines - 1) + 1):length(new2A_v)], ", ", sep = " "))
      cat("\n")
    }

  }

  listname <- mapply(function(l) paste("Dataset ", l), 1:L)
  names(meanx) <- listname
  names(meany) <- listname
  names(normx) <- listname
  names(normy) <- listname
  names(x) <- listname
  names(y) <- listname

  # return objects
  object <- list(
    x = x, y = y, loading.x = what, loading.y = what_v,
    variable.x = new2A, variable.y = new2A_v, meanx = meanx,
    normx = normx, meany = meany, normy = normy, mu1 = mu1, mu2 = mu2
  )
  class(object) <- "meta.scca"

  return(object)
}
