##' @title Plot the results of ispls
##'
##' @description Plot the convergence path graph of the first direction vector w in the integrative sparse partial least squares model or show the regression coefficients.
##'
##' @param x list of "ispls", which is the result of command "ispls".
##' @param type character, "path", "loading" or "heatmap" type, if "path", plot the the convergence path graph of vector w in the integrative sparse partial least squares model, if "loading", plot the the first direction vectors, if "heatmap", show the heatmap of regression coefficients among different datasets.
##'
##' @details See details in \code{\link{ispls}}.
##'
##' @return show the convergence path graph of the first direction vector w or the regression coefficients.
##'
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
##' res_homo_m <- ispls(x = x, y = y, L = L, mu1 = 0.05, mu2 = 0.25,
##'                     eps = 5e-2, trace = FALSE, draw = FALSE)
##' ispls.plot(x = res_homo_m, type = "path")
##' ispls.plot(x = res_homo_m, type = "loading")
##' ispls.plot(x = res_homo_m, type = "heatmap")

ispls.plot <- function(x, type) {
  betahat <- x$betahat
  L <- length(betahat)
  p <- nrow(betahat[[1]])
  loading_trace <- x$loading_trace
  iter <- ncol(loading_trace)

  if(type == "path"){
    for (l in 1:L) {
      matplot(t(loading_trace[((l-1)*p+1):(l*p),]), type = 'l',
              main = paste("Dataset ", l, "\n", "Convergence path of elements in vector w"),
              xlab = "Number of iterations", ylab = "Weight")
    }
  }

  if(type == "loading"){
    for (l in 1:L) {
      plot(x = 1:p, y = as.numeric(loading_trace[((l-1)*p+1):(l*p),iter]),
              main = paste("Dataset ", l, "\n", "The first direction vector"),
              xlab = "Dimension", ylab = "Value", pch = 15)
    }
  }

  if(type == "heatmap"){
    rc <- rainbow(nrow(betahat[[1]]), start = 0, end = .3)
    cc <- rainbow(ncol(betahat[[1]]), start = 0, end = .3)
    for (l in 1:L) {
      heatmap(x = t(betahat[[l]]), RowSideColors = cc, ColSideColors = rc,
              Rowv = NA, Colv = NA, scale = "row", xlab = paste("Dataset ", l, ": X"), ylab = "Y",
              main = expression(paste("Heatmap of coefficient ", beta[PLS])))
    }
  }
}


