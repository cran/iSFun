##' @title Plot the results of ispca
##'
##' @description Plot the convergence path graph or estimated value of the first eigenvector u in the integrative sparse principal component analysis method.
##'
##' @param x list of "ispca", which is the result of command "ispca".
##' @param type character, "path" or "loading" type, if "path", plot the the convergence path graph of the first eigenvector u in the integrative sparse principal component analysis method, if "loading", plot the first eigenvector.
##'
##' @details See details in \code{\link{ispca}}.
##'
##' @return the convergence path graph or the scatter diagrams of the first eigenvector u.
##'
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
##' res_homo_m <- ispca(x = x, L = L, mu1 = 0.5, mu2 = 0.002, trace = FALSE, draw = FALSE)
##' ispca.plot(x = res_homo_m, type = "path")
##' ispca.plot(x = res_homo_m, type = "loading")

ispca.plot <- function(x, type) {
  loadings <- x$eigenvector
  L <- ncol(loadings)
  p <- nrow(loadings)
  loading_trace <- x$loading_trace
  iter <- ncol(loading_trace)

  if(type == "path"){
    for (l in 1:L) {
      matplot(t(loading_trace[((l-1)*p+1):(l*p),]), type = 'l',
              main = paste("Dataset ", l, "\n", "Convergence path of elements in vector u"),
              xlab = "Number of iterations", ylab = "Weight")
    }
  }

  if(type == "loading"){
    for (l in 1:L) {
      plot(x = 1:p, y = as.numeric(loading_trace[((l-1)*p+1):(l*p),iter]),
           main = paste("The first principal component of dataset", l),
           xlab = "Dimension", ylab = "Value", pch = 15)
    }
  }

}


