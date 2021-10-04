##' @title Plot the results of iscca
##'
##' @description Plot the convergence path graph in the integrative sparse canonical correlation analysis method or show the the first pair of canonical vectors.
##'
##' @param x list of "iscca", which is the result of command "iscca".
##' @param type character, "path" or "loading" type, if "path", plot the the convergence path graph of vector u and v in the integrative sparse canonical correlation analysis method, if "loading", show the the first pair of canonical vectors.
##'
##' @details See details in \code{\link{iscca}}.
##'
##' @return the convergence path graph or the scatter diagrams of the first pair of canonical vectors.
##'
##' @import graphics
##' @import stats
##' @importFrom grDevices rainbow
##' @export
##' @examples
##' library(iSFun)
##' data("simData.cca")
##' x <- simData.cca$x
##' y <- simData.cca$y
##' L <- length(x)
##' mu1 <- mu3 <- 0.4
##' mu2 <- mu4 <- 2.5
##'
##' res_homo_m <- iscca(x = x, y = y, L = L, mu1 = mu1, mu2 = mu2, mu3 = mu3,
##'                     mu4 = mu4, eps = 5e-2, maxstep = 100, trace = FALSE, draw = FALSE)
##' iscca.plot(x = res_homo_m, type = "path")
##' iscca.plot(x = res_homo_m, type = "loading")


iscca.plot <- function(x, type) {
  #X <- x$x
  #Y <- x$y
  loading.x <- x$loading.x
  loading.y <- x$loading.y
  L <- ncol(loading.x)
  p <- nrow(loading.x)
  q <- nrow(loading.y)
  loading_trace_u <- x$loading_trace_u
  loading_trace_v <- x$loading_trace_v
  iter_u <- ncol(loading_trace_u)
  iter_v <- ncol(loading_trace_v)

  plot_function <- function(type){
    opar <- par(mfrow = c(1,2))
    on.exit(par(opar))

    if(type == "path"){
      for (l in 1:L) {
        matplot(t(loading_trace_u[((l-1)*p+1):(l*p),]), type = 'l',
                main = paste("Dataset ", l, "\n", "Convergence path of elements in vector u"),
                xlab = "Number of iterations", ylab = "Weight")
        matplot(t(loading_trace_v[((l-1)*q+1):(l*q),]), type = 'l',
                main = paste("Dataset ", l, "\n", "Convergence path of elements in vector v"),
                xlab = "Number of iterations", ylab = "Weight")
      }
    }

    if(type == "loading"){
      for (l in 1:L) {
        plot(x = 1:p, y = as.numeric(loading_trace_u[((l-1)*p+1):(l*p),iter_u]),
             main = paste("The first canonical vector u of dataset", l),
             xlab = "Dimension", ylab = "Value", pch = 15)
        plot(x = 1:p, y = as.numeric(loading_trace_v[((l-1)*p+1):(l*p),iter_v]),
             main = paste("The first canonical vector v of dataset", l),
             xlab = "Dimension", ylab = "Value", pch = 15)
      }
    }
  }
  plot_function(type = type)

}


