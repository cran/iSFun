% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scca.R
\name{scca}
\alias{scca}
\title{Sparse canonical correlation analysis}
\usage{
scca(x, y, mu1, mu2, eps = 1e-04, scale.x = TRUE, scale.y = TRUE,
  maxstep = 50, trace = FALSE)
}
\arguments{
\item{x}{data matrix of explanatory variables}

\item{y}{data matrix of dependent variables.}

\item{mu1}{numeric, sparsity penalty parameter for vector u.}

\item{mu2}{numeric, sparsity penalty parameter for vector v.}

\item{eps}{numeric, the threshold at which the algorithm terminates.}

\item{scale.x}{character, "TRUE" or "FALSE", whether or not to scale the variables x. The default is TRUE.}

\item{scale.y}{character, "TRUE" or "FALSE", whether or not to scale the variables y. The default is TRUE.}

\item{maxstep}{numeric, maximum iteration steps. The default value is 50.}

\item{trace}{character, "TRUE" or "FALSE". If TRUE, prints out its screening results of variables.}
}
\value{
An 'scca' object that contains the list of the following items.
\itemize{
\item{x:}{ data matrix of explanatory variables with centered columns. If scale.x is TRUE, the columns of data matrix are standardized to have mean 0 and standard deviation 1.}
\item{y:}{ data matrix of dependent variables with centered columns. If scale.y is TRUE, the columns of data matrix are standardized to have mean 0 and standard deviation 1.}
\item{loading.x:}{ the estimated canonical vector of variables x.}
\item{loading.y:}{ the estimated canonical vector of variables y.}
\item{variable.x:}{ the screening results of variables x.}
\item{variable.y:}{ the screening results of variables y.}
\item{meanx:}{ column mean of the original dataset x.}
\item{normx:}{ column standard deviation of the original dataset x.}
\item{meany:}{ column mean of the original dataset y.}
\item{normy:}{ column standard deviation of the original dataset y.}
}
}
\description{
This function provides penalty-based sparse canonical correlation analysis to get the first pair of canonical vectors.
}
\examples{
library(iSFun)
data("simData.cca")
x.scca <- do.call(rbind, simData.cca$x)
y.scca <- do.call(rbind, simData.cca$y)
res_scca <- scca(x = x.scca, y = y.scca, mu1 = 0.1, mu2 = 0.1, eps = 1e-3,
                 scale.x = TRUE, scale.y = TRUE, maxstep = 50, trace = FALSE)
}
\seealso{
See Also as \code{\link{iscca}}, \code{\link{meta.scca}}.
}
