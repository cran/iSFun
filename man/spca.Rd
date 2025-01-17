% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spca.R
\name{spca}
\alias{spca}
\title{Sparse principal component analysis}
\usage{
spca(x, mu1, eps = 1e-04, scale.x = TRUE, maxstep = 50,
  trace = FALSE)
}
\arguments{
\item{x}{data matrix of explanatory variables.}

\item{mu1}{numeric, sparsity penalty parameter.}

\item{eps}{numeric, the threshold at which the algorithm terminates.}

\item{scale.x}{character, "TRUE" or "FALSE", whether or not to scale the variables x. The default is TRUE.}

\item{maxstep}{numeric, maximum iteration steps. The default value is 50.}

\item{trace}{character, "TRUE" or "FALSE". If TRUE, prints out its screening results of variables.}
}
\value{
An 'spca' object that contains the list of the following items.
\itemize{
\item{x:}{ data matrix of explanatory variables with centered columns. If scale.x is TRUE, the columns of data matrix are standardized to have mean 0 and standard deviation 1.}
\item{eigenvalue:}{ the estimated first eigenvalue.}
\item{eigenvector:}{ the estimated first eigenvector.}
\item{component:}{ the estimated first principal component.}
\item{variable:}{ the screening results of variables.}
\item{meanx:}{ column mean of the original dataset x.}
\item{normx:}{ column standard deviation of the original dataset x.}
}
}
\description{
This function provides penalty-based integrative sparse principal component analysis to obtain the direction of first principal component of a given dataset with high dimensions.
}
\examples{
library(iSFun)
data("simData.pca")
x.spca <- do.call(rbind, simData.pca$x)
res_spca <- spca(x = x.spca, mu1 = 0.08, eps = 1e-3, scale.x = TRUE,
                 maxstep = 50, trace = FALSE)
}
\seealso{
See Also as \code{\link{ispca}}, \code{\link{meta.spca}}.
}
