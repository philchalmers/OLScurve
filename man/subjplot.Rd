\name{subjplot}
\alias{subjplot}
\alias{subjplot.OLScurve}
\title{Plot individually estimated parameters}
\usage{
  subjplot(object, ...)

  \method{subjplot}{OLScurve} (object, layout = NULL,
    prompt = TRUE, ...)
}
\arguments{
  \item{object}{an object of class \code{OLScurve}}

  \item{layout}{a variable to be passed to \code{xyplot} to
  adjust the graphical layout}

  \item{prompt}{a logical variable indicating whether
  \code{devAskNewPage(ask=TRUE)} should be called}

  \item{...}{additional arguments to be passed}
}
\description{
  A plotting function for displaying the individuals
  trajectories and their modelled functional form. Useful
  for detecting aberrant individual trajectories.
}
\examples{
\dontrun{
data <- t(t(matrix(rnorm(1000),200)) + 1:5)
mod <- OLScurve(~ time, data = data)
subjplot(mod)

##quadratic
data <- t(t(matrix(rnorm(1000),200)) + (0:4)^2)
mod2 <- OLScurve(~ time + I(time^2), data = data)
subjplot(mod2)


}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\keyword{OLS,}
\keyword{growth}

