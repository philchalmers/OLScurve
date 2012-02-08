\name{OLScurve}
\alias{OLScurve}
\alias{plot.OLScurve}
\alias{print.OLScurve}
\title{Ordinary least squares growth curve trajectories}
\usage{
  OLScurve(formula, data, time = 0:(ncol(data) - 1), ...)

  \method{print}{OLScurve} (x, group = NULL, SE = TRUE,
    digits = 3, ...)

  \method{plot}{OLScurve} (x, group = NULL, sep = FALSE,
    ...)
}
\arguments{
  \item{formula}{DESCRIPTION}

  \item{data}{DESCRIPTION}

  \item{time}{DESCRIPTION}

  \item{x}{DESCRIPTION}

  \item{group}{DESCRIPTION}

  \item{SE}{DESCRIPTION}

  \item{digits}{DESCRIPTION}

  \item{sep}{DESCRIPTION}

  \item{...}{DESCRIPTION}
}
\description{
  DESCRIPTION of \code{OLScurve}
}
\examples{
\dontrun{
data <- t(t(matrix(rnorm(1000),200)) + 1:5)
mod <- OLScurve(~ time, data = data)
mod
plot(mod)

}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\keyword{growth}
\keyword{OLS,}

