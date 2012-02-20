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
  \item{formula}{a \code{formula} specifying how the
  functional form of \code{time} should be coded. The
  keyword \code{time} must be used, and additive components
  such as powers, square roots, and exponentials can be
  included}

  \item{data}{a data frame in the wide (one subject per
  row) format containting only the time related variables.
  Can be of class \code{matrix} or \code{data.frame}}

  \item{time}{the relative spacing between time points.
  Default is equal spacing}

  \item{x}{an \code{OLScurve} object}

  \item{group}{a \code{factor} grouping variable used to
  parition the results}

  \item{SE}{logical; print a list containing the standard
  errors?}

  \item{digits}{number of digits to round}

  \item{sep}{logical; should the plots be seperated?}

  \item{...}{additional arguments to be passed}
}
\description{
  The \code{OLScurve} provides a simple way of specifying
  ordinary least squares (OLS) growth curve models in R.
  Individual OLS trajectories are fit to each case and an
  OLScurve object is returned which can be passed to
  several graphical and summary function within the
  package.
}
\details{
  As Bollen and Curran (2006) note, there are a variety of
  advantages to using the case-by-case approach for
  estimating trajectory parameters. First of all, OLS
  estimation is intuitively appealing, making it a good
  pedagogical tool for introducing how to model
  trajectories, and illuminates many essential conditions
  and assumptions necessary for LCMs. Second, prediction of
  the parameters for individual trajectory estimates are
  calculated for each case in the sample, which can lead to
  several diagnostics by statistical and graphical means.
  Also, summary statistics can be computed for these
  estimates (which can also be graphically portrayed) and
  if need be these estimates can be analyzed further by
  other statistical frameworks.

  Unfortunately there are also several limitation to OLS
  estimation for LCMs, namely: overall tests of fit are not
  readily available, the structure of the error variances
  must be unrealistically constrained to estimate a pooled
  standard error, the latent factors cannot be regressed
  without error on other exogenous or time-varying
  variables, and analytic significance tests are often not
  readily available (Bollen & Curran, 2006). However, OLS
  estimation may still be useful in the preliminary stages
  of latent curve modeling for (a) selecting appropriate
  functional forms of growth, (b) examining unconditional
  population homogeneity, (c) observing whether the
  relationship between growth factors are linear, and for
  (d) detecting influential outliers (Carrig et al., 2004).
}
\examples{
\dontrun{
##linear
data <- t(t(matrix(rnorm(1000),200)) + 1:5)
mod1 <- OLScurve(~ time, data = data)
mod1
plot(mod1)

##quadratic
data <- t(t(matrix(rnorm(1000),200)) + (0:4)^2)
mod2 <- OLScurve(~ time + I(time^2), data = data)
mod2
plot(mod2)

##sqrt
data <- t(t(matrix(rnorm(1000),200)) + 20*sqrt(5:1))
mod3 <- OLScurve(~ sqrt(time), data = data)
mod3
plot(mod3)

##exponential
data <- t(t(matrix(rnorm(1000,0,5),200)) + exp(0:4))
mod4 <- OLScurve(~ exp(time), data = data)
mod4
plot(mod4)

##combination
data <- t(t(matrix(rnorm(1000),200)) + 20*sqrt(1:5))
mod5 <- OLScurve(~ time + sqrt(time), data = data)
mod5
plot(mod5)

##two group analysis with linear trajectories
data1 <- t(t(matrix(rnorm(500),100)) + 1:5)
data2 <- t(t(matrix(rnorm(500),100)) + 9:5)
data <- rbind(data1,data2)
group <- c(rep('male',100),rep('female',100))

mod <- OLScurve(~ time, data)
print(mod,group)
plot(mod,group)
}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\references{
  Bollen, K. A. & Curran, P. J. (2006). \emph{Latent Curve
  Models: A Structural Equation Perspective}. John Wiley &
  Sons.

  Carrig, M. M., Wirth, R. J., & Curran, P. J. (2004). A
  SAS Macro for Estimating and Visualizing Individual
  Growth Curves. \emph{Structural Equation Modeling, 11},
  132-149.
}
\seealso{
  \code{\link{parplot}}, \code{\link{subjplot}}
}
\keyword{growth}
\keyword{OLS,}

