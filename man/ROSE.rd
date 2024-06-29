\name{ROSE}
\alias{ROSE}
\title{
Function to compute the region of statistical equivalence (ROSE) for an estimate at a given significance and power level.
}
\description{
ROSE computes the region of statistical equivalence for an estimate at a given significance and power level. The region of statistical equivalence is the smallest region wherein one can statistically significantly bound an estimate using a two one-sided tests (TOST) procedure at a given significance level and a pre-specified power level (Fitzgerald 2024). The command produces the bounds of the region of statistical equivalence, its outer bound (that is, the bound furthest from zero), and the magnitude of this outer bound.
}
\usage{
ROSE(estimate, se, alpha = 0.05, power_target = 0.8, df = NA)
}
\arguments{
  \item{estimate}{
The estimate of interest. Numeric scalar.
  }
  \item{se}{
The standard error of the estimate of interest. Numeric scalar, strictly greater than zero.
  }
  }
  \item{alpha}{
Statistical significance level. Numeric scalar, strictly between 0 and 0.5. Defaults to 0.05.
  }
  \item{power_target}{
Power level. Numeric scalar, strictly between 0.5 and 1. Defaults to 0.8.
  }
  \item{df}{
Degrees of freedom. Numeric scalar, strictly greater than zero (if provided). If not provided, the asymptotically approximate region of statistical equivalence is reported. If provided, the exact region of statistical equivalence is reported.
  }
}
\value{
\item{ROSE}{1x2 data.frame showing the bounds of the (1 - alpha) region of statistical equivalence with power level specified by power_target.}
\item{ROSEOB}{Scalar showing the outer bound of the region of statistical equivalence, which is simply the furthest bound from zero.}
\item{ROSEOB_magnitude}{Scalar showing the magnitude (absolute value) of the outer bound of the region of statistical equivalence.}
}
\references{
Fitzgerald, Jack (2024). "The Need for Equivalence Testing in Economics". Institute for Replication Discussion Paper Series No. 125. https://www.econstor.eu/handle/10419/296190.
}
\author{
Jack Fitzgerald, Vrije Universiteit Amsterdam
}
