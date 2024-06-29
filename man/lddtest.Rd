\name{lddtest}
\alias{lddtest}
\title{
Logarithmic density discontinuity equivalence tests for regression discontinuity designs (Fitzgerald 2024)
}
\description{
lddtest performs logarithmic density discontinuity equivalence tests (Fitzgerald 2024) for regression discontinuity designs. lddtest requires that the user specify a running variable, a cutoff point (where treatment kicks in when the running variable crosses that cutoff), and a maximal acceptable ratio between running variable density estimates to the left and right of the cutoff, beneath which the discontinuity would be deemed to be practically equivalent to zero. Once that maximal acceptable ratio is specified, the command utilizes the rdd package coding of DCdensity (Dimmery 2016), itself a McCrary (2008) density test, to estimate the logarithmic density discontinuity in the running variable at the cutoff. Thereafter, lddtest uses a two one-sided tests approach to assess whether the density discontinuity in the running variable at the cutoff is significantly bounded beneath the maximal acceptable ratio.
}
\usage{
lddtest(runvar, cutpoint, epsilon, alpha = 0.05, power = 0.8, bin = NULL, bw = NULL, verbose = FALSE, plot = TRUE)
}
\arguments{
  \item{runvar}{
The running variable. Should be a column of a dataframe.
  }
  \item{cutpoint}{
The cutoff designating the value of the running variable at which treatment kicks in. Scalar.
  }
  \item{epsilon}{
Maximal acceptable ratio between running variable density estimates to the left and right of the cutoff. Scalar, strictly greater than 1.
  }
  \item{alpha}{
Statistical significance level. Numeric scalar, strictly between zero and one. Defaults to 0.05.
  }
  \item{power}{
Power level. Numeric scalar, strictly between zero and one. Defaults to 0.8.
  }
  \item{bin}{
Bin width. See the help guide to DCdensity for more information.
  }
  \item{bw}{
Bandwidth. See the help guide to DCdensity for more information.
  }
  \item{verbose}{
Indicator for whether you wish to see verbose diagnostics displayed. Boolean, defaults to FALSE.
  }
  \item{plot}{
Indicator for whether you wish for DCdensity to produce a plot. Boolean, defaults to FALSE.
  }
}
\value{
\item{bounds}{2x2 data.frame showing the bounds of the (1 - alpha) ECI and the (1 - alpha) ROSE with the pre-specified power level for the estimated logarithmic density discontinuity.}
\item{test}{1x6 data.frame showing the testing bounds, the logarithmic density discontinuity estimate and its standard error, and test statistics from the equivalence test.}
\item{conclusion}{String detailing the conclusion that can be derived from the equivalence test.}
Running this code automatically prints a citation disclaimer.
}
  \references{
Dimmery, Drew (2016). rdd: Regression discontinuity designs. https://doi.org/10.32614/CRAN.package.rdd.
Fitzgerald, Jack (2024). Manipulation Tests in Regression Discontinuity: The Need for Equivalence Testing. Working paper. https://jack-fitzgerald.github.io/files/RDD_equivalence.pdf.
McCrary, Justin (2008). Manipulation of the running variable in the regression discontinuity design: A density test. Journal of Econometrics 142(2), 698-714.
  }
\author{
Jack Fitzgerald, Vrije Universiteit Amsterdam
}
