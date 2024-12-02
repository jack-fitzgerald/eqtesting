\name{lddtest}
\alias{lddtest}
\title{
Logarithmic density discontinuity equivalence tests for regression discontinuity designs (Fitzgerald 2024)
}
\description{
lddtest performs logarithmic density discontinuity equivalence tests (Fitzgerald 2024) for regression discontinuity designs. lddtest requires that the user specify a running variable, a cutoff point (where treatment kicks in when the running variable crosses that cutoff), and a maximal acceptable ratio between running variable density estimates to the left and right of the cutoff, beneath which the discontinuity would be deemed to be practically equivalent to zero. Once that maximal acceptable ratio is specified, the command utilizes the rdd package coding of DCdensity (Dimmery 2016), itself a McCrary (2008) density test, to estimate the logarithmic density discontinuity in the running variable at the cutoff. Thereafter, lddtest uses a two one-sided tests approach to assess whether the density discontinuity in the running variable at the cutoff is significantly bounded beneath the maximal acceptable ratio.
}
\usage{
lddtest(runvar, data, cutpoint, epsilon, alpha = 0.05, cluster = "", bootstrap = FALSE, breps = 1000, bin = NULL, bw = NULL, verbose = FALSE, plot = TRUE)
}
\arguments{
  \item{runvar}{
The running variable. Required. Must be a single string denoting the column of the data containing the running variable.
  }
  \item{data}{
The dataset. Required. Must be a data.frame object.
  }
  \item{cutpoint}{
The cutoff designating the value of the running variable at which treatment kicks in. Required. Scalar.
  }
  \item{epsilon}{
Maximal acceptable ratio between running variable density estimates to the left and right of the cutoff. Required. Scalar, strictly greater than 1.
  }
  \item{alpha}{
Statistical significance level. Optional, defaults to 0.05. Scalar, strictly between zero and one.
  }
  \item{cluster}{
Clustering variable. Optional, defaults to an empty string. If specified, must be a single string denoting the column of the data containing the clustering variable.
  }
  \item{bootstrap}{
Indicates whether bootstrap inference should be conducted. Optional, defaults to TRUE. Boolean. Automatically set to TRUE if a cluster variable is specified.
  }
  \item{breps}{
Number of bootstrap replications. Optional, defaults to 1000. Scalar. If you specify breps, then you must also set bootstrap to TRUE. 
  }
  \item{bin}{
Bin width. Optional. See the help guide to DCdensity for more information.
  }
  \item{bw}{
Bandwidth. Optional. See the help guide to DCdensity for more information.
  }
  \item{verbose}{
Indicator for whether you wish to see verbose diagnostics displayed. Optional. Boolean, defaults to FALSE.
  }
  \item{plot}{
Indicator for whether you wish for DCdensity to produce a plot. Optional. Boolean, defaults to TRUE.
  }
}
\value{
\item{test}{1x6 data.frame showing the testing bounds, the logarithmic density discontinuity estimate and its standard error, and test statistics from the equivalence test.}
\item{conclusion}{String detailing the conclusion that can be derived from the equivalence test.}
Running this code automatically prints a citation disclaimer.
}
  \references{
Dimmery, Drew (2016). rdd: Regression discontinuity designs. https://doi.org/10.32614/CRAN.package.rdd.
Fitzgerald, Jack (2024). Manipulation Tests in Regression Discontinuity Design: The Need for Equivalence Testing. Institute for Replication Discussion Paper Series, No. 125. https://hdl.handle.net/10419/300277.
McCrary, Justin (2008). Manipulation of the running variable in the regression discontinuity design: A density test. Journal of Econometrics 142(2), 698-714.
  }
\author{
Jack Fitzgerald, Vrije Universiteit Amsterdam
}
