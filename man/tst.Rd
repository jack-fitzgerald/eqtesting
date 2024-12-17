\name{tst}
\alias{tst}
\title{
Function to perform the three-sided testing (TST) procedure (Goeman, Solari, and Stijnen 2010) for a given estimate, and compute both the equivalence confidence interval (ECI) and region of statistical equivalence (ROSE) for that estimate.
}
\description{
tst performs the three-sided testing (TST) procedure (Goeman, Solari, and Stijnen 2010) for a given estimate. tst requires that the user specify an estimate, its standard error, and a region of practical equivalence (ROPE), within which the estimate would be deemed to be practically equivalent to zero. Once the ROPE is specified, the TST procedure combines a two one-sided tests (TOST) procedure testing whether the estimate is bounded in the ROPE with two further one-sided tests assessing whether the estimate is bounded outside the ROPE.
}
\usage{
tst(estimate, se, ROPE, df = NA, alpha = 0.05, power = 0.8)
}
\arguments{
  \item{estimate}{
The estimate of interest. Numeric scalar.
  }
  \item{se}{
The standard error of the estimate of interest. Numeric scalar, strictly greater than zero.
  }
  \item{ROPE}{
The ROPE. Can either be a single numeric scalar (interpreted as the length of a symmetric ROPE around zero) or a vector of two different numeric scalars (interpreted as the bounds of the ROPE; in this case one must be >= 0 and one must be <= 0).
  }
  \item{df}{
Degrees of freedom. Numeric scalar, strictly greater than zero (if provided). If not provided, asymptotically approximate ECIs, ROSEs, and TST results are reported. If provided, exact ECIs, ROSEs, and TST results are reported.
  }
  \item{alpha}{
Statistical significance level. Numeric scalar, strictly between zero and one. Defaults to 0.05.
  }
  \item{beta}{
Power level. Numeric scalar, strictly between zero and one. Defaults to 0.8.
  }
}
\value{
\item{bounds}{2x2 data.frame showing the bounds of the (1 - alpha) ECI and the (1 - alpha) ROSE with the pre-specified power level.}
\item{test}{3x5 data.frame showing the ROPE, test statistics, and p-values for each of the three tests in the TST procedure, as well as which of the three tests is relevant.}
\item{conclusion}{String detailing the }
Running this code automatically prints which test was employed (asymptotically approximate or exact) as well as a citation disclaimer.
}
  \references{
Fitzgerald, J. (2024). "The Need for Equivalence Testing in Economics." Institute for Replication Discusion Paper No. 125. https://www.econstor.eu/handle/10419/296190.
Goeman, J. J., Solari, A., and Stijnen, T. (2024). "Three-sided hypothesis testing: Simultaneous testing of superiority, equivalence and inferiority." Statistics in Medicine 29(20), 2117-2125.
  }
\author{
Jack Fitzgerald, Vrije Universiteit Amsterdam and Tinbergen Institute
}
