\name{medl_example}
\alias{medl_example}
\docType{data}
\title{Simulated marginal summary statistics for one mediator}
\description{Summary-level inputs for a continuous outcome and 20 correlated
synthetic variants. All values were generated for illustration.}
\usage{medl_example}
\format{A named list of four inputs to \code{\link{medl_prepare}}:
\describe{
\item{marginal_beta}{Named vector of 20 marginal outcome effect estimates.}
\item{marginal_se}{Named vector of 20 marginal standard errors.}
\item{weights}{Named vector of 20 signed mediator prediction weights.}
\item{cov_g}{A 20 by 20 reference genotype covariance matrix whose rows and
columns follow the order of the three vectors.}
}}
\details{The reference correlations decay as 0.45 raised to the distance
between variant indices. Joint effects have intercept 0.01, mediator slope
0.15 and independent normal direct-effect deviations with standard deviation
0.025. Marginal estimates include correlated normal sampling error.}
\examples{
names(medl_example)
dim(medl_example$cov_g)
}
\keyword{datasets}
