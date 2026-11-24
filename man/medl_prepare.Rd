\name{medl_prepare}
\alias{medl_prepare}
\title{Reconstruct joint variant effects from marginal summary statistics}
\description{Construct the joint effects and sampling covariance used by
\code{\link{medl_fit}} for a single predicted mediator. Align all inputs to
the same effect alleles and variant order before calling this function.}
\usage{medl_prepare(marginal_beta, marginal_se, weights, cov_g, ridge = 0)}
\arguments{
\item{marginal_beta}{Numeric vector of marginal outcome regression coefficients
(on the log odds scale for binary outcomes).}
\item{marginal_se}{Positive standard errors for the marginal coefficients.}
\item{weights}{Signed prediction weights on the genotype dosage scale, for at
least three variants. The weights must vary to identify the mediator slope.}
\item{cov_g}{Symmetric reference genotype covariance matrix with genotype
variances on the diagonal, using the same dosage scale as the weights.}
\item{ridge}{Nonnegative constant added to the diagonal of the reference
correlation matrix. The resulting matrix must be positive definite.}
}
\value{A list containing the vector of joint effect estimates \code{beta_hat},
the weights \code{W} and the joint sampling covariance matrix \code{Sigma}.}
\examples{
joint <- do.call(medl_prepare, medl_example)
head(joint$beta_hat)
}
