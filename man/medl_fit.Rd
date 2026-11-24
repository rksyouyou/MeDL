\name{medl_fit}
\alias{medl_fit}
\title{Test the mediator slope in a random-effects regression}
\usage{medl_fit(joint, B = 0, seed = NULL, bootstrap_null_interval = NULL)}
\description{Fits a regression of joint genetic effect estimates on signed
prediction weights with a free intercept and an additive direct-effect variance.
Tests the null hypothesis of a zero mediator slope.}
\arguments{
\item{joint}{Output of \code{\link{medl_prepare}}, or a list with aligned
numeric vectors \code{beta_hat} and \code{W} and a positive-definite numeric
sampling covariance matrix \code{Sigma}.}
\item{B}{Number of parametric bootstrap draws. Zero skips the bootstrap.}
\item{seed}{Optional seed for the bootstrap random number generator.}
\item{bootstrap_null_interval}{Optional pair of variance bounds bracketing an
interior minimum of the null negative log likelihood. Changes only the null
parameters used to generate bootstrap draws. Requires positive \code{B}.}
}
\value{A list with \code{p_values} (unadjusted P-values named ML, REML-KR, MoM,
MoM-mKH, LR and, for positive \code{B}, LR_boot), \code{fits} (estimates
and diagnostics by method), \code{null} (intercept-only ML fit), and
\code{bootstrap} (draw count, exceedances, bootstrap P-value and generating
parameters, or \code{NULL}).}
\examples{
joint <- do.call(medl_prepare, medl_example)
# A small bootstrap for illustration.
fit <- medl_fit(joint, B = 99, seed = 42)
fit$p_values
fit$fits[["REML-KR"]][c("theta", "gamma", "tau2")]
}
