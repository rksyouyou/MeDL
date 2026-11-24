# MeDL

MeDL tests whether a genetically predicted mediator is associated with an outcome. It fits a random-effects regression of the estimated joint variant effects on the signed mediator prediction weights, with a free intercept for the mean direct effect. It uses the full sampling covariance of the joint effects.

## Install from GitHub

In R, run:

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("rksyouyou/MeDL", upgrade = "never")
```

Installation compiles the package's native code. Use Rtools on Windows,
Xcode Command Line Tools on macOS, or a C and C++ compiler on Linux.

## Run the simulated example

```r
library(MeDL)
joint <- do.call(medl_prepare, medl_example)
fit <- medl_fit(joint, B = 99, seed = 42)
fit$p_values
fit$fits[["REML-KR"]][c("theta", "gamma", "tau2")]
```

The example contains summary statistics, signed prediction weights and a genotype
covariance matrix for 20 simulated variants. The script
[`data-raw/medl_example.R`](data-raw/medl_example.R) generates every input
from a fixed seed. The 99 bootstrap draws illustrate
the workflow and give a minimum attainable bootstrap P-value of 0.01. See
`?medl_example` for the data-generating parameters and `?medl_fit` for the
test outputs.

## Inputs for one mediator

Supply at least three variants, aligned to the same effect allele and in the same order in every input. Use signed weights that vary across variants. If you prune variants, subset every input before calling MeDL. Supplying consistent variant IDs as vector names and matrix row and column names allows the code to check the order.

- `marginal_beta` and `marginal_se` are numeric vectors of marginal GWAS effect estimates and positive standard errors. Use regression coefficients for continuous outcomes and log odds ratios for binary outcomes.
- `weights` contains mediator prediction coefficients on the same genotype dosage scale.
- `cov_g` is the reference-panel genotype covariance matrix on that scale, including the per-variant genotype variances on the diagonal.
- `ridge` is the amount added to the diagonal of the genotype correlation matrix for reconstruction. The default is zero. The manuscript application used `1 / (sqrt(n_ld) * log(n_ld))`, with `n_ld` equal to the reference-panel sample size. Reconstruction requires a positive-definite matrix after any regularization.

The following calls use objects supplied by the analyst:

```r
joint <- medl_prepare(marginal_beta, marginal_se, weights, cov_g)
result <- medl_fit(joint)
result$p_values
```

Alternatively, supply reconstructed effects and their positive-definite sampling covariance directly:

```r
joint <- list(beta_hat = joint_beta, W = weights, Sigma = joint_covariance)
```

## Inference

`medl_fit()` returns the ML Wald, REML-KR Wald, MoM Wald, MoM-mKH Wald and asymptotic ML likelihood-ratio (LR) tests of a zero mediator slope. REML-KR uses Kenward–Roger standard errors and a t reference with P − 2 degrees of freedom. MoM-mKH applies the modified Knapp–Hartung variance adjustment and uses the same reference degrees of freedom. Here P is the number of retained variants.

Set `B` to a positive integer to include the gene-specific parametric bootstrap LR test. Choose the number of draws according to the precision needed. Its smallest attainable P-value is 1 / (B + 1).

```r
result <- medl_fit(joint, B = n_draws, seed = rng_seed)
result$p_values
result$fits[["REML-KR"]][c("theta", "gamma", "tau2", "gamma_var_KR")]
result$bootstrap
```

The bootstrap generates draws from the fitted intercept-only ML null model and refits both ML models for each draw. Its P-value is (1 + the number of bootstrap LR statistics at least as large as the observed statistic) / (B + 1). The optional `bootstrap_null_interval` argument takes two variance bounds bracketing an interior minimum of the null negative log likelihood. It is intended for sensitivity analysis and changes only the parameters used to generate bootstrap draws. The observed LR statistic and the bootstrap refits retain the standard fitting rule.

The returned `p_values` are unadjusted P-values for one mediator. `fits` contains the estimates and test diagnostics, `null` contains the intercept-only ML fit, and `bootstrap` contains the generating parameters, draw count, exceedances and bootstrap P-value when requested.
