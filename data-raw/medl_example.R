# Run from the package root. All inputs are generated here using base R.
set.seed(20261124)
p <- 20L
ids <- sprintf("sim%02d", seq_len(p))
lag <- abs(outer(seq_len(p), seq_len(p), "-"))
cor_g <- 0.45 ^ lag
genotype_sd <- sqrt(2 * 0.3 * (1 - 0.3))
cov_g <- cor_g * genotype_sd^2
dimnames(cov_g) <- list(ids, ids)

weights <- setNames(runif(p, -0.3, 0.3), ids)
joint_beta <- 0.01 + 0.15 * weights + rnorm(p, sd = 0.025)
marginal_mean <- drop(cov_g %*% joint_beta) / diag(cov_g)
marginal_se <- setNames(runif(p, 0.012, 0.020), ids)
noise <- drop(crossprod(chol(cor_g), rnorm(p)))
marginal_beta <- marginal_mean + marginal_se * noise

medl_example <- list(marginal_beta = marginal_beta, marginal_se = marginal_se,
                     weights = weights, cov_g = cov_g)
save(medl_example, file = "data/medl_example.rda", version = 2)
