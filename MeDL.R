check_dat <- function(dat) {
    y = as.numeric(dat$beta_hat)
    W = as.numeric(dat$W)
    S = as.matrix(dat$Sigma)
    P = length(y)
    if(P <= 2 || length(W) != P || !identical(dim(S), c(P,P)))
        stop("Incompatible dimensions or fewer than three SNPs.")
    if(any(!is.finite(c(y,W,S)))) stop("Non-finite model inputs.")
    if(max(abs(S-t(S))) > 1e-10*max(1,max(abs(S))))
        stop("Sigma is not symmetric.")
    S = (S+t(S))/2
    chol(S) # Fail explicitly rather than silently using a generalized inverse.
    Z = cbind(theta=1, gamma=W)
    if(qr(Z)$rank != 2) stop("Intercept and W are not identifiable: constant W.")
    list(y=y, W=W, Sigma=S, Z=Z, P=P)
}


fit_fixed <- function(d, tau2, gamma=NULL) {
    ## In profile searches, rotate once into the eigenbasis of Sigma.
    ## This is the same GLS likelihood; avoids a matrix inverse at every tau2.
    if(!is.null(d$eigen)) {
        inv = 1/(d$eigen$values + tau2)
        X = d$rotated_Z
        y = d$rotated_y
        if(is.null(gamma)) {
            H = crossprod(X, inv*X)
            bvar = solve(H)
            b = as.numeric(bvar %*% crossprod(X,inv*y))
        } else {
            theta = sum(X[,1]*inv*(y-X[,2]*gamma))/sum(X[,1]^2*inv)
            b = c(theta,gamma)
            bvar = matrix(NA_real_,2,2)
            bvar[1,1] = 1/sum(X[,1]^2*inv)
        }
        residual = y-as.numeric(X%*%b)
        obj = sum(log(d$eigen$values+tau2)) + sum(residual^2*inv)
        return(list(b=b,bvar=bvar,obj=obj,
                    condition=kappa(crossprod(X,inv*X))))
    }
    V = d$Sigma + diag(tau2,d$P)
    R = chol(V)
    Vinv = chol2inv(R)
    if(is.null(gamma)) {
        H = crossprod(d$Z,Vinv%*%d$Z)
        b = as.numeric(solve(H,crossprod(d$Z,Vinv%*%d$y)))
        bvar = solve(H)
    } else {
        ## Fixed slope, but intercept remains free, including under H0.
        theta = sum(Vinv%*%(d$y-d$W*gamma))/sum(Vinv)
        b = c(theta,gamma)
        bvar = matrix(NA_real_,2,2)
        bvar[1,1] = 1/sum(Vinv)
    }
    residual = d$y - as.numeric(d$Z%*%b)
    obj = 2*sum(log(diag(R))) + as.numeric(crossprod(residual,Vinv%*%residual))
    list(b=b, bvar=bvar, obj=obj, condition=kappa(crossprod(d$Z,Vinv%*%d$Z)))
}

method_ML <- function(dat, gamma=NULL) {
    if(!is.null(gamma) && (length(gamma)!=1 || !is.finite(gamma)))
        stop("gamma must be NULL or a finite scalar.")
    d = check_dat(dat)
    d$eigen = eigen(d$Sigma,symmetric=TRUE)
    if(any(d$eigen$values<=0)) stop("Nonpositive eigenvalue of Sigma.")
    d$rotated_Z = crossprod(d$eigen$vectors,d$Z)
    d$rotated_y = as.numeric(crossprod(d$eigen$vectors,d$y))
    ## Scale the search using observed data, never the generating tau2.
    residual = if(is.null(gamma)) lm.fit(d$Z,d$y)$residuals else
        d$y-d$W*gamma-mean(d$y-d$W*gamma)
    scale_tau = max(mean(residual^2),mean(diag(d$Sigma)),1e-12)
    objective = function(x) fit_fixed(d,x*scale_tau,gamma)$obj
    upper = 10
    search_complete = FALSE
    for(attempt in 1:8) {
        ## Grid brackets guard against multiple minima in a profile likelihood.
        grid = c(0,exp(seq(log(1e-8),log(upper),length.out=35)))
        values = vapply(grid,objective,numeric(1))
        candidates = grid
        for(j in 2:(length(grid)-1)) {
            if(values[j]<=values[j-1] && values[j]<=values[j+1]) {
                opt = optimize(objective,c(grid[j-1],grid[j+1]),tol=1e-9)
                candidates = c(candidates,opt$minimum)
            }
        }
        opt = optimize(objective,c(0,upper),tol=1e-9)
        candidates = c(candidates,opt$minimum)
        scores = vapply(candidates,objective,numeric(1))
        best = candidates[which.min(scores)]
        if(best < upper*0.95) {
            search_complete = TRUE
            break
        }
        upper = upper*10
    }
    if(!search_complete) stop("Variance search reached the expansion limit.")
    ## Compare zero explicitly; allow only roundoff-level likelihood ties.
    if(objective(0)<=objective(best)+1e-10) best = 0
    tau2 = best*scale_tau
    fit = fit_fixed(d,tau2,gamma)
    pval = if(is.null(gamma)) pchisq(fit$b[2]^2/fit$bvar[2,2],1,lower.tail=FALSE) else NA_real_
    c(tau2=tau2, theta=fit$b[1], gamma=fit$b[2],
      theta_var=fit$bvar[1,1], gamma_var=fit$bvar[2,2],
      pval=pval, obj=fit$obj, at_zero=as.numeric(tau2==0), condition=fit$condition)
}

## Regression-residual MoM, accounting for BOTH theta and gamma.
## Sigma must be positive definite; Z=[1,W] must have rank 2.
method_MoM <- function(dat) {
    d = check_dat(dat)
    Omega = chol2inv(chol(d$Sigma))
    OZ = Omega%*%d$Z
    A = Omega - OZ%*%solve(crossprod(d$Z,OZ),t(OZ))
    A = (A+t(A))/2
    trace_A = sum(diag(A))
    if(!is.finite(trace_A) || trace_A<=0) stop("Invalid MoM denominator.")
    Q = as.numeric(crossprod(d$y,A%*%d$y))
    tau2 = max(0,(Q-(d$P-2))/trace_A)
    fit = fit_fixed(d,tau2)
    c(tau2=tau2, theta=fit$b[1], gamma=fit$b[2],
      theta_var=fit$bvar[1,1], gamma_var=fit$bvar[2,2],
      pval=pchisq(fit$b[2]^2/fit$bvar[2,2],1,lower.tail=FALSE),
      at_zero=as.numeric(tau2==0), Q=Q, df=d$P-2, condition=fit$condition)
}

## Modified Knapp-Hartung inference, retaining the original MoM estimates.
## Knapp & Hartung (2003), doi:10.1002/sim.1482;
## Rover et al. (2015), doi:10.1186/s12874-015-0091-1.
## Full-covariance adaptation follows the matrix formulation of Jackson &
## Riley (2014), doi:10.1002/sim.5957.
method_MoM_mKH <- function(dat, fit = method_MoM(dat)) {
    d = check_dat(dat)
    V = d$Sigma + diag(unname(fit['tau2']),d$P)
    residual = d$y - as.numeric(d$Z%*%fit[c('theta','gamma')])
    whitened = forwardsolve(t(chol(V)),residual)
    df = d$P-2
    ## This uses fitted V, NOT Sigma from the MoM estimating equation.
    q = sum(whitened^2)/df
    q_adj = max(1,q)
    fit[c('theta_var','gamma_var')] = q_adj*fit[c('theta_var','gamma_var')]
    statistic = unname(fit['gamma']/sqrt(fit['gamma_var']))
    fit['pval'] = 2*pt(abs(statistic),df=df,lower.tail=FALSE)
    ## Q remains the original MoM diagnostic; q is the KH residual scale.
    c(fit,q=q,q_adj=q_adj,statistic=statistic)
}

## REML variance estimate, KR covariance correction, fixed t(P-2) reference.
## V = Sigma + tau2 * I; fixed effects are (theta, gamma).
## Kenward & Roger (1997), Biometrics 53:983-997, doi:10.2307/2533558.
## This uses KR standard errors with a FIXED df, not the full KR test.
method_REML_KR <- function(dat) {
    d = check_dat(dat)
    eig = eigen(d$Sigma,symmetric=TRUE)
    e = eig$values
    if(any(e<=0)) stop('Sigma must be positive definite.')
    Q = crossprod(eig$vectors,d$Z)
    y = as.numeric(crossprod(eig$vectors,d$y))
    fit = .Call('reml_full', y, e, Q, PACKAGE = 'medl')
    tau2 = fit[1]
    iv = 1/(e+tau2)
    A1 = crossprod(Q,Q*iv)
    Phi = chol2inv(chol(A1))
    a = as.numeric(Phi %*% crossprod(Q,y*iv))
    A2 = crossprod(Q,Q*iv^2)
    A3 = crossprod(Q,Q*iv^3)
    ## Expected REML information for the single variance parameter tau2.
    info2 = sum(iv^2)-2*sum(Phi*t(A3)) + sum((Phi%*%A2)*t(Phi%*%A2))
    if(!is.finite(info2) || info2<=0) stop('Invalid REML information for tau2.')
    var_tau2 = 2/info2
    adjusted = Phi + 2*var_tau2*Phi%*%(A3-A2%*%Phi%*%A2)%*%Phi
    adjusted = (adjusted+t(adjusted))/2
    chol(adjusted)
    if(any(diag(adjusted)<diag(Phi)*(1-1e-8))) stop('Unexpected negative KR correction.')
    statistic = a[2]/sqrt(adjusted[2,2])
    df = d$P-2
    ans = c(tau2=tau2,theta=a[1],gamma=a[2],theta_var=Phi[1,1],gamma_var=Phi[2,2],
        theta_var_KR=adjusted[1,1],gamma_var_KR=adjusted[2,2],
        theta_gamma_cov=Phi[1,2],theta_gamma_cov_KR=adjusted[1,2],
        pval=2*pt(abs(statistic),df,lower.tail=FALSE),statistic=statistic,df=df,
        obj=fit[4],at_zero=as.numeric(tau2==0),condition=kappa(A1,exact=TRUE),
        var_tau2=var_tau2)
    if(any(!is.finite(ans)) || ans['pval']<0 || ans['pval']>1 || tau2<0)
        stop('Invalid REML/KR output.')
    ans
}

method_LR <- function(dat, gamma0=0, fit_full=NULL, fit_null=NULL) {
    if(is.null(fit_full)) fit_full = method_ML(dat)
    if(is.null(fit_null)) fit_null = method_ML(dat,gamma=gamma0)
    statistic = unname(fit_null['obj']-fit_full['obj'])
    if(statistic < -1e-7) stop("Full-model likelihood below null: check optimization.")
    statistic = max(0,statistic)
    c(tau2_H0=unname(fit_null['tau2']),tau2_full=unname(fit_full['tau2']),
      theta_H0=unname(fit_null['theta']),theta_full=unname(fit_full['theta']),
      gamma=unname(fit_full['gamma']),statistic=statistic,
      pval=pchisq(statistic,1,lower.tail=FALSE))
}

# When IDs are supplied, require a single order for vectors and both matrix axes.
check_variant_order <- function(inputs, covariance) {
    ids <- Filter(Negate(is.null), c(lapply(inputs, names),
                                     list(rownames(covariance), colnames(covariance))))
    if (length(ids)) {
        first <- ids[[1]]
        if (anyNA(first) || any(!nzchar(first)) || anyDuplicated(first) ||
            !all(vapply(ids, identical, logical(1), first))) {
            stop("Variant IDs must be unique and agree in order across inputs.")
        }
    }
}

check_vector <- function(x, size, label) {
    if (!is.numeric(x) || !is.null(dim(x)) || length(x) != size || any(!is.finite(x))) {
        stop(label, " must be a finite numeric vector with one value per variant.")
    }
}

# Reconstruct joint effects on the allele-dosage scale after harmonization.
medl_prepare <- function(marginal_beta, marginal_se, weights, cov_g, ridge = 0) {
    p <- length(weights)
    check_vector(weights, p, "weights")
    check_vector(marginal_beta, p, "marginal_beta")
    check_vector(marginal_se, p, "marginal_se")
    if (p < 3L || any(marginal_se <= 0)) stop("Need at least three variants and positive standard errors.")
    if (!is.matrix(cov_g) || !is.numeric(cov_g) || !identical(dim(cov_g), c(p,p)) ||
        any(!is.finite(cov_g)) || any(diag(cov_g) <= 0) ||
        !isTRUE(all.equal(cov_g, t(cov_g), check.attributes = FALSE, tolerance = 1e-10))) {
        stop("cov_g must be a symmetric genotype covariance matrix with positive diagonal.")
    }
    check_variant_order(list(marginal_beta, marginal_se, weights), cov_g)
    if (!is.numeric(ridge) || length(ridge) != 1L || !is.finite(ridge) || ridge < 0) {
        stop("ridge must be a nonnegative, finite scalar.")
    }

    genotype_sd <- sqrt(diag(cov_g))
    corr <- cov2cor((cov_g + t(cov_g)) / 2) + diag(ridge, p)
    transform <- chol2inv(chol(corr * outer(genotype_sd, genotype_sd))) %*% diag(diag(cov_g))
    joint_cov <- transform %*% (corr * outer(marginal_se, marginal_se)) %*% t(transform)
    joint_beta <- drop(transform %*% marginal_beta)
    list(beta_hat = joint_beta, W = weights, Sigma = (joint_cov + t(joint_cov)) / 2)
}

# A supplied interval changes only the null parameters used to generate bootstrap draws.
medl_bootstrap <- function(d, null, observed, B, seed, interval) {
    eig <- eigen(d$Sigma, symmetric = TRUE)
    design <- crossprod(eig$vectors, d$Z)
    intercept <- design[, 1L]
    tau2 <- unname(null["tau2"])
    theta <- unname(null["theta"])
    if (!is.null(interval)) {
        if (!is.numeric(interval) || length(interval) != 2L || any(!is.finite(interval)) ||
            interval[1] < 0 || interval[2] <= interval[1]) stop("Invalid bootstrap null interval.")
        y <- drop(crossprod(eig$vectors, d$y))
        intercept_at <- function(v) sum(intercept * y / v) / sum(intercept^2 / v)
        objective <- function(t) {
            v <- eig$values + t
            sum(log(v)) + sum((y - intercept * intercept_at(v))^2 / v)
        }
        local <- optimize(objective, interval, tol = 1e-15)
        if (local$objective >= min(vapply(interval, objective, numeric(1)))) {
            stop("The interval must bracket an interior null-profile minimum.")
        }
        tau2 <- local$minimum
        theta <- intercept_at(eig$values + tau2)
    }

    if (!is.null(seed)) set.seed(seed)
    mean_y <- intercept * theta
    sd_y <- sqrt(eig$values + tau2)
    exceedances <- 0
    remaining <- B
    while (remaining > 0) {
        n <- min(remaining, 1000L)
        y <- matrix(rnorm(d$P * as.double(n)), d$P, n) * sd_y + mean_y
        statistics <- .Call("lr_profiles", y, eig$values, design, PACKAGE = "medl")[, 1L]
        exceedances <- exceedances + sum(statistics >= observed)
        remaining <- remaining - n
    }
    list(B = B, exceedances = exceedances, pval = (1 + exceedances) / (1 + B),
         theta = theta, tau2 = tau2, null_interval = interval, seed = seed)
}

# Test H0: gamma = 0 for a single mediator.
medl_fit <- function(joint, B = 0, seed = NULL, bootstrap_null_interval = NULL) {
    if (!is.list(joint) || !all(c("beta_hat", "W", "Sigma") %in% names(joint))) {
        stop("joint must contain beta_hat, W, and Sigma.")
    }
    if (!is.loaded("reml_full", PACKAGE = "medl")) stop("Load the compiled medl library as described in README.md.")
    check_vector(joint$W, length(joint$W), "W")
    check_vector(joint$beta_hat, length(joint$W), "beta_hat")
    if (!is.matrix(joint$Sigma) || !is.numeric(joint$Sigma)) stop("Sigma must be a numeric matrix.")
    check_variant_order(list(joint$W, joint$beta_hat), joint$Sigma)
    d <- check_dat(joint)
    if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 0 ||
        B != floor(B) || B > .Machine$integer.max) stop("B must be a nonnegative integer.")
    if (B == 0 && !is.null(bootstrap_null_interval)) stop("A null interval requires B > 0.")

    ml <- method_ML(joint)
    null <- method_ML(joint, gamma = 0)
    reml <- method_REML_KR(joint)
    mom <- method_MoM(joint)
    mkh <- method_MoM_mKH(joint, mom)
    lr <- method_LR(joint, fit_full = ml, fit_null = null)
    fits <- list(ML = ml, `REML-KR` = reml, MoM = mom, `MoM-mKH` = mkh, LR = lr)
    p_values <- vapply(fits, function(fit) unname(fit["pval"]), numeric(1))
    bootstrap <- NULL
    if (B > 0) {
        bootstrap <- medl_bootstrap(d, null, unname(lr["statistic"]), as.double(B),
                                    seed, bootstrap_null_interval)
        p_values["LR_boot"] <- bootstrap$pval
    }
    list(p_values = p_values, fits = fits, null = null, bootstrap = bootstrap)
}
