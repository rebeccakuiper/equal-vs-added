# estimate_betas() with the CORRECT latent error variances for the probit and logit models.
# In sim/functions.R, these two were swapped (probit = pi^2/3, logit = 1), whereas the data
# are generated with pnorm() (probit: standard-normal latent error, variance 1) and with the
# logistic CDF (logit: logistic latent error, variance pi^2/3). Consequently, the probit
# studies had a (much) larger and the logit studies a (much) smaller latent R2 than intended
# (e.g., for R2 = .15: .37 and .05). Apart from these two values, the function is identical.
estimate_betas <- function(true_R2, ratios, mu, sigma, n, model, intercept, type) {
  ratio_study <- ratios[1:length(mu)]
  switch(type,
    "linear" = {var.e <- 1-true_R2},  # var(y) = 1
    "probit" = {var.e <- 1},          # CORRECTED: latent error ~ N(0, 1)
    "logit"  = {var.e <- pi^2/3}      # CORRECTED: latent error ~ logistic(0, 1), variance pi^2/3
  )
  fun <- function(x){
    b <- ratio_study * x
    drop(t(b) %*% sigma %*% b) / (drop(t(b) %*% sigma %*% b) + var.e) - true_R2
  }
  x <- uniroot(fun, lower = 0, upper = 100)$root
  betas <- ratio_study * x
  sample <- as.data.frame(mvrnorm(n = n, mu, sigma, empirical = FALSE))
  sample <- as.data.frame(scale(sample))
  switch(type,
     "linear" = {
       epsilon <- rnorm(n, sd = sqrt(var.e))
       y <- intercept + as.matrix(sample) %*% matrix(betas, ncol = 1) + epsilon
     },
     "probit" = {
       sample$ystar <- intercept + as.matrix(sample) %*% matrix(betas, ncol = 1)
       sample$prob <- pnorm(sample$ystar)
       sample$y <- rbinom(n, size = 1, prob = sample$prob)
     },
     "logit" = {
       sample$ystar <- intercept + as.matrix(sample) %*% matrix(betas, ncol = 1)
       sample$prob <- exp(sample$ystar) / (1 + exp(sample$ystar))
       sample$y <- rbinom(n, size = 1, prob = sample$prob)
     })
  fit <- eval(parse(text = model))
  betahat <- coef(fit); VCovMatrix <- vcov(fit)
  names(betahat) <- sprintf('beta%d', 0:(length(betahat)-1))
  list(ests = betahat, vcovmat = VCovMatrix)
}
