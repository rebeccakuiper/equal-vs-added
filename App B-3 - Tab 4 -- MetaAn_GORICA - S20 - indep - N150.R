###############################################################################
# Data: S = 20 independently generated datasets, n_s = 150 per study; H_123 TRUE.
# Used in: Appendix B.3, Table 4: the overall GORICA weights for H_123 for FE, RE, Equal,
#   and Added (columns), versus the unconstrained and the complement (rows);
#   these are the rows of the 'Final GORICA weights' tables printed by run_showcase().
# The cumulative plots (MA_FE_..., MA_RE_...) are not in the paper, but they are
#   saved by this script.
###############################################################################

###############################################################################
# Meta-analysis (FE and RE) + GORICA  vs.  equal / added evSyn
#
# Same data as in the evSyn showcase script
#   'Sec 5-2-3 & App D - Tab 5 -- evSyn - S20 - indep - N150.R'
# so that the results can be compared one-to-one.
#
# Steps:
#  1. Generate the study data exactly as in the evSyn script (same seed).
#  2. Meta-analyze the S studies (multivariate; the 3 regression slopes) with
#     a) a fixed-effect (FE; common-effect) model and
#     b) a random-effects (RE) model (between-study covariance matrix: 'UN',
#        with fallback to 'DIAG' if 'UN' cannot be estimated),
#     using metafor::rma.mv().
#  3. Evaluate H_123 with the GORICA on the pooled estimates and their
#     covariance matrix (vs. unconstrained and vs. complement).
#  4. Do the same cumulatively (studies 1..s, s = 1, ..., S) such that the
#     plots can be compared to the evSyn plots.
#  5. Run equal and added evSyn (as in the evSyn script) and put all final
#     GORICA weights in one table.
#
# Note on the intercept: H_123 only concerns the slopes. Because the intercept
# is unconstrained under every hypothesis, the GORICA based on the slopes only
# gives the same log-likelihood differences and penalties as the GORICA based
# on intercept + slopes (as used in the evSyn scripts). The meta-analyses are
# therefore done on the 3 slopes only (this keeps the RE model estimable).
###############################################################################

if (!require("restriktor")) install.packages("restriktor")
library("restriktor")
if (!require("metafor")) install.packages("metafor")
library("metafor")
if (!require("ggplot2")) install.packages("ggplot2")
library("ggplot2")
library(MASS) # mvrnorm


# Specs situation (identical to the evSyn script)
S <- 20 # 2 # 4 # 20
n_s <- 150 # 25, 50, 75, 100, 150, 200, 300, 500, 1,000
f2 <- 0.15 # medium; f2 = R2 / (1 - R2)
R2 <- f2 / (1 + f2)
rho <- 0.1
opposite <- FALSE # TRUE: H_123 is 'x1 < x2 < x3' (false theory)
indep_samples <- TRUE # TRUE: S independent random samples (as in the '... - indep' scripts)

name_folder <- paste0('images_S', S, "_N", n_s, if (opposite) " - opposite" else "")
dir.create(file.path(name_folder), showWarnings = FALSE)
suffix <- paste0("_S", S, if (opposite) "_Opp" else "", if (indep_samples) "_indep" else "",
                 if (n_s != 85) paste0("_N", n_s) else "")


###############################################################################
# 1. Generate data
###############################################################################

n.coef <- 3
var.e <- 1-R2
sigma <- matrix(rho, n.coef, n.coef); diag(sigma) <- 1
ratio_study <- c(3, 2, 1)
fun <- function(x){
  (t(ratio_study*x) %*% sigma %*% ratio_study*x) / (t(ratio_study*x) %*% sigma %*% ratio_study*x + var.e) - R2
}
x <- uniroot(fun, lower=0, upper=100)$root
betas <- ratio_study*x
mu <- rep(0, n.coef)
intercept <- 0

# One dataset, as in the evSyn scripts ('empirical = TRUE' when exact = TRUE)
gen_one <- function(exact) {
  if(n_s %% 2 == 0){ # even
    e <- mvrnorm(n = n_s/2, 0, var.e, empirical=exact)
    e_even <- c(e, -e)
    epsilon <- e_even*sqrt(var.e/var(e_even))
    x <- as.data.frame(mvrnorm(n = n_s/2, mu, sigma, empirical=exact))
    x_even <- rbind(x, x)
    sample <- x_even*sqrt(diag(sigma/var(x_even)))
  }else{ # odd
    e <- mvrnorm(n = (n_s-1)/2, 0, var.e, empirical=exact)
    e_odd <- c(e , 0, -e)
    epsilon <- e_odd*sqrt(var.e/var(e_odd))
    x <- as.data.frame(mvrnorm(n = (n_s-1)/2, mu, sigma, empirical=exact))
    x_odd <- rbind(x, rep(0, n.coef), x)
    sample <- x_odd*sqrt(diag(sigma/var(x_odd)))
  }
  sample$y <- intercept + as.matrix(sample, ncol = n.coef) %*% matrix(betas, ncol = 1) + epsilon
  names(sample) <- c("x1", "x2", "x3", "y")
  lm(y ~ 1 + x1 + x2 + x3, data = sample)
}

if (indep_samples) {
  # S independent random samples (as in the '... - indep' scripts)
  Param_studies <- list()
  CovMx_studies <- list()
  set.seed(123)
  for(s in 1:S){
    x <- as.data.frame(mvrnorm(n = n_s, mu, sigma)) # predictors
    epsilon <- rnorm(n = n_s, mean = 0, sd = sqrt(var.e)) # error terms
    sample <- x
    sample$y <- intercept + as.matrix(sample, ncol = n.coef) %*% matrix(betas, ncol = 1) + epsilon
    names(sample) <- c("x1", "x2", "x3", "y")
    fit.lm <- lm(y ~ 1 + x1 + x2 + x3, data = sample)
    Param_studies[[s]] <- coef(fit.lm)
    CovMx_studies[[s]] <- vcov(fit.lm)
  }
} else {
  # S times the same dataset
  set.seed(123)
  fit.lm <- gen_one(exact = TRUE)
  est <- coef(fit.lm)
  VCOV <- vcov(fit.lm)
  Param_studies <- rep(list(est), S)
  CovMx_studies <- rep(list(VCOV), S)
}

# Informative hypothesis
H_123 <- if (opposite) 'x1 < x2 < x3' else 'x1 > x2 > x3'


###############################################################################
# 2.-4. Helper functions: meta-analysis + GORICA
###############################################################################

slopes <- c("x1", "x2", "x3")

# Multivariate meta-analysis of the slopes of the studies in 'Param' / 'CovMx'.
# Returns the pooled estimates and their covariance matrix (named x1, x2, x3).
metaan <- function(Param, CovMx, model = c("FE", "RE")) {
  model <- match.arg(model)
  k <- length(Param)
  dat <- data.frame(study = rep(1:k, each = length(slopes)),
                    coef  = factor(rep(slopes, k), levels = slopes),
                    yi    = unlist(lapply(Param, function(p) p[slopes])))
  V <- bldiag(lapply(CovMx, function(v) v[slopes, slopes]))
  if (model == "FE" || k < 2) {
    # Note: with 1 study, the RE model cannot be estimated, then RE = FE.
    res <- rma.mv(yi, V, mods = ~ 0 + coef, data = dat, method = "ML")
    struct <- "none"
  } else {
    res <- tryCatch({
      r <- rma.mv(yi, V, mods = ~ 0 + coef, random = ~ coef | study,
                  struct = "UN", data = dat, method = "ML")
      struct <- "UN"; r
    }, error = function(e) NULL)
    if (is.null(res)) {
      res <- rma.mv(yi, V, mods = ~ 0 + coef, random = ~ coef | study,
                    struct = "DIAG", data = dat, method = "ML")
      struct <- "DIAG"
    }
  }
  est <- coef(res); names(est) <- slopes
  vc  <- vcov(res); dimnames(vc) <- list(slopes, slopes)
  list(fit = res, est = est, VCOV = vc, struct = struct)
}

gorica_ma <- function(ma, comparison) {
  goric(ma$est, VCOV = ma$VCOV,
        hypotheses = list(H_123 = H_123),
        comparison = comparison, # default in the case of one specified hypothesis is 'complement'
        type = 'gorica')         # default is 'gorica' when est & VCOV are the input
}

# Cumulative: meta-analysis + GORICA on studies 1..s, for s = 1, ..., S
cumulative_ma <- function(Param, CovMx, model, comparison) {
  W <- t(sapply(1:length(Param), function(s) {
    ma <- metaan(Param[1:s], CovMx[1:s], model = model)
    gorica_ma(ma, comparison)$result$gorica.weights
  }))
  colnames(W) <- if (comparison == "complement") c("H_123", "Complement") else c("H_123", "Unconstrained")
  rownames(W) <- paste0("Studies 1-", 1:length(Param))
  W
}

plot_cumulative <- function(W, title) {
  df <- data.frame(Study = rep(1:nrow(W), ncol(W)),
                   Hypothesis = factor(rep(colnames(W), each = nrow(W)), levels = colnames(W)),
                   Weight = as.vector(W))
  ggplot(df, aes(x = Study, y = Weight, colour = Hypothesis, linetype = Hypothesis)) +
    geom_line() + geom_point() +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = 1:nrow(W)) +
    labs(x = "Number of studies meta-analyzed (cumulative)", y = "GORICA weight",
         title = title) +
    theme_minimal()
}

save_plot <- function(p, fileName) {
  ggsave(filename = fileName, path = name_folder, plot = p,
         width = 8, height = 6, dpi = 300)
}


###############################################################################
# 5. Analyses per show case
###############################################################################

run_showcase <- function(Param, CovMx, label) {
  cat("\n\n================", label, "================\n")
  out <- NULL
  for (model in c("FE", "RE")) {
    ma <- metaan(Param, CovMx, model = model)
    cat("\n---", model, "meta-analysis (RE struct:", ma$struct, ") ---\n")
    print(ma$fit, digits = 3)
    for (comparison in c("unconstrained", "complement")) {
      g <- gorica_ma(ma, comparison)
      cat("\nGORICA on", model, "estimates, H_123 vs", comparison, ":\n")
      print(g)
      out <- rbind(out, data.frame(Method = paste0("GORICA on ", model),
                                   Failsafe = comparison,
                                   w_H123 = g$result$gorica.weights[1],
                                   w_failsafe = g$result$gorica.weights[2]))
      # cumulative + plot
      W <- cumulative_ma(Param, CovMx, model, comparison)
      fs <- if (comparison == "complement") "compl" else "unc"
      p <- plot_cumulative(W, paste0("GORICA on ", model, " meta-analysis: H_123 vs ", comparison,
                                     " (", label, ")"))
      print(p)
      save_plot(p, paste0("MA_", model, "_123_", fs, label_file(label), suffix, ".png"))
    }
  }
  # equal and added evSyn (as in the evSyn scripts), for comparison
  for (type_ev in c("equal", "added")) {
    for (comparison in c("unconstrained", "complement")) {
      evS <- evSyn(object = Param, VCOV = CovMx,
                   hypotheses = list(H_123 = H_123),
                   type_ev = type_ev, comparison = comparison)
      fw <- evS$Cumulative_GORICA_weights[nrow(evS$Cumulative_GORICA_weights), ]
      out <- rbind(out, data.frame(Method = paste0(type_ev, " evSyn"),
                                   Failsafe = comparison,
                                   w_H123 = fw[1], w_failsafe = fw[2]))
    }
  }
  rownames(out) <- NULL
  out$ratio <- out$w_H123 / out$w_failsafe
  cat("\n--- Final GORICA weights (", label, ") ---\n")
  print(out, digits = 3)
  invisible(out)
}
label_file <- function(label) if (grepl("^1 outlier", label)) "_out" else ""


## Show case 1a: no outliers ##
# >>> PAPER: Appendix B.3, Table 4, rows 'Outlier = no'
res_1a <- run_showcase(Param_studies, CovMx_studies, "no outlier")


