###############################################################################
# App B-4 -- MetaAn_GORICA vs evSyn - Berkey (FE vs RE).R
# >>> PAPER: Appendix B, subsection 'When Do FE and RE Differ?':
#   data of Berkey et al. (1998) (dat.berkey1998 in metafor): five studies comparing
#   surgical and non-surgical treatment of periodontal disease on two outcomes,
#   attachment level (AL) and probing depth (PD); effect size: (raw) mean difference.
#   Hypothesis (Set 2 in 'Example_GORICAonMeta-an.R' on
#   https://github.com/rebeccakuiper/GORICA_on_MetaAn; Kuiper, 2022):
#     H2: |theta_AL| > 0.2 & |theta_PD| > 0.2 (both effects at least small),
#   evaluated against its complement.
# Output (printed): the numbers reported in the text:
#   - heterogeneity (Q = 128.2, df = 8);
#   - FE and RE estimates and standard errors;
#   - GORICA weights based on the FE and RE meta-analyzed estimates;
#   - study-specific estimates and GORICA weights;
#   - equal- and added-evidence GORICA evidence synthesis (GORICA and log-likelihood weights);
#   - 'RE-like' equal-evidence approach: the estimated between-study covariance matrix
#     is added to the covariance matrix of each study.
###############################################################################

library(metafor)
library(restriktor)

dat <- dat.berkey1998
pars <- c("theta_AL", "theta_PD")

# within-study covariance matrices (block-diagonal V; per trial in the order PD, AL as in the data)
V <- bldiag(lapply(split(dat[, c("v1i", "v2i")], dat$trial), as.matrix))

# 1. Multivariate meta-analysis: fixed-effect (FE) and random-effects (RE)
#    (method = "ML", as in Kuiper (2022), to mimic Berkey et al. (1998))
fe <- rma.mv(yi, V, mods = ~ outcome - 1, data = dat, method = "ML")
re <- rma.mv(yi, V, mods = ~ outcome - 1, random = ~ outcome | trial, struct = "UN",
             data = dat, method = "ML")
cat("Heterogeneity: Q =", round(fe$QE, 1), ", df =", fe$k - fe$p, ", p =", signif(fe$QEp, 3), "\n")

# estimates and covariance matrices (coefficients in the order outcomeAL, outcomePD)
get_est  <- function(m) setNames(coef(m), pars)
get_vcov <- function(m) { v <- vcov(m); dimnames(v) <- list(pars, pars); v }
print(round(rbind(est_FE = get_est(fe), est_RE = get_est(re),
                  SE_FE = sqrt(diag(get_vcov(fe))), SE_RE = sqrt(diag(get_vcov(re)))), 3))

H2 <- "abs(theta_AL) > 0.2; abs(theta_PD) > 0.2"

# 2. GORICA on the FE and RE meta-analyzed estimates (versus the complement)
gorica_FE <- goric(get_est(fe), VCOV = get_vcov(fe), hypotheses = list(H2 = H2),
                   comparison = "complement", type = "gorica")
gorica_RE <- goric(get_est(re), VCOV = get_vcov(re), hypotheses = list(H2 = H2),
                   comparison = "complement", type = "gorica")
gorica_FE
gorica_RE

# 3. GORICA evidence synthesis (equal and added evidence)
est_studies  <- lapply(split(dat, dat$trial), function(d) setNames(d$yi[match(c("AL", "PD"), d$outcome)], pars))
vcov_studies <- lapply(split(dat, dat$trial), function(d) {
  v <- as.matrix(d[, c("v1i", "v2i")]); dimnames(v) <- list(d$outcome, d$outcome)
  v <- v[c("AL", "PD"), c("AL", "PD")]; dimnames(v) <- list(pars, pars); v })
cat("Study-specific estimates:\n"); print(round(do.call(rbind, est_studies), 2))

evS_equal <- evSyn(est_studies, VCOV = vcov_studies, hypotheses = list(H2 = H2),
                   type_ev = "equal", comparison = "complement")
evS_added <- evSyn(est_studies, VCOV = vcov_studies, hypotheses = list(H2 = H2),
                   type_ev = "added", comparison = "complement")
evS_equal
evS_added
cat("Study-specific GORICA weights:\n"); print(round(evS_equal$GORICA_weight_m, 2))
cat("Number of studies with more support for H2 than for its complement:",
    sum(evS_equal$GORICA_weight_m[, 1] > evS_equal$GORICA_weight_m[, 2]), "of", length(est_studies), "\n")

# 4. 'RE-like' equal-evidence approach: add the estimated between-study covariance
#    matrix (Tau) of the RE model to the covariance matrix of each study
Tau <- re$G                                  # rows/columns in the order AL, PD
dimnames(Tau) <- list(pars, pars)
evS_equal_RElike <- evSyn(est_studies, VCOV = lapply(vcov_studies, function(v) v + Tau),
                          hypotheses = list(H2 = H2), type_ev = "equal", comparison = "complement")
evS_equal_RElike

# Summary (numbers in the text): overall GORICA weight and log-likelihood (LL) weight of H2
w_goric <- function(x) c(GORICA_weight_H2 = x$result$gorica.weights[1], LL_weight_H2 = x$result$loglik.weights[1])
w_evSyn <- function(x) c(GORICA_weight_H2 = tail(x$Cumulative_GORICA_weights[, 1], 1),
                         LL_weight_H2 = tail(x$Cumulative_LL_weights[, 1], 1))
print(round(rbind(GORICA_FE    = w_goric(gorica_FE),
                  GORICA_RE    = w_goric(gorica_RE),
                  Equal        = w_evSyn(evS_equal),
                  Added        = w_evSyn(evS_added),
                  Equal_RElike = w_evSyn(evS_equal_RElike)), 3))
