###############################################################################
# Data: S = 20 independently generated datasets, n_s = 150 per study; H_123 FALSE ('opposite': x1 < x2 < x3; failsafe true).
# Used in:
#   - Appendix D, Table 5 (unconstrained as failsafe)
#   - Section 5.2.3 (text: overall weights for n_s = 150)
# Each evSyn() call below is preceded by a '# >>> PAPER:' line stating where its
# result is used.
# The plots for the unconstrained as failsafe are not in the paper, but they are
# saved by this script.
###############################################################################

if (!require("restriktor")) install.packages("restriktor")
library("restriktor")
#?goric
#?evSyn

# If you want to use restriktor from github:
#if (!require("devtools")) install.packages("devtools")
#library(devtools) 
#install_github("LeonardV/restriktor")
#library(restriktor) # for goric function
##?goric
##?evSyn

if (!require("metafor")) install.packages("metafor")
library("metafor")

if (!require("ggplot2")) install.packages("ggplot2")
library("ggplot2")
# ggsave

###


# Instead of having S times the same dataset (as in other showcases / R scripts),
# this R script draws S independent random samples from the same population.


# Show case #

# Aggregate S studies 
# regression with 3 predictors, 
# each study has N observations (n_s)


# Specs situation
S <- 20 # 2 # 4 # 20
n_s <- 150 # 25, 50, 75, 100, 150, 200, 300, 500, 1,000
#n_tot <- S*n_s
f2 <- 0.15 # medium; f2 = R2 / (1 - R2)
R2 <- f2 / (1 + f2) # medium: R2 = f2 / (1 + f2) with f2 = 0.15 is approx R2 = 0.1304348
rho <- 0.1


#name_folder <- paste0('images_S', S, "_N", n_s)
name_folder <- paste0('images_S', S, "_N", n_s, " - opposite")
dir.create(file.path(name_folder), showWarnings = FALSE)



# generate data (linear model)
n.coef <- 3 # + 1 for intercept
var.e <- 1-R2
sigma <- matrix(rho, n.coef, n.coef) # correlation matrix of predictors
diag(sigma) <- 1
ratio_study <- c(3, 2, 1) # ratio_study <- c(4, 2, 1) # ratio_study <- c(1.77, 1.33, 1.00) # all length n.coef
# obtain beta values based on sample R2 and ratios
fun <- function(x){ # Function that must be minimized (or maximized) to obtain sample R2 value.
  (t(ratio_study*x) %*% sigma %*% ratio_study*x) / (t(ratio_study*x) %*% sigma %*% ratio_study*x + var.e) - R2
}
x <- uniroot(fun, lower=0, upper=100)$root 
# x <- 0.09448432
# R2:
#(t(ratio_study*x) %*% sigma %*% ratio_study*x) / (t(ratio_study*x) %*% sigma %*% ratio_study*x + var.e)
#
betas <- ratio_study*x # beta estimates for the sample.
#betas
#betas[1]-betas[2]; betas[2]-betas[3]
#
# sample residuals ('epsilon') and predictors ('sample')
# Using exact parameters:
# One approach is to use a perfectly symmetrical noise. The noise cancels itself so the estimated parameters are exactly the input parameters, yet the residuals appear normally distributed.
library(MASS) # mvrnorm
#set.seed(123)
mu <- rep(0, n.coef)
#
#resid <- matrix(NA, nrow = n_s, ncol = S)
Param_studies <- list()
CovMx_studies <- list()
set.seed(123)
for(s in 1:S){
  # s = 1
  #
  # Independent random sample from the population (no mirrored residuals and
  # no rescaling), such that the estimates vary over studies as they would
  # in practice.
  x <- as.data.frame(mvrnorm(n = n_s, mu, sigma)) # predictors
  epsilon <- rnorm(n = n_s, mean = 0, sd = sqrt(var.e)) # error terms
  sample <- x
  #
  # Create y (outcome) based on determined betas (regr coeff) and sampled residuals ('epsilon') and predictors ('sample')
  intercept <- 0
  sample$y <- intercept + as.matrix(sample, ncol = n.coef) %*% matrix(betas, ncol = 1) + epsilon   # sample outcome variables
  names(sample) <- c("x1", "x2", "x3", "y")
  #var(sample)
  #colMeans(sample)
  #
  ## I already assumed and created standardized data, otherwise:
  #sample <- as.data.frame(scale(sample)) # not needed when ones on diag of sigma
  ##var(sample)
  ##colMeans(sample)
  #
  # Run regression model (with intercept which will be approx 0):
  fit.lm <- lm(y ~ 1 + x1 + x2 + x3, data = sample)
  #coef(fit.lm)[-1]; betas
  #((summary(fit.lm)$sigma)^2)*(n_s-(n.coef+1))/(n_s-1); var.e

  Param_studies[[s]] <- coef(fit.lm)
  CovMx_studies[[s]] <- vcov(fit.lm)
}

###


# Informative hypothesis
H_123 <- 'x1 < x2 < x3'


## H_123 vs its Compl
# GORICA
set.seed(123)
est <- coef(fit.lm)
VCOV <- vcov(fit.lm)
gorica_lm <- goric(est, VCOV = VCOV, 
                   hypotheses = list(H_123 = H_123),
                   comparison = "complement",
                   type = 'gorica') 
gorica_lm
gorica_lm$ratio.lw


## Combined:


## Show case 1a: no outliers ##

# # Make lists of param. estimates and their cov. matrices from each study
# rep_est_text <- paste0(rep(" est", (S)), collapse=",")
# est_text <- noquote(paste("list(", rep_est_text, ")"))
# Param_studies <- eval(parse(text = est_text))
# #Param_studies
# rep_vcov_text <- paste0(rep(" VCOV", (S)), collapse=",")
# vcov_text <- noquote(paste("list(", rep_vcov_text, ")"))
# CovMx_studies <- eval(parse(text = vcov_text))
#
#
H_123
# vs unc
#
# >>> PAPER: no outlier, added-evidence: 
# Appendix D, Table 5 (unconstrained as failsafe; row: Indep., false, n_s = 150, outlier = no; column: Added, S = 20)
evS_added_123_unc <- evSyn(object = Param_studies, VCOV = CovMx_studies, 
                           hypotheses = list(H_123 = H_123),
                   type_ev = "added", # default
                   comparison = "unconstrained",
                   order_studies = "ascending") 
evS_added_123_unc$Cumulative_GORICA_weights
evS_added_123_unc$Cumulative_GORICA
plot(evS_added_123_unc)
#png(paste0(name_folder, "/evS_added_123_unc_S", S, "_Opp.png"))
#plot(evS_added_123_unc)
#dev.off()
plot_evS_added_123_unc_indep_N150 <- plot(evS_added_123_unc)
fileName <- paste0("evS_added_123_unc_S", S, "_Opp_indep_N150.png")
ggsave(
  filename = fileName,  # file name
  path = name_folder,
  plot = plot_evS_added_123_unc_indep_N150,     # the plot object
  width = 8,                                    # width in inches
  height = 6,                                   # height in inches
  dpi = 300                                     # resolution
)
#
# >>> PAPER: no outlier, equal-evidence: 
# Appendix D, Table 5 (unconstrained as failsafe; row: Indep., false, n_s = 150, outlier = no; column: Equal, S = 20)
evS_equal_123_unc <- evSyn(object = Param_studies, VCOV = CovMx_studies, 
                           hypotheses = list(H_123 = H_123),
                   type_ev = "equal", 
                   comparison = "unconstrained",
                   order_studies = "ascending") 
evS_equal_123_unc$Cumulative_GORICA_weights
evS_equal_123_unc$Cumulative_GORICA
plot(evS_equal_123_unc)
#png(paste0(name_folder, "/evS_equal_123_unc_S", S, "_Opp.png"))
#plot(evS_equal_123_unc)
#dev.off()
plot_evS_equal_123_unc_indep_N150 <- plot(evS_equal_123_unc)
fileName <- paste0("evS_equal_123_unc_S", S, "_Opp_indep_N150.png")
ggsave(
  filename = fileName,  # file name
  path = name_folder,
  plot = plot_evS_equal_123_unc_indep_N150,     # the plot object
  width = 8,                                    # width in inches
  height = 6,                                   # height in inches
  dpi = 300                                     # resolution
)
#
#
H_123
# vs its compl
#
# >>> PAPER: no outlier, added-evidence: 
# Section 5.2.3 (text: overall weights for n_s = 150) [complement as failsafe]
evS_added_123 <- evSyn(object = Param_studies, VCOV = CovMx_studies, 
                       hypotheses = list(H_123 = H_123),
                   type_ev = "added", 
                   comparison = "complement",
                   order_studies = "ascending") 
evS_added_123$Cumulative_GORICA_weights
evS_added_123$Cumulative_GORICA
plot(evS_added_123)
#png(paste0(name_folder, "/evS_added_123_S", S, "_Opp.png"))
#plot(evS_added_123)
#dev.off()
plot_evS_added_123_compl_indep_N150 <- plot(evS_added_123)
fileName <- paste0("evS_added_123_compl_S", S, "_Opp_indep_N150.png")
ggsave(
  filename = fileName,  # file name
  path = name_folder,
  plot = plot_evS_added_123_compl_indep_N150,   # the plot object
  width = 8,                                    # width in inches
  height = 6,                                   # height in inches
  dpi = 300                                     # resolution
)
#
# >>> PAPER: no outlier, equal-evidence: 
# Section 5.2.3 (text: overall weights for n_s = 150) [complement as failsafe]
evS_equal_123 <- evSyn(object = Param_studies, VCOV = CovMx_studies, 
                   hypotheses = list(H_123 = H_123),
                   type_ev = "equal", 
                   comparison = "complement",
                   order_studies = "ascending") 
evS_equal_123$Cumulative_GORICA_weights
evS_equal_123$Cumulative_GORICA
plot(evS_equal_123) 
#png(paste0(name_folder, "/evS_equal_123_S", S, "_Opp.png"))
#plot(evS_equal_123)
#dev.off()
plot_evS_equal_123_compl_indep_N150 <- plot(evS_equal_123)
fileName <- paste0("evS_equal_123_compl_S", S, "_Opp_indep_N150.png")
ggsave(
  filename = fileName,  # file name
  path = name_folder,
  plot = plot_evS_equal_123_compl_indep_N150,   # the plot object
  width = 8,                                    # width in inches
  height = 6,                                   # height in inches
  dpi = 300                                     # resolution
)

