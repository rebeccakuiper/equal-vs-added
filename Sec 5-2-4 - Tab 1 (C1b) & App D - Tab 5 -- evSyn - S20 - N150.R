###############################################################################
# Data: S = 20 times the same dataset, n_s = 150 per study; H_123 TRUE (x1 > x2 > x3).
# Also: 1 outlier study (estimates of the last study reversed).
# Used in:
#   - Appendix D, Table 5 (unconstrained as failsafe)
#   - Section 5.2.4 (text: overall weights for n_s = 150)
#   - Table 1, Case C1b (Section 5.2.4, Leave One Study Out)
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


name_folder <- paste0('images_S', S, "_N", n_s)
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
set.seed(123)
mu <- rep(0, n.coef)
if(n_s %% 2 == 0){ # even
  e <- mvrnorm(n = n_s/2, 0, var.e, empirical=TRUE) # sample error terms with exact that res. var.
  #epsilon <- c(e, -e)
  #var(e); var(-e)
  e_even <- c(e, -e)
  #var(e_even)
  epsilon <- e_even*sqrt(var.e/var(e_even))
  #var(epsilon); var.e
  x <- as.data.frame(mvrnorm(n = n_s/2, mu, sigma, empirical=TRUE)) # Create predictor variables with exact that sigma
  #var(x)
  x_even <- rbind(x, x)
  #var(x_even)
  sample <- x_even*sqrt(diag(sigma/var(x_even)))
  #var(sample); sigma
}else{ # odd
  e <- mvrnorm(n = (n_s-1)/2, 0, var.e, empirical=TRUE) # sample error terms with exact that res. var.
  e_odd <- c(e , 0, -e)
  epsilon <- e_odd*sqrt(var.e/var(e_odd))
  x <- as.data.frame(mvrnorm(n = (n_s-1)/2, mu, sigma, empirical=TRUE)) # Create predictor variables with exact that sigma
  x_odd <- rbind(x, rep(0, n.coef), x)
  sample <- x_odd*sqrt(diag(sigma/var(x_odd)))
  #var(epsilon); var.e
  #mean(epsilon)
  #var(sample); sigma
}
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


###


# Informative hypothesis
H_123 <- 'x1 > x2 > x3'


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
#The order-restricted hypothesis ‘H_123’ has 3.12 times more support than its complement.
gorica_lm$ratio.lw
# LLs not equal (ratio is approx. 1.35548)


## Combined:


## Show case 1a: no outliers ##

# Make lists of param. estimates and their cov. matrices from each study
rep_est_text <- paste0(rep(" est", (S)), collapse=",")
est_text <- noquote(paste("list(", rep_est_text, ")"))
Param_studies <- eval(parse(text = est_text))
#Param_studies
rep_vcov_text <- paste0(rep(" VCOV", (S)), collapse=",")
vcov_text <- noquote(paste("list(", rep_vcov_text, ")"))
CovMx_studies <- eval(parse(text = vcov_text))
#
#
H_123
# vs unc
#
# >>> PAPER: no outlier, added-evidence: 
# Appendix D, Table 5 (unconstrained as failsafe; row: Same, true, n_s = 150, outlier = no; column: Added, S = 20); 
# the study-specific weight (column S = 1) is also taken from this output
evS_added_123_unc <- evSyn(object = Param_studies, VCOV = CovMx_studies, 
                           hypotheses = list(H_123 = H_123),
                   type_ev = "added", # default
                   comparison = "unconstrained") 
evS_added_123_unc$Cumulative_GORICA_weights
evS_added_123_unc$Cumulative_GORICA
plot(evS_added_123_unc)
#png(paste0(name_folder, "/evS_added_123_unc_S", S, "_N150.png"))
#plot(evS_added_123_unc)
#dev.off()
plot_evS_added_123_unc <- plot(evS_added_123_unc)
fileName <- paste0("evS_added_123_unc_S", S, "_N150.png")
ggsave(
  filename = fileName,  # file name
  path = name_folder,
  plot = plot_evS_added_123_unc,                # the plot object
  width = 8,                                    # width in inches
  height = 6,                                   # height in inches
  dpi = 300                                     # resolution
)
#
# >>> PAPER: no outlier, equal-evidence: 
# Appendix D, Table 5 (unconstrained as failsafe; row: Same, true, n_s = 150, outlier = no; column: Equal, S = 20); 
# the study-specific weight (column S = 1) is also taken from this output
evS_equal_123_unc <- evSyn(object = Param_studies, VCOV = CovMx_studies, 
                           hypotheses = list(H_123 = H_123),
                   type_ev = "equal", 
                   comparison = "unconstrained") 
evS_equal_123_unc$Cumulative_GORICA_weights
evS_equal_123_unc$Cumulative_GORICA
plot(evS_equal_123_unc)
#png(paste0(name_folder, "/evS_equal_123_unc_S", S, "_N150.png"))
#plot(evS_equal_123_unc)
#dev.off()
plot_evS_equal_123_unc <- plot(evS_equal_123_unc)
fileName <- paste0("evS_equal_123_unc_S", S, "_N150.png")
ggsave(
  filename = fileName,  # file name
  path = name_folder,
  plot = plot_evS_equal_123_unc,                # the plot object
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
# not shown in the paper (H_123 true, n_s = 150; only used for the outlier case, Section 5.2.4) [complement as failsafe]
evS_added_123 <- evSyn(object = Param_studies, VCOV = CovMx_studies, 
                       hypotheses = list(H_123 = H_123),
                   type_ev = "added", 
                   comparison = "complement") 
evS_added_123$Cumulative_GORICA_weights
evS_added_123$Cumulative_GORICA
plot(evS_added_123)
#png(paste0(name_folder, "/evS_added_123_S", S, "_N150.png"))
#plot(evS_added_123)
#dev.off()
plot_evS_added_123_compl <- plot(evS_added_123)
fileName <- paste0("evS_added_123_compl_S", S, "_N150.png")
ggsave(
  filename = fileName,  # file name
  path = name_folder,
  plot = plot_evS_added_123_compl,              # the plot object
  width = 8,                                    # width in inches
  height = 6,                                   # height in inches
  dpi = 300                                     # resolution
)
#
# >>> PAPER: no outlier, equal-evidence: 
# not shown in the paper (H_123 true, n_s = 150; only used for the outlier case, Section 5.2.4) [complement as failsafe]
evS_equal_123 <- evSyn(object = Param_studies, VCOV = CovMx_studies, 
                   hypotheses = list(H_123 = H_123),
                   type_ev = "equal", 
                   comparison = "complement") 
evS_equal_123$Cumulative_GORICA_weights
evS_equal_123$Cumulative_GORICA
plot(evS_equal_123) 
#png(paste0(name_folder, "/evS_equal_123_S", S, "_N150.png"))
#plot(evS_equal_123)
#dev.off()
plot_evS_equal_123_compl <- plot(evS_equal_123)
fileName <- paste0("evS_equal_123_compl_S", S, "_N150.png")
ggsave(
  filename = fileName,  # file name
  path = name_folder,
  plot = plot_evS_equal_123_compl,              # the plot object
  width = 8,                                    # width in inches
  height = 6,                                   # height in inches
  dpi = 300                                     # resolution
)


## Show case 1b: 1 outlier (last study) ##

# Adjust estimates for last study, so for Study S 
# (Note that VCOV stays the same)
Param_studies_out <- Param_studies
Param_studies_out[[S]][2:4] <- rev(Param_studies_out[[S]][2:4])
#Param_studies_out
#
#
H_123
# vs unc
#
# >>> PAPER: with outlier study, added-evidence: 
# Appendix D, Table 5 (unconstrained as failsafe; row: Same, true, n_s = 150, outlier = yes; column: Added, S = 20)
evS_added_123_unc_out <- evSyn(object = Param_studies_out, VCOV = CovMx_studies, 
                           hypotheses = list(H_123 = H_123),
                           type_ev = "added", 
                           comparison = "unconstrained") 
evS_added_123_unc_out$Cumulative_GORICA_weights
evS_added_123_unc_out$Cumulative_GORICA
plot(evS_added_123_unc_out)
#png(paste0(name_folder, "/evS_added_123_unc_out_S", S, "_N150.png"))
#plot(evS_added_123_unc_out)
#dev.off()
plot_evS_added_123_unc_out <- plot(evS_added_123_unc_out)
fileName <- paste0("evS_added_123_unc_out_S", S, "_N150.png")
ggsave(
  filename = fileName,  # file name
  path = name_folder,
  plot = plot_evS_added_123_unc_out,            # the plot object
  width = 8,                                    # width in inches
  height = 6,                                   # height in inches
  dpi = 300                                     # resolution
)
#
# >>> PAPER: with outlier study, equal-evidence: 
# Appendix D, Table 5 (unconstrained as failsafe; row: Same, true, n_s = 150, outlier = yes; column: Equal, S = 20)
evS_equal_123_unc_out <- evSyn(object = Param_studies_out, VCOV = CovMx_studies, 
                           hypotheses = list(H_123 = H_123),
                           type_ev = "equal", 
                           comparison = "unconstrained") 
evS_equal_123_unc_out$Cumulative_GORICA_weights
evS_equal_123_unc_out$Cumulative_GORICA
plot(evS_equal_123_unc_out)
#png(paste0(name_folder, "/evS_equal_123_unc_out_S", S, "_N150.png"))
#plot(evS_equal_123_unc_out)
#dev.off()
plot_evS_equal_123_unc_out <- plot(evS_equal_123_unc_out)
fileName <- paste0("evS_equal_123_unc_out_S", S, "_N150.png")
ggsave(
  filename = fileName,  # file name
  path = name_folder,
  plot = plot_evS_equal_123_unc_out,            # the plot object
  width = 8,                                    # width in inches
  height = 6,                                   # height in inches
  dpi = 300                                     # resolution
)
#
# Check
# >>> PAPER: with outlier study, equal-evidence: 
# not used in the paper (check, with the studies in ascending order)
evS_equal_123_unc_out_asc <- evSyn(object = Param_studies_out, VCOV = CovMx_studies, 
                               hypotheses = list(H_123 = H_123),
                               type_ev = "equal", 
                               comparison = "unconstrained",
                               order_studies = "ascending") 
summary(evS_equal_123_unc_out_asc)
#evS_equal_123_unc_out_asc$Cumulative_GORICA_weights
#evS_equal_123_unc_out_asc$Cumulative_GORICA
plot(evS_equal_123_unc_out_asc)
# Cum. GORICA values voor H123 nu -3.752 & -14.521
# Moet zijn:
# Cum. LL moet zijn (4.709 + 5.384) = 10.093
# Cum PT weight (2.833+2.833)/2=2.833
# -2 * 10.094 + 2 * 2.833 = -14.522 -- klopt dus!
#
# Check op added dan ook, daar gaat het ook goed:
# >>> PAPER: with outlier study, added-evidence: 
# not used in the paper (check, with the studies in ascending order)
evS_added_123_unc_out_asc <- evSyn(object = Param_studies_out, VCOV = CovMx_studies, 
                                   hypotheses = list(H_123 = H_123),
                                   type_ev = "added", 
                                   comparison = "unconstrained",
                                   order_studies = "ascending") 
summary(evS_added_123_unc_out_asc)
# Cum. GORICA values voor H123 nu idd -3.752 & (-3.752 + -5.102) = -8.854
#
#
H_123
# vs its compl
#
# >>> PAPER: with outlier study, added-evidence: 
# Section 5.2.4 (text: overall weights for n_s = 150) [complement as failsafe]
evS_added_123_out <- evSyn(object = Param_studies_out, VCOV = CovMx_studies, 
                       hypotheses = list(H_123 = H_123),
                       type_ev = "added", 
                       comparison = "complement") 
evS_added_123_out$Cumulative_GORICA_weights
evS_added_123_out$Cumulative_GORICA
plot(evS_added_123_out)
#png(paste0(name_folder, "/evS_added_123_out_S", S, "_N150.png"))
#plot(evS_added_123_out)
#dev.off()
plot_evS_added_123_compl_out <- plot(evS_added_123_out)
fileName <- paste0("evS_added_123_compl_out_S", S, "_N150.png")
ggsave(
  filename = fileName,  # file name
  path = name_folder,
  plot = plot_evS_added_123_compl_out,          # the plot object
  width = 8,                                    # width in inches
  height = 6,                                   # height in inches
  dpi = 300                                     # resolution
)
#
# >>> PAPER: with outlier study, equal-evidence: 
# Section 5.2.4 (text: overall weights for n_s = 150) [complement as failsafe]
evS_equal_123_out <- evSyn(object = Param_studies_out, VCOV = CovMx_studies, 
                       hypotheses = list(H_123 = H_123),
                       type_ev = "equal", 
                       comparison = "complement") 
evS_equal_123_out$Cumulative_GORICA_weights
evS_equal_123_out$Cumulative_GORICA
plot(evS_equal_123_out) 
#png(paste0(name_folder, "/evS_equal_123_out_S", S, "_N150.png"))
#plot(evS_equal_123_out)
#dev.off()
plot_evS_equal_123_compl_out <- plot(evS_equal_123_out)
fileName <- paste0("evS_equal_123_compl_out_S", S, "_N150.png")
ggsave(
  filename = fileName,  # file name
  path = name_folder,
  plot = plot_evS_equal_123_compl_out,          # the plot object
  width = 8,                                    # width in inches
  height = 6,                                   # height in inches
  dpi = 300                                     # resolution
)


#Bij outlier cases ook leave1studyout doen:

# >>> PAPER: checked; not reported in the paper (no meaningful influence of the outlier study in this case)
leave1studyout(evS_added_123_unc_out)
#Leave Study 18 out: -139.1849     -97.28485
#Leave Study 19 out: -139.1849     -97.28485
#Leave Study 20 out: -141.6182     -97.28485

# >>> PAPER: Table 1, Case C1b (Section 5.2.4, Leave One Study Out)
leave1studyout(evS_equal_123_unc_out)
# Leave Study 18 out: -241.1849     -241.2849
# Leave Study 19 out: -241.1849     -241.2849
# Leave Study 20 out: -243.6182     -241.2849
# Leave Study 18 out: 0.4875138     0.5124862
# Leave Study 19 out: 0.4875138     0.5124862
# Leave Study 20 out: 0.7625420     0.2374580
# Leave Study 18 out: "Unconstrained"
# Leave Study 19 out: "Unconstrained"
# Leave Study 20 out: "H_123"
0.4875138/0.5124862 # 0.9512721
0.7625420/0.2374580 # 3.211271
0.5124862/0.4875138 # 1.051224

# >>> PAPER: checked; not reported in the paper (no meaningful influence of the outlier study in this case)
leave1studyout(evS_added_123_out)
# Leave Study 18 out: -139.1849  -99.00192
# Leave Study 19 out: -139.1849  -99.00192
# Leave Study 20 out: -141.6182  -98.39361

# >>> PAPER: checked; not reported in the paper (no meaningful influence of the outlier study in this case)
leave1studyout(evS_equal_123_out)
# Leave Study 18 out: -241.1849  -231.0019
# Leave Study 19 out: -241.1849  -231.0019
# Leave Study 20 out: -243.6182  -230.3936
# Leave Study 18 out: 0.9938888 0.006111152
# Leave Study 19 out: 0.9938888 0.006111152
# Leave Study 20 out: 0.9986580 0.001341952
