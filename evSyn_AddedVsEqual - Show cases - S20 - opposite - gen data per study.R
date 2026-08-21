
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

if (!require("ggplot2")) install.packages("metafor")
library("ggplot2")
# ggsave

###


# NU NIET S keer zelfde data, maar S keer gesampeld.


# Show case #

# Aggregate S studies 
# regression with 3 predictors, 
# each study has N observations (n_s)


# Specs situation
S <- 20 # 2 # 4 # 20
n_s <- 85 # 25, 50, 75, 100, 150, 200, 300, 500, 1,000
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
  if(n_s %% 2 == 0){ # even
    e <- mvrnorm(n = n_s/2, 0, var.e) #, empirical=TRUE) # sample error terms with exact that res. var.
    #epsilon <- c(e, -e)
    #var(e); var(-e)
    e_even <- c(e, -e)
    #var(e_even)
    epsilon <- e_even*sqrt(var.e/var(e_even))
    #var(epsilon); var.e
    x <- as.data.frame(mvrnorm(n = n_s/2, mu, sigma)) #, empirical=TRUE)) # Create predictor variables with exact that sigma
    #var(x)
    x_even <- rbind(x, x)
    #var(x_even)
    sample <- x_even*sqrt(diag(sigma/var(x_even))) 
    #var(sample); sigma
  }else{ # odd
    e <- mvrnorm(n = (n_s-1)/2, 0, var.e) #, empirical=TRUE) # sample error terms with exact that res. var.
    e_odd <- c(e , 0, -e)
    epsilon <- e_odd*sqrt(var.e/var(e_odd))
    x <- as.data.frame(mvrnorm(n = (n_s-1)/2, mu, sigma)) #, empirical=TRUE)) # Create predictor variables with exact that sigma
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
#The order-restricted hypothesis ‘H_123’ has less support 
#(namely,  0.200 / 0.800 < 1 times more support) than its complement. 
#That is, the complement has 4.00 times more support than ‘H_123’
gorica_lm$ratio.lw
# LLs not equal (ratio is approx. 9.015397)


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
plot_evS_added_123_unc_gendata <- plot(evS_added_123_unc)
fileName <- paste0("evS_added_123_unc_S", S, "_Opp_gendata.png")
ggsave(
  filename = fileName,  # file name
  path = name_folder,
  plot = plot_evS_added_123_unc_gendata,             # the plot object
  width = 8,                                    # width in inches
  height = 6,                                   # height in inches
  dpi = 300                                     # resolution
)
#
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
plot_evS_equal_123_unc_gendata <- plot(evS_equal_123_unc)
fileName <- paste0("evS_equal_123_unc_S", S, "_Opp_gendata.png")
ggsave(
  filename = fileName,  # file name
  path = name_folder,
  plot = plot_evS_equal_123_unc_gendata,             # the plot object
  width = 8,                                    # width in inches
  height = 6,                                   # height in inches
  dpi = 300                                     # resolution
)
#
#
H_123
# vs its compl
#
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
plot_evS_added_123_compl_gendata <- plot(evS_added_123)
fileName <- paste0("evS_added_123_compl_S", S, "_Opp_gendata.png")
ggsave(
  filename = fileName,  # file name
  path = name_folder,
  plot = plot_evS_added_123_compl_gendata,             # the plot object
  width = 8,                                    # width in inches
  height = 6,                                   # height in inches
  dpi = 300                                     # resolution
)
#
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
plot_evS_equal_123_compl_gendata <- plot(evS_equal_123)
fileName <- paste0("evS_equal_123_compl_S", S, "_Opp_gendata.png")
ggsave(
  filename = fileName,  # file name
  path = name_folder,
  plot = plot_evS_equal_123_compl_gendata,             # the plot object
  width = 8,                                    # width in inches
  height = 6,                                   # height in inches
  dpi = 300                                     # resolution
)


###################################################################################

# Not done yet!!!

# ## Show case 1b: 1 outlier (last study) ##
# 
# # Pas estimates aan (NB VCOV blijft gelijk):
# Param_studies_out <- Param_studies
# Param_studies_out[[S]][2:4] <- rev(Param_studies_out[[S]][2:4])
# #Param_studies_out
# #
# #
# H_123
# # vs unc
# #
# evS_added_123_unc_out <- evSyn(object = Param_studies_out, VCOV = CovMx_studies, 
#                            hypotheses = list(H_123 = H_123),
#                            type_ev = "added", 
#                            comparison = "unconstrained") 
# evS_added_123_unc_out$Cumulative_GORICA_weights
# evS_added_123_unc_out$Cumulative_GORICA
# plot(evS_added_123_unc_out)
# #png(paste0(name_folder, "/evS_added_123_unc_out_S", S, "_Opp.png"))
# #plot(evS_added_123_unc_out)
# #dev.off()
# plot_evS_added_123_unc_out <- plot(evS_added_123_unc_out)
# fileName <- paste0("evS_added_123_unc_out_S", S, "_Opp.png")
# ggsave(
#   filename = fileName,  # file name
#   path = name_folder,
#   plot = plot_evS_added_123_unc_out,             # the plot object
#   width = 8,                                    # width in inches
#   height = 6,                                   # height in inches
#   dpi = 300                                     # resolution
# )
# #
# evS_equal_123_unc_out <- evSyn(object = Param_studies_out, VCOV = CovMx_studies, 
#                            hypotheses = list(H_123 = H_123),
#                            type_ev = "equal", 
#                            comparison = "unconstrained") 
# evS_equal_123_unc_out$Cumulative_GORICA_weights
# evS_equal_123_unc_out$Cumulative_GORICA
# plot(evS_equal_123_unc_out)
# #png(paste0(name_folder, "/evS_equal_123_unc_out_S", S, "_Opp.png"))
# #plot(evS_equal_123_unc_out)
# #dev.off()
# plot_evS_equal_123_unc_out <- plot(evS_equal_123_unc_out)
# fileName <- paste0("evS_equal_123_unc_out_S", S, "_Opp.png")
# ggsave(
#   filename = fileName,  # file name
#   path = name_folder,
#   plot = plot_evS_equal_123_unc_out,             # the plot object
#   width = 8,                                    # width in inches
#   height = 6,                                   # height in inches
#   dpi = 300                                     # resolution
# )
# #
# # Check
# evS_equal_123_unc_out_asc <- evSyn(object = Param_studies_out, VCOV = CovMx_studies, 
#                                hypotheses = list(H_123 = H_123),
#                                type_ev = "equal", 
#                                comparison = "unconstrained",
#                                order_studies = "ascending") 
# summary(evS_equal_123_unc_out_asc)
# #evS_equal_123_unc_out_asc$Cumulative_GORICA_weights
# #evS_equal_123_unc_out_asc$Cumulative_GORICA
# plot(evS_equal_123_unc_out_asc)
# # Cum. GORICA values voor H123 nu -3.752 & -14.521
# # Moet zijn:
# # Cum. LL moet zijn (4.709 + 5.384) = 10.093
# # Cum PT weight (2.833+2.833)/2=2.833
# # -2 * 10.094 + 2 * 2.833 = -14.522 -- klopt dus!
# #
# # Check op added dan ook, daar gaat het ook goed:
# evS_added_123_unc_out_asc <- evSyn(object = Param_studies_out, VCOV = CovMx_studies, 
#                                    hypotheses = list(H_123 = H_123),
#                                    type_ev = "added", 
#                                    comparison = "unconstrained",
#                                    order_studies = "ascending") 
# summary(evS_added_123_unc_out_asc)
# # Cum. GORICA values voor H123 nu idd -3.752 & (-3.752 + -5.102) = -8.854
# #
# #
# H_123
# # vs its compl
# #
# evS_added_123_out <- evSyn(object = Param_studies_out, VCOV = CovMx_studies, 
#                        hypotheses = list(H_123 = H_123),
#                        type_ev = "added", 
#                        comparison = "complement") 
# evS_added_123_out$Cumulative_GORICA_weights
# evS_added_123_out$Cumulative_GORICA
# plot(evS_added_123_out)
# #png(paste0(name_folder, "/evS_added_123_out_S", S, "_Opp.png"))
# #plot(evS_added_123_out)
# #dev.off()
# plot_evS_added_123_compl_out <- plot(evS_added_123_out)
# fileName <- paste0("evS_added_123_compl_out_S", S, "_Opp.png")
# ggsave(
#   filename = fileName,  # file name
#   path = name_folder,
#   plot = plot_evS_added_123_compl_out,             # the plot object
#   width = 8,                                    # width in inches
#   height = 6,                                   # height in inches
#   dpi = 300                                     # resolution
# )
# #
# evS_equal_123_out <- evSyn(object = Param_studies_out, VCOV = CovMx_studies, 
#                        hypotheses = list(H_123 = H_123),
#                        type_ev = "equal", 
#                        comparison = "complement") 
# evS_equal_123_out$Cumulative_GORICA_weights
# evS_equal_123_out$Cumulative_GORICA
# plot(evS_equal_123_out) 
# #png(paste0(name_folder, "/evS_equal_123_out_S", S, "_Opp.png"))
# #plot(evS_equal_123_out)
# #dev.off()
# plot_evS_equal_123_compl_out <- plot(evS_equal_123_out)
# fileName <- paste0("evS_equal_123_compl_out_S", S, "_Opp.png")
# ggsave(
#   filename = fileName,  # file name
#   path = name_folder,
#   plot = plot_evS_equal_123_compl_out,             # the plot object
#   width = 8,                                    # width in inches
#   height = 6,                                   # height in inches
#   dpi = 300                                     # resolution
# )
# 
