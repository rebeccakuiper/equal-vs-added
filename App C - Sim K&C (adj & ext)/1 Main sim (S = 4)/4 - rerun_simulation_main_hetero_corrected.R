###############################################################################
# 4 - rerun_simulation_main_hetero_corrected.R
# Re-runs the 'hetero' model condition of the main simulation of the Added paper
# (study 1: probit, studies 2-3: linear, study 4: logit) with the CORRECT latent error
# variances for the probit (1) and logit (pi^2/3) models (see
# rerun_functions_hetero_corrected.R), saving both the Added and the Equal weights.
# The same row-specific seeds as in '1 - rerun_simulation_main.R' are used, so the results
# are paired with the re-run with the original (swapped) variances.
#
# By default only N = 25, 50, 75, 100, 150 (set N_levels to all 9 levels for the full set;
# runtime roughly 0.1 sec per iteration per core, more for large N).
#
# Working directory: the folder of this script.
proj <- file.path("..", "Orig code K&C", "Sim - Added - Final")
N_levels <- c(25, 50, 75, 100, 150)          # all: c(25,50,75,100,150,200,300,500,1000)
###############################################################################
library(MASS); library(Matrix); library(matrixStats); library(restriktor)
library(foreach); library(doSNOW)

load(file.path(proj, "sim", "conditions.RData"))
source(file.path(proj, "sim", "functions.R"))          # original functions
source("rerun_functions.R")                            # GoricEvSynEssentials() for current restriktor
source("rerun_functions_hetero_corrected.R")           # overrides estimate_betas(): corrected variances
dir.create("results", showWarnings = FALSE)

rows <- which(conditions$models == "hetero" & conditions$N %in% N_levels)
run_row <- function(rownum){
  set.seed(6164900 + rownum)                           # same seeds as in 1 - rerun_simulation_main.R
  tryCatch(analysis(conditions$sets[rownum], conditions$models[rownum],
                    conditions$true_R2[rownum], conditions$ratios[rownum],
                    conditions$rho[rownum], conditions$N[rownum], specs, rownum),
           error = function(e) list(c(rownum, rep(0, 4)), c(rownum, rep(0, 4))))
}
pad <- function(x) c(x, rep(NA, 5 - length(x)))

cl <- makeCluster(max(1, parallel::detectCores() - 1))
registerDoSNOW(cl)
clusterExport(cl, c("analysis", "print_model", "estimate_betas", "GoricEvSynEssentials",
                    "conditions", "specs", "run_row", "pad"))
pb <- txtProgressBar(min = 1, max = length(rows), style = 3)
opts <- list(progress = function(n) setTxtProgressBar(pb, n))
out <- foreach(rownum = rows, .packages = c("MASS", "Matrix", "matrixStats", "restriktor"),
               .options.snow = opts) %dopar% { r <- run_row(rownum); list(pad(r[[1]]), pad(r[[2]])) }
stopCluster(cl)

for (k in 1:2) {
  typeEv <- c("added", "equal")[k]
  res <- as.data.frame(do.call(rbind, lapply(out, `[[`, k)))
  colnames(res) <- c("rownum", sprintf("H%d", 1:4))
  res <- res[order(res$rownum), ]
  final_res <- cbind(conditions[res$rownum, 2:ncol(conditions)], res[, 2:ncol(res)])
  rownames(final_res) <- NULL
  saveRDS(final_res, file.path("results", paste0("results_", typeEv, "_hetero_corrected.RData")))
}
