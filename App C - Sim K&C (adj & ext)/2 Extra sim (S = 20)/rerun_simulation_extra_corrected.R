###############################################################################
# rerun_simulation_extra_corrected.R
# Re-runs the extra simulation (based on the review; 20 studies, all 'hetero': probit,
# linear, and logit studies) with the CORRECT latent error variances for the probit (1)
# and logit (pi^2/3) models (in sim/functions.R these were swapped), saving both the
# Added and the Equal weights. Uses row-specific seeds (6164900 + row number), so the
# results are reproducible (also in parallel).
# Output: results/results_{added,equal}_extra_corrected.RData (same format as the
# original sim/results/results_*.RData), to be processed with process_results_extra_Equal.R
# (set corrected <- TRUE there).
# Runtime: roughly 0.2-0.3 sec per iteration per core (21,600 iterations).
#
# Working directory: the folder of this script.
proj  <- file.path("..", "Orig code K&C", "Sim - extra (review)")
fdir1 <- file.path("..", "1 Main sim (S = 4)")
###############################################################################
library(MASS); library(Matrix); library(matrixStats); library(restriktor)
library(foreach); library(doSNOW)

load(file.path(proj, "sim", "conditions-2024-09-22.RData"))   # conditions and specs
source(file.path(proj, "sim", "functions.R"))                 # original functions
source(file.path(fdir1, "rerun_functions.R"))                 # GoricEvSynEssentials() for current restriktor
source(file.path(fdir1, "rerun_functions_hetero_corrected.R"))# corrected estimate_betas()
dir.create("results", showWarnings = FALSE)

run_row <- function(rownum){
  set.seed(6164900 + rownum)
  tryCatch(analysis(conditions$sets[rownum], conditions$n_studies[rownum], conditions$models[rownum],
                    conditions$true_R2[rownum], conditions$ratios[rownum], conditions$rho[rownum],
                    conditions$N[rownum], conditions$inverse_ratio[rownum], conditions$outlier_study[rownum],
                    specs, rownum),
           error = function(e) list(c(rownum, rep(0, 4)), c(rownum, rep(0, 4))))
}
pad <- function(x) c(x, rep(NA, 5 - length(x)))

cl <- makeCluster(max(1, parallel::detectCores() - 1))
registerDoSNOW(cl)
clusterExport(cl, c("analysis", "print_model", "estimate_betas", "GoricEvSynEssentials",
                    "get_n_study_proportion", "conditions", "specs", "run_row", "pad"))
pb <- txtProgressBar(min = 1, max = nrow(conditions), style = 3)
opts <- list(progress = function(n) setTxtProgressBar(pb, n))
out <- foreach(rownum = 1:nrow(conditions), .packages = c("MASS", "Matrix", "matrixStats", "restriktor"),
               .options.snow = opts) %dopar% { r <- run_row(rownum); list(pad(r[[1]]), pad(r[[2]])) }
stopCluster(cl)

for (k in 1:2) {
  typeEv <- c("added", "equal")[k]
  res <- as.data.frame(do.call(rbind, lapply(out, `[[`, k)))
  colnames(res) <- c("rownum", sprintf("H%d", 1:4))
  res <- res[order(res$rownum), ]
  final_res <- cbind(conditions[res$rownum, 2:ncol(conditions)], res[, 2:ncol(res)])
  rownames(final_res) <- NULL
  saveRDS(final_res, file.path("results", paste0("results_", typeEv, "_extra_corrected.RData")))
}
