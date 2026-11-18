###############################################################################
# 1 - rerun_simulation_main.R
# The original simulation of the Added-evidence paper (folder 'Sim - Added - Final')
# computed both the added- and the equal-evidence weights (see sim/simulation.R and
# GoricEvSynEssentials() in sim/functions.R), but only the added-evidence results were
# kept (sim/results/results_added_*.txt and results_added.RData; there is no
# results_equal folder). Therefore, this script re-runs the simulation with exactly
# the same design and code (sim/conditions.RData and sim/functions.R), keeping BOTH
# the added- and equal-evidence weights, which are thus based on the same data sets.
#
# Changes w.r.t. the original sim/simulation.R:
#  - GoricEvSynEssentials() is replaced by the version in rerun_functions.R, since
#    the current restriktor needs the hypotheses as a list in goric(); the
#    computations (LL, PT, added and equal weights) are the same.
#  - Each row (iteration) gets its own seed (6164900 + row number), such that the
#    results are reproducible, also when run in parallel (in the original run, the
#    seed was not passed to the parallel workers). Hence, the generated data differ
#    from the original run, but the Added results can be compared to the original ones
#    (see 2 - process_results_main_Equal.R).
#
# Runtime: about 0.075 sec per iteration per core (388,800 iterations), so roughly
#   1-1.5 hours with 8 cores. The results are included:
#   results/results_added_rerun.RData and results/results_equal_rerun.RData; these cover
#   N = 25, 50, 75, 100, 150 (rows 1-216000 of conditions), not N = 200-1000.
#   NB: the hetero condition here uses the original (swapped) probit/logit error variances;
#   see 4 - rerun_simulation_main_hetero_corrected.R for the corrected version.
#
# Working directory: the folder of this script.
proj <- file.path("..", "Orig code K&C", "Sim - Added - Final")
###############################################################################
library(MASS); library(Matrix); library(matrixStats); library(restriktor)
library(foreach); library(doSNOW)

load(file.path(proj, "sim", "conditions.RData"))    # conditions and specs
source(file.path(proj, "sim", "functions.R"))       # original simulation functions
source("rerun_functions.R")                         # GoricEvSynEssentials() for current restriktor
dir.create("results", showWarnings = FALSE)

run_row <- function(rownum){
  set.seed(6164900 + rownum)                        # row-specific seed
  tryCatch(analysis(conditions$sets[rownum], conditions$models[rownum],
                    conditions$true_R2[rownum], conditions$ratios[rownum],
                    conditions$rho[rownum], conditions$N[rownum], specs, rownum),
           error = function(e) list(c(rownum, rep(0, 4)), c(rownum, rep(0, 4))))
}
pad <- function(x) c(x, rep(NA, 5 - length(x)))     # rownum + max. 4 weights

cl <- makeCluster(max(1, parallel::detectCores() - 1))
registerDoSNOW(cl)
clusterExport(cl, c("analysis", "print_model", "estimate_betas", "GoricEvSynEssentials",
                    "conditions", "specs", "run_row", "pad"))
pb <- txtProgressBar(min = 1, max = nrow(conditions), style = 3)
opts <- list(progress = function(n) setTxtProgressBar(pb, n))
t0 <- Sys.time()
out <- foreach(rownum = 1:nrow(conditions), .packages = c("MASS", "Matrix", "matrixStats", "restriktor"),
               .options.snow = opts) %dopar% {
  r <- run_row(rownum); list(pad(r[[1]]), pad(r[[2]]))
}
stopCluster(cl)
print(Sys.time() - t0)

# Combine output + conditions (as in the original simulation.R)
for (k in 1:2) {
  typeEv <- c("added", "equal")[k]
  res <- as.data.frame(do.call(rbind, lapply(out, `[[`, k)))
  colnames(res) <- c("rownum", sprintf("H%d", 1:4))
  res <- res[order(res$rownum), ]
  final_res <- cbind(conditions[, 2:ncol(conditions)], res[, 2:ncol(res)])
  saveRDS(final_res, file.path("results", paste0("results_", typeEv, "_rerun.RData")))
}
