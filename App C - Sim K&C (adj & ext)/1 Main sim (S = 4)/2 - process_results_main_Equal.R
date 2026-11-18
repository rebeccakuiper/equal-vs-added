###############################################################################
# 2 - process_results_main_Equal.R
# Does for the EQUAL-evidence results of the main simulation of the Added-evidence
# paper what analysis/process_results.R does for the ADDED-evidence results:
# performance per iteration (correct = true (most parsimonious) hypothesis selected,
# support = relative support selected vs second best, Weight_Ht, Weight_Hs, grace)
# and per unique condition (avg_per_condition).
# The Added results of the re-run are processed in the same way, and compared with
# the original Added results (as a check of the re-run).
#
# Run 1 - rerun_simulation_main.R first (or use the included results/*.RData).
# Working directory: the folder of this script.
proj <- file.path("..", "Orig code K&C", "Sim - Added - Final")
###############################################################################
source(file.path(proj, "sim", "functions.R"))   # contains performance()
library(kit)
dir.create("data", showWarnings = FALSE)

# configurations, Ht per set, etc.: exactly as in the original process_results.R
n_conditions <- 6
true_hyps <- c(set1 = 2, set2 = 2, set3 = 4, set4 = 1, set5 = 1, set6 = 3, set7 = 1, set8 = 1, set9 = 2)
grace_hyps <- list(set1 = c(2), set2 = c(2, 4), set3 = c(4), set4 = c(1), set5 = c(1, 3), set6 = c(3),
                   set7 = c(1), set8 = c(1), set9 = c(2))
true_hyps_ratio4 <- c(set1 = 2, set2 = 4, set3 = 3, set4 = 1, set5 = 3, set6 = 2, set7 = 1, set8 = 2, set9 = 2)
grace_hyps_ratio4 <- list(set1 = c(2), set2 = c(4), set3 = c(3, 4), set4 = c(1), set5 = c(3), set6 = c(2, 3),
                          set7 = c(1), set8 = c(2), set9 = c(2))
# Info hypo sets (set nr in sim (set nr in article)):
# set 1 (1): H10: beta3 = 0, H11: beta3 > 0, H13: beta3 < 0
# set 2 (2): H20: b1=b2=b3, H21: b1>b2>b3, H22: b1<b2,b1>b3, H2u
# set 3 (3): H30: b1=b2=b3, H31: b1<b2>b3, H32: b1<b2<b3, H3u
# set 4 (4): H11, H13;  set 5 (5): H21, H22, H2u;  set 6 (6): H31, H32, H3u
# set 8 (7): H21 vs its complement;  set 9 (8): H31 vs its complement

for (typeEv in c("added", "equal")) {
  output <- readRDS(file.path("results", paste0("results_", typeEv, "_rerun.RData")))
  rownames(output) <- 1:nrow(output)
  hypothesis_data <- output[, (n_conditions + 1):ncol(output)]
  noncon <- which(rowSums(hypothesis_data, na.rm = TRUE) < 0.999)
  cat(typeEv, ": number of non-converged iterations:", length(noncon), "\n")

  metrics <- matrix(NA, nrow(output), 5)
  colnames(metrics) <- c("correct", "support", "Weight_Ht", "Weight_Hs", "grace")
  for (row in 1:nrow(output)) metrics[row, ] <- performance(output[row, ])   # takes a few minutes
  res <- cbind(output, metrics)
  saveRDS(res, file.path("data", paste0("data_to_analyse_main_", typeEv, ".RData")))

  # performance per unique condition (median weights and THR), as in the original
  niter <- 100
  Sequence <- seq(1, nrow(res), niter)
  perf_per_condition <- matrix(NA, length(Sequence), 8)
  colnames(perf_per_condition) <- c(sprintf("H%d", 1:4), "support", "Weight_Ht", "Weight_Hs", "THR")
  for (i in seq_along(Sequence)) {
    Subset <- res[Sequence[i]:(Sequence[i] + niter - 1), ]
    perf_per_condition[i, ] <- c(apply(Subset[, c(sprintf("H%d", 1:4), "support", "Weight_Ht", "Weight_Hs")], 2, median),
                                 sum(Subset$correct) / nrow(Subset))
  }
  res_per_condition <- cbind(unique(output[, 1:n_conditions]), perf_per_condition)
  saveRDS(res_per_condition, file.path("data", paste0("avg_per_condition_main_", typeEv, ".RData")))
}

#### ---- Check of the re-run: Added (re-run) versus Added (original) ----
orig  <- readRDS(file.path(proj, "analysis", "data", "data_to_analyse.RData"))
rerun <- readRDS(file.path("data", "data_to_analyse_main_added.RData"))
key <- function(d) do.call(paste, d[, 1:n_conditions])
orig <- orig[key(orig) %in% key(rerun), ]   # same conditions (only relevant if the re-run is incomplete)
cmp <- function(f) merge(aggregate(correct ~ get(f), data = orig, FUN = mean),
                         aggregate(correct ~ get(f), data = rerun, FUN = mean), by = "get(f)",
                         suffixes = c("_original", "_rerun"))
cat("\nTHR Added: original run versus re-run\n")
cat("overall:", round(mean(orig$correct), 3), "vs", round(mean(rerun$correct), 3), "\n")
for (f in c("sets", "N", "true_R2", "rho", "models", "ratios")) { x <- cmp(f); names(x)[1] <- f; print(x, digits = 3) }
