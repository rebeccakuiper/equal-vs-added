###############################################################################
# process_results_extra_Equal.R
# Does for the EQUAL-evidence results of the extra simulation (based on the review;
# 20 studies, possibly with an outlier study) what analysis/process_results.R does
# for the ADDED-evidence results: performance per iteration (correct = true
# hypothesis selected, support, weight true hypothesis (Weight_Ht), weight
# selected hypothesis (Weight_Hs), grace).
# For checking purposes, the Added results are processed in the same way (this
# reproduces analysis/data/data_to_analyse_2024_09_22.RData).
#
# Working directory: the folder of this script
#   (RStudio: Session > Set Working Directory > To Source File Location).
# The original code of Kuiper & Clapper is in the folder 'Orig code K&C':
proj <- file.path("..", "Orig code K&C", "Sim - extra (review)")
# corrected = FALSE: original results (probit/logit error variances swapped);
# corrected = TRUE : results of rerun_simulation_extra_corrected.R (correct variances)
corrected <- FALSE
suffix <- if (corrected) "_corrected" else ""
###############################################################################

source(file.path(proj, "sim", "functions.R"))  # contains performance()
library(kit)                                    # kit::topn(), used in performance()
dir.create("data", showWarnings = FALSE)

# configurations (as in the original process_results.R)
n_conditions <- 9
# Ht (= true hypothesis) per set; set 5: H1, H2, Hu; set 7: H1 vs its complement.
# performance() changes Ht to Hu (set 5) / complement (set 7) when inverse_ratio == 1
true_hyps <- c(set1 = -1, set2 = -1, set3 = -1, set4 = -1, set5 = 1, set6 = -1, set7 = 1)
grace_hyps <- list(set1 = c(-1), set2 = c(-1), set3 = c(-1),
                   set4 = c(-1), set5 = c(1, 3), set6 = c(-1), set7 = c(1))

dat_list <- list()
for (typeEv in c("added", "equal")) {
  output <- if (corrected) readRDS(file.path("results", paste0("results_", typeEv, "_extra_corrected.RData"))) else
                           readRDS(file.path(proj, "sim", "results", paste0("results_", typeEv, ".RData")))
  rownames(output) <- 1:nrow(output)

  # non-convergence (weights all 0): the same iteration in both (one iteration)
  hypothesis_data <- output[, (n_conditions + 1):ncol(output)]
  noncon <- which(rowSums(hypothesis_data, na.rm = TRUE) < 0.1)
  cat(typeEv, ": non-converged iterations:", noncon, "\n")
  if (length(noncon) > 0) output <- output[-noncon, ]

  metrics <- matrix(NA, nrow(output), 5)
  colnames(metrics) <- c("correct", "support", "Weight_Ht", "Weight_Hs", "grace")
  for (row in 1:nrow(output)) metrics[row, ] <- performance(output[row, ])
  res <- cbind(output, metrics)
  dat_list[[typeEv]] <- res
  saveRDS(res, file.path("data", paste0("data_to_analyse_extra_", typeEv, suffix, ".RData")))
}

# check: the Added results are identical to the ones in the original project folder
if (!corrected) {
orig <- readRDS(file.path(proj, "analysis", "data", "data_to_analyse_2024_09_22.RData"))
cat("Added identical to data_to_analyse_2024_09_22.RData:",
    isTRUE(all.equal(orig, dat_list$added, check.attributes = FALSE)), "\n")
}
