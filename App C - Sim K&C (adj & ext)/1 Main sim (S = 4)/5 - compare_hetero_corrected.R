###############################################################################
# 5 - compare_hetero_corrected.R
# Processes the results of '4 - rerun_simulation_main_hetero_corrected.R' (hetero condition
# with the correct probit/logit error variances) as in process_results.R, and compares them
# with the hetero results with the original (swapped) variances from the re-run
# (same rows and seeds), for both Added and Equal.
# Run 2 - process_results_main_Equal.R and 4 - ... first.
# Working directory: the folder of this script.
proj <- file.path("..", "Orig code K&C", "Sim - Added - Final")
###############################################################################
source(file.path(proj, "sim", "functions.R"))   # performance()
library(kit); library(knitr)
dir.create("images", showWarnings = FALSE)
n_conditions <- 6
true_hyps <- c(set1 = 2, set2 = 2, set3 = 4, set4 = 1, set5 = 1, set6 = 3, set7 = 1, set8 = 1, set9 = 2)
grace_hyps <- list(set1 = c(2), set2 = c(2, 4), set3 = c(4), set4 = c(1), set5 = c(1, 3), set6 = c(3),
                   set7 = c(1), set8 = c(1), set9 = c(2))
true_hyps_ratio4 <- c(set1 = 2, set2 = 4, set3 = 3, set4 = 1, set5 = 3, set6 = 2, set7 = 1, set8 = 2, set9 = 2)
grace_hyps_ratio4 <- list(set1 = c(2), set2 = c(4), set3 = c(3, 4), set4 = c(1), set5 = c(3), set6 = c(2, 3),
                          set7 = c(1), set8 = c(2), set9 = c(2))

key <- function(d) do.call(paste, d[, 1:n_conditions])
D <- list()
for (typeEv in c("added", "equal")) {
  output <- readRDS(file.path("results", paste0("results_", typeEv, "_hetero_corrected.RData")))
  metrics <- t(sapply(1:nrow(output), function(r) performance(output[r, ])))
  colnames(metrics) <- c("correct", "support", "Weight_Ht", "Weight_Hs", "grace")
  corr <- cbind(output, metrics)
  saveRDS(corr, file.path("data", paste0("data_to_analyse_main_hetero_corrected_", typeEv, ".RData")))
  orig <- readRDS(file.path("data", paste0("data_to_analyse_main_", typeEv, ".RData")))
  orig <- orig[orig$models == "hetero" & key(orig) %in% key(corr), ]
  D[[typeEv]] <- rbind(cbind(orig[, names(corr)], variance = "original (swapped)"),
                       cbind(corr, variance = "corrected"))
}

#### ---- Tables: THR and median weight true hypothesis, original vs corrected ----
summ <- function(f) {
  out <- NULL
  for (typeEv in names(D)) {
    x <- do.call(rbind, lapply(split(D[[typeEv]], list(D[[typeEv]][[f]], D[[typeEv]]$variance)), function(d)
      data.frame(level = d[[f]][1], variance = d$variance[1], THR = round(mean(d$correct), 3),
                 medWt = round(median(d$Weight_Ht), 3), conf_wrong = round(mean(d$Weight_Ht < .05), 3))))
    out <- rbind(out, cbind(evidence = typeEv, x))
  }
  out <- out[order(out$evidence, out$level, out$variance), ]; rownames(out) <- NULL
  names(out)[names(out) == "level"] <- f
  out
}
for (f in c("sets", "N", "true_R2", "ratios")) { cat("\n### hetero: per", f, "\n"); print(kable(summ(f))) }
cat("\n### hetero overall\n")
for (typeEv in names(D)) print(kable(aggregate(cbind(THR = correct) ~ variance, D[[typeEv]], mean), caption = typeEv))

#### ---- Plot: THR per set, Added/Equal x original/corrected ----
allsets <- c(1, 4, 2, 5, 8, 3, 6, 9); labelnames <- c(1, 4, 2, 5, 7, 3, 6, 8)
thr <- function(d, v) sapply(allsets, function(s) mean(d$correct[d$sets == s & d$variance == v]))
pdf(file.path("images", "THR_sets_hetero_original_vs_corrected.pdf"), width = 8, height = 6)
x <- seq_along(allsets)
plot(x - .2, thr(D$added, "original (swapped)"), pch = 0, ylim = c(0, 1), xaxt = "n", xlab = "sets", ylab = "THR",
     xlim = c(.6, length(x) + .4), main = "hetero condition (N <= 150): original vs corrected error variances")
points(x - .07, thr(D$added, "corrected"), pch = 15)
points(x + .07, thr(D$equal, "original (swapped)"), pch = 1, col = "darkorange3")
points(x + .2, thr(D$equal, "corrected"), pch = 19, col = "darkorange3")
axis(1, at = x, labels = labelnames)
legend("bottomleft", legend = c("Added, original", "Added, corrected", "Equal, original", "Equal, corrected"),
       pch = c(0, 15, 1, 19), col = c("black", "black", "darkorange3", "darkorange3"), bty = "n", cex = .8)
dev.off()
