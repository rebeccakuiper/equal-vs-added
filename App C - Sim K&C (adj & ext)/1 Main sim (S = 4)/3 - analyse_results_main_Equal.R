###############################################################################
# 3 - analyse_results_main_Equal.R
# Does for the EQUAL-evidence results of the main simulation what
# analysis/analyse_results.R and analysis/post_analysis.R do for the ADDED-evidence
# results (tables and the plots in 'images Added article'), and compares Equal with
# Added (re-run; same data sets).
#
# Run 2 - process_results_main_Equal.R first.
# Working directory: the folder of this script.
proj <- file.path("..", "Orig code K&C", "Sim - Added - Final")
###############################################################################
source(file.path(proj, "analysis", "functions_analysis.R"))   # original functions
source(file.path("..", "functions_analysis_Equal.R"))         # extra functions
library(knitr)
dir.create("images", showWarnings = FALSE)

numerical   <- c("true_R2", "rho", "N")      # as in the original script
categorical <- c("sets", "models", "ratios")

D <- list(added = readRDS(file.path("data", "data_to_analyse_main_added.RData")),
          equal = readRDS(file.path("data", "data_to_analyse_main_equal.RData")))

# Is the true hypothesis the failsafe (Hu or complement) or an informative hypothesis?
# Hu true: sets 3, 6 with ratio1-3, sets 2, 5 with ratio4; complement true: set 9 (and set 8 with ratio4)
Ht_failsafe <- function(d) (d$sets %in% c(3, 6) & d$ratios != "ratio4") | (d$sets %in% c(2, 5) & d$ratios == "ratio4") |
                           (d$sets == 9) | (d$sets == 8 & d$ratios == "ratio4")
for (k in names(D)) D[[k]]$Ht_is_failsafe <- Ht_failsafe(D[[k]])

#### ---- Overall performance (as in analyse_results.R) ----
good_practice_subset <- function(dat){
  gp <- dat[!(dat$true_R2 == 0.02 & dat$sets %in% c(3, 6) & dat$ratios != "ratio4"), ]
  gp <- gp[!(gp$sets %in% c(2, 5) & gp$true_R2 == 0.02 & gp$ratios == "ratio4"), ]
  gp[!(gp$sets %in% c(1, 2, 3) & gp$true_R2 == 0.02), ]
}
overall <- sapply(D, function(dat) c(
  overall     = performance_condition(dat)[c("THR", "Weight_ht_p50", "Support_p50")],
  good_pract  = performance_condition(good_practice_subset(dat))[c("THR", "Weight_ht_p50", "Support_p50")],
  Ht_informative = performance_condition(dat[!dat$Ht_is_failsafe, ])[c("THR", "Weight_ht_p50")],
  Ht_failsafe    = performance_condition(dat[dat$Ht_is_failsafe, ])[c("THR", "Weight_ht_p50")],
  conf_wrong  = mean(dat$Weight_Ht < .05)))   # 'confidently wrong': weight true hypothesis < .05
cat("\n### Overall performance: Added vs Equal\n"); print(kable(round(overall, 3)))

for (k in names(D)) {
  dat <- D[[k]]
  cat("\n\n######## ", toupper(k), " evidence ########\n")
  cat("### overall\n");             print(kable(t(performance_condition(dat))))
  cat("### true selection cases\n");  print(kable(t(performance_condition(dat[dat$correct == 1, ]))))
  cat("### false selection cases\n"); print(kable(t(performance_condition(dat[dat$correct == 0, ]))))
  for (f in c(categorical, numerical)) { cat("### per", f, "\n"); print(kable(marginal_effect(dat, f))) }
  cat("### THR: sets x N\n"); print(kable(xtabs(correct ~ sets + N, data = aggregate(correct ~ sets + N, dat, mean)), digits = 2))
  cat("### THR: sets x true_R2\n"); print(kable(xtabs(correct ~ sets + true_R2, data = aggregate(correct ~ sets + true_R2, dat, mean)), digits = 2))
  cat("### THR: true hypothesis is failsafe (Hu/complement) or not, x N\n")
  print(kable(xtabs(correct ~ Ht_is_failsafe + N, data = aggregate(correct ~ Ht_is_failsafe + N, dat, mean)), digits = 2))
}

#### ---- Plots as in the Added article, for Equal (and for the Added re-run) ----
for (k in names(D)) {
  dat <- D[[k]]
  # marginal effects of N, R2, rho (cf. samplesize_marginal.pdf etc.);
  # marginal_effect_Eq() = marginal_effect() with a y-axis (right panel) that also shows the
  # 'all' and 'false' cases (with marginal_effect(), these fall off the plot for Equal)
  pdf(file.path("images", paste0("samplesize_marginal_", k, ".pdf")), height = 5, width = 7)
  marginal_effect_Eq(dat, "N", plot_out = TRUE, main = toupper(k)); dev.off()
  pdf(file.path("images", paste0("trueR2_marginal_", k, ".pdf")), height = 5, width = 7)
  marginal_effect_Eq(dat, "true_R2", plot_out = TRUE, main = toupper(k)); dev.off()
  pdf(file.path("images", paste0("rho_marginal_", k, ".pdf")), height = 5, width = 7)
  marginal_effect_Eq(dat, "rho", plot_out = TRUE, MARSlegendpos = "topright", main = toupper(k)); dev.off()
  # THR per set (cf. THR_sets.pdf)
  pdf(file.path("images", paste0("THR_sets_", k, ".pdf"))); THR_sets_plot(dat); title(toupper(k)); dev.off()
  # HSR (cf. HSR_selection_criteria.pdf; post_analysis.R)
  pdf(file.path("images", paste0("HSR_selection_criteria_", k, ".pdf")), width = 7, height = 5)
  HSR_plot(dat); mtext(toupper(k), side = 3, outer = TRUE, line = -1.5); dev.off()
  # weight of the selected false hypothesis (cf. effect_on_hf.pdf)
  pdf(file.path("images", paste0("effect_on_hf_", k, ".pdf")), height = 5, width = 7)
  effect_on_hf(dat); dev.off()
}
# Added re-run with the original marginal_effect() (to compare with the article figures)
pdf(file.path("images", "samplesize_marginal_added_original_function.pdf"), height = 5, width = 7)
marginal_effect(D$added, "N", plot_out = TRUE); dev.off()

#### ---- Comparison plots Added vs Equal ----
pdf(file.path("images", "THR_sets_added_vs_equal.pdf"), width = 8, height = 6)
THR_sets_plot(D$added, D$equal); dev.off()
for (f in c("N", "true_R2", "rho")) {
  pdf(file.path("images", paste0("compare_", f, "_added_vs_equal.pdf")), height = 5, width = 9)
  compare_marginal(D$added, D$equal, f, main = "all sets"); dev.off()
}
# per set (article numbering: 1, 4, 2, 5, 7, 3, 6, 8) and per type of true hypothesis
setlab <- c("1" = 1, "4" = 4, "2" = 2, "5" = 5, "8" = 7, "3" = 3, "6" = 6, "9" = 8)
pdf(file.path("images", "compare_N_added_vs_equal_per_set.pdf"), height = 5, width = 9)
for (s in names(setlab)) for (hf in c(FALSE, TRUE)) {
  sa <- D$added[D$added$sets == as.numeric(s) & D$added$Ht_is_failsafe == hf, ]
  if (nrow(sa) == 0) next
  se <- D$equal[D$equal$sets == as.numeric(s) & D$equal$Ht_is_failsafe == hf, ]
  compare_marginal(sa, se, "N", main = paste0("Set ", setlab[s], if (hf) " (true: failsafe)" else " (true: informative hyp.)"))
}
dev.off()
