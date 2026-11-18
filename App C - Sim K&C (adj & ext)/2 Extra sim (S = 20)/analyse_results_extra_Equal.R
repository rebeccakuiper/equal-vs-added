###############################################################################
# analyse_results_extra_Equal.R
# Does for the EQUAL-evidence results of the extra simulation (20 studies, based on
# the review) what analysis/analyse_results.R does for the ADDED-evidence results
# (tables and the N_*.pdf and R2_*.pdf plots), and adds plots comparing Added and
# Equal. Additionally, it also inspects the situation where the informative
# hypothesis is NOT true (inverse_ratio == 1), which was not inspected in Kuiper & Clapper.
#
# Run process_results_extra_Equal.R first.
# Working directory: the folder of this script.
proj      <- file.path("..", "Orig code K&C", "Sim - extra (review)")
proj_main <- file.path("..", "1 Main sim (S = 4)")  # for the 4-study comparison (re-run)
# corrected = FALSE: original results (probit/logit error variances swapped);
# corrected = TRUE : results with the correct variances (rerun_simulation_extra_corrected.R and,
#                    for the 4-study comparison, 4 - rerun_simulation_main_hetero_corrected.R)
corrected <- FALSE
suffix <- if (corrected) "_corrected" else ""
###############################################################################

source(file.path(proj, "analysis", "functions_analysis.R"))  # original functions
source(file.path("..", "functions_analysis_Equal.R"))        # extra functions
library(knitr)
img <- paste0("images", suffix)
dir.create(img, showWarnings = FALSE)

numerical   <- c("true_R2", "rho", "N", "inverse_ratio")      # as in the original script
categorical <- c("sets", "models", "ratios", "outlier_study")

dat_type <- list(added = readRDS(file.path("data", paste0("data_to_analyse_extra_added", suffix, ".RData"))),
                 equal = readRDS(file.path("data", paste0("data_to_analyse_extra_equal", suffix, ".RData"))))

# Subsets as in the original script (inverse_ratio == 0: informative hypothesis true)
# None = same R2 and ratio as the other studies; "inverse" = last study has the inverse
# ratio and R2 = 0.02; "extreme" = last study has the inverse ratio and R2 = 0.25.
# Set 5 (in article: Set 5): H1: b1>b2>b3, H2: b1<b2,b1>b3, Hu.  Set 7 (article: Set 4): H1 vs its complement.
subsets <- function(dat, inv = 0) list(
  no_outlier  = dat[dat$inverse_ratio == inv & dat$outlier_study == "None", ],
  inv_outlier = dat[dat$inverse_ratio == inv & dat$outlier_study == "inverse", ],
  ext_outlier = dat[dat$inverse_ratio == inv & dat$outlier_study == "extreme", ])

# 4-study comparison data (from the main simulation; sets 5 and 8, hetero, ratio1/2,
# rho = .25, N = 25, 75, 200), as in the original script; here from the re-run of the
# main simulation (which contains both Added and Equal, see folder 1), if available.
dat4 <- list()
for (typeEv in c("added", "equal")) {
  f <- file.path(proj_main, "data", paste0(if (corrected) "data_to_analyse_main_hetero_corrected_" else "data_to_analyse_main_", typeEv, ".RData"))
  if (file.exists(f)) {
    d <- readRDS(f)
    dat4[[typeEv]] <- d[d$sets %in% c(5, 8) & d$models == "hetero" & d$ratios %in% c("ratio1", "ratio2") &
                          d$rho == 0.25 & d$N %in% c(25, 75, 200), ]
  }
}

safe_plot <- function(expr) tryCatch(expr, error = function(e) {plot.new(); text(.5, .5, conditionMessage(e))})

#### ---- Tables: overall performance and per N / R2, per evidence type ----
for (typeEv in c("added", "equal")) {
  for (inv in 0:1) {
    sb <- subsets(dat_type[[typeEv]], inv)
    cat("\n\n######## ", toupper(typeEv), " evidence; ",
        if (inv == 0) "informative hypothesis TRUE (inverse_ratio = 0)" else
          "informative hypothesis NOT true (inverse_ratio = 1; Hu resp. complement true)", " ########\n")
    if (inv == 0 && !is.null(dat4[[typeEv]])) {
      cat("### Simulation results 4 studies \n"); print(kable(t(performance_condition(dat4[[typeEv]]))))
    }
    for (k in names(sb)) {
      cat("### 20 studies:", k, "\n"); print(kable(t(performance_condition(sb[[k]]))))
      cat("# per set\n");  print(kable(marginal_effect(sb[[k]], "sets")))
      cat("# per N\n");    print(kable(marginal_effect(sb[[k]], "N")))
      cat("# per R2\n");   print(kable(marginal_effect(sb[[k]], "true_R2")))
    }
  }
}

#### ---- Plots as in the original script (now for Equal; Added re-made for comparison) ----
# marginal_effect_Eq() = marginal_effect() with a y-axis (right panel) that also shows the
# all and false cases (see functions_analysis_Equal.R)
for (typeEv in c("added", "equal")) {
  for (inv in 0:1) {
    sb <- subsets(dat_type[[typeEv]], inv)
    tag <- paste0("_", typeEv, if (inv == 1) "_Hfalse" else "")
    for (k in names(sb)) {
      pdf(file.path(img, paste0("N_", k, tag, ".pdf")), height = 5, width = 7)
      safe_plot(marginal_effect_Eq(sb[[k]], "N", plot_out = TRUE, main = paste0(toupper(typeEv), " - ", k, if (inv == 1) " (H not true)" else " (H true)"))); dev.off()
      pdf(file.path(img, paste0("R2_", k, tag, ".pdf")), height = 5, width = 7)
      safe_plot(marginal_effect_Eq(sb[[k]], "true_R2", plot_out = TRUE, main = paste0(toupper(typeEv), " - ", k, if (inv == 1) " (H not true)" else " (H true)"))); dev.off()
    }
    if (inv == 0 && !is.null(dat4[[typeEv]])) {
      pdf(file.path(img, paste0("N_4_studies", tag, ".pdf")), height = 5, width = 7)
      safe_plot(marginal_effect_Eq(dat4[[typeEv]], "N", plot_out = TRUE, main = paste(toupper(typeEv), "- 4 studies"))); dev.off()
      pdf(file.path(img, paste0("R2_4_studies", tag, ".pdf")), height = 5, width = 7)
      safe_plot(marginal_effect_Eq(dat4[[typeEv]], "true_R2", plot_out = TRUE, main = paste(toupper(typeEv), "- 4 studies"))); dev.off()
    }
  }
}

#### ---- Comparison plots Added vs Equal (THR and median weight true hypothesis) ----
for (inv in 0:1) {
  sa <- subsets(dat_type$added, inv); se <- subsets(dat_type$equal, inv)
  lab <- if (inv == 0) "H true" else "H not true"
  for (k in names(sa)) for (cond in c("N", "true_R2")) {
    pdf(file.path(img, paste0("compare_", cond, "_", k, if (inv == 1) "_Hfalse" else "", ".pdf")), height = 5, width = 9)
    compare_marginal(sa[[k]], se[[k]], cond, main = paste0(k, ", ", lab)); dev.off()
  }
}
if (length(dat4) == 2) for (cond in c("N", "true_R2")) {
  pdf(file.path(img, paste0("compare_", cond, "_4_studies.pdf")), height = 5, width = 9)
  compare_marginal(dat4$added, dat4$equal, cond, main = "4 studies, H true"); dev.off()
}

#### ---- Compact overview table: THR Added vs Equal per set, N, R2 (no outlier) ----
overview <- NULL
for (inv in 0:1) for (s in c(5, 7)) for (N in c(25, 75, 200)) for (R2 in c(0.02, 0.15, 0.25)) {
  f <- function(d) d[d$inverse_ratio == inv & d$outlier_study == "None" & d$sets == s & d$N == N & d$true_R2 == R2, ]
  overview <- rbind(overview, data.frame(H = if (inv == 0) "true" else "not true", set = s, N = N, R2 = R2,
    THR_added = mean(f(dat_type$added)$correct), THR_equal = mean(f(dat_type$equal)$correct),
    medWt_added = round(median(f(dat_type$added)$Weight_Ht), 3), medWt_equal = round(median(f(dat_type$equal)$Weight_Ht), 3)))
}
print(kable(overview))
write.csv(overview, paste0("THR_overview_extra_Added_vs_Equal", suffix, ".csv"), row.names = FALSE)
