###############################################################################
# App C - Fig 8 -- Equal vs Added in sim K&C.R
# >>> PAPER (Equal-vs-Added paper): Appendix C 'Equal-Evidence in the Simulations of
#     Kuiper and Clapper (2025)', Figure 8 (images/sim_KC_THR.png).
# Uses (with the CORRECT probit/logit error variances in the hetero studies):
#  - main simulation (S = 4 studies; N = 25-150): re-run of folder 1, where the hetero
#    rows are replaced by the corrected re-run (4 - rerun_simulation_main_hetero_corrected.R);
#  - extra simulation (S = 20 studies): corrected re-run of folder 2
#    (rerun_simulation_extra_corrected.R; processed with corrected <- TRUE).
# Panels: (a) S = 4: THR per set (article numbering), Added vs Equal;
#         (b)-(c) S = 20, no outlier study: THR per N for Set 5 (H1, H2, Hu) and Set 4
#         (H1 vs its complement), when H1 is true (b) and when H1 is not true (c).
# Working directory: the folder of this script.
f1 <- file.path("..", "1 Main sim (S = 4)", "data")
f2 <- file.path("..", "2 Extra sim (S = 20)", "data")
###############################################################################
cols <- c(Added = "black", Equal = "darkorange3")

# main simulation: hetero rows replaced by the corrected ones
main <- list()
for (typeEv in c("added", "equal")) {
  d  <- readRDS(file.path(f1, paste0("data_to_analyse_main_", typeEv, ".RData")))
  hc <- readRDS(file.path(f1, paste0("data_to_analyse_main_hetero_corrected_", typeEv, ".RData")))
  d  <- d[d$models != "hetero", ]
  main[[typeEv]] <- rbind(d[, names(hc)], hc)
}
# extra simulation (corrected)
extra <- list(added = readRDS(file.path(f2, "data_to_analyse_extra_added_corrected.RData")),
              equal = readRDS(file.path(f2, "data_to_analyse_extra_equal_corrected.RData")))

allsets <- c(1, 4, 2, 5, 8, 3, 6, 9); labelnames <- c(1, 4, 2, 5, 7, 3, 6, 8)
thr_set <- function(d) sapply(allsets, function(s) mean(d$correct[d$sets == s]))

# Figure 8: drawn twice, as png (folder "images") and as pdf (vector graphics; for the journal),
# the latter in the folder "Figures (pdf)" in the main folder of the R files (two levels up)
pdf_folder <- file.path("..", "..", "Figures (pdf)")
dir.create(pdf_folder, showWarnings = FALSE)
fig8 <- function() {
  par(mfrow = c(1, 3), mar = c(4.6, 4.4, 3.1, 1.1))
  # (a)
  x <- seq_along(allsets)
  plot(x - .12, thr_set(main$added), pch = 15, col = cols["Added"], ylim = c(0, 1), xlim = c(.6, 8.4),
       xaxt = "n", xlab = "Set", ylab = "THR", main = "(a) S = 4 (N = 25-150)", font.main = 1, cex = 1.4)
  points(x + .12, thr_set(main$equal), pch = 19, col = cols["Equal"], cex = 1.4)
  axis(1, at = x, labels = labelnames)
  abline(v = c(2.5, 5.5), lty = 3, col = "gray60")
  legend("bottomleft", legend = names(cols), pch = c(15, 19), col = cols, bty = "n")
  # (b), (c)
  for (inv in 0:1) {
    lv <- c(25, 75, 200)
    plot(NA, xlim = c(.8, 3.2), ylim = c(0, 1), xaxt = "n", xlab = "N", ylab = "THR",
         main = if (inv == 0) expression("(b) S = 20, " * H[1] * " true") else expression("(c) S = 20, " * H[1] * " not true"))
    axis(1, at = 1:3, labels = lv)
    for (typeEv in c("added", "equal")) for (s in c(5, 7)) {
      d <- extra[[typeEv]]; d <- d[d$sets == s & d$inverse_ratio == inv & d$outlier_study == "None", ]
      y <- sapply(lv, function(n) mean(d$correct[d$N == n]))
      lines(1:3, y, type = "b", col = cols[if (typeEv == "added") "Added" else "Equal"],
            pch = if (s == 5) 15 else 17, lty = if (s == 5) 1 else 2)
    }
    if (inv == 0) legend("right", legend = c("Added, Set 5", "Added, Set 4", "Equal, Set 5", "Equal, Set 4"),
                         col = rep(cols, each = 2), pch = c(15, 17, 15, 17), lty = c(1, 2, 1, 2), bty = "n", cex = .85)
  }
}
png(file.path("images", "sim_KC_THR.png"), width = 10, height = 4.2, units = "in", res = 300)
fig8()
dev.off()
pdf(file.path(pdf_folder, "Fig8_sim_KC_THR.pdf"), width = 10, height = 4.2)
fig8()
dev.off()

# numbers used in the text
cat("S = 4: overall THR Added/Equal:", round(mean(main$added$correct), 2), round(mean(main$equal$correct), 2), "\n")
print(round(rbind(Added = thr_set(main$added), Equal = thr_set(main$equal)), 2))
Hf <- function(d) (d$sets %in% c(3, 6) & d$ratios != "ratio4") | (d$sets %in% c(2, 5) & d$ratios == "ratio4") |
                  (d$sets == 9) | (d$sets == 8 & d$ratios == "ratio4")
for (k in names(main)) { d <- main[[k]]; h <- Hf(d)
  cat(k, ": THR true=informative", round(mean(d$correct[!h]), 2), "; true=failsafe", round(mean(d$correct[h]), 2),
      "; conf. wrong", round(mean(d$Weight_Ht < .05), 3), "\n") }
for (k in names(extra)) for (inv in 0:1) for (s in c(5, 7)) {
  d <- extra[[k]]; d <- d[d$sets == s & d$inverse_ratio == inv, ]
  cat(k, "S=20 set", s, "inv", inv, ": THR", round(mean(d$correct), 2), " conf.wrong", round(mean(d$Weight_Ht < .05), 2),
      " per N:", round(tapply(d$correct, d$N, mean), 2), "\n") }

# proportion of replications with a weight for the true hypothesis smaller than .05, .25, .45
th <- c(.05, .25, .45); pb <- function(w) round(sapply(th, function(t) mean(w < t)), 3)
for (k in names(main)) { d <- main[[k]]; h <- Hf(d)
  cat("S = 4,", k, ": all", pb(d$Weight_Ht), "; true = informative", pb(d$Weight_Ht[!h]), "; true = failsafe", pb(d$Weight_Ht[h]), "\n") }
for (k in names(extra)) for (inv in 0:1) for (s in c(5, 7)) {
  d <- extra[[k]]; d <- d[d$sets == s & d$inverse_ratio == inv & d$outlier_study == "None", ]
  cat("S = 20,", k, ", set", s, if (inv == 0) "(H1 true)" else "(H1 not true)", ":", pb(d$Weight_Ht), "\n") }

# >>> PAPER: Appendix C (text): S = 20, without an outlier study (as in Figure 8b-c):
#     THR overall and per n_s
for (k in names(extra)) for (inv in 0:1) for (s in c(5, 7)) {
  d <- extra[[k]]; d <- d[d$sets == s & d$inverse_ratio == inv & d$outlier_study == "None", ]
  cat(k, "S=20 set", s, if (inv == 0) "(H1 true)" else "(H1 not true)", ", no outlier: THR", round(mean(d$correct), 2),
      " per N:", round(tapply(d$correct, d$N, mean), 2), "\n") }

# >>> PAPER: Section 6 (text on Figure 7, comparison with Kuiper & Clapper): added-evidence
#     approach in the main simulation (S = 4; re-run, n_s <= 150): median weight of the
#     selected hypothesis in the false and in the true selection cases, per n_s
d <- main$added
print(round(rbind(median_Ws_false = tapply(d$Weight_Hs[d$correct == 0], d$N[d$correct == 0], median),
                  median_Ws_true  = tapply(d$Weight_Hs[d$correct == 1], d$N[d$correct == 1], median)), 2))
