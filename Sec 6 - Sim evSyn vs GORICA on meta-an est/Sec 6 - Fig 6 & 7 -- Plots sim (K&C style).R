###############################################################################
# Sec 6 - Fig 6 & 7 -- Plots sim (K&C style).R
#
# Plots, in the style of the plots in the Added-evidence paper (Kuiper & Clapper;
# see 'images Added article'), for the simulation of the Equal-vs-Added paper
# (Section 6 & Appendix D; script 'Section 6 - Tab 2 & App D - Tab 6 -- Simulation -
# evSyn vs MetaAn_GORICA.R'), which saves the weights per replication in
# 'sim_results.rds' (list with elements 'S20_n25', ..., 'S4_n600'; columns
# 'H|failsafe|method' with the overall weight of H_123).
#
# Design: S = 20 and S = 4 studies; n_s = 25, 50, 85, 150, 300, 600; H_123 true or
# false; failsafe = complement or unconstrained; methods: Equal, Added, and GORICA
# on fixed-effect (FE) and random-effects (RE) meta-analyzed estimates; 100 replications.
# True hypothesis: H_123 if H_123 is true (also versus the unconstrained, as in the
# paper), the failsafe if H_123 is false.
#
# Performance measures as in the Added paper:
#   correct   = 1 if the true hypothesis has the highest weight (THR = mean of correct);
#   Weight_Ht = weight of the true hypothesis; Weight_Hs = weight of the selected one;
#   support   = relative support of the selected hypothesis versus the other one.
#
# Working directory: folder of this script (with sim_results.rds in it).
###############################################################################

source("functions_analysis_Equal.R")  # HSR_plot(), marginal_effect_Eq(), ...
# performance_condition() as in analysis/functions_analysis.R of the Added paper:
performance_condition <- function(condition){
  niters <- nrow(condition)
  THR <- sum(condition$correct) / niters
  p_ht <- quantile(condition$Weight_Ht, c(.05, .5, .95)); names(p_ht) <- sprintf("Weight_ht_p%d", c(5, 50, 95))
  p_hs <- quantile(condition$Weight_Hs, c(.05, .5, .95)); names(p_hs) <- sprintf("Weight_hs_p%d", c(5, 50, 95))
  p_s  <- quantile(condition$support,   c(.05, .5, .95)); names(p_s)  <- sprintf("Support_p%d",   c(5, 50, 95))
  round(c(THR = THR, p_ht, p_s, p_hs, niters = niters), 3)
}
dir.create("images", showWarnings = FALSE)
# PDF versions (vector graphics; for the journal) of Figures 6 and 7 are saved in the folder
# "Figures (pdf)" in the main folder of the R files (one level up from this script's folder)
pdf_folder <- file.path("..", "Figures (pdf)")
dir.create(pdf_folder, showWarnings = FALSE)
# opens a png device (in "images") or, for a file name ending in .pdf, a pdf device (in pdf_folder)
open_device <- function(file, width, height) {
  if (grepl("\\.pdf$", file)) pdf(file.path(pdf_folder, file), width = width, height = height)
  else png(file.path("images", file), width = width, height = height, units = "in", res = 300)
}

#### ---- Data in long format ----
allres <- readRDS("sim_results.rds")
methods <- c(FE = "FE", RE = "RE", Equal = "equal", Added = "added")   # order as in the tables
dl <- NULL
for (k in names(allres)) {
  m <- allres[[k]]
  S <- as.numeric(sub("_n.*", "", sub("^S", "", k))); N <- as.numeric(sub(".*_n", "", k))
  for (H in c("true", "false")) for (cmp in c("complement", "unconstrained")) for (meth in names(methods)) {
    w <- m[, paste(H, cmp, methods[meth], sep = "|")]          # weight of H_123
    wt <- if (H == "true") w else 1 - w                         # weight of the true hypothesis
    ws <- pmax(w, 1 - w)                                        # weight of the selected hypothesis
    dl <- rbind(dl, data.frame(S = S, N = N, H = H, failsafe = cmp, method = meth,
                               Weight_Ht = wt, Weight_Hs = ws, correct = as.numeric(wt > .5),
                               support = pmin(ws / (1 - ws), 1e12)))  # capped at 1e12 (weights of 1)
  }
}
dl$method <- factor(dl$method, levels = names(methods))
saveRDS(dl, "data_to_analyse_sim_section6.rds")

# THR / mean weight table (cf. Table 2 (Section 6) and Table 6 (Appendix D) of the paper)
tab <-do.call(rbind, lapply(split(dl, list(dl$S, dl$failsafe, dl$H, dl$N, dl$method), drop = TRUE), function(d)
  data.frame(S = d$S[1], failsafe = d$failsafe[1], H = d$H[1], N = d$N[1], method = d$method[1],
             THR = mean(d$correct), mean_Wt = round(mean(d$Weight_Ht), 3), median_Wt = round(median(d$Weight_Ht), 3),
             p5_Wt = round(quantile(d$Weight_Ht, .05), 3), p95_Wt = round(quantile(d$Weight_Ht, .95), 3),
             conf_wrong = mean(d$Weight_Ht < .05),     # 'confidently wrong': weight true hypothesis < .05
             wt_below_25 = mean(d$Weight_Ht < .25),    # weight true hypothesis < .25
             wt_below_45 = mean(d$Weight_Ht < .45))))  # weight true hypothesis < .45
tab <- tab[order(tab$S, tab$failsafe, tab$H, tab$method, tab$N), ]
rownames(tab) <- NULL
write.csv(tab, "summary_sim_section6_per_condition.csv", row.names = FALSE)

# colors per method (no red/green: these are used for the false/true selection cases in Figure 7)
cols <- c(FE = "blue", RE = "purple", Equal = "darkorange2", Added = "black")
pchs <- c(FE = 1, RE = 2, Equal = 0, Added = 5)

#### ---- 1. 'samplesize_marginal' style, per method (THR & median weight; relative support) ----
# One pdf per S and failsafe; per page one method and scenario (H_123 true / false).
for (S in c(20, 4)) for (cmp in c("complement", "unconstrained")) {
  pdf(file.path("images", sprintf("samplesize_marginal_S%d_%s.pdf", S, cmp)), height = 5, width = 7)
  for (H in c("true", "false")) for (meth in names(methods)) {
    d <- dl[dl$S == S & dl$failsafe == cmp & dl$H == H & dl$method == meth, ]
    marginal_effect_Eq(d, "N", plot_out = TRUE,
                       main = sprintf("%s; S = %d; H_123 %s vs %s", meth, S, H, cmp))
  }
  dev.off()
}

#### ---- 2. Comparison of the methods: THR and median weight true hypothesis vs n ----
compare_methods <- function(d, main){
  lv <- sort(unique(d$N)); x <- seq_along(lv)
  par(mfrow = c(1, 2), mar = c(5.1, 4.4, 4.1, 1.1), oma = c(0, 0, 1.5, 0))
  plot(NA, xlim = c(0.8, length(x) + 0.2), ylim = c(0, 1), xaxt = "n", xlab = "n", ylab = "THR")
  axis(1, at = x, labels = lv)
  for (meth in names(methods)) {
    thr <- sapply(lv, function(n) mean(d$correct[d$N == n & d$method == meth]))
    lines(x, thr, type = "b", col = cols[meth], pch = pchs[meth], lty = if (meth %in% c("Equal", "Added")) 1 else 2)
  }
  legend("bottomright", legend = names(methods), col = cols, pch = pchs, lty = c(1, 1, 2, 2), cex = 0.75, bg = "white")
  plot(NA, xlim = c(0.8, length(x) + 0.2), ylim = c(0, 1), xaxt = "n", xlab = "n",
       ylab = expression("median w"["t"]^"S" * " (5th-95th perc.)"))
  axis(1, at = x, labels = lv)
  off <- c(Equal = -0.15, Added = -0.05, FE = 0.05, RE = 0.15)
  for (meth in names(methods)) {
    q <- sapply(lv, function(n) quantile(d$Weight_Ht[d$N == n & d$method == meth], c(.05, .5, .95)))
    lines(x + off[meth], q[2, ], type = "b", col = cols[meth], pch = pchs[meth], lty = if (meth %in% c("Equal", "Added")) 1 else 2)
    suppressWarnings(arrows(x + off[meth], q[1, ], x + off[meth], q[3, ], angle = 90, code = 3, length = 0.02, col = cols[meth]))  # (zero-length bars give a warning)
  }
  mtext(main, side = 3, outer = TRUE, line = 0)
}
for (S in c(20, 4)) for (cmp in c("complement", "unconstrained")) {
  pdf(file.path("images", sprintf("compare_methods_S%d_%s.pdf", S, cmp)), height = 5, width = 9)
  for (H in c("true", "false")) {
    compare_methods(dl[dl$S == S & dl$failsafe == cmp & dl$H == H, ],
                    sprintf("S = %d; H_123 %s; failsafe: %s", S, H, cmp))
  }
  dev.off()
}

#### ---- 3. 'THR_sets' style: THR per scenario (S, failsafe, H) for each method ----
scen <- expand.grid(H = c("true", "false"), failsafe = c("complement", "unconstrained"), S = c(4, 20), stringsAsFactors = FALSE)
scen$label1 <- with(scen, paste0("S = ", S, "; ", ifelse(failsafe == "complement", "compl.", "unc.")))
scen$label2 <- with(scen, paste0("H_123 ", H))
pdf(file.path("images", "THR_scenarios.pdf"), height = 5, width = 9)
par(mar = c(6.1, 4.1, 2.1, 1.1))
plot(NA, xlim = c(0.5, nrow(scen) + 0.5), ylim = c(0, 1), xaxt = "n", xlab = "", ylab = "THR (averaged over n)")
off <- c(Equal = -0.24, Added = -0.08, FE = 0.08, RE = 0.24)
for (i in 1:nrow(scen)) for (meth in names(methods)) {
  d <- dl[dl$S == scen$S[i] & dl$failsafe == scen$failsafe[i] & dl$H == scen$H[i] & dl$method == meth, ]
  points(i + off[meth], mean(d$correct), pch = 15, col = cols[meth], cex = 1.8)
}
axis(1, at = 1:nrow(scen), labels = scen$label1, cex.axis = 0.75)
axis(1, at = 1:nrow(scen), labels = scen$label2, cex.axis = 0.75, line = 1.2, tick = FALSE)
abline(v = 4.5, lty = 3)
legend("bottomleft", legend = names(methods), col = cols, pch = 15, bty = "n", horiz = TRUE)
dev.off()

#### ---- 4. 'effect_on_hf' style: weight of the selected but false hypothesis (false selection cases) ----
pdf(file.path("images", "effect_on_hf.pdf"), height = 5, width = 9)
for (S in c(20, 4)) for (cmp in c("complement", "unconstrained")) {
  par(mfrow = c(1, 4), oma = c(0, 0, 3, 0), mar = c(5.1, 4.1, 2.6, 0.6))
  for (meth in names(methods)) {
    d <- dl[dl$S == S & dl$failsafe == cmp & dl$method == meth, ]
    lv <- sort(unique(d$N)); x <- seq_along(lv)
    q <- sapply(lv, function(n) { w <- d$Weight_Hs[d$N == n & d$correct == 0]
      if (length(w) == 0) rep(NA, 3) else quantile(w, c(.05, .5, .95)) })
    pf <- sapply(lv, function(n) round(100 * sum(d$N == n & d$correct == 0) / nrow(d), 1))
    plot(x, q[2, ], type = "b", col = "red", ylim = c(0.5, 1), xaxt = "n", xlab = "n",
         ylab = if (meth == "Equal") expression("Median " * " w "["f"]^"S") else "", main = meth)
    points(x, q[1, ], pch = "-", col = "red", cex = 1.2); points(x, q[3, ], pch = "+", col = "red")
    axis(1, at = x, labels = lv); axis(3, at = x, labels = pf, cex.axis = 0.7, col.axis = "red")
  }
  mtext(sprintf("S = %d; failsafe: %s (H_123 true and false together); top: %% of iterations that are false selection cases", S, cmp),
        side = 3, outer = TRUE, line = 1, cex = 0.8)
}
dev.off()

#### ---- 5. 'HSR_selection_criteria' style, per method (S = 20 and S = 4; both failsafes) ----
pdf(file.path("images", "HSR_selection_criteria.pdf"), width = 7, height = 5)
for (S in c(20, 4)) for (meth in names(methods)) {
  HSR_plot(dl[dl$S == S & dl$method == meth, ])
  mtext(sprintf("%s; S = %d (all n, both scenarios and failsafes)", meth, S), side = 3, outer = TRUE, line = -1.5)
}
dev.off()

print(tab[tab$failsafe == "complement", ])

###############################################################################
#### ---- 6. >>> PAPER: Section 6, Figures 6 and 7 (complement as failsafe) ----
# In the style of Figure 2 of Kuiper & Clapper (2025).
# The results are plotted at the actual n_s values (x-axis to scale up to n_s = 150; the parts
# between 150 and 300 and between 300 and 600 are skipped, denoted by axis breaks); the numbers
# are placed close to their tick mark (pushed apart only as far as needed to avoid overlap) and
# linked to it by a thin line.
# Figure 6: median overall weight of the true hypothesis and its 5th and 95th percentile,
#   per method (FE, RE, Equal, Added) and n_s; rows: S = 20, S = 4; columns: H_123 true,
#   H_123 false (cf. left plot in Figure 2 of Kuiper & Clapper). Two versions:
#   images/sim_weight_true_compl_whiskers.png: percentiles as whiskers;
#   images/sim_weight_true_compl_pm.png: percentiles as '-' and '+' (as in Kuiper & Clapper).
# Figure 7: median overall weight of the SELECTED hypothesis for all cases (method color, o),
#   for the true selection cases (green, t), and for the false selection cases (red, f),
#   per method and n_s (H_123 true and false together); the red numbers on top denote the
#   percentage of replications that are false selection cases (cf. right plot in Figure 2 of
#   Kuiper & Clapper, which displays the relative support instead of the weight). Versions:
#   images/sim_weight_false_compl_pctl_false.png: 5th (-) and 95th (+) percentiles for the false cases;
#   images/sim_weight_false_compl_pctl_both.png: idem for the false and the true cases;
#   images/sim_weight_false_compl_whiskers_false.png / _whiskers_both.png: idem, with whiskers.
# The versions used in the paper (Figures 6 and 7) are sim_weight_true_compl_whiskers.png and
# sim_weight_false_compl_whiskers_both.png; at the end of this script, these are also saved as
# images/sim_weight_true_compl.png and images/sim_weight_false_compl.png (the names used in the TeX file).
###############################################################################
cmp <- "complement"
lv_all <- sort(unique(dl$N))                                 # 25, 50, 85, 150, 300, 600
# x-axis to scale up to n_s = 150; the parts between 150 and 300 and between 300 and 600
# are skipped (axis breaks)
xmap <- function(n) ifelse(n <= 150, n, ifelse(n <= 300, 150 + 60 * (n - 150) / 150, 210 + 60 * (n - 300) / 300))
xbreaks <- c((xmap(150) + xmap(300)) / 2, (xmap(300) + xmap(600)) / 2)
xr <- c(xmap(min(lv_all)) - 15, xmap(max(lv_all)) + 15)
offs <- c(FE = -10.5, RE = -3.5, Equal = 3.5, Added = 10.5)  # offsets in units of n_s
# axis breaks: a small 'square wave' on the frame between 150 and 300 and between 300 and 600 (bottom and top)
axis_break <- function(side = 1) {
  usr <- par("usr"); y <- if (side == 1) usr[3] else usr[4]
  w <- 0.012 * diff(usr[1:2]); h <- 0.012 * diff(usr[3:4])
  for (xbreak in xbreaks) {
    rect(xbreak - 1.6 * w, y - 1.5 * h, xbreak + 1.6 * w, y + 1.5 * h, col = "white", border = NA, xpd = NA)
    lines(c(xbreak - 1.6 * w, xbreak - w, xbreak - w, xbreak, xbreak, xbreak + w, xbreak + w, xbreak + 1.6 * w),
          c(y, y, y + h, y + h, y - h, y - h, y, y), xpd = NA)
  }
}
# whiskers (5th-95th percentile) with an open space around the median symbol
whisk <- function(x, lo, mid, hi, col) {
  usr <- par("usr"); g <- 0.85 * strheight("M"); cw <- 0.008 * diff(usr[1:2])
  for (i in seq_along(x)) {
    if (any(is.na(c(lo[i], mid[i], hi[i])))) next
    # (a whisker part that falls within the open space around the symbol is not drawn)
    if (lo[i] < mid[i] - g) segments(c(x[i], x[i] - cw), c(lo[i], lo[i]), c(x[i], x[i] + cw), c(mid[i] - g, lo[i]), col = col)
    if (hi[i] > mid[i] + g) segments(c(x[i], x[i] - cw), c(mid[i] + g, hi[i]), c(x[i], x[i] + cw), c(hi[i], hi[i]), col = col)
  }
}
# label positions: as close as possible to the actual n_s, but pushed apart (only as far
# as needed) such that the labels do not overlap
spread_labels <- function(x, labs, cex) {
  usr <- par("usr")
  w <- strwidth(as.character(labs), cex = cex)
  gap <- 0.5 * (w[-1] + w[-length(w)]) + 0.035 * diff(usr[1:2])   # minimal distance between centres
  xl <- x
  for (i in 2:length(x)) xl[i] <- max(xl[i], xl[i - 1] + gap[i - 1])
  ub <- usr[2] + 0.02 * diff(usr[1:2])                               # do not run far beyond the plot
  if (xl[length(xl)] > ub) {
    xl[length(xl)] <- ub
    for (i in (length(x) - 1):1) xl[i] <- min(xl[i], xl[i + 1] - gap[i])
  }
  xl
}
# x-axis: tick marks at the actual n_s; numbers close to them (without overlap), linked to their tick mark
axis_ns <- function(lv, cex.lab = 1) {
  xp <- xmap(lv)
  axis(1, at = xp, labels = FALSE, tcl = -0.3)
  usr <- par("usr"); pin <- par("pin"); plt_h <- diff(usr[3:4]) / pin[2]   # y-units per inch
  y0 <- usr[3] - 0.3 * par("csi") * plt_h      # end of the tick marks
  y1 <- usr[3] - 0.75 * par("csi") * plt_h     # start of the labels
  xl <- spread_labels(xp, lv, cex.lab)
  segments(xp, y0, xl, y1, xpd = NA, col = "gray40", lwd = 0.8)
  text(xl, y1, labels = lv, pos = 1, offset = 0.25, xpd = NA, cex = cex.lab)
  axis_break(1); axis_break(3)
}
# idem for the top axis (Figure 7: percentage of false selection cases, in red)
axis_ns_top <- function(lv, labs, cex.lab = 0.75, col = "red") {
  xp <- xmap(lv)
  axis(3, at = xp, labels = FALSE, tcl = -0.3)
  usr <- par("usr"); pin <- par("pin"); plt_h <- diff(usr[3:4]) / pin[2]
  y0 <- usr[4] + 0.3 * par("csi") * plt_h
  y1 <- usr[4] + 0.75 * par("csi") * plt_h
  xl <- spread_labels(xp, labs, cex.lab)
  segments(xp, y0, xl, y1, xpd = NA, col = "gray40", lwd = 0.8)
  text(xl, y1, labels = labs, pos = 3, offset = 0.2, xpd = NA, cex = cex.lab, col = col)
  axis_break(3)
}
xlab_ns <- expression(n[s])

fig6 <- function(file, type = c("whiskers", "pm")) {
  type <- match.arg(type)
  open_device(file, width = 9, height = 7.5)
  par(mfrow = c(2, 2), mar = c(5.1, 4.6, 2.6, 1.1), oma = c(2.5, 0, 0, 0), mgp = c(3.4, 1, 0))
  for (S in c(20, 4)) for (H in c("true", "false")) {
    d <- dl[dl$S == S & dl$failsafe == cmp & dl$H == H, ]
    lv <- sort(unique(d$N))
    plot(NA, xlim = xr, ylim = c(0, 1), xaxt = "n", xlab = xlab_ns,
         ylab = expression("median w"["t"]^"S"),
         main = if (H == "true") bquote(S == .(S) * ";  " * H[123] ~ "true (true hypothesis: " * H[123] * ")") else
                                 bquote(S == .(S) * ";  " * H[123] ~ "false (true hypothesis: complement)"),
         font.main = 1, cex.main = 0.95)
    axis_ns(lv)
    abline(h = .5, lty = 3, col = "gray40")   # dotted line: weight of .5
    for (meth in names(methods)) {
      q <- sapply(lv, function(n) quantile(d$Weight_Ht[d$N == n & d$method == meth], c(.05, .5, .95)))
      x <- xmap(lv) + offs[meth]
      if (type == "whiskers") {
        whisk(x, q[1, ], q[2, ], q[3, ], cols[meth])
      } else {
        points(x, q[1, ], pch = "-", col = cols[meth], cex = 1.3)
        points(x, q[3, ], pch = "+", col = cols[meth], cex = 1.0)
      }
      lines(x, q[2, ], type = "b", col = cols[meth], pch = pchs[meth], lty = 2, bg = "white")
    }
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  if (type == "whiskers") {
    legend("bottom", legend = c(paste("median,", names(methods)), "whiskers: 5th-95th percentile"),
           col = c(cols, "gray30"), pch = c(pchs, NA), lty = c(2, 2, 2, 2, 1), ncol = 5, bty = "n", cex = 0.8)
  } else {
    legend("bottom", legend = c(paste("median,", names(methods)), "95th percentile", "5th percentile"),
           col = c(cols, "gray30", "gray30"), pch = c(pchs, 3, 45), lty = c(2, 2, 2, 2, NA, NA),
           ncol = 6, bty = "n", cex = 0.8)
  }
  dev.off()
}
fig6("sim_weight_true_compl_whiskers.png", "whiskers")
fig6("sim_weight_true_compl_pm.png", "pm")

fig7 <- function(file, pctl = c("false", "both"), whisk = FALSE) {
  pctl <- match.arg(pctl)
  open_device(file, width = 10, height = 6.5)
  par(mfrow = c(2, 4), mar = c(5.1, 4.6, 4.6, 0.6), oma = c(0, 0, 0, 0), mgp = c(3.4, 1, 0))
  for (S in c(20, 4)) for (meth in names(methods)) {
    d <- dl[dl$S == S & dl$failsafe == cmp & dl$method == meth, ]
    lv <- sort(unique(d$N))
    qq <- function(w) if (length(w) == 0) rep(NA, 3) else quantile(w, c(.05, .5, .95))
    q_all <- sapply(lv, function(n) qq(d$Weight_Hs[d$N == n]))
    q_t   <- sapply(lv, function(n) qq(d$Weight_Hs[d$N == n & d$correct == 1]))
    q_f   <- sapply(lv, function(n) qq(d$Weight_Hs[d$N == n & d$correct == 0]))
    pf    <- sapply(lv, function(n) round(100 * mean(d$correct[d$N == n] == 0), 1))
    xp <- xmap(lv)
    plot(xp, q_all[2, ], type = "b", col = cols[meth], pch = 1, lty = 2, ylim = c(0.5, 1), xaxt = "n",
         xlim = xr, xlab = xlab_ns, ylab = if (meth == "FE") expression("median w"["s"]^"S") else "")
    lines(xp - 6, q_t[2, ], type = "b", col = "olivedrab4", pch = "t")
    lines(xp + 6, q_f[2, ], type = "b", col = "red", pch = "f")
    if (whisk) {
      whisk(xp + 6, q_f[1, ], q_f[2, ], q_f[3, ], "red")
      if (pctl == "both") whisk(xp - 6, q_t[1, ], q_t[2, ], q_t[3, ], "olivedrab4")
    } else {
      points(xp + 6, q_f[1, ], pch = "-", col = "red", cex = 1.3)
      points(xp + 6, q_f[3, ], pch = "+", col = "red", cex = 1.0)
      if (pctl == "both") {
        points(xp - 6, q_t[1, ], pch = "-", col = "olivedrab4", cex = 1.3)
        points(xp - 6, q_t[3, ], pch = "+", col = "olivedrab4", cex = 1.0)
      }
    }
    axis_ns(lv, cex.lab = 0.8)
    axis_ns_top(lv, pf)
    mtext(paste0(meth, ";  S = ", S), side = 3, line = 2.6, cex = 0.75, font = 2)
    if (S == 20 && meth == "FE")
      if (whisk) legend("bottomleft", legend = c("all cases", "true cases", "false cases", "whiskers: 5th-95th perc."),
             col = c(cols[meth], "olivedrab4", "red", "gray30"), pch = c("o", "t", "f", NA), lty = c(2, 1, 1, 1),
             cex = 0.75, bg = "white") else
      legend("bottomleft", legend = c("all cases", "true cases", "false cases", "95th percentile", "5th percentile"),
             col = c(cols[meth], "olivedrab4", "red", "gray30", "gray30"), pch = c("o", "t", "f", "+", "-"),
             cex = 0.75, bg = "white")
  }
  dev.off()
}
fig7("sim_weight_false_compl_pctl_false.png", "false")
fig7("sim_weight_false_compl_pctl_both.png", "both")
fig7("sim_weight_false_compl_whiskers_false.png", "false", whisk = TRUE)
fig7("sim_weight_false_compl_whiskers_both.png", "both", whisk = TRUE)

# versions used in the paper (names as in the TeX file)
file.copy(file.path("images", "sim_weight_true_compl_whiskers.png"), file.path("images", "sim_weight_true_compl.png"), overwrite = TRUE)
file.copy(file.path("images", "sim_weight_false_compl_whiskers_both.png"), file.path("images", "sim_weight_false_compl.png"), overwrite = TRUE)
# PDF versions of Figures 6 and 7 (saved in pdf_folder)
fig6("Fig6_sim_weight_true_compl.pdf", "whiskers")
fig7("Fig7_sim_weight_false_compl.pdf", "both", whisk = TRUE)

# >>> PAPER: Section 6 (text on Figure 7): weight of the selected hypothesis (w_s) for all
#     replications, for the true selection cases, and for the false selection cases (median,
#     5th and 95th percentile), and the percentage of false selection cases (complement as
#     failsafe; H_123 true and false together)
q3 <- function(w) if (length(w) == 0) c(NA, NA, NA) else round(quantile(w, c(.05, .5, .95)), 3)
sel <- do.call(rbind, lapply(split(dl[dl$failsafe == "complement", ], list(dl$S[dl$failsafe == "complement"],
              dl$N[dl$failsafe == "complement"], dl$method[dl$failsafe == "complement"]), drop = TRUE), function(d) {
  a <- q3(d$Weight_Hs); t <- q3(d$Weight_Hs[d$correct == 1]); f <- q3(d$Weight_Hs[d$correct == 0])
  data.frame(S = d$S[1], N = d$N[1], method = d$method[1], perc_false = round(100 * mean(d$correct == 0), 1),
             p5_Ws_all = a[1], median_Ws_all = a[2], p95_Ws_all = a[3],
             p5_Ws_true = t[1], median_Ws_true = t[2], p95_Ws_true = t[3],
             p5_Ws_false = f[1], median_Ws_false = f[2], p95_Ws_false = f[3]) }))
sel <- sel[order(-sel$S, sel$method, sel$N), ]
rownames(sel) <- NULL
print(sel)
write.csv(sel, "summary_sim_section6_selected_hypothesis.csv", row.names = FALSE)
