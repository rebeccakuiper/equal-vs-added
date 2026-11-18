###############################################################################
# functions_analysis_Equal.R
# Extra functions for summarizing and plotting the EQUAL-evidence results (and
# comparing them with the ADDED-evidence results) of the simulations of the
# Added-evidence paper (Kuiper & Clapper), in the same style as the original
# plots of that paper (see 'images Added article').
#
# These functions are used next to the original functions in
#   analysis/functions_analysis.R  (performance_condition(), marginal_effect(), ...)
# which are sourced from the original code of Kuiper & Clapper and used unchanged.
#
# New functions:
#   effect_on_hf()          : the 'effect_on_hf.pdf' plot (the original code for
#                             this plot is not available): median (and 5th/95th percentile)
#                             weight of the selected but false hypothesis, in the
#                             false selection cases, per level of N, R2, and rho.
#   THR_sets_plot()         : the 'THR_sets.pdf' plot (THR per set), optionally
#                             with Added and Equal next to each other.
#   HSR_plot()              : the 'HSR_selection_criteria.pdf' plot (post_analysis.R).
#   compare_marginal()      : THR and median weight of the true hypothesis per level
#                             of a design factor, for Added and Equal in one plot.
###############################################################################

# colours / symbols for the comparison plots
col_added <- "black"
col_equal <- "darkorange3"

effect_on_hf <- function(dat, factors = c("N", "true_R2", "rho")){
  # dat: data_to_analyse (performance per iteration) of ONE evidence type
  nf <- length(factors)
  par(mfrow = c(1, nf), oma = c(0, 0, 2, 0))
  for (j in seq_along(factors)) {
    cond <- factors[j]
    lvls <- names(table(dat[[cond]]))
    tab <- t(sapply(lvls, function(lvl){
      d <- dat[dat[[cond]] == lvl & dat$correct == 0, ]
      c(quantile(d$Weight_Hs, c(.05, .5, .95)), p_iters = round(nrow(d) / nrow(dat) * 100, 1))
    }))
    xlab <- switch(cond, "true_R2" = expression(R^2), "rho" = expression(rho), "N" = "N", cond)
    par(mar = c(5.1, if (j == 1) 4.4 else 2.1, 2.1, 1.1))
    x <- as.numeric(lvls)
    plot(x, tab[, 2], type = "b", col = "red", ylim = c(0, 1), xaxt = "n",
         xlab = xlab, ylab = if (j == 1) expression("Median " * " w "["f"]^"S") else "")
    axis(1, at = x, labels = lvls)
    points(x, tab[, 1], col = "red", pch = "-", cex = 1.1)
    points(x, tab[, 3], col = "red", pch = "+", cex = 1.1)
    axis(3, at = x, labels = tab[, "p_iters"])
    if (j == 1) legend("bottomleft", legend = c("5th percentile", "95th percentile"),
                       col = "red", pch = c("-", "+"), bty = "n", cex = 0.8)
  }
  mtext("% of total iterations", side = 3, outer = TRUE, line = 0.5)
  invisible(NULL)
}

THR_per_set <- function(dat, allsets = c(1, 4, 2, 5, 8, 3, 6, 9)){
  vapply(allsets, function(set){ cr <- dat[dat$sets == set, "correct"]; sum(cr) / length(cr) }, numeric(1))
}

THR_sets_plot <- function(dat, dat2 = NULL, labels2 = c("Added", "Equal"),
                          allsets = c(1, 4, 2, 5, 8, 3, 6, 9),
                          labelnames = c(1, 4, 2, 5, 7, 3, 6, 8)){
  # As in analyse_results.R: THR per set, grouped by colour (black: no Hu/compl.;
  # gray50: Hu; gray: complement / Hu true). With dat2: Added (squares) and
  # Equal (circles) next to each other.
  Colors <- c(rep("black", 2), rep("gray50", 3), rep("gray", 3))
  THRs <- THR_per_set(dat, allsets)
  xs <- seq_along(allsets)
  if (is.null(dat2)) {
    plot(xs, THRs, pch = 15, xlab = "sets", xaxt = "n", cex = 2.2, ylim = c(0, 1), ylab = "THR", col = Colors)
  } else {
    THRs2 <- THR_per_set(dat2, allsets)
    plot(xs - 0.15, THRs, pch = 15, xlab = "sets", xaxt = "n", cex = 2.0, ylim = c(0, 1), xlim = c(0.6, length(xs) + 0.4),
         ylab = "THR", col = Colors)
    points(xs + 0.15, THRs2, pch = 19, cex = 2.0, col = Colors)
    legend("bottomleft", legend = labels2, pch = c(15, 19), bty = "n")
  }
  for (k in xs) axis(1, at = k, labels = labelnames[k])
  invisible(rbind(THRs, if (!is.null(dat2)) THR_per_set(dat2, allsets)))
}

HSR_table <- function(dat, Q = 1:200){
  support <- dat$support; s_c <- dat$support[dat$correct == 1]; s_f <- dat$support[dat$correct == 0]
  res <- t(vapply(Q, function(q) c(sum(support > q) / length(support),
                                   sum(s_c > q) / max(1, length(s_c)),
                                   sum(s_f > q) / max(1, length(s_f))), numeric(3)))
  res <- cbind(Q, round(res, 3)); colnames(res) <- c("Q", "HSR", "HSR_correct", "HSR_incorrect")
  res
}

HSR_plot <- function(dat, Q = 1:200){
  # As in post_analysis.R
  res <- HSR_table(dat, Q)
  par(mfrow = c(1, 2))
  plot(res[, "Q"], res[, "HSR"], ylab = "Hypothesis Selection Rate", xlab = "q", ylim = c(0, 1),
       type = "l", xaxt = "n", col = "blue")
  axis(1, at = c(1, 50, 100, 150, 200))
  lines(res[, "Q"], res[, "HSR_correct"], col = "olivedrab4", lty = 4)
  lines(res[, "Q"], res[, "HSR_incorrect"], col = "red", lty = 2)
  legend("topright", legend = c("All cases", "True selection cases", "False selection cases"),
         col = c("blue", "olivedrab4", "red"), lty = c(1, 4, 2), cex = 0.7)
  r <- 1:30
  plot(res[r, "Q"], res[r, "HSR"], ylab = "", xlab = "q", ylim = c(0, 1), type = "l", xaxt = "n", col = "blue")
  axis(1, at = c(1, 5, 10, 15, 20, 25, 30))
  lines(res[r, "Q"], res[r, "HSR_correct"], col = "olivedrab4", lty = 4)
  lines(res[r, "Q"], res[r, "HSR_incorrect"], col = "red", lty = 2)
  invisible(res)
}

compare_marginal <- function(dat_added, dat_equal, condition, main = NULL){
  # Left panel: THR (squares) per level, Added (black, solid) vs Equal (orange, dashed).
  # Right panel: median weight of the true hypothesis with 5th/95th percentiles.
  lvls <- names(table(dat_added[[condition]]))
  perf <- function(d) t(sapply(lvls, function(l) performance_condition(d[d[[condition]] == l, ])))
  pa <- perf(dat_added); pe <- perf(dat_equal)
  xlab <- switch(condition, "true_R2" = expression(R^2), "rho" = expression(rho), "N" = "n", condition)
  x <- seq_along(lvls)
  par(mfrow = c(1, 2), mar = c(5.1, 4.4, 4.1, 1.1))
  plot(x, pa[, "THR"], type = "b", pch = 15, col = col_added, ylim = c(0, 1), xaxt = "n",
       xlab = xlab, ylab = "THR", main = main)
  lines(x, pe[, "THR"], type = "b", pch = 0, col = col_equal, lty = 2)
  axis(1, at = x, labels = lvls)
  legend("bottomright", legend = c("Added", "Equal"), col = c(col_added, col_equal), pch = c(15, 0), lty = c(1, 2), cex = 0.8)
  plot(x, pa[, "Weight_ht_p50"], type = "b", pch = 16, col = col_added, ylim = c(0, 1), xaxt = "n",
       xlab = xlab, ylab = expression("median w"["t"]^"S" * " (5th-95th perc.)"), main = main)
  suppressWarnings(arrows(x - 0.05, pa[, "Weight_ht_p5"], x - 0.05, pa[, "Weight_ht_p95"], angle = 90, code = 3, length = 0.03, col = col_added))  # (zero-length bars give a warning)
  lines(x + 0.05, pe[, "Weight_ht_p50"], type = "b", pch = 1, col = col_equal, lty = 2)
  suppressWarnings(arrows(x + 0.05, pe[, "Weight_ht_p5"], x + 0.05, pe[, "Weight_ht_p95"], angle = 90, code = 3, length = 0.03, col = col_equal, lty = 2))  # (zero-length bars give a warning)
  axis(1, at = x, labels = lvls)
  invisible(list(added = pa, equal = pe))
}

marginal_effect_Eq <- function(dat, condition, plot_out = FALSE, MARSlegendpos = "topleft", main = NULL){
  # Identical to marginal_effect() in analysis/functions_analysis.R, except for the
  # y-axis of the right panel (median relative aggregated support): there, the
  # y-axis was based on the true selection cases only, so the 'all cases' and
  # 'false cases' lines fall off the plot when the false selection cases have a
  # (much) higher support, which happens for Equal. Here, the y-axis covers all
  # three lines (on a log scale when the support exceeds 100).
  condition_name <- condition
  lvls <- names(table(dat[[condition_name]]))
  pc <- function(d) if (nrow(d) == 0) {x <- performance_condition(dat[1, ]); x[] <- NA; x["niters"] <- 0; x} else performance_condition(d)
  dat_per_level <- t(sapply(lvls, function(lvl) pc(dat[dat[[condition_name]] == lvl, ])))
  if (!plot_out) return(dat_per_level)
  xlab <- switch(condition, "true_R2" = expression(R^2), "rho" = expression(rho), "n")
  dat_per_level_true  <- t(sapply(lvls, function(lvl) pc(dat[dat[[condition_name]] == lvl & dat$correct == 1, ])))
  dat_per_level_false <- t(sapply(lvls, function(lvl) {
    perf <- pc(dat[dat[[condition_name]] == lvl & dat$correct == 0, ])
    c(perf, p_iters = round(unname(perf["niters"]) / nrow(dat), 3) * 100)}))
  par(mfrow = c(1, 2), oma = c(0, 0, if (is.null(main)) 0 else 1.5, 0))
  par(mar = c(5.1, 4.4, 4.1, 2.1))
  plot(lvls, dat_per_level[, "THR"], xlab = xlab, ylab = expression("THR & median " * " w"["t"]^"S"),
       xaxt = "n", type = "b", ylim = c(0, 1), pch = 0)
  axis(1, at = lvls, labels = lvls)
  lines(lvls, dat_per_level[, "Weight_ht_p50"], col = "blue", type = "b")
  points(lvls, dat_per_level[, "Weight_ht_p5"], col = "blue", pch = "-", cex = 1.1)
  points(lvls, dat_per_level[, "Weight_ht_p95"], col = "blue", pch = "+", cex = 1.1)
  legend("right", legend = c("THR", expression("median w"["t"]^"S"), "95th percentile", "5th percentile"),
         col = c("black", "blue", "blue", "blue"), cex = 0.70, pch = c(0, 1, 3, charToRaw("-")))
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  allsup <- c(dat_per_level[, "Support_p50"], dat_per_level_true[, "Support_p50"], dat_per_level_false[, "Support_p50"])
  mx <- max(allsup, na.rm = TRUE)
  lg <- if (mx > 100) "y" else ""
  plot(lvls, dat_per_level[, "Support_p50"], xlab = xlab, log = lg,
       ylab = paste0("Median Relative Aggregated Support", if (lg == "y") " (log)" else ""),
       xaxt = "n", type = "b", ylim = if (lg == "y") c(1, mx * 2) else c(1, mx + 5), col = "blue")
  axis(1, at = lvls, labels = lvls)
  lines(lvls, dat_per_level_true[, "Support_p50"], col = "olivedrab4", type = "b", pch = "t")
  lines(lvls, dat_per_level_false[, "Support_p50"], col = "red", type = "b", pch = "f")
  legend(MARSlegendpos, legend = c("all cases", "true cases", "false cases"),
         col = c("blue", "olivedrab4", "red"), cex = 0.70, pch = c(1, charToRaw("t"), charToRaw("f")))
  axis(3, at = rownames(dat_per_level_false), labels = dat_per_level_false[, "p_iters"], col.axis = "red", cex.axis = 0.8)
  if (!is.null(main)) mtext(main, side = 3, outer = TRUE, line = 0)
  invisible(dat_per_level)
}
