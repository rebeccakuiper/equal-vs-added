# Version of GoricEvSynEssentials() for current restriktor versions (>= 0.6):
# goric() now needs the hypotheses as a list (hypotheses = list(H1 = ..., ...)).
# Apart from that, the computation is identical to the one in sim/functions.R.
GoricEvSynEssentials <- function(S, Param_studies, CovMx_studies, SameHypo, NrHypos, Hypo_studies,
                                 Safeguard = "unconstrained", Name_studies = 1:S) {
  NrHypos_incl <- if (Safeguard == "none") NrHypos else NrHypos + 1
  LL <- PT <- matrix(NA, nrow = S, ncol = NrHypos_incl)
  hyps <- as.list(unname(Hypo_studies)); names(hyps) <- paste0("H", 1:NrHypos)
  for (s in 1:S) {
    est <- Param_studies[[s]]; cov <- CovMx_studies[[s]]
    if (length(cov) == 1) cov <- matrix(cov)
    res_goric <- restriktor::goric(est, VCOV = cov, hypotheses = hyps,
                                   type = "gorica", comparison = Safeguard)
    LL[s, ] <- res_goric$result[, 2]
    PT[s, ] <- res_goric$result[, 3]
  }
  sumLL <- colSums(LL); sumPT <- colSums(PT)
  IC_added <- -2 * sumLL + 2 * sumPT         # added evidence
  IC_equal <- -2 * sumLL + 2 * sumPT / S     # equal evidence
  w <- function(ic) { d <- exp(-0.5 * (ic - min(ic))); d / sum(d) }
  # same output structure as before: the final weights are in row S+1
  list(added = rbind(matrix(NA, S, NrHypos_incl), w(IC_added)),
       equal = rbind(matrix(NA, S, NrHypos_incl), w(IC_equal)))
}
