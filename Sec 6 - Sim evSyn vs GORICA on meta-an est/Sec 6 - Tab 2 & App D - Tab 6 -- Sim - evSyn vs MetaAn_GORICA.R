###############################################################################
# Used in: Section 6 (Simulation), Table 2 (complement as failsafe; S = 20 and S = 4)
#   and Appendix D, Table 6 (unconstrained as failsafe; S = 20 and S = 4): mean
#   overall GORICA weight for the true hypothesis ('mean_') and true hypothesis rate
#   ('THR_') per S and n_s, for H_123 true/false (column H) and the
#   complement/unconstrained as failsafe (cmp);
#   see sim_summary.csv / the printed 'out' at the end of this script.
# The text in Section 6 and in Appendix D (subsection Simulation) is based on these tables.
###############################################################################
###############################################################################
# Simulation (Section 6): equal & added evSyn vs GORICA on FE & RE meta-analysis
# S = 20 and S = 4 independent studies per replication; R replications per S and n_s.
# Output: sim_results.rds (weights per replication) and sim_summary.csv
# (mean weight of the TRUE hypothesis and true hypothesis rate (weight > .5)).
# Note: mclapply runs in parallel on Linux/Mac; on Windows set mc.cores = 1
# (or use parallel::parLapply).
###############################################################################
suppressMessages({library(restriktor); library(metafor); library(MASS); library(parallel)})
f2 <- 0.15; R2 <- f2/(1+f2); rho <- 0.1; var.e <- 1-R2
sigma <- matrix(rho,3,3); diag(sigma) <- 1; ratio_study <- c(3,2,1)
fun <- function(x) (t(ratio_study*x) %*% sigma %*% ratio_study*x)/(t(ratio_study*x) %*% sigma %*% ratio_study*x + var.e) - R2
betas <- ratio_study*uniroot(fun, c(0,100))$root
slopes <- c("x1","x2","x3")
gen_study <- function(n){ X <- mvrnorm(n, rep(0,3), sigma); y <- X %*% betas + rnorm(n,0,sqrt(var.e))
  d <- data.frame(X,y); names(d) <- c(slopes,"y"); f <- lm(y ~ x1+x2+x3, d); list(est=coef(f), V=vcov(f))}
metaan <- function(P,V,model){
  dat <- data.frame(study=rep(seq_along(P),each=3), coef=factor(rep(slopes,length(P)),levels=slopes), yi=unlist(lapply(P,function(p)p[slopes])))
  VV <- bldiag(lapply(V,function(v)v[slopes,slopes]))
  res <- if(model=="FE") rma.mv(yi,VV,mods=~0+coef,data=dat,method="ML") else
    tryCatch(rma.mv(yi,VV,mods=~0+coef,random=~coef|study,struct="UN",data=dat,method="ML"),
      error=function(e) rma.mv(yi,VV,mods=~0+coef,random=~coef|study,struct="DIAG",data=dat,method="ML"))
  est <- coef(res); names(est) <- slopes; vc <- vcov(res); dimnames(vc) <- list(slopes,slopes); list(est=est,VCOV=vc)}
w_from <- function(LL, PT, S, type){ IC <- if(type=="added") -2*LL + 2*PT else -2*LL + 2*PT/S; w <- exp(-.5*(IC-min(IC))); w/sum(w)}
one_rep <- function(r, n, S=20){
  # seed per replication (r), study-specific sample size (n), and number of studies (S)
  set.seed(1000*n + r + if (S == 20) 0 else 1e6 * S)
  st <- replicate(S, gen_study(n), simplify=FALSE); P <- lapply(st,`[[`,"est"); V <- lapply(st,`[[`,"V")
  FE <- metaan(P,V,"FE"); RE <- metaan(P,V,"RE")
  out <- c()
  for(Hn in c("true","false")){
    H <- list(H_123 = if(Hn=="true") 'x1 > x2 > x3' else 'x1 < x2 < x3')
    for(cmp in c("complement","unconstrained")){
      ev <- suppressMessages(evSyn(P, VCOV=V, hypotheses=H, type_ev="added", comparison=cmp))
      LL <- colSums(ev$LL_m); PT <- colSums(ev$PT_m)
      res <- c(equal = w_from(LL,PT,S,"equal")[1], added = w_from(LL,PT,S,"added")[1],
        FE = goric(FE$est, VCOV=FE$VCOV, hypotheses=H, comparison=cmp, type="gorica")$result$gorica.weights[1],
        RE = goric(RE$est, VCOV=RE$VCOV, hypotheses=H, comparison=cmp, type="gorica")$result$gorica.weights[1])
      names(res) <- paste(Hn, cmp, c("equal","added","FE","RE"), sep="|")
      out <- c(out, res)
    }
  }
  out
}
# sanity check of w_from against evSyn equal
set.seed(5); st <- replicate(20, gen_study(85), simplify=FALSE); P <- lapply(st,`[[`,"est"); V <- lapply(st,`[[`,"V")
H <- list(H_123='x1 > x2 > x3')
ev <- suppressMessages(evSyn(P, VCOV=V, hypotheses=H, type_ev="added", comparison="complement"))
eq <- suppressMessages(evSyn(P, VCOV=V, hypotheses=H, type_ev="equal", comparison="complement"))
stopifnot(isTRUE(all.equal(unname(w_from(colSums(ev$LL_m),colSums(ev$PT_m),20,"equal")), unname(eq$Cumulative_GORICA_weights[21,]))))
stopifnot(isTRUE(all.equal(unname(w_from(colSums(ev$LL_m),colSums(ev$PT_m),20,"added")), unname(ev$Cumulative_GORICA_weights[21,]))))
cat("sanity ok\n")
R <- 100; ns <- c(25, 50, 85, 150, 300, 600); Ss <- c(20, 4)
allres <- list()
for(S in Ss) for(n in ns){
  t0 <- Sys.time()
  m <- do.call(rbind, mclapply(1:R, one_rep, n=n, S=S, mc.cores = if (.Platform$OS.type == "windows") 1 else max(1, detectCores()-1)))
  allres[[paste0("S", S, "_n", n)]] <- m
  cat("S =", S, "n =", n, "done in", round(difftime(Sys.time(),t0,units="mins"),1), "min\n")
  saveRDS(allres, "sim_results.rds")
}

# >>> PAPER: Section 6, Table 2 (complement) and Appendix D, Table 6 (unconstrained)
# Summary: mean weight of the true hypothesis and true hypothesis rate
out <- NULL
for (k in names(allres)) { m <- allres[[k]]
  S <- sub("_n.*", "", sub("^S", "", k)); n <- sub(".*_n", "", k)
  for (H in c("true","false")) for (cmp in c("complement","unconstrained")) {
    r <- c(S=S, n=n, H=H, cmp=cmp, R=nrow(m))
    for (meth in c("FE","RE","equal","added")) {
      w <- m[, paste(H, cmp, meth, sep="|")]
      wt <- if (H=="true") w else 1-w   # weight of the TRUE hypothesis
      r[paste0("mean_",meth)] <- round(mean(wt), 3); r[paste0("THR_",meth)] <- mean(wt > .5)
    }
    out <- rbind(out, r)
  }}
out <- as.data.frame(out); write.csv(out, "sim_summary.csv", row.names=FALSE)
print(out)
