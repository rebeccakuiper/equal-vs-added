# Extra simulation (not in the paper): performance of Equal, Added, and the GORICA on FE
# meta-analyzed estimates for various numbers of studies, sample sizes, and effect sizes.
# Conditions: S in {2,4,10,20}, n_s in {50,85,150}, effect f2 in {.15 (medium), .35 (large)};
# H_123 true / false (opposite), evaluated vs complement. 100 reps. Methods: Equal, Added, FE-GORICA.
suppressMessages({library(restriktor); library(MASS); library(parallel)})
rho <- .1; sigma <- matrix(rho,3,3); diag(sigma) <- 1; ratio <- c(3,2,1); slopes <- c("x1","x2","x3")
betas_for <- function(f2){ R2 <- f2/(1+f2); ve <- 1-R2
  fun <- function(x) (t(ratio*x)%*%sigma%*%ratio*x)/(t(ratio*x)%*%sigma%*%ratio*x+ve)-R2
  list(b = ratio*uniroot(fun,c(0,100))$root, ve = ve) }
gen <- function(n, pb){ X <- mvrnorm(n, rep(0,3), sigma); y <- X%*%pb$b + rnorm(n,0,sqrt(pb$ve))
  d <- data.frame(X,y); names(d) <- c(slopes,"y"); f <- lm(y~x1+x2+x3,d); list(est=coef(f)[slopes], V=vcov(f)[slopes,slopes]) }
w_from <- function(LL,PT,S,type){ IC <- if(type=="added") -2*LL+2*PT else -2*LL+2*PT/S; w <- exp(-.5*(IC-min(IC))); w/sum(w) }
one <- function(r, S, n, f2){
  set.seed(r + 1000*n + 10*S + round(1e5*f2)); pb <- betas_for(f2)
  st <- replicate(S, gen(n,pb), simplify=FALSE); P <- lapply(st,`[[`,"est"); V <- lapply(st,`[[`,"V")
  W <- Reduce(`+`, lapply(V, solve)); Vfe <- solve(W); th <- drop(Vfe %*% Reduce(`+`, Map(function(p,v) solve(v)%*%p, P, V)))
  names(th) <- slopes; dimnames(Vfe) <- list(slopes,slopes)
  out <- c()
  for(Hn in c("true","false")){
    H <- list(H_123 = if(Hn=="true") 'x1 > x2 > x3' else 'x1 < x2 < x3')
    ev <- suppressMessages(evSyn(P, VCOV=V, hypotheses=H, type_ev="added", comparison="complement"))
    LL <- colSums(ev$LL_m); PT <- colSums(ev$PT_m)
    fe <- goric(th, VCOV=Vfe, hypotheses=H, comparison="complement", type="gorica")$result$gorica.weights[1]
    w <- c(equal=w_from(LL,PT,S,"equal")[1], added=w_from(LL,PT,S,"added")[1], FE=fe)
    wt <- if(Hn=="true") w else 1-w     # weight of the TRUE hypothesis
    names(wt) <- paste(Hn, c("equal","added","FE"), sep="|"); out <- c(out, wt)
  }
  out }
res <- NULL
for(f2 in c(.15,.35)) for(S in c(2,4,10,20)) for(n in c(50,85,150)){
  m <- do.call(rbind, mclapply(1:100, one, S=S, n=n, f2=f2, mc.cores=max(1,detectCores()-1)))
  row <- data.frame(f2=f2, S=S, n=n)
  for(k in colnames(m)){ row[[paste0("mean|",k)]] <- mean(m[,k]); row[[paste0("THR|",k)]] <- mean(m[,k]>.5)
                         row[[paste0("wrong95|",k)]] <- mean(m[,k]<.05) }
  res <- rbind(res,row); saveRDS(res,"sim2.rds"); cat(f2,S,n,"done\n")
}
