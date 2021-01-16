# Survey Pepe-Mori --------------------------------------------------------

# CURRENTLY DOES NOT MATCH compCIF FROM PENTILE 2006
svypepemori <- function(design) {
  ajw <- survfit(Surv(ftime, factor(fstatus)) ~ x1, se = F,
                 weights = weights(design),
                 data = model.frame(design))
  
  # Overall CI (for event 1)
  # Ni <- ajw$n
  Ni <- as.vector(
    tapply(weights(design),
           model.frame(design)[,'x1'], sum))
  
  N <- sum(Ni)
  t <- sort(ajw$time)
  ntimes <- length(t)
  ngroups <- length(Ni)
  
  Fi <- array(dim = c(ntimes, ngroups))
  for(i in 1:ngroups) {
    fn <- stepfun(ajw[i,]$time, c(0, ajw[i,]$pstate[,2]))
    Fi[,i] <- fn(t)
  }
  F0 <- rowSums(t(t(Fi)*Ni/N))
  
  # Weight function
  kmc <- survfit(Surv(ftime, fstatus==0) ~ x1, se = F,
                 weights = weights(design),
                 data = model.frame(design))
  
  C <- array(dim = c(ntimes, ngroups))
  for(i in 1:ngroups) {
    fn <- stepfun(kmc[i]$time, c(1, kmc[i]$surv))
    C[,i] <- fn(t)
  }
  W <- N*apply(C, 1, prod)/rowSums(t(t(C)*Ni))
  
  # Test statistic
  Z <- c(1, W[-ntimes])*(Fi - F0)*diff(c(0,t))
  Z <- t(sqrt(Ni)*t(Z))
  
  # return(Z)
  return(list(t, c(1, W[-ntimes]), Fi, F0, diff(c(0,t))))
}


