# svygray<-function(time, status, group, design,rho=0,gamma=0){
  # fSurv.cause1 <- eval(substitute(Surv(time, status==1)), model.frame(design))
  # fSurv.crisks <- eval(substitute(Surv(time, factor(status))), model.frame(design))
  # 
  # f.cause1 <- as.formula(paste0('fSurv.cause1 ~ ', substitute(group)))
  # f.crisks <- as.formula(paste0('fSurv.crisks ~ ', substitute(group)))
svygray<-function(f.crisks, design,rho=0,gamma=0, ecode = '1'){
  terms <- all.vars(f.crisks)
  
  expr <- paste0('Surv(', terms[1], ',', terms[2], '==', ecode, ')')
  fSurv.cause1 <- eval(parse(text=expr), model.frame(design))
  f.cause1 <- as.formula(paste0('fSurv.cause1 ~ ', terms[3]))
  
  null.f.cause1<-update(f.cause1,.~1)
  S<-svykm(null.f.cause1,design,se=FALSE)
  epsilon<-min(diff(sort(unique(S$time))))/10
  w<-approxfun(S$time+epsilon,S$surv^rho*(1-S$surv)^gamma,
               method="constant",rule=2)
  
  # Remove fns environments
  environment(f.cause1)<-environment()
  environment(f.crisks)<-environment()
  
  # fg <- finegray(f.crisks, data=model.frame(design),
  fg <<- finegray(Surv(ftime, factor(fstatus)) ~ ., data=model.frame(design),
                 weights=weights(design,"sampling"), etype = ecode)
  fg.dsgn <- svydesign(id = ~1, strata = ~strat, weights = ~fgwt*`(weights)`, data = fg)
  # fg.dsgn <- svydesign(id = ~1, weights = ~fgwt*`(weights)`, data = fg)
  # fg.dsgn <- update(design, weights = ~fgwt*`(weights)`, data = fg)
  coxmodel <- coxph(Surv(fgstart, fgstop, fgstatus) ~ x1,
                    weight=fgwt, data=fg, iter.max=0)
  
  # coxmodel <- coxph(f.cause1,data=model.frame(design),
  #                    weights=weights(design,"sampling"),
  #                    iter.max=0)
  
  x<-model.matrix(coxmodel)
  detail <- coxph.detail(coxmodel, riskmat=TRUE)	
  
  # cinc <- survfit(f.crisks, data=model.frame(design),
  #                 weights = weights(design,"sampling"),
  #                 se = F)
  # F1 <- cinc$pstate[,2]
  # OS <- cinc$pstate[,1]
  # Y<-t(detail$riskmat*(1-F1)/OS)
  
  Y<-t(detail$riskmat)
  dLambda<-detail$hazard
  E<-as.matrix(detail$means)
  N<-coxmodel$y[,"status"]
  
  # times<-coxmodel$y[,"time"]
  times<-coxmodel$y[,"stop"]
  U<-matrix(nrow=nrow(x),ncol=ncol(x))
  index<-match(times[N==1],detail$time)
  ZmEdN<- matrix(0,nrow=nrow(x),ncol=ncol(x))
  ZmEdN[N==1,]<-x[N==1,,drop=FALSE]-E[index,]
  for(p in 1:ncol(x)){
    ZmE <- -outer(E[,p], x[,p], "-")  ##times are rows, people are columns
    U[,p]<- ZmEdN[,p]*w(times)-colSums(w(detail$time)*ZmE*dLambda*Y)
  }
  # means <- svytotal(U,design)
  means <- svytotal(U,fg.dsgn)
  zstat<-coef(means)/SE(means)
  chisqstat<-coef(means)%*%solve(vcov(means),coef(means))
  
  rval<-list(cbind(score=coef(means),
                   se=SE(means),
                   z=coef(means)/SE(means),
                   p=2*pnorm(-abs(coef(means)/SE(means)))),
             c(chisq=chisqstat,p=pchisq(chisqstat, df=ncol(x), lower.tail=FALSE)))
  rval
}

