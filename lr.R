

library(tidyverse)
library(survey)
library(survival)


source('new_logrank.R')
# source('survey/R/logrank.R')
assignInNamespace("svylogrank", "svylogrank", ns = "survey")


data(pbc, package="survival")
pbc$randomized <- with(pbc, !is.na(trt) & trt>0)
biasmodel<-glm(randomized~age*edema,data=pbc)
pbc$randprob<-fitted(biasmodel)

pbcfg <- finegray(Surv(time,factor(status))~., data = pbc)

dpbc <- svydesign(id=~1, prob=~randprob, strata=~edema, data=pbcfg)

svylogrank(Surv(time, status==1)~trt,design=dpbc)


finegray(Surv(time,factor(status))~., data = model.frame(dpbc))


# Treat time to death and plasma cell malignancy as competing risks
etime <- with(mgus2, ifelse(pstat==0, futime, ptime))
event <- with(mgus2, ifelse(pstat==0, 2*death, 1))
event <- factor(event, 0:2, labels=c("censor", "pcm", "death"))

# FG model for PCM
pdata <- finegray(Surv(etime, event) ~ ., data=mgus2)
fgfit <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age + sex,
               weight=fgwt, data=pdata)







