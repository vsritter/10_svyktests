


rm(list = ls()); gc()

library(tidyverse)
library(survey)


# list.files("./", full.names = T) %>%
#   {grep("new_", ., value = T)} %>%
#   sapply(source,.GlobalEnv)

# assignInNamespace("onestage", onestage_fast, ns = "survey")

# source('new_svykm.R')
# assignInNamespace("svykm", "svykm", ns = "survey")

data(pbc, package="survival")
pbc$randomized <- with(pbc, !is.na(trt) & trt>0)
biasmodel<-glm(randomized~age*edema,data=pbc)
pbc$randprob<-fitted(biasmodel)

dpbc<-svydesign(id=~1, prob=~randprob, strata=~edema, data=subset(pbc,randomized))


aj <- survfit(Surv(time,factor(status))~1,
               data=model.frame(dpbc))

ajw <- survfit(Surv(time,factor(status))~1,
              weights = weights(dpbc,"sampling"),
              data=model.frame(dpbc))

# s1 <- svykm(Surv(time, status==2)~1, design=dpbc, se=F)

str(summary(aj))



plot(aj)
lines(ajw, col = 2)
lines(s1$time, 1-s1$surv, col = 3, type = 's')






plot(s1)




dpbc



