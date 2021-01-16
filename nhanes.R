rm(list = ls()); gc()

library(tidyverse)
library(survey)

source('svygray.R')
source('compCIF.R')

load('../../Application/nhanes_v4/data/nhanes_cvd_risk.rda')

nhanes <- nhanes_cvd_risk %>% 
  mutate(race = case_when(
    race == 3 ~ "3: Non-Hispanic White",
    race == 4 ~ "4: Non-Hispanic Black",
    TRUE ~ "Other"),
    
    pce_risk_cat = case_when(
      pce_risk_cat <= 2  ~ "<7.5: Low risk",
      is.na(pce_risk_cat) ~ NA_character_,
      TRUE ~ ">7.5: High risk"),
    
    gender = case_when(
      gender == 1  ~ "1: Male",
      TRUE ~ "2: Female"))

dt_etm <- nhanes %>%
  filter(time > 0) %>%
  # mutate(event = event + 1) %>% 
  select(id, to = event, exit = time, clust, strat, wt,
         eligible, gender, race, pce_risk_cat) %>%
  mutate(id = as.numeric(as.factor(id)),
         entry = 0,
         exit = as.numeric(exit),
         from = "1",
         # to = ifelse(to == 1, "cens", as.character(to)),
         to2 = (to == "2"))


dt_etm <- dt_etm %>% 
  # mutate(to = ifelse(to == 'cens', '0', to)) %>% 
  rename(ftime = exit,
         fstatus = to,
         x1 = race)

dt_all <- mutate(dt_etm, w_all = ifelse(eligible == 1, wt, 0.00000001))
# dt_mw <- mutate(dt_etm, w_mw = ifelse(eligible == 1 & gender == "1: Male" & race == "3: Non-Hispanic White", wt, 0.00000001))
# dt_mb <- mutate(dt_etm, w_mb = ifelse(eligible == 1 & gender == "1: Male" & race == "4: Non-Hispanic Black", wt, 0.00000001))
# dt_mo <- mutate(dt_etm, w_mo = ifelse(eligible == 1 & gender == "1: Male" & race == "Other", wt, 0.00000001))
# dt_fb <- mutate(dt_etm, w_fb = ifelse(eligible == 1 & gender == "2: Female" & race == "4: Non-Hispanic Black", wt, 0.00000001))
# dt_fo <- mutate(dt_etm, w_fo = ifelse(eligible == 1 & gender == "2: Female" & race == "Other", wt, 0.00000001))
# dt_fw <- mutate(dt_etm, w_fw = ifelse(eligible == 1 & gender == "2: Female" & race == "3: Non-Hispanic White", wt, 0.00000001))

dsgn_all <- svydesign(id = ~clust, strata = ~strat, weights = ~w_all, data = dt_all, nest = TRUE)
# dsgn_mw <- svydesign(id = ~clust, strata = ~strat, weights = ~w_mw, data = dt_mw, nest = TRUE)
# dsgn_mb <- svydesign(id = ~clust, strata = ~strat, weights = ~w_mb, data = dt_mb, nest = TRUE)
# dsgn_mo <- svydesign(id = ~clust, strata = ~strat, weights = ~w_mo, data = dt_mo, nest = TRUE)
# dsgn_fb <- svydesign(id = ~clust, strata = ~strat, weights = ~w_fb, data = dt_fb, nest = TRUE)
# dsgn_fo <- svydesign(id = ~clust, strata = ~strat, weights = ~w_fo, data = dt_fo, nest = TRUE)
# dsgn_fw <- svydesign(id = ~clust, strata = ~strat, weights = ~w_fw, data = dt_fw, nest = TRUE)




test <- with(dt_all %>% 
               filter(eligible == 1) %>% 
               filter(gender == "1: Male") %>%
               filter(x1 != "Other") %>% 
               mutate(fstatus = ifelse(fstatus == '3', '2', fstatus)),
             cmprsk::cuminc(ftime, fstatus, x1, cencode = '0'))
test$Tests
# plot(test)


library(survival)
km <- with(dt_all %>% 
             filter(eligible == 1) %>% 
             filter(gender == "1: Male") %>%
             filter(x1 != "Other") %>% 
             mutate(fstatus = ifelse(fstatus == '3', '2', fstatus)),
           survfit(Surv(ftime, factor(fstatus)) ~ x1, weights = wt))
plot(km)
plot(km[1:2,2], conf.int = T)


dt <- dt_all %>% 
  # filter(eligible == 1) %>% 
  # filter(gender == "1: Male") %>%
  # filter(x1 != "Other") %>% 
  mutate(fstatus = ifelse(fstatus == '3', '2', fstatus),
         w = ifelse((eligible == 1) & (gender == "1: Male") & (x1 != "Other"), wt, 0.0000001))
dsgn <- svydesign(id = ~clust, strata = ~strat, data = dt, weights = ~wt, nest = TRUE)
svygray(Surv(ftime, factor(fstatus)) ~ x1, design = dsgn)
# 0.9027254

# coxmodel <- coxph(Surv(fgstart, fgstop, fgstatus) ~ x1,
#                   weight=fgwt, data=fg)
# summary(coxmodel)

dt <- dt_all %>% 
  filter(eligible == 1) %>%
  filter(gender == "1: Male") %>%
  filter(x1 != "Other") %>%
  mutate(fstatus = ifelse(fstatus == '3', '2', fstatus),
         w = ifelse((eligible == 1) & (gender == "1: Male") & (x1 != "Other"), wt, 0.0000001)) %>% 
  mutate(ftime = ftime + seq(0.0001, 0.0009, length.out = length(ftime)))
dsgn <- svydesign(id = ~clust, strata = ~strat, data = dt, weights = ~wt, nest = TRUE)

z1i2 <- with(dsgn$variables,
             compCIF(ftime, fstatus, x1, weights = wt,
                     return_si = T))
means2 <- svytotal(c(z1i2[[1]], 0, 0), dsgn)
svyzstat2 <- coef(means2)/SE(means2)
svyvar2 <- SE(means2)^2
svychisq2 <- svyzstat2^2
svyp2 <- 1-pchisq(svychisq2, 1)
svyp2
# 0.4889472


pp <- with(dsgn$variables,
     compCIF(ftime, fstatus, x1))
1-pchisq(pp$s^2/pp$v, 1)







