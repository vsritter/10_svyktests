rm(list = ls()); gc()

library(tidyverse)
library(survey)
library(doParallel)

source('compCIF.R')
source('aux_code/svypepemori.R')

# On Windows, set cores to be 1
if (.Platform$OS.type == "windows") {
  cores = 1
} else {
  cores = detectCores()
}

set.seed(64753)

N_pop <- 100
# beta1 <- cbind(0)
beta1 <- cbind(0)
beta2 <- cbind(1)

x1 <- sample(rep(0:1, N_pop*c(.6, .4))) #rbinom(N_pop, 1, .4)
X <- cbind(x1)

z1 <- rnorm(N_pop)
Z <- cbind(z1)

ph <- c(0.8, 0.2)
psu_size <- 5

dat <- fastcmprsk::simulateTwoCauseFineGrayModel(
  N_pop, beta1, beta2, X, u.min = 1, u.max = 2, p = 0.7)

dat <- bind_cols(dat, X, Z) %>%
  mutate(strat = sample(rep(1:length(ph), N_pop*ph), replace = F),
         w = 1/runif(N_pop,.4,.6))

# dat <- arrange(dat, ftime) %>% mutate(w = rep(2:1, each = 5))
# dat2 <- rbind(dat, dat[1:5,]) %>% arrange(ftime)
# dat;dat2
# 
# dat <- sample_n(dat, nrow(dat))
# dat2 <- sample_n(dat2, nrow(dat2))

# time = with(dat, ftime)
# cens = with(dat, fstatus)
# group = with(dat, x1)
# weights = dat$w
# table(time,cens,group)*rep(2:1, each = 5)

with(dat, compCIF(ftime, fstatus, x1))

with(dat, compCIF(ftime, fstatus, x1, weights = w))



# dat$w <- 1/runif(nrow(dat),.06,.09)
# t1 <- with(dat, compCIF(ftime, fstatus, x1, weights = w))
# dsgn <- survey::svydesign(id = ~1, weights = ~w, data = dat)
# t2 <- svypepemori(dsgn)
# km <- survfit(Surv(ftime, factor(fstatus))~x1, data = dat, weights = w)

# t1[[1]] # dt
# t2[[4]]

# t1[[5]] # f1
# t2[[2]]

library(cmprsk)
with(dat, plot(cuminc(ftime, fstatus, x1)))

lines(km)
lines(as.numeric(names(t1[[7]])), t1[[7]], type = 's', col=2)
lines(t2[[1]], t2[[3]][,1], type = 's', col=3)



RNGkind("L'Ecuyer-CMRG")
set.seed(245753)

N_pop <- 1000
pbmcapply::pbmclapply(1:5000, mc.set.seed = T, function(x){
  x1 <- sample(rep(0:1, N_pop*c(.6, .4))) #rbinom(N_pop, 1, .4)
  X <- cbind(x1)
  
  z1 <- rnorm(N_pop)
  Z <- cbind(z1)
  
  dat <- fastcmprsk::simulateTwoCauseFineGrayModel(
    N_pop, beta1, beta2, X, u.min = 1, u.max = 2, p = 0.3)
  
  dat <- bind_cols(dat, X, Z) %>%
    mutate(strat = sample(rep(1:length(ph), N_pop*ph), replace = F))
  
  with(dat, compCIF(ftime, fstatus, x1))
}) %>%
  bind_rows() %>% 
  with(mean(1-pchisq(chisquare, 1) < 0.05, na.rm = T))















data('follic', package = 'randomForestSRC')
fit <- with(follic, compCIF(time, status, 1*(age>65)))
# fit <- with(follic, cmprsk::cuminc(time, status, 1*(age>65))$Tests)

dat <- follic %>% 
  mutate(ftime = time,
         fstatus = status,
         x1 = 1*(age>65))

pepemori(dat) %>% colSums()











# N_pop <- 2000
# beta1 <- cbind(0)
# beta2 <- cbind(1)
# 
# x1 <- sample(rep(0:1, N_pop*c(.6, .4))) #rbinom(N_pop, 1, .4)
# X <- cbind(x1)
# 
# z1 <- rnorm(N_pop)
# Z <- cbind(z1)
# 
# ph <- c(0.8, 0.2)
# # ph <- c(0.5, 0.5)
# psu_size <- 5
# 
# B <- lapply(1:1000, function(b){
#   dat <- fastcmprsk::simulateTwoCauseFineGrayModel(
#     N_pop, beta1, beta2, X, u.min = 1, u.max = 2, p = 0.3)
# 
#   dat <- bind_cols(dat, X, Z) %>%
#     mutate(strat = sample(rep(1:length(ph), N_pop*ph), replace = F))
# 
#   Z <- pepemori(dat)
#   z1 <- t(colSums(Z))
# 
#   data.frame(z1)
# }) %>% bind_rows()
# 
# mean(B[,1])
# sd(B[,1])
# 
# chisqstat <- mean(B[,1])/sd(B[,1])
# 1-pchisq(chisqstat^2, 1)





