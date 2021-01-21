rm(list = ls()); gc()

library(tidyverse)
library(survey)
library(doParallel)

source('compCIF.R')
source('aux_code/my_TwoCauseFineGray.R')

# On Windows, set cores to be 1
if (.Platform$OS.type == "windows") {
  cores = 1
} else {
  cores = detectCores()
}

# find seed for this task
RNGkind("L'Ecuyer-CMRG")
# set a 'master' seed
set.seed(3419564)
# pbmcapply::pbmclapply(1:4, mc.set.seed = T, function(x){
#   runif(3)
# })

# -------------------------------------------------------------------------

N_pop <- 10000
# pwbeta1 <- list(cbind(0), cbind(1))
beta2 <- cbind(1)

z1 <- rnorm(N_pop)
Z <- cbind(z1)

ph <- c(0.5, 0.25, 0.15, 0.05, 0.05)
# ph <- c(0.8, 0.2)
# ph <- c(0.5, 0.5)
psu_size <- 5

x1 <- rep(0:1, N_pop*c(.6, .4)) #rbinom(N_pop, 1, .4)
# x1 <- c(rep(0:1, N_pop*ph[1]*c(.9, .1)),
#         rep(0:1, N_pop*ph[2]*c(.7, .3)),
#         rep(0:1, N_pop*ph[3]*c(.6, .4)),
#         rep(0:1, N_pop*ph[4]*c(.5, .5)),
#         rep(0:1, N_pop*ph[5]*c(.5, .5)))
X <- cbind(x1)

n <- 1000
f <- rep(1/length(ph), length(ph))

# dat <- fastcmprsk::simulateTwoCauseFineGrayModel(
#   N_pop, beta1, beta2, X, u.min = 1, u.max = 2, p = 0.3) %>%
#   bind_cols(X, Z)
# 
# dat_srs <-  sample_n(dat, n)
# 
# dat_s <- dat_srs %>%
#   mutate(w = N_pop/n)
# 
# dsgn <- survey::svydesign(
#   id = ~1, weights = ~w, data = dat_s)
# 
# z1 <- with(dsgn$variables,
#            compCIF(ftime, fstatus, x1, weights = w,
#                    return_si = T))

# z1[[8]] %>% sum()
# 
# # with(rbind(dsgn$variables, dsgn$variables, dsgn$variables, dsgn$variables,
# #            dsgn$variables),
# #      compCIF(ftime, fstatus, x1))
# # 
# # tau <- min(tapply(dat_s$ftime, dat_s$x1, max))
# 
# # unweighted
# z1i <- with(dsgn$variables,
#             compCIF(ftime, fstatus, x1,
#                     return_si = T))
# z1i[[8]] %>% sum()
# 
# means <- svytotal(z1i, subset(dsgn, lead(ftime) <= tau))
# 
# z1
# means
# 
# 0.1141651*5



# B <- pbmcapply::pbmclapply(1:500, mc.set.seed = T, function(b){
B <- mclapply(1:500, mc.set.seed = T, function(r){
  # dat <- fastcmprsk::simulateTwoCauseFineGrayModel(
  #   N_pop, beta1, beta2, X, u.min = 1, u.max = 2, p = 0.7) %>% 
  #   bind_cols(X, Z) %>% 
  #   mutate(strat = sample(rep(1:length(ph), N_pop*ph), replace = F))
  # 
  # dat_srs <-  sample_n(dat, n)
  # 
  # dat_s <- dat %>%
  #   group_by(strat) %>%
  #   do({
  #     nh = n*f[.$strat[1]]
  #     Nh = N_pop*ph[.$strat[1]]
  #     sample_n(., nh) %>%
  #       mutate(prob = nh/Nh,
  #              w = 1/prob)
  #   })
  # 
  # dsgn <- survey::svydesign(
  #   id = ~1, strata = ~strat,
  #   weights = ~w, data = dat_s)
  
  # dat_s <- dat_srs %>% 
  #   mutate(w = N_pop/n)
  # 
  # dsgn <- survey::svydesign(
  #     id = ~1, weights = ~w, data = dat_s)
  
  # b <- list(rep(0, 5), rep(0, 5))
  # b <- list(rep(.3, 5), rep(-1, 5))
  # b <- list(seq(.25, .5, length.out = 5),
  #           seq(-.9, -1.1, length.out = 5))
  b <- list(c(.3,.3,.3,0,0),
            c(-1,-1,-1,0,0))
  b <- list(c(.3,.3,0,0,0),
            c(.3,.3,0,0,0))
  
  x1 <- rep(0:1, N_pop*ph[1]*c(.6, .4)); X <- cbind(x1)
  dat1 <- my_simulateTwoCauseFineGrayModel(
    N_pop*ph[1], pwbeta1 = lapply(b, '[', 1), beta2, ucp = .7, X, u.min = 1, u.max = 2, p = 0.7) %>% 
    bind_cols(X)
  
  x1 <- rep(0:1, N_pop*ph[2]*c(.6, .4)); X <- cbind(x1)
  dat2 <- my_simulateTwoCauseFineGrayModel(
    N_pop*ph[2], pwbeta1 = lapply(b, '[', 2), beta2, ucp = .7, X, u.min = 1, u.max = 2, p = 0.7) %>% 
    bind_cols(X)
  
  x1 <- rep(0:1, N_pop*ph[3]*c(.6, .4)); X <- cbind(x1)
  dat3 <- my_simulateTwoCauseFineGrayModel(
    N_pop*ph[3], pwbeta1 = lapply(b, '[', 3), beta2, ucp = .7, X, u.min = 1, u.max = 2, p = 0.7) %>% 
    bind_cols(X)
  
  x1 <- rep(0:1, N_pop*ph[4]*c(.6, .4)); X <- cbind(x1)
  dat4 <- my_simulateTwoCauseFineGrayModel(
    N_pop*ph[4], pwbeta1 = lapply(b, '[', 4), beta2, ucp = .7, X, u.min = 1, u.max = 2, p = 0.7) %>% 
    bind_cols(X)
  
  x1 <- rep(0:1, N_pop*ph[5]*c(.6, .4)); X <- cbind(x1)
  dat5 <- my_simulateTwoCauseFineGrayModel(
    N_pop*ph[5], pwbeta1 = lapply(b, '[', 5), beta2, ucp = .7, X, u.min = 1, u.max = 2, p = 0.7) %>% 
    bind_cols(X)
  
  dat <- bind_rows(dat1, dat2, dat3, dat4, dat5) %>% 
    mutate(strat = rep(1:length(ph), N_pop*ph))#,#[rank(z1)],
  # clust = rep(1:(N_pop/psu_size), each = psu_size)[rank(z1)])
  
  dat_srs <-  sample_n(dat, n)
  
  # s <- dat %>% 
  #   select(strat, clust) %>% 
  #   filter(!duplicated(.)) %>% 
  #   group_by(strat) %>% 
  #   do({
  #     nh = n*f[.$strat[1]]
  #     Nh = N_pop*ph[.$strat[1]]
  #     sample_n(., nh/psu_size) %>%
  #       mutate(w = Nh/nh)
  #   }) %>% 
  #   ungroup()
  # 
  # dat_s <- dat %>% 
  #   filter(paste(strat, clust) %in% paste(s$strat, s$clust)) %>%
  #   left_join(s, by = c('strat', 'clust'))
  # 
  # dsgn <- survey::svydesign(
  #   id = ~clust, strata = ~strat, nest = T,
  #   weights = ~w, data = dat_s)
  
  dat_s <- dat %>%
    group_by(strat) %>%
    do({
      nh = n*f[.$strat[1]]
      Nh = N_pop*ph[.$strat[1]]
      sample_n(., nh) %>%
        mutate(prob = nh/Nh,
               w = 1/prob)
    })
  
  dsgn <- survey::svydesign(
    id = ~1, strata = ~strat,
    weights = ~w, data = dat_s)
  
  # dat_s <- dat_srs %>% 
  #   mutate(w = N_pop/n)
  
  # dsgn <- survey::svydesign(
  #   id = ~1, weights = ~w, data = dat_s)
  
  z1 <- tryCatch({
    with(dat_srs, compCIF(ftime, fstatus, x1))
  }, error = function(e) NULL)
  
  z2 <- tryCatch({
    with(dat_s, compCIF(ftime, fstatus, x1))
  }, error = function(e) NULL)
  
  svy.pepe <- function(){
    tau <- min(tapply(dat_s$ftime, dat_s$x1, max))
  
    # unweighted
    z1i <- with(dsgn$variables,
                compCIF(ftime, fstatus, x1,
                        return_si = T))
    # means <- svytotal(z1i[[1]], dsgn2)
    means <- svytotal(z1i[[1]], subset(dsgn, lead(ftime) <= tau))
    # svyzstat <- coef(means)/SE(means)
    # svyvar <- SE(means)^2
    # svychisq <- svyzstat^2
    # svyp <- 1-pchisq(svychisq, 1)
    # 
    # out <- data.frame(svyzstat, svyvar, svychisq, svyp)
    out <- data.frame(svy.s = coef(means),
                      w.s = z1i[[3]],
                      svy.v = SE(means)^2)
    
    # weighted
    z1i2 <- with(dsgn$variables,
                 compCIF(ftime, fstatus, x1, weights = w,
                         return_si = T))
    # means2 <- svytotal(z1i2[[1]], dsgn2)
    means2 <- svytotal(z1i2[[1]], subset(dsgn, lead(ftime) <= tau))
    # svyzstat2 <- coef(means2)/SE(means2)
    # svyvar2 <- SE(means2)^2
    # svychisq2 <- svyzstat2^2
    # svyp2 <- 1-pchisq(svychisq2, 1)
    # 
    # out2 <- data.frame(svyzstat2, svyvar2, svychisq2, svyp2)
    out2 <- data.frame(svy2.s = coef(means2),
                       w2.s = z1i2[[3]],
                       svy2.v = SE(means2)^2)
    cbind(out, out2)
  }
  
  svyz1 <- tryCatch(svy.pepe(), error = function(e) NULL)
  
  if(is.null(z1) | is.null(svyz1))
    return(data.frame(B = r))
  
  data.frame(B = r, z1, z2, svyz1)
}) %>% bind_rows()
# B

# srs
mean(1-pchisq(B$s^2/B$v, 1) < 0.05, na.rm = T)

# s
mean(1-pchisq(B$s.1^2/B$v.1, 1) < 0.05, na.rm = T)

# svy
mean(1-pchisq(B$svy.s^2/B$svy.v, 1) < 0.05, na.rm = T)

# weighted
mean(1-pchisq(B$svy2.s^2/B$svy2.v, 1) < 0.05, na.rm = T)

# mean(1-pchisq(B$w.s^2/B$svy.v, 1) < 0.05, na.rm = T)
# mean(1-pchisq(B$w2.s^2/B$svy2.v/2, 1) < 0.05, na.rm = T)



head(B)

colMeans(B, na.rm = T) %>% round(5)
colSums(is.na(B), na.rm = T)




