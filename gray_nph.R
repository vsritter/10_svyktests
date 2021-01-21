rm(list = ls()); gc()

library(tidyverse)
library(survey)
library(doParallel)

source('svygray.R')
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
# set.seed(357345)
set.seed(3419564)
# pbmcapply::pbmclapply(1:4, mc.set.seed = T, function(x){
#   runif(3)
# })

# -------------------------------------------------------------------------

N_pop <- 10000
pwbeta1 <- list(cbind(.3), cbind(-1))
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

# dat <- my_simulateTwoCauseFineGrayModel(
#   N_pop, pwbeta1, beta2, X, ucp = .5, u.min = 1, u.max = 2, p = 0.7) %>%
#   bind_cols(X, Z) %>%
#   mutate(strat = rep(1:length(ph), N_pop*ph)[rank(z1)],
#          clust = rep(1:(N_pop/psu_size), each = psu_size)[rank(z1)])

# dat <- fastcmprsk::simulateTwoCauseFineGrayModel(
#   N_pop, beta1=pwbeta1[[1]], beta2, X, u.min = 1, u.max = 2, p = 0.7) %>%
#   bind_cols(X, Z) %>%
#   mutate(strat = rep(1:length(ph), N_pop*ph)[rank(z1)],
#          clust = rep(1:(N_pop/psu_size), each = psu_size)[rank(z1)])

# library(cmprsk)
# table(dat$fstatus)
# table(dat$fstatus) %>% prop.table() %>% round(2)
# with(mutate(dat, x1=paste0('x',x1)),
#      plot(cuminc(ftime, fstatus, x1),
#           col=c(1,1,2,2),
#           lty=c(1,2,1,2)))

# with(dat, cuminc(ftime, fstatus, x1))


# dat_srs <-  sample_n(dat, n)
# 
# s <- dat %>%
#   select(strat, clust) %>%
#   filter(!duplicated(.)) %>%
#   group_by(strat) %>%
#   do({
#     nh = n*f[.$strat[1]]
#     Nh = N_pop*ph[.$strat[1]]
#     sample_n(., nh/psu_size)
#   }) %>%
#   ungroup()
# 
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
# 
# svylogrank(Surv(ftime, fstatus==1) ~ x1, design = dsgn)
# svygray(Surv(ftime, factor(fstatus)) ~ x1, design = dsgn)

# km <- survfit(Surv(ftime, fstatus==1) ~ x1, data = dsgn$variables,
#               weights = w)
# survival:::survmean(km, rmean = max(km$time))

# B <- pbmcapply::pbmclapply(1:5, mc.set.seed = T, function(r){
# B <- mclapply(1:500, mc.set.seed = T, function(r){
  # dat <- fastcmprsk::simulateTwoCauseFineGrayModel(
  #   N_pop, beta1, beta2, X, u.min = 1, u.max = 2, p = 0.7) %>% 
  #   bind_cols(X, Z) %>% 
  #   mutate(strat = rep(1:length(ph), N_pop*ph),#[rank(z1)],
  #          clust = rep(1:(N_pop/psu_size), each = psu_size)[rank(z1)])
  
  # b <- list(rep(0, 5), rep(0, 5))
  # b <- list(rep(.3, 5), rep(-1, 5))
  # b <- list(seq(0, .5, length.out = 5),
  #           seq(-.8, -1.2, length.out = 5))
  # b <- list(seq(.25, .5, length.out = 5),
  #           seq(-.9, -1.1, length.out = 5))
  b <- list(c(.3,.3,.3,0,0),
            c(-1,-1,-1,0,0))
  
  x1 <- rep(0:1, N_pop*ph[1]*c(.6, .4)); X <- cbind(x1)
  dat1 <- my_simulateTwoCauseFineGrayModel(
    N_pop*ph[1], pwbeta1 = lapply(b, '[', 1), beta2, ucp = .7, X, u.min = 1, u.max = 2, p = 0.7) %>% 
    bind_cols(X) %>% 
    mutate(clust = rep(1:(N_pop*ph[1]/psu_size), each = psu_size)[rank(ftime)])
  
  x1 <- rep(0:1, N_pop*ph[2]*c(.6, .4)); X <- cbind(x1)
  dat2 <- my_simulateTwoCauseFineGrayModel(
    N_pop*ph[2], pwbeta1 = lapply(b, '[', 2), beta2, ucp = .7, X, u.min = 1, u.max = 2, p = 0.7) %>% 
    bind_cols(X) %>% 
    mutate(clust = rep(1:(N_pop*ph[2]/psu_size), each = psu_size)[rank(ftime)])
  
  x1 <- rep(0:1, N_pop*ph[3]*c(.6, .4)); X <- cbind(x1)
  dat3 <- my_simulateTwoCauseFineGrayModel(
    N_pop*ph[3], pwbeta1 = lapply(b, '[', 3), beta2, ucp = .7, X, u.min = 1, u.max = 2, p = 0.7) %>% 
    bind_cols(X) %>% 
    mutate(clust = rep(1:(N_pop*ph[3]/psu_size), each = psu_size)[rank(ftime)])
  
  x1 <- rep(0:1, N_pop*ph[4]*c(.6, .4)); X <- cbind(x1)
  dat4 <- my_simulateTwoCauseFineGrayModel(
    N_pop*ph[4], pwbeta1 = lapply(b, '[', 4), beta2, ucp = .7, X, u.min = 1, u.max = 2, p = 0.7) %>% 
    bind_cols(X) %>% 
    mutate(clust = rep(1:(N_pop*ph[4]/psu_size), each = psu_size)[rank(ftime)])
  
  x1 <- rep(0:1, N_pop*ph[5]*c(.6, .4)); X <- cbind(x1)
  dat5 <- my_simulateTwoCauseFineGrayModel(
    N_pop*ph[5], pwbeta1 = lapply(b, '[', 5), beta2, ucp = .7, X, u.min = 1, u.max = 2, p = 0.7) %>% 
    bind_cols(X) %>% 
    mutate(clust = rep(1:(N_pop*ph[5]/psu_size), each = psu_size)[rank(ftime)])
  
  dat <- bind_rows(dat1, dat2, dat3, dat4, dat5) %>% 
    mutate(strat = rep(1:length(ph), N_pop*ph))
  
  # table(dat$fstatus) %>% prop.table()
  
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
      c <- sample(sample(Nh/psu_size, nh/psu_size))
      
      filter(., clust %in% c) %>%
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
  
  gray.srs <- tryCatch({
    test <- with(dat_srs,
         cmprsk::cuminc(ftime, fstatus, paste(x1), cencode = '0'))
    g1 <- test$Tests[1,1:2]
    g2 <- test$Tests[2,1:2]
    
    names(g1) <- paste0(names(g1), '.g1.srs')
    names(g2) <- paste0(names(g2), '.g2.srs')
    
    c(g1, g2)
  }, error = function(e) NULL)
  
  gray <- tryCatch({
    test <- with(dat_s,
                 cmprsk::cuminc(ftime, fstatus, paste(x1), cencode = '0'))
    g1 <- test$Tests[1,1:2]
    g2 <- test$Tests[2,1:2]
    
    names(g1) <- paste0(names(g1), '.g1.s')
    names(g2) <- paste0(names(g2), '.g2.s')
    
    c(g1, g2)
  }, error = function(e) NULL)
  
  svy.gray <- tryCatch({
    g1 <- svygray(Surv(ftime, factor(fstatus)) ~ x1, design = dsgn)[[2]]
    # g1 <- c(0,0)
    g2 <- g1
    # g2 <- svygray(Surv(ftime, factor(fstatus)) ~ x1, design = dsgn, ecode = '2')[[2]]
    # g1 <- svylogrank(Surv(ftime, fstatus==1) ~ x1, design = dsgn)[[2]]
    # g2 <- svylogrank(Surv(ftime, fstatus==2) ~ x1, design = dsgn)[[2]]
    
    names(g1) <- c('stat.g1.svy', 'pv.g1.svy')
    names(g2) <- c('stat.g2.svy', 'pv.g2.svy')
    
    c(g1, g2)
  }, error = function(e) NULL)
  
  if(is.null(gray.srs) & is.null(gray) & is.null(svy.gray))
    return(data.frame(B = r))
  
  data.frame(B = r, t(gray.srs), t(gray), t(svy.gray))
}) %>% bind_rows()
# B

# size
mean(B$pv.g1.srs < 0.05, na.rm = T)
mean(B$pv.g1.s < 0.05, na.rm = T)
mean(B$pv.g1.svy < 0.05, na.rm = T)


# power
# mean(B$pv.1 < 0.05, na.rm = T)
# mean(B$p.1 < 0.05, na.rm = T)



