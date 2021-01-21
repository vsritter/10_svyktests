
# -------------------------------------------------------------------------

B <- mclapply(1:ENV_REP, mc.cores = cores, mc.set.seed = T, function(r){
  x1 <- rep(0:1, N_pop*ph[1]*c(.5, .5)); X <- cbind(x1)
  dat1 <- my_simulateTwoCauseFineGrayModel(
    N_pop*ph[1], pwbeta1 = lapply(b, '[', 1), beta2, ucp = .7, X,
    u.min = u1, u.max = u2, p = 0.7) %>% 
    bind_cols(X) %>% 
    mutate(clust = rep(1:(N_pop*ph[1]/psu_size),
                       each = psu_size)[rank(ftime)])
  
  x1 <- rep(0:1, N_pop*ph[2]*c(.5, .5)); X <- cbind(x1)
  dat2 <- my_simulateTwoCauseFineGrayModel(
    N_pop*ph[2], pwbeta1 = lapply(b, '[', 2), beta2, ucp = .7, X,
    u.min = u1, u.max = u2, p = 0.7) %>% 
    bind_cols(X) %>% 
    mutate(clust = rep(1:(N_pop*ph[2]/psu_size),
                       each = psu_size)[rank(ftime)])
  
  x1 <- rep(0:1, N_pop*ph[3]*c(.5, .5)); X <- cbind(x1)
  dat3 <- my_simulateTwoCauseFineGrayModel(
    N_pop*ph[3], pwbeta1 = lapply(b, '[', 3), beta2, ucp = .7, X,
    u.min = u1, u.max = u2, p = 0.7) %>% 
    bind_cols(X) %>% 
    mutate(clust = rep(1:(N_pop*ph[3]/psu_size),
                       each = psu_size)[rank(ftime)])
  
  x1 <- rep(0:1, N_pop*ph[4]*c(.5, .5)); X <- cbind(x1)
  dat4 <- my_simulateTwoCauseFineGrayModel(
    N_pop*ph[4], pwbeta1 = lapply(b, '[', 4), beta2, ucp = .7, X,
    u.min = u1, u.max = u2, p = 0.7) %>% 
    bind_cols(X) %>% 
    mutate(clust = rep(1:(N_pop*ph[4]/psu_size),
                       each = psu_size)[rank(ftime)])
  
  x1 <- rep(0:1, N_pop*ph[5]*c(.5, .5)); X <- cbind(x1)
  dat5 <- my_simulateTwoCauseFineGrayModel(
    N_pop*ph[5], pwbeta1 = lapply(b, '[', 5), beta2, ucp = .7, X,
    u.min = u1, u.max = u2, p = 0.7) %>% 
    bind_cols(X) %>% 
    mutate(clust = rep(1:(N_pop*ph[5]/psu_size),
                       each = psu_size)[rank(ftime)])
  
  dat <- bind_rows(dat1, dat2, dat3, dat4, dat5) %>% 
    mutate(strat = rep(1:length(ph), N_pop*ph))
  
  # table(dat$fstatus) %>% prop.table()
  
  dat_srs <-  sample_n(dat, n)
  
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

out <- c(sum(!is.na(B$pv.g1.svy)),
         mean(B$pv.g1.srs < 0.05, na.rm = T),
         mean(B$pv.g1.s < 0.05, na.rm = T),
         mean(B$pv.g1.svy < 0.05, na.rm = T))

names(out) <- c('n_rep', 'gray_srs', 'gray_s', 'gray_svy')
