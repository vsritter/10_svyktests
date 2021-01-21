
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
    means <- svytotal(z1i[[1]], subset(dsgn, rank(ftime) <= z1i[[5]]))
    out <- data.frame(svy.s = coef(means),
                      w.s = z1i[[3]],
                      svy.v = SE(means)^2)
    
    # weighted
    z1i2 <- with(dsgn$variables,
                 compCIF(ftime, fstatus, x1, weights = w,
                         return_si = T))
    means2 <- svytotal(z1i2[[1]], subset(dsgn, rank(ftime) <= z1i2[[5]]))
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

out <- c(
  # n
  nrow(complete.cases(B)),
  # srs
  mean(1-pchisq(B$s^2/B$v, 1) < 0.05, na.rm = T),
  # s
  mean(1-pchisq(B$s.1^2/B$v.1, 1) < 0.05, na.rm = T),
  # svy
  mean(1-pchisq(B$svy.s^2/B$svy.v, 1) < 0.05, na.rm = T),
  # weighted svy
  mean(1-pchisq(B$svy2.s^2/B$svy2.v, 1) < 0.05, na.rm = T))

names(out) <- c('n', 'pepe_srs', 'pepe_s', 'pepe_svy', 'pepe_w_svy')
