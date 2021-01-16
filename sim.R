

library(tidyverse)
library(survival)


set.seed(643372)
nobs <- 10000
# beta1 <- cbind(c(0,0))
# beta2 <- cbind(c(1,-.5))
beta1 <- cbind(0)
beta2 <- cbind(1)

x1 <- rep(0:1, each = nobs/2)
# x2 <- rnorm(nobs, mean = -.5*(1-x1))
# X <- cbind(x1, x2)
X <- cbind(x1)

z1 <- rnorm(nobs)
# z2 <- rnorm(nobs)
z2 <- rnorm(nobs, mean = -.5*(1-x1))
Z <- cbind(z1, z2)
# Z <- cbind(z1)

# ph <- c(0.50, 0.25, 0.15, 0.05, 0.05)
ph <- c(0.8, 0.2)
Nc <- 5

dat <- fastcmprsk::simulateTwoCauseFineGrayModel(
  nobs, beta1, beta2, X, u.min = 1, u.max = 2, p = 0.3)

dat <- bind_cols(dat, X, Z) %>% 
  arrange(x1) %>% 
  mutate(strat = rep(1:length(ph), nobs*ph)) #%>% 
  # arrange(x1 + z1, x2 + z2) %>% 
  # mutate(clust = rep(1:(nobs/Nc), each = Nc))

# fit <- fastcmprsk::fastCrr(fastcmprsk::Crisk(ftime, fstatus) ~ x1+strat, variance = F, data = dat)
# summary(fit)

km <- survfit(Surv(ftime, factor(fstatus)) ~ x1 + strat, se=F, data = dat)
# survminer::ggcompetingrisks(km, multiple_panels = FALSE)

ci_fit <- with(dat %>%
                 mutate(fstatus = paste0('f=', fstatus),
                        x1 = paste0('x1=', x1),
                        strat = paste0('st=', strat)), {
  cmprsk::cuminc(ftime, fstatus, paste(x1), cencode = 'f=0')
})

# dat %>%
#   mutate(fstatus = paste0('f-', fstatus),
#          x1 = paste0('x1-', x1),
#          strat = paste0('st-', strat),
#          rhs = paste(x1, strat)) %>% 
#   with(factor(rhs))

plot(ci_fit, col=rep(1:4, 2), lty = rep(1:2, each=4))



# ci_fit$Tests
km[1,]$strata
lines(km[1,], lwd=2, type = 's')

km$p0


str(km)


table(dat$fstatus)

ci_fit <- with(dat %>%
                 mutate(fstatus = paste0('f-', fstatus),
                        x1 = paste0('x1-', x1),
                        strat = paste0('st-', strat)) %>% 
                 group_by(fstatus) %>% 
                 sample_n(1500), {
                          cmprsk::cuminc(ftime, fstatus, paste(x1), cencode = 'f-0')
                        })

plot(ci_fit, col=rep(1:4, 2), lty = rep(1:2, each=4))
ci_fit$Tests

plot(ci_fit, col=rep(1:4, 2), lty = rep(1:2, each=4))


str(summary(km))

km$time
km$pstate %>% head()

km$p0
nobs








