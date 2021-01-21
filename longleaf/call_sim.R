rm(list = ls()); gc()

library(tidyverse)
library(survey)
library(doParallel)

source('svygray.R')
source('compCIF.R')
source('aux_code/my_TwoCauseFineGray.R')

# On Windows, set cores to be 1
if (.Platform$OS.type == 'windows') {
  cores = 1
} else {
  if (.Platform$GUI == 'RStudio') {
    cores = detectCores()
  } else {
    cores = 24
  }
}

# -------------------------------------------------------------------------

ENV_SIMID <- 11
ENV_SEED <- 42867464
ENV_CENS <- 20
ENV_S_SIZE <- 1000
ENV_PSU_SIZE <- 5
ENV_EFFECT <- 2

ENV_REP <- 5
ENV_TEST <- 'pepe'

# ENV_SIMID <- Sys.getenv('ENV_SIMID') %>% as.numeric()
# ENV_SEED <- Sys.getenv('ENV_SEED') %>% as.numeric()
# ENV_CENS <- Sys.getenv('ENV_CENS') %>% as.numeric()
# ENV_S_SIZE <- Sys.getenv('ENV_S_SIZE') %>% as.numeric()
# ENV_PSU_SIZE <- Sys.getenv('ENV_PSU_SIZE') %>% as.numeric()
# ENV_EFFECT <- Sys.getenv('ENV_EFFECT') %>% as.numeric()
# 
# ENV_REP <- Sys.getenv('ENV_REP') %>% as.numeric()
# ENV_TEST <- Sys.getenv('ENV_TEST')

# -------------------------------------------------------------------------

TASK_ID <- Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.integer()

# set a 'master' seed
RNGkind("L'Ecuyer-CMRG")
set.seed(ENV_SEED)

seed <- .Random.seed
for (i in 1:TASK_ID) {
  seed <- parallel::nextRNGStream(seed)
}
.GlobalEnv$.Random.seed <- seed

# -------------------------------------------------------------------------

N_pop <- 10000
beta2 <- cbind(1)

ph <- c(0.5, 0.25, 0.15, 0.05, 0.05)
psu_size <- ENV_PSU_SIZE

if(ENV_CENS == 20){
  u1 <- 1; u2 <- 2
} else {
  u1 <- 0; u2 <- 1.5
}

n <- ENV_S_SIZE
f <- rep(1/length(ph), length(ph))

if(ENV_EFFECT == 0){
  b <- list(c(0,0,0,0,0),
            c(0,0,0,0,0))
}

if(ENV_EFFECT == 1){
  b <- list(c(.3,.25,.1,.05,0),
            c(.3,.25,.1,.05,0))
}

if(ENV_EFFECT == 2){
  b <- list(c(.3,.3,.3,0,0),
            c(-1,-1,-1,0,0))
}

# -------------------------------------------------------------------------

JOB_ID <- Sys.getenv('SLURM_ARRAY_JOB_ID') %>% as.integer()

sim <- c(ENV_SIMID, ENV_TEST, TASK_ID)
names(sim) <- c('ENV_SIMID', 'ENV_TEST', 'TASK_ID')

if (ENV_TEST == 'pepe'){
  # Pepe
  source('./longleaf/pepe_nph.R')
}

if (ENV_TEST == 'gray'){
  # Gray
  source('./longleaf/gray_nph.R')
}

sim <- t(c(sim, out))
write_csv(as.data.frame(sim),
          paste0('./longleaf/slurm_out/',
                 toupper(ENV_TEST), '_',
                 sprintf("%02d", ENV_SIMID), '_JOB_',
                 JOB_ID, '.csv'))
