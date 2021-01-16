my_simulateTwoCauseFineGrayModel <-
  function (nobs, pwbeta1, beta2, ucp, X = NULL,
            u.min = 0, u.max, p = 0.5, returnX = FALSE) {
    
    beta1 <- pwbeta1[[1]]
    beta1_2 <- pwbeta1[[2]]
    
    if (length(beta1) != length(beta2)) 
      stop("Dimension of beta1 and beta2 should be the same")
    ncovs <- length(beta1)
    if (is.null(X)) {
      X <- matrix(rnorm(nobs * ncovs), nrow = nobs)
      returnX <- TRUE
    }
    
    ftime <- numeric(nobs)
    c.ind <- 1 + rbinom(nobs, 1, prob = (1 - p)^exp(X %*% beta1))
    # eta1 <- X[c.ind == 1, ] %*% beta1
    # eta2 <- X[c.ind == 2, ] %*% beta2
    eta1 <- X %*% beta1
    eta2 <- X %*% beta2
    u1 <- runif(length(eta1))
    t1 <- -log(1 - (1 - (1 - u1 * (1 - (1 - p)^exp(eta1)))^(1/exp(eta1)))/p)
    
    c.ind_2 <- 1 + rbinom(nobs, 1, prob = (1 - p)^exp(X %*% beta1_2))
    eta1_2 <- X %*% beta1_2
    u1_2 <- runif(length(eta1_2))
    t1_2 <- -log(1 - (1 - (1 - u1_2 * (1 - (1 - p)^exp(eta1_2)))^(1/exp(eta1_2)))/p)

    t1 <- ifelse(u1 < ucp,
                 -log(1 - (1 - (1 - u1 * (1 - (1 - p)^exp(eta1)))^(1/exp(eta1)))/p),
                 -log(1 - (1 - (1 - u1 * (1 - (1 - p)^exp(eta1_2)))^(1/exp(eta1_2)))/p))
    
    # c.ind[c.ind == 1] <- ifelse(t1[c.ind == 1] < tcp, c.ind[c.ind == 1], c.ind_2[c.ind == 1])
    # t1 <- ifelse(t1 < tcp, t1, t1_2)
    
    t2 <- rexp(length(eta2), rate = exp(eta2))
    ci <- runif(nobs, min = u.min, max = u.max)
    
    ftime[c.ind == 1] <- t1[c.ind == 1]
    ftime[c.ind == 2] <- t2[c.ind == 2]
    ftime <- pmin(ftime, ci)
    fstatus <- ifelse(ftime == ci, 0, 1)
    fstatus <- fstatus * c.ind
    out <- list()
    out$ftime <- ftime
    out$fstatus <- fstatus
    if (returnX) 
      out$X <- X
    return(out)
  }