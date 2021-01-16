draw_sample <- function(pop, n, type, u_cens){
  pars <- pop$pars
  dt_pop <- pop$data

  # -----------------------------------------------------------------------------
  # Draw sample
  # -----------------------------------------------------------------------------
  if(type == "clust") {
    n_strata <- nrow(pars)

    pars <- pars %>%
      mutate(nh = n/n_strata,
             n_clust_sample = nh/Nc)

    # Select <nc> or more clusters from each strata
    draw <- pars %>%
      group_by(strata) %>%
      do(c = sample.int(.$n_clust, size = .$n_clust_sample, replace = F) + .$ini_clust)

    # Filter census data
    dt_s <- dt_pop %>%
      filter(clust %in% unlist(draw$c)) %>%
      mutate(id = as.numeric(as.factor(id)))

    # Sampling weights
    s <- dt_s %>%
      # select(id, clust, strata, domain) %>%
      filter(!duplicated(.)) %>%
      mutate(prob = with(pars, rep(n_clust_sample/n_clust, times = Nc * n_clust_sample)),
             w = 1/prob)

    # Define design
    dsgn <- survey::svydesign(id = ~clust,
                              strata = ~strata,
                              weights = ~w, data = s)
  }

  if (type == "strat") {
    n_strata <- nrow(pars)

    pars <- pars %>%
      mutate(nh = n/n_strata)

    draw <- pars %>%
      group_by(strata) %>%
      do(id = sample.int(.$Nh, size = .$nh, replace = F) + .$ini_id)

    # Filter census data
    dt_s <- dt_pop %>%
      filter(id %in% unlist(draw$id)) %>%
      mutate(id = as.numeric(as.factor(id)),
             clust = NULL)

    # Sampling weights
    s <- dt_s %>%
      # select(id, strata, domain) %>%
      filter(!duplicated(.)) %>%
      mutate(prob = with(pars, rep(nh/Nh, times = nh)),
             w = 1/prob)

    # Define design
    dsgn <- survey::svydesign(id = ~1,
                              strata = ~strata,
                              weights = ~w, data = s)
  }

  if (type == "srs") {
    N <- sum(pars$Nh)

    draw <- pars %>%
      do(id = sample.int(N, size = n, replace = F))

    # Filter census data
    dt_s <- dt_pop %>%
      filter(id %in% unlist(draw$id)) %>%
      mutate(id = as.numeric(as.factor(id)),
             strata = NULL, clust = NULL)

    # Sampling weights
    s <- dt_s %>%
      # select(id, domain) %>%
      filter(!duplicated(.)) %>%
      mutate(prob = n/N, w = 1/prob)

    # Define design
    dsgn <- survey::svydesign(id = ~1,
                              weights = ~w, data = s)
  }

  # filtering domain
  dsgn <- subset(dsgn, domain == 1)

  dt_s <- dt_s %>%
    filter(domain == 1) %>%
    mutate(id = as.numeric(as.factor(id)))

  return(list(dt_s = dt_s, dsgn = dsgn))
}
