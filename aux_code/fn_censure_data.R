censure_data <- function(data, pars) {
  data %>%
    left_join(select(pars, strata, u1, u2), by = c("strata")) %>%
    group_by(strata, id) %>%
    mutate(c = runif(1, min = u1, max = u2)) %>%
    ungroup() %>%
    mutate(
      to = ifelse(exit > c, "cens", to),
      # d = ifelse(exit > c, 1, 0),
      entry = ifelse(entry > c, NA, entry),
      exit = ifelse(exit > c, c, exit)) %>%
    filter(!is.na(entry)) %>%
    select(-u1, -u2, -c)
}

censure_data_d <- function(data, pars) {
  data %>%
    left_join(select(pars, strata, u1, u2), by = "strata") %>%
    group_by(strata, id) %>%
    mutate(c = runif(1, min = u1, max = u2)) %>%
    ungroup() %>%
    mutate(
      # to = ifelse(exit > c, "cens", to),
      d = ifelse(exit > c, 1, 0),
      entry = ifelse(entry > c, NA, entry),
      exit = ifelse(exit > c, c, exit)) %>%
    filter(!is.na(entry)) %>%
    select(-u1, -u2, -c)
}
