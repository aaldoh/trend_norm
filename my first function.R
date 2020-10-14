test.bfrr <- function (sample_mean = 0, sample_se = 0.1, sample_df = 99, model = "normal", 
          ..., criterion = 3, rr_interval = NA, precision = 0.05) 
{
  if (!(tolower(model) %in% c("normal", "norm", "n", "uniform", 
                              "unif", "u"))) {
    stop("The model type '", model, "' is not supported")
  }
  params <- list(...)
  model <- toupper(substr(model, 1, 1))
  x <- list(null = likelihood(sample_mean, sample_se, sample_df), 
            theory = likelihood(sample_mean, sample_se, sample_df, 
                                model, ...))
  x$BF <- x$theory/x$null
  if (x$BF >= criterion) {
    conclusion <- "H1"
  }
  else if (x$BF <= 1/criterion) {
    conclusion <- "H0"
  }
  else {
    conclusion <- "no"
  }
  if (model == "N") {
    params$mean <- default(params$mean, 0)
    params$sd <- default(params$sd, 1)
    params$tail <- default(params$tail, 2)
    sd_intv <- default(rr_interval$sd, c(0, 2 * params$sd))
    rr_sd <- seq(sd_intv[1], sd_intv[2], precision)
    rr_sd <- c(rr_sd[rr_sd > 0], params$sd)
    mean_intv <- default(rr_interval$mean, c(0, 2 * params$mean))
    rr_mean <- c(seq(mean_intv[1], mean_intv[2], precision), 
                 params$mean)
    rr_params <- expand.grid(M = rr_mean, SD = rr_sd)
    rr_BF <- purrr::map2_dbl(rr_params$M, rr_params$SD, function(rr_m, 
                                                                 rr_s) {
      likelihood(sample_mean, sample_se, sample_df, model, 
                 mean = rr_m, sd = rr_s, tail = params$tail)/x$null
    })
    H1 <- sprintf("%sN(%g, %g)", ifelse(identical(tail, 1), 
                                        "H", ""), params$mean, params$sd)
  }
  else if (model == "U") {
    lower <- default(params$lower, 0)
    upper <- default(params$upper, 1)
    lower_intv <- default(rr_interval$lower, c(lower, sample_mean - 
                                                 precision))
    upper_intv <- default(rr_interval$upper, c(upper, sample_mean + 
                                                 precision))
    rr_lower <- seq(lower_intv[1], lower_intv[2], precision)
    rr_upper <- seq(upper_intv[1], upper_intv[2], precision)
    rr_params <- expand.grid(lower = rr_lower, upper = rr_upper)
    rr_BF <- purrr::map2_dbl(rr_params$lower, rr_params$upper, 
                             function(rr_l, rr_u) {
                               likelihood(sample_mean, sample_se, sample_df, 
                                          model, lower = rr_l, upper = rr_u)/x$null
                             })
    H1 <- sprintf("U(%g, %g)", lower, upper)
  }
  dat <- rr_params
  dat$BF <- rr_BF
  dat <- dplyr::distinct(dat)
  dat$support <- dplyr::case_when(dat$BF >= criterion ~ "H1", 
                                  dat$BF <= 1/criterion ~ "H0", TRUE ~ "no")
  if (model == "N") {
    rrM <- dplyr::filter(dat, SD == params$sd) %>% dplyr::arrange(M) %>% 
      dplyr::mutate(i = support != dplyr::lag(support), 
                    i = ifelse(dplyr::row_number() == 1, TRUE, i), 
                    grp = cumsum(i))
    mgrp <- dplyr::filter(rrM, M == params$mean)$grp[1]
    rrSD <- dplyr::filter(dat, M == params$mean) %>% dplyr::arrange(SD) %>% 
      dplyr::mutate(i = support != dplyr::lag(support), 
                    i = ifelse(dplyr::row_number() == 1, TRUE, i), 
                    grp = cumsum(i))
    sdgrp <- dplyr::filter(rrSD, SD == params$sd)$grp[1]
    rr <- list(mean = range(dplyr::filter(rrM, grp == mgrp)$M), 
               sd = range(dplyr::filter(rrSD, grp == sdgrp)$SD))
  }
  else if (model == "U") {
    rr <- list(lower = range(dplyr::filter(dat, upper == 
                                             upper, support == conclusion)$lower), upper = range(dplyr::filter(dat, 
                                                                                                               lower == lower, support == conclusion)$upper))
  }
  
    precision <- dat[1,2] - dat[1, 1]
    half <- ifelse(identical(params$tail, 1), "H", "")
    
    RR <- sapply(names(rr), function(n) {
      rng <- rr[[n]]
      if (unique(rng) %>% length() > 1) {
        dplyr::case_when(
          n == "mean" ~ paste0(half, "N([", toString(rng), "], ", params$sd, ")"),
          n == "sd" ~ paste0(half, "N(", params$mean, ", [", toString(rng), "])"),
          n == "lower" ~ paste0("U([", toString(rng), "], ", params$upper, ")"),
          n == "upper" ~ paste0("U(", params$lower, ", [", toString(rng), "])")
        )
      } else {
        "no"
      }
      
    })
    
    test.RR <- RR[RR != "no"] %>% paste(collapse = "; ")
  
  ret <- list(theory = x$theory, null = x$null, BF = x$BF, RR = rr, 
              RR_apa = test.RR, conclusion = conclusion, sample_mean = sample_mean, 
              sample_se = sample_se, sample_df = sample_df, model = model, 
              params = params, criterion = criterion, H0_model = sprintf("T(%g)", 
                                                                         sample_df), H1_model = H1)
  class(ret) <- c("bfrr", "list")
  invisible(ret)
}
