# Load packages ---------------------------
library(psych)
library(tidyverse)
library(papaja)

options(mc.cores = parallel::detectCores()) ## Run chains in parallel
export <- haven::read_spss("data/01_EXPORTED.sav")

# Data cleaning ---------------------------
export[, c(2:5, 15, 17)] <- haven::as_factor(export[, c(2:5, 15, 17)])

raw <- export %>%
  mutate(att_checkbi = ifelse(att_check == 3, 1, 0),
         condition = factor(condition, labels = c("Dynamic", "Static", "No norm")),
         conditionbi = na_if(condition, "No norm") %>% droplevels(),
         cons_proj = select(., cons_now_perc_1, cons_next_perc_1, cons_six_perc_1) %>% rowMeans(., na.rm = TRUE),
         genderbi = na_if(gender, "Other") %>% droplevels(),
         conformity_3 = 8 - conformity_3,
         conformity_mean = rowMeans(select(., starts_with("conformity_")), na.rm = TRUE)) %>%
  rename_with(., .fn = ~ str_remove(., "_1"), .cols = c(interest_1:politics_1, meat_cons_1))

# excluding vegetarians and attention check fails
complete <- raw %>% 
  filter(att_check == 1 | veg != 1)

# excluding unneeded fields
clean <- complete %>%
  select(-RecordedDate, -veg, -att_check, -att_checkbi)

# Data Overview ---------------------------
describeBy(clean, clean$condition)

# outliers
mahalfiltered = mahalanobis(clean[ , -c(1:4, 14, 17:22, 23, 25)],
                            colMeans(clean[ , -c(1:4, 14, 17:22,23, 25)], na.rm = T),
                            cov(clean[ , -c(1:4, 14, 17:22,23, 25)]))
cutoff = qchisq(1-.001, ncol(clean[ , -c(1:4, 14, 17:22,23, 25)]))
ncol(clean[ , -c(1:4, 14, 17:22,23, 25)]) #df
no_out = subset(clean, mahalfiltered < cutoff)

# participants
data_desc <- c(total_n = nrow(raw),
               clean_n = nrow(clean))

# age
age_desc <- clean %>%
  summarise(min_age = min(age),
            max_age = max(age),
            m_age = printnum(mean(age)),
            sd_age = printnum(sd(age)))

# gender
gender_freq <- round(100 * prop.table(table(clean$gender)), digits = 2)

# outcomes
out_plots <- clean %>%
  pivot_longer(cols = interest:expectation, names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = condition, y = value)) +
  facet_wrap(~ variable) +
  geom_violin(trim = FALSE) + 
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), geom = "pointrange", color = "black") +
  ggtitle("Distribution of outcome variables by condition")

# Randomization check ---------------------------

# age
agebycond <- clean %>%
  group_by(condition) %>%
  summarise(age_m = mean(age),
            age_sd = sd(age), .groups = "rowwise") %>%
  printnum()

age_check <- paste(agebycond$age_m, agebycond$age_sd, sep = ' $\\pm$ ')
age_stat <- apa_print(aov(age ~ condition, clean))

# political position
polbycond <- clean %>%
  group_by(condition) %>%
  summarise(pol_m = mean(politics), pol_sd = sd(politics), .groups = "rowwise") %>%
  printnum()

pol_check <- paste(polbycond$pol_m, polbycond$pol_sd, sep = ' $\\pm$ ')
pol_stat <- apa_print(aov(politics ~ condition, clean))

# gender
gender_check <- clean %>%
  count(condition, gender) %>%
  group_by(condition) %>%
  mutate(rel.freq = paste0(round(100 * n / sum(n), 2), "\\%")) %>%
  summarise(Gender = paste0(gender, " (", rel.freq, ")", collapse = "\n"), .groups = "rowwise") %>%
  pull(Gender)

gender_stat <- apa_print(chisq.test(clean$condition, clean$gender), n = nrow(clean))

# nation
nation_check <- clean %>%
  count(condition, country) %>%
  group_by(condition) %>%
  mutate(rel.freq = round(100 * n / sum(n), 2)) %>%
  summarise(Nationality = paste0(country, " (", rel.freq, "\\%)", collapse = "\n"), .groups = "rowwise") %>%
  pull(Nationality)

nation_stat <- apa_print(chisq.test(clean$condition, clean$country), n = nrow(clean))

# table
random_check <- tibble(Item = c("Age (years)", "Gender (\\%)", "Political position", "Home country (\\%)"),
                       Dynamic = c(age_check[1], gender_check[1], pol_check[1], nation_check[1]),
                       Static = c(age_check[2], gender_check[2], pol_check[2], nation_check[2]),
                       None = c(age_check[3], gender_check[3], pol_check[3], nation_check[3]),
                       Test = c(age_stat$statistic$condition, gender_stat$statistic, pol_stat$statistic$condition, nation_stat$statistic)) %>%
  mutate_all(linebreak)

# Correlation matrix ---------------------------
# reliability
cron <- apply(matrix(6:14, ncol = 3), 2, function(x) printnum(cronbach(clean[x])$alpha))

# correlation matrix
measure_sum <- clean %>%
  summarise(m_interest = mean(interest),
            m_attitude = mean(attitude),
            m_expect = mean(expectation),
            m_intent = mean(intention),
            sd_interest = sd(interest),
            sd_attitude = sd(attitude),
            sd_expect = sd(expectation),
            sd_intent = sd(intention)) %>%
  printnum()

mcor <- clean %>% select("interest", "attitude", "expectation", "intention", "cons_proj", "conformity_mean") %>% corstars()

measure_sum <- tibble(Measure = c("1. Interest", "2. Attitudes", "3. Expectations", "4. Intentions", "5. Projected consumption", "6. Conformity"),
                      Mean = unlist(measure_sum[, 1:6]),
                      SD = unlist(measure_sum[, 7:12]),
                      Alpha = c("-", cron, "-", "-"))

cortable <- cbind(measure_sum, mcor) %>% as_tibble()

# H1 ---------------------------
outcomes_desc <- clean %>%
  group_by(condition) %>%
  summarise(n = n(),
            m_interest = mean(interest),
            m_attitude = mean(attitude),
            m_expect = mean(expectation),
            m_intent = mean(intention),
            sd_interest = sd(interest),
            sd_attitude = sd(attitude),
            sd_expect = sd(expectation),
            sd_intent = sd(intention), .groups = "rowwise")

# effect sizes
H1.effect <- list(
  DYST = printnum(d.ind.t(interest_desc[[1, 3]], interest_desc[[2, 3]], interest_desc[[1, 4]], interest_desc[[2, 4]], interest_desc[[1, 2]], interest_desc[[2, 2]], a = .05)),
  DYNO = printnum(d.ind.t(interest_desc[[1, 3]], interest_desc[[3, 3]], interest_desc[[1, 4]], interest_desc[[3, 4]], interest_desc[[1, 2]], interest_desc[[3, 2]], a = .05)),
  STNO = printnum(d.ind.t(interest_desc[[2, 3]], interest_desc[[3, 3]], interest_desc[[2, 4]], interest_desc[[3, 4]], interest_desc[[2, 2]], interest_desc[[3, 2]], a = .05)),
  DYST.Bf = Bf(sd = ls_interest$SE[1], obtained = ls_interest$estimate[1], dfdata = ls_interest$df[1], meanoftheory = 0,
               sdtheory = 0.69, dftheory = 10^10, tail = 1),
  rr      = bfrr(sample_mean = ls_interest$estimate[1], sample_se = ls_interest$SE[1], sample_df = ls_interest$df[1], model = "normal", mean = 0, sd = 0.69, tail = 1, criterion = 5,
                 rr_interval = list(mean = c(-2, 2), sd = c(0, 2)), precision = 0.05))

demoreg.out$table <- demoreg.out$table %>%
  mutate(predictor = c("Intercept", "Condition", "Gender", "Political position", "Condition $\\times$ Gender", "Condition $\\times$ Political position"),
         rownames = NULL)

# H2 ---------------------------

# Future consumption
future_cons <- clean %>%
  group_by(condition) %>%
  summarise(n = n(),
            Mean = mean(PERCEPTCHANGE, na.rm = TRUE),
            SD = sd(PERCEPTCHANGE), .groups = "rowwise")

ls_change <- lm(PERCEPTCHANGE ~ condition, clean) %>%
  emmeans(., "condition", contr = contrasts) %>%
  .$contrasts %>%
  c(summary(.), apa_print(.))

ls_change$table <- ls_change$table %>%
  mutate(contrast = c("Dynamic, static", "Dynamic, control", "Static, control", "Dynamic, both", "Norms, control"),
         rownames = NULL)

H2change.effect <- list(
  DYST = printnum(d.ind.t(future_cons[[1, 3]], future_cons[[2, 3]], future_cons[[1, 4]], future_cons[[2, 4]], future_cons[[1, 2]], future_cons[[2, 2]], a = .05)),
  DYNO = printnum(d.ind.t(future_cons[[1, 3]], future_cons[[3, 3]], future_cons[[1, 4]], future_cons[[3, 4]], future_cons[[1, 2]], future_cons[[3, 2]], a = .05)),
  STNO = printnum(d.ind.t(future_cons[[2, 3]], future_cons[[3, 3]], future_cons[[2, 4]], future_cons[[3, 4]], future_cons[[2, 2]], future_cons[[3, 2]], a = .05)),
  DYST.Bf = Bf(sd = ls_change$SE[1], obtained = ls_change$estimate[1], dfdata = ls_change$df[1], meanoftheory = 0, sdtheory = 0.40, dftheory = 10^10, tail = 1),
  rr      = bfrr(sample_mean = ls_change$estimate[1], sample_se = ls_change$SE[1], sample_df = ls_change$df[1], model = "normal", mean = 0, sd = 0.40, tail = 1, criterion = 5,
                 rr_interval = list(mean = c(-2, 2), sd = c(0, 2)), precision = 0.05))

# Preconformity

# summary table
preconformity <- clean %>%
  dplyr::group_by(condition) %>%
  dplyr::summarise(n = n(),
                   Mean = mean(PRECONFORMITY),
                   SD = sd(PRECONFORMITY), .groups = "rowwise")

# contrasts
ls_preconformity <- lm(PRECONFORMITY ~ condition, clean) %>%
  emmeans(., "condition", contr = contrasts) %>%
  .$contrasts %>%
  c(summary(.), apa_print(.))

ls_preconformity$table <- ls_preconformity$table %>%
  mutate(contrast = c("Dynamic, static", "Dynamic, control", "Static, control", "Dynamic, both", "Norms, control"),
         rownames = NULL)

H2preconformity.effect <- list(
  DYST = printnum(d.ind.t(preconformity[[1, 3]], preconformity[[2, 3]], preconformity[[1, 4]], preconformity[[2, 4]], preconformity[[1, 2]], preconformity[[2, 2]], a = .05)),
  DYNO = printnum(d.ind.t(preconformity[[1, 3]], preconformity[[3, 3]], preconformity[[1, 4]], preconformity[[3, 4]], preconformity[[1, 2]], preconformity[[3, 2]], a = .05)),
  STNO = printnum(d.ind.t(preconformity[[2, 3]], preconformity[[3, 3]], preconformity[[2, 4]], preconformity[[3, 4]], preconformity[[2, 2]], preconformity[[3, 2]], a = .05)),
  DYST.Bf = Bf(sd = ls_preconformity$SE[1], obtained = ls_preconformity$estimate[1], dfdata = ls_preconformity$df[1], meanoftheory = 0, sdtheory = 0.40, dftheory = 10^10, tail = 1),
  rr      = bfrr(sample_mean = ls_preconformity$estimate[1], sample_se = ls_preconformity$SE[1], sample_df = ls_preconformity$df[1], model = "normal", mean = 0, sd = 0.40, tail = 1, criterion = 5,
                 rr_interval = list(mean = c(-2, 2), sd = c(0, 2)), precision = 0.05))



# Secondary analyses ---------------------------

# Perception change by condition on numerican 1-100 scale
PERCEPTNUM_cond <- clean %>%
  group_by(condition) %>%
  summarise(n = n(),
            Mean = mean(PERCEPTNUM),
            SD = sd(PERCEPTNUM), .groups = "rowwise")

TOSTtwo.sci(m1 = PERCEPTNUM_cond[[1, 3]], m2 = PERCEPTNUM_cond[[2, 3]], sd1 = PERCEPTNUM_cond[[1, 4]], sd2 = PERCEPTNUM_cond[[2, 4]], n1 = PERCEPTNUM_cond[[1, 2]], n2 = PERCEPTNUM_cond[[2, 2]], low_eqbound = -5, high_eqbound = 5)

# Perception change by condition on Likert scale
PERCEPTSCALE_cond <- clean %>%
  group_by(condition) %>%
  summarise(n = n(),
            Mean = mean(PERCEPTSCALE),
            SD = sd(PERCEPTSCALE), .groups = "rowwise")

# construal of consumption
construal_cond <- clean %>%
  group_by(condition) %>%
  summarise(n = n(),
            Mean = mean(CONSTRUAL_1),
            SD = sd(CONSTRUAL_1), .groups = "rowwise") %>%
  mutate_if(is.numeric, round, 2)

# H3 ---------------------------
## Does dynamic norm (versus static or no norm) information increase participants’ positive attitude, intentions, and expectations to reduce their meat consumption?
simple.mod <- '
interest       ~ conditionbi
attitude  ~ conditionbi
intention ~ conditionbi
expectation    ~ conditionbi'

### Fitting models
simpleuninf.fit <- bsem(model = simple.mod, data = clean, target = "stan", seed = 2019) # uninformative
simplelowinf.fit <- bsem(model = simple.mod, data = clean, target = "stan", seed = 2019, dp = dpriors(beta = "normal(0.5,0.447)")) # informative
simplehighinf.fit <- bsem(model = simple.mod, data = clean, target = "stan", seed = 2019, dp = dpriors(beta = "normal(0.5, 0.316)")) # informative

### Summary outputs
h3.out <- list(lowinf = simplelowinf.fit,
               highinf = simplehighinf.fit,
               uninf = simpleuninf.fit) %>%
  lapply(., bsem.summary)

h3.table <- cbind("Parameter" = rbind("Interest", "Attitude", "Intention", "Expectation"),
                  "Mean (SD)" = paste0(h3.out$uninf[, 2], " (", h3.out$uninf[, 3], ")"),
                  "95% PPI" = paste(h3.out$uninf[, 4], h3.out$uninf[, 5], sep = ", "),
                  "neff" = h3.out$uninf[, 7],
                  "PSRF" = round(blavInspect(simpleuninf.fit, "psrf"), 3),
                  "Prior" = h3.out$uninf[, 8],
                  "Mean (SD)" = paste0(h3.out$lowinf[, 2], " (", h3.out$lowinf[, 3], ")"),
                  "95% PPI" = paste(h3.out$lowinf[, 4], h3.out$lowinf[, 5], sep = ", "),
                  "neff" = h3.out$lowinf[, 7],
                  "PSRF" = round(blavInspect(simplelowinf.fit, "psrf"), 3),
                  "Prior" = h3.out$lowinf[, 8],
                  "Mean (SD)" = paste0(h3.out$highinf[, 2], " (", h3.out$highinf[, 3], ")"),
                  "95% PPI" = paste(h3.out$highinf[, 4], h3.out$highinf[, 5], sep = ", "),
                  "neff" = h3.out$highinf[, 7],
                  "PSRF" = round(blavInspect(simplehighinf.fit, "psrf"), 3),
                  "Prior" = h3.out$highinf[, 8])

### Estimate bias
h3.estimates <- lapply(h3.out, function(x) as.numeric(x[, 2][1:4]))
h3.bias <- c(lowinf_uninf = round(100 * ((h3.estimates$lowinf - h3.estimates$uninf) / h3.estimates$uninf), 2),
             highinf_uninf = round(100 * ((h3.estimates$highinf - h3.estimates$uninf) / h3.estimates$uninf), 2))

posterior1.2.3 <- list(interest = bind_rows("uninformative prior"      = enframe(as.matrix(blavInspect(simpleuninf.fit, what = "mcmc"))[, 'bet_sign[1]']),
                                            "informative prior (low)"  = enframe(as.matrix(blavInspect(simplelowinf.fit, what = "mcmc"))[, 'bet_sign[1]']),
                                            "informative prior (high)" = enframe(as.matrix(blavInspect(simplehighinf.fit, what = "mcmc"))[, 'bet_sign[1]']),
                                            .id = "id1"),
                       attitude = bind_rows("uninformative prior"      = enframe(as.matrix(blavInspect(simpleuninf.fit, what = "mcmc"))[, 'bet_sign[2]']),
                                            "informative prior (low)"  = enframe(as.matrix(blavInspect(simplelowinf.fit, what = "mcmc"))[, 'bet_sign[2]']),
                                            "informative prior (high)" = enframe(as.matrix(blavInspect(simplehighinf.fit, what = "mcmc"))[, 'bet_sign[2]']),
                                            .id = "id1"),
                       intention = bind_rows("uninformative prior"      = enframe(as.matrix(blavInspect(simpleuninf.fit, what = "mcmc"))[, 'bet_sign[3]']),
                                             "informative prior (low)"  = enframe(as.matrix(blavInspect(simplelowinf.fit, what = "mcmc"))[, 'bet_sign[3]']),
                                             "informative prior (high)" = enframe(as.matrix(blavInspect(simplehighinf.fit, what = "mcmc"))[, 'bet_sign[3]']),
                                             .id = "id1"),
                       expectation = bind_rows("uninformative prior"      = enframe(as.matrix(blavInspect(simpleuninf.fit, what = "mcmc"))[, 'bet_sign[4]']),
                                               "informative prior (low)"  = enframe(as.matrix(blavInspect(simplelowinf.fit, what = "mcmc"))[, 'bet_sign[4]']),
                                               "informative prior (high)" = enframe(as.matrix(blavInspect(simplehighinf.fit, what = "mcmc"))[, 'bet_sign[4]']),
                                               .id = "id1"))

prior1.2.3 <- bind_rows("uninformative prior" = enframe(rnorm(10000, mean = 0, sd = sqrt(1 / 1e-2))),
                        "informative prior (low)" = enframe(rnorm(10000, mean = 0.5, sd = sqrt(0.2))),
                        "informative prior (high)" = enframe(rnorm(10000, mean = 0.5, sd = sqrt(0.01))),
                        .id = "id1") # here we sample a large number of values from the prior distributions to be able to plot them.

h3priors.posterior <- list(interest = bind_rows("posterior" = posterior1.2.3$interest, "prior" = prior1.2.3, .id = "id2"),
                           attitude = bind_rows("posterior" = posterior1.2.3$attitude, "prior" = prior1.2.3, .id = "id2"),
                           intention = bind_rows("posterior" = posterior1.2.3$intention, "prior" = prior1.2.3, .id = "id2"),
                           expectation = bind_rows("posterior" = posterior1.2.3$expectation, "prior" = prior1.2.3, .id = "id2"))

h3.plots <- lapply(seq_along(h3priors.posterior), function(x) bsem.plot(h3priors.posterior, names(h3priors.posterior)[x]))
ggarrange(plotlist = h3.plots, common.legend = TRUE, legend = "bottom", labels = c("A", "B", "C", "D"))
