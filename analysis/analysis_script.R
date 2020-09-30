# Load packages ---------------------------
library(psych)
library(tidyverse)
library(papaja)
library(lavaan)
library(rlang)

options(mc.cores = parallel::detectCores()) ## Run chains in parallel
export <- haven::read_spss("data/01_EXPORTED.sav")

# Data cleaning ---------------------------
export[, c(2:5, 15, 17)] <- haven::as_factor(export[, c(2:5, 15, 17)])

raw <- export %>%
  mutate(att_checkbi = ifelse(att_check == 3, 1, 0) %>% as_factor(.),
         condition = as_factor(export$condition) %>% set_attrs(., labels = c("dynamic" = 1, "static" = 2, "none" = 3)),
         cons_proj = select(., cons_now_perc_1, cons_next_perc_1, cons_six_perc_1) %>% rowMeans(., na.rm = TRUE),
         genderbi = na_if(gender, "Other") %>% droplevels() %>% recode_factor(., "Male" = -1, "Female" = 1),
         women = recode_factor(gender, "Male" = -1, "Female" = 1, "Other" = 0),
         conformity_3 = 8 - conformity_3,
         conformity_mean = rowMeans(select(., starts_with("conformity_")), na.rm = TRUE)) %>%
  rename_with(., .fn = ~ str_remove(., "_1"), .cols = c(interest_1:politics_1, meat_cons_1)) %>%
  cbind(., psych::dummy.code(.$condition))

# excluding vegetarians and attention check fails
complete <- raw %>% 
  filter(att_check == 1 | veg != 1)

# excluding unneeded fields
clean <- complete %>%
  select(-RecordedDate, -veg, -att_check, -att_checkbi)

# outliers
mahalfiltered = mahalanobis(clean[,c(5:8, 24, 26)], colMeans(clean[,c(5:8, 24, 26)], na.rm = T), cov(clean[,c(5:8, 24, 26)]))
cutoff = qchisq(1-.001, ncol(clean[,c(5:8, 24, 26)]))
ncol(clean[,c(5:8, 24, 26)]) #df
no_out = subset(clean, mahalfiltered < cutoff)

# Data Overview ---------------------------
describeBy(clean, clean$condition) # check distribution and normality

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
measure_sum <- clean %>%
  summarise(m_interest = mean(interest),
            m_attitude = mean(attitude),
            m_expect = mean(expectation),
            m_intent = mean(intention),
            m_cons = mean(meat_cons),
            m_cons_proj = mean(cons_proj),
            m_conformity = mean(conformity_mean),
            sd_interest = sd(interest),
            sd_attitude = sd(attitude),
            sd_expect = sd(expectation),
            sd_intent = sd(intention),
            sd_cons = sd(meat_cons),
            sd_cons_proj = sd(cons_proj),
            sd_conformity = sd(conformity_mean)) %>%
  printnum()

out_plots <- no_out %>%
  pivot_longer(cols = interest:expectation, names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = condition, y = value)) +
  facet_wrap(~ variable) +
  geom_violin(trim = FALSE) + 
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), geom = "pointrange", color = "black") +
  ggtitle("Distribution of outcome variables by condition")

# reliability
cron <- clean %>% select(conformity_1:conformity_6) %>% psych::alpha() # scale reliability

# correlation matrix
mcor = cor(clean[,c(5:8, 15, 23, 25)])
symnum(mcor)

measure.tib <- tibble(Measure = c("1. Interest", "2. Attitude", "3. Expectation", "4. Intention", "5. Own consumption", "6. Projected consumption", "7. Conformity"),
                      Mean = unlist(measure_sum[, 1:7]),
                      SD = unlist(measure_sum[, 8:14]))

cortable <- cbind(measure.tib, printnum(mcor)) %>% as_tibble()

# Randomization check ---------------------------
age_stat <- apa_print(aov(age ~ condition, clean)) # age
pol_stat <- apa_print(aov(politics ~ condition, clean)) # political position
gender_stat <- apa_print(chisq.test(clean$condition, clean$gender), n = nrow(clean)) # gender
nation_stat <- apa_print(chisq.test(clean$condition, clean$country), n = nrow(clean)) # nation

# H1 ---------------------------
## Does communicating a trending minority norm increase interest over and above communicating a minority norm only?
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

h1.mod <- '
interest       ~ static + none
attitude  ~ static + none
intention ~ static + none
expectation    ~ static + none'

h1.fit <- sem(model = h1.mod, data = clean)
h1.out <-  summary(h1.fit)

# H2 ---------------------------
##Will participants in the trending minority norm condition be more likely (than minority norm only) to expect a decrease in meat consumption by British people? 
cons_desc <- clean %>%
  group_by(condition) %>%
  summarise(n = n(),
            m_current = mean(cons_now_perc),
            m_next = mean(cons_next_perc),
            m_six = mean(cons_six_perc),
            m_composite = mean(cons_proj),
            sd_current = sd(cons_now_perc),
            sd_next = sd(cons_next_perc),
            sd_six = sd(cons_six_perc),
            sd_composite = sd(cons_proj), .groups = "rowwise")

h2.test <- apa_print(aov(cons_proj ~ condition, clean)) 

# H3 ---------------------------
##Does the perceived current and future popularity of sustainable eating behaviours correlate with interest, attitudes, expectations, and intentions to limit own meat consumption? 

h3.mod <- '
interest       ~ cons_proj
attitude  ~ cons_proj
intention ~ cons_proj
expectation    ~ cons_proj'

h3.fit <- sem(model = h3.mod, data = clean)
h3.out <-  summary(h3.fit)

# H4 ---------------------------
##Is projected meat consumption a mediator of the effect of trending minority norms vs. minority only on meat consumption outcomes?

h4.mod <- '
interest       ~ a*cons_proj + b*static
attitude  ~ c*cons_proj + d*static
intention ~ e*cons_proj + f*static
expectation    ~ g*cons_proj + h*static
cons_proj ~ k*static

ka := k*a
kc := k*c
ke := k*e
kg := k*g
totint := ka + b
totatt := kc + d
totintent := ke + f
totexp := kg + h'

h4.fit <- sem(h4.mod,data=clean, se="bootstrap", test="bootstrap", bootstrap = 5000, meanstructure=TRUE)
h4.out <- summary(h4.fit, standardized=TRUE)
h4.pam <- parameterEstimates(h4.fit)

# H5 ---------------------------
##How do demographic factors such as age, gender, and political position predict primary dependent variables relating to meat consumption? 

h5.mod <- '
interest       ~ static + cons_proj + age + genderbi + politics
attitude  ~ static + cons_proj + age + genderbi + politics
intention ~ static + cons_proj + age + genderbi + politics
expectation    ~ static + cons_proj + age + genderbi + politics'

h5.fit <- sem(model = h5.mod, data = clean)
h5.out <- summary(h5.fit)
