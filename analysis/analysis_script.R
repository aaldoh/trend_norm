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
mcor = cor(clean[,c(5:8, 15, 24, 26)])
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
interest       ~ conditionbi
attitude  ~ conditionbi
intention ~ conditionbi
expectation    ~ conditionbi'

h1.fit <- sem(model = simple.mod, data = clean, seed = 2019)
h1.out <- summary(h1.fit)
