# author: Alaa AlDoh
# contact: az.aldoh@gmail.com/a.aldoh@sussex.ac.uk

# Load packages ---------------------------
list.of.packages <- c("svglite", "ggplot2", "papaja", "lavaan", "tidyverse", "knitr", "kableExtra", "codebook", "psych", "rlang", "bfrr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

set.seed(1234)
options(mc.cores = parallel::detectCores(), ## Run chains in parallel
        knitr.kable.NA = "") ## Don't show NAs in tables

export <- haven::read_spss("../data/01_EXPORTED.sav")

# Data cleaning ---------------------------
export[, c(2:5, 15, 17)] <- haven::as_factor(export[, c(2:5, 15, 17)])

raw <- export %>%
  mutate(att_checkbi = ifelse(att_check == 3, 1, 0) %>% as_factor(.),
         condition = as_factor(export$condition) %>% set_attrs(., labels = c("dynamic" = 1, "static" = 2, "none" = 3)) %>% relevel(., "dynamic"),
         cons_proj = select(., cons_now_perc_1, cons_next_perc_1, cons_six_perc_1) %>% rowMeans(., na.rm = TRUE),
         genderbi = na_if(gender, "Other") %>% droplevels() %>% recode_factor(., "Male" = -1, "Female" = 1),
         conformity_3 = 8 - conformity_3,
         conformity_mean = rowMeans(select(., starts_with("conformity_")), na.rm = TRUE)) %>%
  rename_with(., .fn = ~ str_remove(., "_1"), .cols = c(interest_1:politics_1, meat_cons_1))

# excluding vegetarians and attention check fails
#complete <- raw %>% filter(att_check == 1)
complete <- raw %>% filter(!is.na(conformity_1))
noveg <- complete %>% filter(veg != "Yes")

# excluding unneeded fields
clean <- noveg %>%
  select(-RecordedDate, -veg, -att_check) %>%
  mutate(condition = relevel(condition, "dynamic"),
         cons_proj_c = cons_proj - mean(cons_proj), #centering projected cons
         conformity_mean_c = conformity_mean - mean(conformity_mean)) %>% #centering conformity
  cbind(., psych::dummy.code(.$condition))

# outliers
mahalfiltered = mahalanobis(clean[,c(5:13, 17:22)], colMeans(clean[,c(5:13, 17:22)], na.rm = T), cov(clean[,c(5:13, 17:22)]))
cutoff = qchisq(1-.001, ncol(clean[,c(5:13, 17:22)]))
no_out = subset(clean, mahalfiltered < cutoff)

# correlation
mcor = cor(clean[,c(5:8, 15, 24, 26)], use = "pairwise.complete.obs")
symnum(mcor)

##assumption set up
random = rchisq(nrow(clean), 7)
fake = lm(random ~., data = clean[ , -1])
standardized = rstudent(fake)
fitted = scale(fake$fitted.values)

assumptions <- list(additivity = cor(mcor[ , -1]) %>% symnum(),
                    homogeneity = plot(fitted, standardized) + abline(0,0) + abline(v = 0), ##homog and s
                    normality = hist(standardized), ##multivariate normality
                    linearity = qqnorm(standardized)) ##multivariate linearity

# Data Overview ---------------------------
describeBy(clean, clean$condition) # check distribution and normality

# participants
data_desc <- c(total_n = nrow(raw),
               incomp_n = nrow(raw) - nrow(complete),
               veg_n = nrow(complete) - nrow(noveg),
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
            m_intent = mean(intention),
            m_expect = mean(expectation),
            m_cons = mean(meat_cons),
            m_cons_proj = mean(cons_proj),
            m_conformity = mean(conformity_mean, na.rm = T),
            sd_interest = sd(interest),
            sd_attitude = sd(attitude),
            sd_intent = sd(intention),
            sd_expect = sd(expectation),
            sd_cons = sd(meat_cons),
            sd_cons_proj = sd(cons_proj),
            sd_conformity = sd(conformity_mean, na.rm = T))

out_plots <- clean %>%
  pivot_longer(cols = interest:expectation, names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = condition, y = value)) +
  facet_wrap(~ variable, labeller = as_labeller(c('attitude'="Attitude", 'expectation'="Expectation", 'intention'="Intention", 'interest'="Interest"))) +
  geom_violin(trim = FALSE) + 
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), geom = "pointrange", color = "black") +
  ggtitle("Figure 1. Distribution of outcome variables by condition") +
  xlab("Norm Condition") +
  ylab("Value (%)") +
  scale_x_discrete(limits = c("dynamic", "static", "none"), labels=c("dynamic" = "Trending", "static" = "Minority","none" = "None")) +
  papaja::theme_apa()

scat_plots <- clean %>%
  pivot_longer(cols = interest:expectation, names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = cons_proj, y = value)) +
  facet_wrap(~ variable, labeller = as_labeller(c('attitude'="Attitude", 'expectation'="Expectation", 'intention'="Intention", 'interest'="Interest"))) +
  geom_smooth(aes(color = condition), method = loess, se = F) +
  stat_summary(fun.data = "mean_cl_boot", geom = "point", size = 0.5, alpha = 0.4) +
  ggtitle("Figure 1. Relationship between projected consumption and outcome variables") +
  xlab("Projected consumption (% limiting meat eating)") +
  ylab("Value of outcome (%)") +
  papaja::theme_apa()

ggsave(file="out_dist.svg", plot=out_plots)
ggsave(file="corr_scat.svg", plot=scat_plots)

# reliability
cron <- clean %>% select(conformity_1:conformity_6) %>% psych::alpha() # scale reliability

# correlation matrix
measure.tib <- tibble(Measure = c("1. Interest", "2. Attitude", "3. Intention", "4. Expectation", "5. Own consumption", "6. Projected consumption", "7. Conformity"),
                      Mean = unlist(measure_sum[, 1:7]),
                      SD = unlist(measure_sum[, 8:14])) %>%
  cbind(., mcor) %>% select(-conformity_mean) %>% remove_rownames() 

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
            m_intent = mean(intention),
            m_expect = mean(expectation),
            sd_interest = sd(interest),
            sd_attitude = sd(attitude),
            sd_intent = sd(intention),
            sd_expect = sd(expectation), .groups = "rowwise") %>% printnum()

outcomes_tab <- paste(unlist(outcomes_desc[3:6]), unlist(outcomes_desc[7:10]), sep = ' $\\pm$ ')

outcomes.tib <- tibble(Measure = c("1. Interest", "2. Attitude", "3. Intention", "4. Expectation"),
                      Trending = unlist(outcomes_tab[c(3,6,9,12)]),
                      Minority = unlist(outcomes_tab[c(1,4,7,10)]),
                      None     = unlist(outcomes_tab[c(2,5,8,11)]))

h1.mod <- '
interest       ~ static + none
attitude  ~ static + none
intention ~ static + none
expectation    ~ static + none'

h1.fit <- sem(model = h1.mod, data = clean)
h1.out <-  summary(h1.fit, standardized = T, ci = T, fit.measures = T, rsq = T)

h1.effect <- sapply(1:8, function(x) bfrr(-1*h1.out$PE[x, 5],h1.out$PE[x,6], sample_df = h1.out$FIT["ntotal"] - 1, model = "normal", mean = 0, sd = 5, tail = 1, criterion = 3,
                                          rr_interval = list(mean = c(-15, 15), sd = c(0, 15)), precision = 0.05))[-14,] # effect sizes

h1.rr <- sapply(1:8, function(x) paste0("HN[", toString(h1.effect[,x]$RR$sd), "]"))

h1.table <- cbind(h1.out$PE[1:8,c(1, 3, 5, 12, 6:10)], unlist(h1.effect[3,]), h1.rr, unlist(h1.effect[5,])) %>% .[with(., order(rhs, decreasing = TRUE)), ] %>% select(., -2) %>%
  mutate(pvalue = printp(pvalue))

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
interest  ~ cons_proj
attitude  ~ cons_proj
intention ~ cons_proj
expectation ~ cons_proj'

h3.fit <- sem(model = h3.mod, data = clean)
h3.out <-  summary(h3.fit, standardized = T, ci = T, fit.measures = T, rsq = T)

h3.effect <- sapply(1:4, function(x) bfrr(h3.out$PE[x, 5],h3.out$PE[x,6], sample_df = h3.out$FIT["ntotal"] - 1, model = "normal", mean = 0, sd = 5, tail = 1, criterion = 3,
                                          rr_interval = list(mean = c(-15, 15), sd = c(0, 15)), precision = 0.05))[-14,] # effect sizes

h3.rr <- sapply(1:4, function(x) paste0("HN[", toString(h3.effect[,x]$RR$sd), "]"))

h3.table <- cbind(h3.out$PE[1:4, c(1, 3, 5, 11, 6:10)], unlist(h3.effect[3,]), h3.rr, unlist(h3.effect[5,])) %>% select(., -2) %>%
  mutate(pvalue = printp(pvalue))

# H4 ---------------------------
##Is projected meat consumption a mediator of the effect of trending minority norms vs. minority only on meat consumption outcomes?

h4.mod <- '
interest  ~ a*cons_proj + b*static
attitude  ~ c*cons_proj + d*static
intention ~ e*cons_proj + f*static
expectation ~ g*cons_proj + h*static
cons_proj ~ t*static

ta := t*a
tc := t*c
te := t*e
tg := t*g
totint := ta + b
totatt := tc + d
totintent := te + f
totexp := tg + h'

h4.fit <- sem(h4.mod,data=clean, se="bootstrap", test="bootstrap", bootstrap = 5000, meanstructure=TRUE)
h4.out <- summary(h4.fit, standardized = T, ci = T, fit.measures = T)
h4.pam <- parameterEstimates(h4.fit)

h4.table <- h4.out$PE[c(1:9, 28:35),c(1, 3, 6, 13, 7:11)] %>% .[with(., order(rhs, decreasing = FALSE)), ]
h4.table <- rbind("", h4.table[1:4,], "", h4.table[5:17,] ) %>%
  mutate_at(vars(3:9), ~as.numeric(as.character(.))) %>%
  mutate(pvalue = printp(pvalue),
         lhs = c("~ Projected consumption", "Interest", "Attitude", "Intention", "Expectation", "~ Condition", "Interest", "Attitude", "Intention", "Expectation", "Projected consumption", 
                 "Interest", "Attitude", "Intention", "Expectation", "Interest", "Attitude", "Intention", "Expectation")) %>% 
  select(., -2) # effect sizes

# H5 ---------------------------
##How do demographic factors such as age, gender, and political position predict primary dependent variables relating to meat consumption? 

h5.mod <- '
interest  ~ static + cons_proj + age + genderbi + politics
attitude  ~ static + cons_proj + age + genderbi + politics
intention ~ static + cons_proj + age + genderbi + politics
expectation ~ static + cons_proj + age + genderbi + politics'

h5.fit <- sem(model = h5.mod, data = clean)
h5.out <- summary(h5.fit, standardized = T, ci = T, fit.measures = T, rsq = T)

h5.table <- h5.out$PE[1:20,c(3, 5, 12, 6:10)] %>%
  mutate(pvalue = printp(pvalue),
         rhs = rep(c("Condition$^a$", "Projected consumption", "Age", "Gender$^b$", "Politics$^c$"), 4))

# Exploratory ---------------------------
##moderated mediation - CONDITIONAL PROCESS ANALYSIS
clean$cons_conformity <- clean$cons_proj_c*clean$conformity_mean_c # create the interaction (product) term

exp.mod <- '
cons_proj_c ~ a*static
interest ~ b1*cons_proj_c + b2*conformity_mean_c + b3*cons_conformity + c*static

# indirect effect when conformity = 0
ab0 := a*b1 
total0 := ab0 + b1 + c

# indirect effect when negexp = -1.0sd (sd=0.8422133)
ablow := a + b1 + b3*-0.8422133
totallow := ablow + c
abhigh := a + b1 + b3*0.8422133
totalhigh := abhigh + c'

exp.out <- sem(exp.mod, data=clean, meanstructure=TRUE)
summary(exp.out, fit.measures=TRUE,  rsq=TRUE)


