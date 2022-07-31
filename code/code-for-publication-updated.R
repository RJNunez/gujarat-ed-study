# -- Set up 
source("code/init.R")

# -- Loading municipality data
load("data/rdas/municipality-data.rda")

##########################################################################################################
################## -- MARGINAL ANALYSIS ------------------------------------------------- ################
##########################################################################################################
# -- Dates
dates <- tibble(date = seq(make_date(2019, 01, 01), make_date(2021, 12, 31), by = "days")) %>%
  mutate(week = ifelse(week(date) == 53, 52, week(date)),
         year = year(date)) %>%
  group_by(year, week) %>%
  filter(date == first(date)) %>%
  ungroup()

# -- Wrangling mortality data
counts <- dat %>%
  mutate(year = year(date),
         week = ifelse(week(date) == 53, 52, week(date))) %>%
  dplyr::select(-date) %>%
  left_join(dates, by = c("week", "year")) %>%
  group_by(date, village) %>%
  summarize(outcome = n()) %>%
  ungroup() %>%
  mutate(year = year(date)) %>%
  arrange(date) %>%
  filter(date >= "2019-01-01",
         date <= "2021-04-23")

# -- All municipalities
all_municipalities <- sort(unique(counts$village))

# -- Municipalities to exclude
p           <- 0.40
N_train     <- sum(unique(counts$date) <= "2020-02-28")
N_test      <- sum(unique(counts$date) > "2020-02-28")
village_out <- counts %>%
  mutate(flag = ifelse(date <= "2020-02-28", "train", "test")) %>%
  group_by(village, flag) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = flag, values_from = n) %>%
  mutate(p_train = train / N_train,
         p_test  = test / N_test) %>%
  filter(p_train < p | p_test  < p) %>%
  pull(village) %>%
  unique(.)

# -- Viz municipalities to take out
counts %>%
  filter(village %in% village_out) %>%
  ggplot(aes(date, outcome)) +
  facet_wrap(~village) +
  geom_vline(xintercept = make_date(2020,03,01),
             color      = "gray",
             linetype   = 2) +
  geom_point(alpha = 0.50) +
  scale_x_date(date_labels = "%b") +
  labs(x = "Date",
       y = "Number of deaths",
       color = "")

# -- Getting rid of bad quality data
counts <- filter(counts, !village %in% village_out)

# -- Marginal municipalities
marginal_municipalities <- sort(unique(counts$village))

# -- Training index: Jan 2019 to Feb 2020
train_idx <- which(counts$date %in% seq(make_date(2019, 01, 01), make_date(2020, 02, 28), by = "days"))

# -- Design matrix for municipality-specific intercepts
x_i <- model.matrix(outcome ~ village, data = counts)[,-1]

# -- Design matrix for seasonal component
yd  <- noleap_yday(counts$date)
x_s <- cbind(1, fourier_trend(yd, k = 2))
colnames(x_s)[1] <- "Intercept"

# -- Design matrix for model fit
x <- cbind(x_s, x_i)

# -- Outcome
y <- counts$outcome

# -- Putting everything together
X <- mutate(as_tibble(x), y = y)

# -- Fitting model and getting fitted values
fit <- glm(y ~ . - 1, family = "quasipoisson", data = X[train_idx,])

# -- Extracting municipality-specific intercepts from the model fit
intercepts <- tibble(intercept = exp(c(coef(fit)[1], coef(fit)[1] + coef(fit)[-c(1:5)])),
                     village   = sort(unique(counts$village)))

# -- Design matrix for municipality-specific intercepts
v            <- factor(counts$village)
contrasts(v) <- contr.sum(length(levels(v)), contrasts = TRUE)
x_i          <- model.matrix(~v)[, -1]

# -- Design matrix for seasonal component
yd  <- noleap_yday(counts$date)
x_s <- cbind(1, fourier_trend(yd, k = 2))
colnames(x_s)[1] <- "Intercept"

# -- Design matrix for model fit
x <- cbind(x_s, x_i)

# -- Outcome
y <- counts$outcome

# -- Putting everything together
X <- mutate(as_tibble(x), y = y)

# -- Fitting model and getting fitted values
fit <- glm(y ~ . - 1, family = "quasipoisson", data = X[train_idx,])

# -- Sorting municipalities by the estimated intercept
levs <- intercepts %>%
  arrange(desc(intercept)) %>%
  pull(village)

##
intercepts %>%
  ggplot(aes(intercept)) +
  geom_density(fill  = "black",
               alpha = 0.50) +
  geom_point(aes(intercept, 0), 
             shape = "|",
             size  = 3) +
  labs(x = expression(gamma[i]),
       y = "Density")

# -- Getting fitted values and ses
preds_po <- predict(fit, newdata = X, se.fit = TRUE, type = "link")

# -- Putting everything together
fit_res <- counts %>%
  mutate(expected        = exp(preds_po$fit),
         log_expected_se = preds_po$se.fit,
         expected_lwr    = exp(log(expected) - 1.96 * log_expected_se),
         expected_upr    = exp(log(expected) + 1.96 * log_expected_se),
         type            = "glm-po") %>%
  left_join(intercepts, by = "village") %>%
  mutate(village = factor(village, levels = levs))

# -- Subsetting data only for Covid times
test_counts    <- filter(fit_res, date >= "2020-03-01")
municipalities <- sort(unique(test_counts$village))

# -- Fitting second model
res_municipalities <- map_df(municipalities, function(x) {
  
  # -- Subsetting data
  tmp_dat <- filter(test_counts, village == x)
  
  if(nrow(tmp_dat) <= 26) { 
    # -- Number of knots to use
    nknots <- nrow(tmp_dat) 
    
    # -- Fitting second model
    tmp_fit <- glm(outcome ~ ns(date, df = nknots - 2, intercept = FALSE) + offset(log(expected)), 
                   family = "quasipoisson",
                   data   = tmp_dat)
  } else {
    # -- Number of knots to use
    nknots <- 24
    
    # -- Fitting second model
    tmp_fit <- glm(outcome ~ ns(date, df = nknots + 1, intercept = FALSE) + offset(log(expected)), 
                   family = "quasipoisson",
                   data   = tmp_dat)
  }
  
  # -- Getting fitted values and se
  preds <- predict(tmp_fit, se.fit = TRUE, type = "response")
  
  # -- Percent change and SE 
  mu          <- tmp_dat$expected
  lambda      <- preds$fit
  dispersion  <- summary(tmp_fit)$dispersion
  cova        <- summary(tmp_fit)$cov.unscaled * dispersion
  lambda_vari <- lambda^2 * diag(model.matrix(tmp_fit) %*% cova %*% t(model.matrix(tmp_fit)))
  mu_vari     <- mu^2 * tmp_dat$log_expected_se^2
  fitted      <- lambda / mu - 1
  se          <- sqrt((lambda_vari / mu^2) + (lambda^2 * mu_vari / mu^4))
  
  # -- Excess deaths
  ed     <- lambda - mu
  ed_se  <- sqrt(lambda * dispersion + mu_vari) 
  ced    <- cumsum(ed)
  ced_se <- sqrt(cumsum(ed_se^2))
  
  # -- Putting everything together
  tmp_dat <- tmp_dat %>%
    mutate(fitted = fitted,
           se     = se,
           ed     = ed, 
           ed_se  = ed_se,
           ced    = ced,
           ced_se = ced_se)
})

# -- Marginal percent change from average
pc_marginal <- res_municipalities %>%
  mutate(expected_se = expected * log_expected_se) %>%
  group_by(date) %>%
  mutate(sum_age_expected       = sum(expected),
         sum_age_gamma_expected = sum(fitted * expected)) %>%
  ungroup() %>%
  mutate(gamma_derivative = expected / sum_age_expected,
         mu_derivative    = (fitted / sum_age_expected) + (sum_age_gamma_expected / (sum_age_expected^2))) %>%
  mutate(tmp = se^2 * gamma_derivative^2 + expected_se^2 * mu_derivative^2) %>%
  group_by(date) %>%
  summarize(fitted = sum(fitted * expected) / sum(expected),
            observed    = sum(outcome),
            expected    = sum(expected),
            expected_se = sqrt(sum(expected_se^2)),
            se          = sqrt(sum(tmp))) %>%
  ungroup()

# -- Cumulative excess deaths from the overall analysis
ced_overall <- res_municipalities %>%
  group_by(date) %>%
  summarize(ed        = sum(ed), 
            ed_se     = sqrt(sum(ed_se^2)),
            intercept = sum(intercept),
            expected  = sum(expected)) %>%
  ungroup() %>%
  mutate(ced         = cumsum(ed),
         ced_se      = sqrt(cumsum(ed_se^2)), 
         lwr         = ced - 1.96*ced_se,
         upr         = ced + 1.96*ced_se,
         ced_adj     = ced / intercept,
         ced_adj_lwr = ced / intercept - 1.96*ced_se / intercept,
         ced_adj_upr = ced / intercept + 1.96*ced_se / intercept)

##
ced_overall %>%
  ggplot(aes(date, ced)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              alpha = 0.30) +
  geom_line(size  = 1, 
            color = "white") +
  geom_line() +
  geom_label(aes(x = make_date(2020,09,01), y = 20000, label = "21,300 [95% CI: 20,700 to 22,000]"),
             fontface = "bold") +
  scale_x_date(date_labels = "%b") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Date",
       y = "Cumulative excess deaths")
##########################################################################################################
################## -- END MARGINAL ANALYSIS --------------------------------------------- ################
##########################################################################################################

##########################################################################################################
################## -- GENDER ANALYSIS --------------------------------------------------- ################
##########################################################################################################
# -- Dates
dates <- tibble(date = seq(make_date(2019, 01, 01), make_date(2021, 12, 31), by = "days")) %>%
  mutate(week = ifelse(week(date) == 53, 52, week(date)),
         year = year(date)) %>%
  group_by(year, week) %>%
  filter(date == first(date)) %>%
  ungroup()

# -- Wrangling mortality data
counts <- dat %>%
  mutate(year = year(date),
         week = ifelse(week(date) == 53, 52, week(date))) %>%
  dplyr::select(-date) %>%
  left_join(dates, by = c("week", "year")) %>%
  group_by(date, village, gender) %>%
  summarize(outcome = n()) %>%
  ungroup() %>%
  mutate(year = year(date)) %>%
  arrange(date) %>%
  filter(date >= "2019-01-01",
         date <= "2021-04-23")

# -- Municipalities to exclude
p           <- 0.40
N_train     <- sum(unique(counts$date) <= "2020-02-28")
N_test      <- sum(unique(counts$date) > "2020-02-28")
village_out <- counts %>%
  mutate(flag = ifelse(date <= "2020-02-28", "train", "test")) %>%
  filter(gender != "other") %>%
  group_by(village, gender, flag) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  group_by(village, flag) %>%
  summarize(n = min(n)) %>%
  ungroup() %>%
  pivot_wider(names_from = flag, values_from = n) %>%
  mutate(p_train = train / N_train,
         p_test  = test / N_test) %>%
  filter(p_train < p | p_test  < p) %>%
  pull(village) %>%
  unique(.)
# Specifically, we excluded municiaplities with at least 60% missing data in either the train or test set
  
# -- Adding municipalities with no data
village_out <- c(village_out, "Boriavi", "Halvad", "Idar", "Okha", "Pardi", "Zalod")

# -- Viz municipalities to take out
counts %>%
  filter(village %in% village_out) %>%
  ggplot(aes(date, outcome, color = gender)) +
  geom_vline(xintercept = make_date(2020,03,01),
             color      = "gray",
             linetype   = 2) +
  facet_wrap(~village) +
  geom_point(alpha = 0.50) +
  scale_color_manual(values = as.vector(my_palette[c(2,3,9)])) +
  scale_x_date(date_labels = "%b") +
  labs(x = "Date",
       y = "Number of deaths",
       color = "") +
  theme(legend.title     = element_blank(),
        legend.direction = "vertical",
        legend.position  = c(0.88, 0.08))
# width = 10
# height = 6

# -- Subsetting data
counts <- filter(counts, gender != "other", !village %in% village_out)

# -- Fitting mean model for each municipality
res_sex <- map_df(c("male", "female"), function(i) {
  
  # -- Subsetting data
  tmp_counts <- filter(counts, gender == i)
  
  # -- Training index: Jan 2019 to Feb 2020
  train_idx <- which(tmp_counts$date %in% seq(make_date(2019, 01, 01), make_date(2020, 02, 28), by = "days"))
  
  # -- Design matrix for municipality-specific intercepts
  x_i   <- model.matrix(outcome ~ village, data = tmp_counts)[,-1]
  
  # -- Design matrix for seasonal component
  yd  <- noleap_yday(tmp_counts$date)
  x_s <- cbind(1, fourier_trend(yd, k = 2))
  colnames(x_s)[1] <- "Intercept"
  
  # -- Design matrix for model fit
  x <- cbind(x_s, x_i)
  
  # -- Outcome
  y <- tmp_counts$outcome
  
  # -- Putting everything together
  X <- mutate(as_tibble(x), y = y)
  
  # -- Fitting model and getting fitted values
  fit <- glm(y ~ . - 1, family = "quasipoisson", data = X[train_idx,])
  
  # -- Extracting municipality-specific intercepts from the model fit
  intercepts <- tibble(intercept = exp(c(coef(fit)[1], coef(fit)[1] + coef(fit)[-c(1:5)])),
                       village   = sort(unique(counts$village)))
  
  # -- Getting fitted values and ses
  preds_po <- predict(fit, newdata = X, se.fit = TRUE, type = "link")
  
  # -- Putting everything together
  fit_res <- tmp_counts %>%
    mutate(expected        = exp(preds_po$fit),
           log_expected_se = preds_po$se.fit,
           expected_lwr    = exp(log(expected) - 1.96 * log_expected_se),
           expected_upr    = exp(log(expected) + 1.96 * log_expected_se)) %>%
    left_join(intercepts, by = "village")
})

# -- Ordering for villages by intercept value
levs <- res_sex %>%
  select(village, gender, intercept) %>%
  unique() %>%
  group_by(village) %>%
  summarize(avg = sum(intercept),
            n   = n()) %>%
  ungroup() %>%
  arrange(desc(avg)) %>%
  pull(village)

# -- Adding ordering the data
res_sex <- mutate(res_sex, village = factor(village, levels = levs))

# -- Subsetting data only for Covid times
municipalities    <- sort(unique(res_sex$village))
test_counts_males <- filter(res_sex, date >= "2020-03-01", gender == "male")
test_counts_females <- filter(res_sex, date >= "2020-03-01", gender == "female")

# -- Fitting second model to male data
res_males <- map_df(municipalities, function(x) {
  
  # -- Subsetting data
  tmp_dat <- filter(test_counts_males, village == x)
  
  # -- Number of knots to use
  nknots <- 24
  
  # -- Fitting second model
  tmp_fit <- glm(outcome ~ ns(date, df = nknots, intercept = FALSE) + offset(log(expected)), 
                 family = "quasipoisson",
                 data   = tmp_dat)
  
  # -- Getting fitted values and se
  preds <- predict(tmp_fit, se.fit = TRUE, type = "response")
  
  # -- Percent change and SE 
  mu          <- tmp_dat$expected
  lambda      <- preds$fit
  dispersion  <- summary(tmp_fit)$dispersion
  cova        <- summary(tmp_fit)$cov.unscaled * dispersion
  lambda_vari <- lambda^2 * diag(model.matrix(tmp_fit) %*% cova %*% t(model.matrix(tmp_fit)))
  mu_vari     <- mu^2 * tmp_dat$log_expected_se^2
  fitted      <- lambda / mu - 1
  se          <- sqrt((lambda_vari / mu^2) + (lambda^2 * mu_vari / mu^4))
  
  # -- Excess deaths
  ed     <- lambda - mu
  ed_se  <- sqrt(lambda * dispersion + mu_vari) 
  ced    <- cumsum(ed)
  ced_se <- sqrt(cumsum(ed_se^2))
  
  # -- Putting everything together
  tmp_dat <- tmp_dat %>%
    mutate(fitted = fitted,
           se     = se,
           ed     = ed, 
           ed_se  = ed_se,
           ced    = ced,
           ced_se = ced_se,
           gender = "male")
})

# -- Fitting second model to female data
res_females <- map_df(municipalities, function(x) {
  
  # -- Subsetting data
  tmp_dat <- filter(test_counts_females, village == x)
  
  # -- Number of knots to use
  nknots <- 24
  
  # -- Fitting second model
  tmp_fit <- glm(outcome ~ ns(date, df = nknots, intercept = FALSE) + offset(log(expected)), 
                 family = "quasipoisson",
                 data   = tmp_dat)
  
  # -- Getting fitted values and se
  preds <- predict(tmp_fit, se.fit = TRUE, type = "response")
  
  # -- Percent change and SE 
  mu          <- tmp_dat$expected
  lambda      <- preds$fit
  dispersion  <- summary(tmp_fit)$dispersion
  cova        <- summary(tmp_fit)$cov.unscaled * dispersion
  lambda_vari <- lambda^2 * diag(model.matrix(tmp_fit) %*% cova %*% t(model.matrix(tmp_fit)))
  mu_vari     <- mu^2 * tmp_dat$log_expected_se^2
  fitted      <- lambda / mu - 1
  se          <- sqrt((lambda_vari / mu^2) + (lambda^2 * mu_vari / mu^4))
  
  # -- Excess deaths
  ed     <- lambda - mu
  ed_se  <- sqrt(lambda * dispersion + mu_vari) 
  ced    <- cumsum(ed)
  ced_se <- sqrt(cumsum(ed_se^2))
  
  # -- Putting everything together
  tmp_dat <- tmp_dat %>%
    mutate(fitted = fitted,
           se     = se,
           ed     = ed, 
           ed_se  = ed_se,
           ced    = ced,
           ced_se = ced_se,
           gender = "female")
})

# -- Putting data together
res_sex <- bind_rows(res_males, res_females)

# -- Marginal percent Change from average
pc_marginal_sex <- res_sex %>%
  mutate(expected_se = expected * log_expected_se) %>%
  group_by(date, gender) %>%
  mutate(sum_age_expected       = sum(expected),
         sum_age_gamma_expected = sum(fitted * expected)) %>%
  ungroup() %>%
  mutate(gamma_derivative = expected / sum_age_expected,
         mu_derivative    = (fitted / sum_age_expected) + (sum_age_gamma_expected / (sum_age_expected^2))) %>%
  mutate(tmp = se^2 * gamma_derivative^2 + expected_se^2 * mu_derivative^2) %>%
  group_by(date, gender) %>%
  summarize(fitted = sum(fitted * expected) / sum(expected),
            observed    = sum(outcome),
            expected    = sum(expected),
            expected_se = sqrt(sum(expected_se^2)),
            se          = sqrt(sum(tmp))) %>%
  ungroup()

# -- Cumulative excess deaths from the gender analysis
ced_gender <- res_sex %>%
  group_by(date, gender) %>%
  summarize(ed        = sum(ed), 
            ed_se     = sqrt(sum(ed_se^2)),
            intercept = sum(intercept),
            expected  = sum(expected)) %>%
  ungroup() %>%
  group_by(gender) %>%
  mutate(ced         = cumsum(ed),
         ced_se      = sqrt(cumsum(ed_se^2)), 
         lwr         = ced - 1.96*ced_se,
         upr         = ced + 1.96*ced_se,
         ced_adj     = ced / intercept,
         ced_adj_lwr = ced / intercept - 1.96*ced_se / intercept,
         ced_adj_upr = ced / intercept + 1.96*ced_se / intercept) %>%
  ungroup()
##########################################################################################################
################## -- END GENDER ANALYSIS ----------------------------------------------- ################
##########################################################################################################

##########################################################################################################
################## -- AGE ANALYSIS ------------------------------------------------------ ################
##########################################################################################################
# -- Dates
dates <- tibble(date = seq(make_date(2019, 01, 01), make_date(2021, 12, 31), by = "days")) %>%
  mutate(week = ifelse(week(date) == 53, 52, week(date)),
         year = year(date)) %>%
  group_by(year, week) %>%
  filter(date == first(date)) %>%
  ungroup()

# -- Wrangling mortality data
counts <- dat %>%
  mutate(year     = year(date),
         week     = ifelse(week(date) == 53, 52, week(date)),
         agegroup = cut(age, breaks = c(0, 20, 40, 65, Inf))) %>%
  dplyr::select(-date) %>%
  left_join(dates, by = c("week", "year")) %>%
  group_by(date, village, agegroup) %>%
  summarize(outcome = n()) %>%
  ungroup() %>%
  mutate(year = year(date)) %>%
  arrange(date) %>%
  filter(date >= "2019-01-01",
         date <= "2021-04-23")

# -- Municipalities to exclude
p           <- 0.40
N_train     <- sum(unique(counts$date) <= "2020-02-28")
N_test      <- sum(unique(counts$date) > "2020-02-28")

# -- For (0,20] group
village_out_20 <- counts %>%
  mutate(flag = ifelse(date <= "2020-02-28", "train", "test")) %>%
  filter(agegroup == "(0,20]") %>%
  group_by(village, flag) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = flag, values_from = n) %>%
  mutate(p_train = train / N_train,
         p_test  = test / N_test) %>%
  filter(p_train < p | p_test  < p) %>%
  pull(village) %>%
  unique(.)

# -- For (20,40] group
village_out_40 <- counts %>%
  mutate(flag = ifelse(date <= "2020-02-28", "train", "test")) %>%
  filter(agegroup == "(20,40]") %>%
  group_by(village, flag) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = flag, values_from = n) %>%
  mutate(p_train = train / N_train,
         p_test  = test / N_test) %>%
  filter(p_train < p | p_test  < p) %>%
  pull(village) %>%
  unique(.)

# -- For (40,65] group
village_out_65 <- counts %>%
  mutate(flag = ifelse(date <= "2020-02-28", "train", "test")) %>%
  filter(agegroup == "(40,65]") %>%
  group_by(village, flag) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = flag, values_from = n) %>%
  mutate(p_train = train / N_train,
         p_test  = test / N_test) %>%
  filter(p_train < p | p_test  < p) %>%
  pull(village) %>%
  unique(.)

# -- For (65,Inf] group
village_out_Inf <- counts %>%
  mutate(flag = ifelse(date <= "2020-02-28", "train", "test")) %>%
  filter(agegroup == "(65,Inf]") %>%
  group_by(village, flag) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = flag, values_from = n) %>%
  mutate(p_train = train / N_train,
         p_test  = test / N_test) %>%
  filter(p_train < p | p_test  < p) %>%
  pull(village) %>%
  unique(.)

# -- Fitting mean model for each municipality
res_age <- map_df(c("(0,20]", "(20,40]", "(40,65]", "(65,Inf]"), function(i) {
  
  # -- Subsetting data
  if(i == "(0,20]") {
    tmp_counts <- filter(counts, agegroup == i, !village %in% village_out_20)
  } else if(i == "(20,40]") {
    tmp_counts <- filter(counts, agegroup == i, !village %in% village_out_40)
  } else if(i == "(40,65]") {
    tmp_counts <- filter(counts, agegroup == i, !village %in% village_out_65)
  } else if(i == "(65,Inf]") {
    tmp_counts <- filter(counts, agegroup == i, !village %in% village_out_Inf)
  }
  
  # -- Training index: Jan 2019 to Feb 2020
  train_idx <- which(tmp_counts$date %in% seq(make_date(2019, 01, 01), make_date(2020, 02, 28), by = "days"))
  
  # -- Design matrix for municipality-specific intercepts
  x_i <- model.matrix(outcome ~ village, data = tmp_counts)[,-1]
  
  # -- Design matrix for seasonal component
  yd  <- noleap_yday(tmp_counts$date)
  x_s <- cbind(1, fourier_trend(yd, k = 2))
  colnames(x_s)[1] <- "Intercept"
  
  # -- Design matrix for model fit
  x <- cbind(x_s, x_i)
  
  # -- Outcome
  y <- tmp_counts$outcome
  
  # -- Putting everything together
  X <- mutate(as_tibble(x), y = y)
  
  # -- Fitting model and getting fitted values
  fit <- glm(y ~ . - 1, family = "quasipoisson", data = X[train_idx,])
  
  # -- Extracting municipality-specific intercepts from the model fit
  intercepts <- tibble(intercept = exp(c(coef(fit)[1], coef(fit)[1] + coef(fit)[-c(1:5)])),
                       village   = sort(unique(tmp_counts$village)))
  
  # -- Getting fitted values and ses
  preds_po <- predict(fit, newdata = X, se.fit = TRUE, type = "link")
  
  # -- Putting everything together
  fit_res <- tmp_counts %>%
    mutate(expected        = exp(preds_po$fit),
           log_expected_se = preds_po$se.fit,
           expected_lwr    = exp(log(expected) - 1.96 * log_expected_se),
           expected_upr    = exp(log(expected) + 1.96 * log_expected_se)) %>%
    left_join(intercepts, by = "village")
})

# -- Ordering for villages by intercept value
levs <- res_age %>%
  select(village, agegroup, intercept) %>%
  unique() %>%
  group_by(village) %>%
  summarize(avg = sum(intercept),
            n   = n()) %>%
  ungroup() %>%
  arrange(desc(avg)) %>%
  pull(village)

# -- Adding ordering the data
res_age <- mutate(res_age, village = factor(village, levels = levs))

# -- Subsetting data only for Covid times
municipalities  <- sort(unique(res_age$village))
test_counts_20  <- filter(res_age, date >= "2020-03-01", agegroup == "(0,20]")
test_counts_40  <- filter(res_age, date >= "2020-03-01", agegroup == "(20,40]")
test_counts_65  <- filter(res_age, date >= "2020-03-01", agegroup == "(40,65]")
test_counts_Inf <- filter(res_age, date >= "2020-03-01", agegroup == "(65,Inf]")

# -- Fitting second model to agegroup data
res_20 <- map_df(setdiff(municipalities, village_out_20), function(x) {
  
  # -- Subsetting data
  tmp_dat <- filter(test_counts_20, village == x)
  
  # -- Number of knots to use
  nknots <- 22
  
  # -- Fitting second model
  tmp_fit <- glm(outcome ~ ns(date, df = nknots, intercept = FALSE) + offset(log(expected)), 
                 family = "quasipoisson",
                 data   = tmp_dat)

  # -- Getting fitted values and se
  preds <- predict(tmp_fit, se.fit = TRUE, type = "response")
  
  # -- Percent change and SE 
  mu          <- tmp_dat$expected
  lambda      <- preds$fit
  dispersion  <- summary(tmp_fit)$dispersion
  cova        <- summary(tmp_fit)$cov.unscaled * dispersion
  lambda_vari <- lambda^2 * diag(model.matrix(tmp_fit) %*% cova %*% t(model.matrix(tmp_fit)))
  mu_vari     <- mu^2 * tmp_dat$log_expected_se^2
  fitted      <- lambda / mu - 1
  se          <- sqrt((lambda_vari / mu^2) + (lambda^2 * mu_vari / mu^4))
  
  # -- Excess deaths
  ed     <- lambda - mu
  ed_se  <- sqrt(lambda * dispersion + mu_vari) 
  ced    <- cumsum(ed)
  ced_se <- sqrt(cumsum(ed_se^2))
  
  # -- Putting everything together
  tmp_dat <- tmp_dat %>%
    mutate(fitted = fitted,
           se     = se,
           ed     = ed, 
           ed_se  = ed_se,
           ced    = ced,
           ced_se = ced_se,
           agegroup = "(0,20]")
})

# -- Fitting second model to agegroup data
res_40 <- map_df(setdiff(municipalities, village_out_40), function(x) {
  
  # -- Subsetting data
  tmp_dat <- filter(test_counts_40, village == x)
  
  # -- Number of knots to use
  nknots <- 22
  
  # -- Fitting second model
  tmp_fit <- glm(outcome ~ ns(date, df = nknots, intercept = FALSE) + offset(log(expected)), 
                 family = "quasipoisson",
                 data   = tmp_dat)
  
  # -- Getting fitted values and se
  preds <- predict(tmp_fit, se.fit = TRUE, type = "response")
  
  # -- Percent change and SE 
  mu          <- tmp_dat$expected
  lambda      <- preds$fit
  dispersion  <- summary(tmp_fit)$dispersion
  cova        <- summary(tmp_fit)$cov.unscaled * dispersion
  lambda_vari <- lambda^2 * diag(model.matrix(tmp_fit) %*% cova %*% t(model.matrix(tmp_fit)))
  mu_vari     <- mu^2 * tmp_dat$log_expected_se^2
  fitted      <- lambda / mu - 1
  se          <- sqrt((lambda_vari / mu^2) + (lambda^2 * mu_vari / mu^4))
  
  # -- Excess deaths
  ed     <- lambda - mu
  ed_se  <- sqrt(lambda * dispersion + mu_vari) 
  ced    <- cumsum(ed)
  ced_se <- sqrt(cumsum(ed_se^2))
  
  # -- Putting everything together
  tmp_dat <- tmp_dat %>%
    mutate(fitted = fitted,
           se     = se,
           ed     = ed, 
           ed_se  = ed_se,
           ced    = ced,
           ced_se = ced_se,
           agegroup = "(20,40]")
})

# -- Fitting second model to agegroup data
res_65 <- map_df(setdiff(municipalities, village_out_65), function(x) {

  # -- Subsetting data
  tmp_dat <- filter(test_counts_65, village == x)
  
  # -- Number of knots to use
  nknots <- 22
  
  # -- Fitting second model
  tmp_fit <- glm(outcome ~ ns(date, df = nknots, intercept = FALSE) + offset(log(expected)), 
                 family = "quasipoisson",
                 data   = tmp_dat)
  
  # -- Getting fitted values and se
  preds <- predict(tmp_fit, se.fit = TRUE, type = "response")
  
  # -- Percent change and SE 
  mu          <- tmp_dat$expected
  lambda      <- preds$fit
  dispersion  <- summary(tmp_fit)$dispersion
  cova        <- summary(tmp_fit)$cov.unscaled * dispersion
  lambda_vari <- lambda^2 * diag(model.matrix(tmp_fit) %*% cova %*% t(model.matrix(tmp_fit)))
  mu_vari     <- mu^2 * tmp_dat$log_expected_se^2
  fitted      <- lambda / mu - 1
  se          <- sqrt((lambda_vari / mu^2) + (lambda^2 * mu_vari / mu^4))
  
  # -- Excess deaths
  ed     <- lambda - mu
  ed_se  <- sqrt(lambda * dispersion + mu_vari) 
  ced    <- cumsum(ed)
  ced_se <- sqrt(cumsum(ed_se^2))
  
  # -- Putting everything together
  tmp_dat <- tmp_dat %>%
    mutate(fitted = fitted,
           se     = se,
           ed     = ed, 
           ed_se  = ed_se,
           ced    = ced,
           ced_se = ced_se,
           agegroup = "(40,65]")
})

# -- Fitting second model to agegroup data
res_Inf <- map_df(setdiff(municipalities, village_out_Inf), function(x) {
  
  # -- Subsetting data
  tmp_dat <- filter(test_counts_Inf, village == x)
  
  # -- Number of knots to use
  nknots <- 22
  
  # -- Fitting second model
  tmp_fit <- glm(outcome ~ ns(date, df = nknots + 1, intercept = FALSE) + offset(log(expected)), 
                 family = "quasipoisson",
                 data   = tmp_dat)
  
  # -- Getting fitted values and se
  preds <- predict(tmp_fit, se.fit = TRUE, type = "response")
  
  # -- Percent change and SE 
  mu          <- tmp_dat$expected
  lambda      <- preds$fit
  dispersion  <- summary(tmp_fit)$dispersion
  cova        <- summary(tmp_fit)$cov.unscaled * dispersion
  lambda_vari <- lambda^2 * diag(model.matrix(tmp_fit) %*% cova %*% t(model.matrix(tmp_fit)))
  mu_vari     <- mu^2 * tmp_dat$log_expected_se^2
  fitted      <- lambda / mu - 1
  se          <- sqrt((lambda_vari / mu^2) + (lambda^2 * mu_vari / mu^4))
  
  # -- Excess deaths
  ed     <- lambda - mu
  ed_se  <- sqrt(lambda * dispersion + mu_vari) 
  ced    <- cumsum(ed)
  ced_se <- sqrt(cumsum(ed_se^2))
  
  # -- Putting everything together
  tmp_dat <- tmp_dat %>%
    mutate(fitted = fitted,
           se     = se,
           ed     = ed, 
           ed_se  = ed_se,
           ced    = ced,
           ced_se = ced_se,
           agegroup = "(65,Inf]")
})

# -- Putting data together
res_age <- bind_rows(res_20, res_40, res_65, res_Inf)

# -- Marginal percent change from average
pc_marginal_age <- res_age %>%
  mutate(expected_se = expected * log_expected_se) %>%
  group_by(date, agegroup) %>%
  mutate(sum_age_expected       = sum(expected, na.rm = TRUE),
         sum_age_gamma_expected = sum(fitted * expected, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(gamma_derivative = expected / sum_age_expected,
         mu_derivative    = (fitted / sum_age_expected) + (sum_age_gamma_expected / (sum_age_expected^2))) %>%
  mutate(tmp = se^2 * gamma_derivative^2 + expected_se^2 * mu_derivative^2) %>%
  group_by(date, agegroup) %>%
  summarize(fitted = sum(fitted * expected, na.rm = TRUE) / sum(expected, na.rm = TRUE),
            observed    = sum(outcome, na.rm = TRUE),
            expected    = sum(expected, na.rm = TRUE),
            expected_se = sqrt(sum(expected_se^2, na.rm = TRUE)),
            se          = sqrt(sum(tmp, na.rm = TRUE))) %>%
  ungroup()

# -- Cumulative excess deaths from the age analysis
ced_age <- res_age %>%
  group_by(date, agegroup) %>%
  summarize(ed        = sum(ed, na.rm = TRUE), 
            ed_se     = sqrt(sum(ed_se^2, na.rm = TRUE)),
            intercept = sum(intercept, na.rm = TRUE),
            expected  = sum(expected, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(agegroup) %>%
  mutate(ced         = cumsum(ed),
         ced_se      = sqrt(cumsum(ed_se^2)), 
         lwr         = ced - 1.96*ced_se,
         upr         = ced + 1.96*ced_se,
         ced_adj     = ced / intercept,
         ced_adj_lwr = ced / intercept - 1.96*ced_se / intercept,
         ced_adj_upr = ced / intercept + 1.96*ced_se / intercept) %>%
  ungroup()
##########################################################################################################
################## -- END AGE ANALYSIS -------------------------------------------------- ################
##########################################################################################################

##########################################################################################################
################## -- FIGURE 1 ---------------------------------------------------------- ################
##########################################################################################################
# -- Dates 
dates <- tibble(date = seq(make_date(2019, 01, 01), make_date(2021, 12, 31), by = "days")) %>%
  mutate(week = ifelse(week(date) == 53, 52, week(date)),
         year = year(date)) %>%
  group_by(year, week) %>%
  filter(date == first(date)) %>%
  ungroup()

# -- Wrangling mortality data
counts <- dat %>%
  filter(!village %in% c("Jetpur", "Modasa")) %>%
  mutate(year = year(date),
         week = ifelse(week(date) == 53, 52, week(date))) %>%
  dplyr::select(-date) %>%
  left_join(dates, by = c("week", "year")) %>%
  group_by(date) %>%
  summarize(outcome = n()) %>%
  ungroup() %>%
  mutate(year = year(date)) %>%
  arrange(date) %>%
  filter(date >= "2019-01-01",
         date <= "2021-04-23")

# -- Training index: Jan 2019 to Feb 2020
train_idx <- which(counts$date %in% seq(make_date(2019, 01, 01), make_date(2020, 02, 28), by = "days"))

# -- Design matrix for seasonal component
yd <- noleap_yday(counts$date)
x  <- cbind(1, fourier_trend(yd, k = 2))
colnames(x)[1] <- "Intercept"

# -- Outcome
y <- counts$outcome

# -- Putting everything together
X <- mutate(as_tibble(x), y = y)

# -- Fitting model and getting fitted values
fit <- glm(y ~ . - 1, family = "quasipoisson", data = X[train_idx,])

# -- Getting fitted values and ses
preds_po <- predict(fit, newdata = X, se.fit = TRUE, type = "link")

# -- Putting everything together
fit_res <- counts %>%
  mutate(expected        = exp(preds_po$fit),
         log_expected_se = preds_po$se.fit,
         expected_lwr    = exp(log(expected) - 1.96 * log_expected_se),
         expected_upr    = exp(log(expected) + 1.96 * log_expected_se),
         type            = "glm-po")

# -- Index for test dates
test_counts <- filter(fit_res, date >= "2020-03-01")

# -- Fitting second model
nknots <- 24
tmp_fit <- glm(outcome ~ ns(date, df = nknots + 1) + offset(log(expected)), 
               family = "quasipoisson",
               data   = test_counts)

# -- Getting fitted values and se
preds <- predict(tmp_fit, se.fit = TRUE, type = "response")

# -- Percent change and SE 
mu          <- test_counts$expected
lambda      <- preds$fit
dispersion  <- summary(tmp_fit)$dispersion
cova        <- summary(tmp_fit)$cov.unscaled * dispersion
lambda_vari <- lambda^2 * diag(model.matrix(tmp_fit) %*% cova %*% t(model.matrix(tmp_fit)))
mu_vari     <- mu^2 * test_counts$log_expected_se^2
fitted      <- lambda / mu - 1
se          <- sqrt((lambda_vari / mu^2) + (lambda^2 * mu_vari / mu^4))

# -- Excess deaths
ed     <- lambda - mu
ed_se  <- sqrt(lambda * dispersion + mu_vari) 
ced    <- cumsum(ed)
ced_se <- sqrt(cumsum(ed_se^2))

# -- Putting everything together
test_counts <- test_counts %>%
  mutate(fitted = fitted,
         se     = se,
         ed     = ed, 
         ed_se  = ed_se,
         ced    = ced,
         ced_se = ced_se)

# -- Figure 1
ggplot() +
  geom_vline(xintercept = ymd("2020-03-01"),
             linetype   = 2,
             color      = "gray") +
  geom_point(aes(date, outcome),
             color = my_palette[["black"]],
             alpha = 0.50,
             data  = fit_res) +
  geom_line(aes(date, expected),
            color = "white",
            size  = 1,
            data  = fit_res) +
  geom_line(aes(date, expected, color = "Baseline"),
            data  = fit_res) +
  geom_line(aes(date, expected * (fitted + 1)),
            color = "white",
            size  = 1,
            data  = test_counts) +
  geom_line(aes(date, expected * (fitted + 1), color = "Smooth observed"),
            data  = test_counts) +
  ggtext::geom_richtext(
    data = tibble(
      x = make_date(2020, 06, 15), 
      y = 4000,
      species = c("Onset of Covid-19"),
      lab = c("<b style='font-family:Helvetica;font-size:14pt;color:#525252'>Onset of Covid-19</b><br><i style='color:darkgrey;'>March 2020</i>"),
      angle = c(0)
    ),
    aes(x, y, label = lab, angle = angle), 
    size        = 4, 
    fill        = NA,
    color       = my_palette[["blue"]],
    lineheight  = 0.30,
    label.color = NA
  ) +
  scale_color_manual(values = as.vector(my_palette[c(1,8)])) +
  scale_x_date(date_labels = "%b") +
  scale_y_continuous(labels = scales::comma) +
  coord_cartesian(ylim = c(700, 5500)) +
  labs(x     = "Date",
       y     = "Number of deaths",
       color = "") +
  theme(legend.position = c(0.50, 0.90),
        legend.text = element_text(face = "bold"))
##########################################################################################################
################## -- END FIGURE 1 ---------------------------------------------------------- ############
##########################################################################################################

##########################################################################################################
################## -- FIGURE 2 ---------------------------------------------------------- ################
##########################################################################################################
# -- Percent change results from all three analyses
pc_all <- pc_marginal %>%
  mutate(demo = "marginal") %>%
  select_at(c("date", "demo", "fitted", "observed", "expected", "expected_se", "se")) %>%
  bind_rows(rename(pc_marginal_age, demo = agegroup)) %>% 
  bind_rows(rename(pc_marginal_sex, demo = gender)) %>%
  mutate(demo = factor(demo, levels = c("marginal", "male", "female", "(65,Inf]", "(40,65]", "(20,40]", "(0,20]")))

# -- Percent change for all the municipalities combined for Figure 2
txt <- pc_all %>%
  group_by(demo) %>%
  filter(date == "2021-04-23") %>%
  ungroup() %>%
  mutate(label = as.character(demo),
         label = case_when(demo == "marginal" ~ "Overall",
                           demo == "male" ~ "Male",
                           demo == "female" ~ "Female",
                           demo == "(0,20]" ~ "Less than 20",
                           demo == "(20,40]" ~ "20 to 40",
                           demo == "(40,65]" ~ "40 to 65",
                           demo == "(65,Inf]" ~ "65 and over"),
         txt   = paste0(" ", label,": ", 100 * round(fitted, 2), "%\n[CI 95%: ", 100 * round(fitted - 2 * se, 2), "%, ", 100 * round(fitted + 2 * se, 2), "%]")) %>%
  select(demo, label, txt) %>%
  mutate(demo = factor(demo, levels = c("marginal", "male", "female", "(65,Inf]", "(40,65]", "(20,40]", "(0,20]")))

# -- Data for figure 2
viz_dat <- pc_all %>%
  left_join(txt, by = "demo") %>%
  mutate(group = case_when(demo == "marginal" ~ "marginal",
                           demo %in% c("male", "female") ~ "gender",
                           TRUE ~ "age"),
         group = factor(group, levels = c("marginal", "gender", "age")))

# -- Figure 2: Percent change for maginal, gender, and age analyses
viz_dat %>%
  ggplot(aes(date, fitted, color = demo, fill = demo)) +
  geom_hline(yintercept = 0, 
             linetype   = 2,
             color      = "gray") +
  facet_wrap(~group, ncol = 1,
             labeller = as_labeller(c("gender"   = "Gender analysis",
                                      "age"      = "Age analysis",
                                      "marginal" = "Overall analysis"))) +
  geom_line(size  = 1, 
            color = "white") +
  geom_line(show.legend = FALSE) +
  geom_point(shape = 21,
             color = "#FBFCFC",
             data = filter(viz_dat, date == "2021-04-23"),
             show.legend = FALSE) +
  geom_dl(aes(label = txt),
          method = list("last.qp",
                        fontface = "bold",
                        cex      = 0.65)) +
  coord_cartesian(xlim = c(make_date(2020,03,01), make_date(2021,04,30) + months(2))) +
  scale_x_date(date_labels = "%b",
               breaks      = seq(make_date(2020,04,01), make_date(2021,04,01), by = "3 months")) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_manual(values = as.vector(my_palette[c(9, 2, 3, 1, 8, 9, 4)]),
                     labels = c("Overall", "Male", "Female", "65 and over", "40 to 65", "20 to 40", "Less than 20")) +
  scale_fill_manual(values = as.vector(my_palette[c(9, 2, 3, 1, 8, 9, 4)]),
                    labels = c("Overall", "Male", "Female", "65 and over", "40 to 65", "20 to 40", "Less than 20")) +
  labs(x = "Date",
       y = "Percent change in mortality from average") +
  theme(strip.text       = element_text(face = "bold", color = "black", hjust = 0, vjust = 0, size = 12),
        strip.background = element_rect(fill = "#FBFCFC", color = "#FBFCFC"),
        legend.direction = "horizontal",
        legend.title     = element_blank(),
        legend.position  = c(0.50, 0.55))
##########################################################################################################
################## -- END FIGURE 2 ---------------------------------------------------------- ############
##########################################################################################################

##########################################################################################################
################## -- FIGURE 3 ---------------------------------------------------------- ################
##########################################################################################################
# -- Figure 3: Marginal percent change for all municipalities
res_municipalities %>%
  mutate(fitted_new = ifelse(fitted <= 0, 0, fitted),
         fitted_new = ifelse(fitted >= 6, 6, fitted_new)) %>%
  ggplot(aes(date, reorder(village, intercept), fill = fitted_new)) +
  geom_tile(color = "black",
            size  = 0.10) +
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "Reds"),
                       breaks = 0:6,
                       labels = c(expression(phantom(x) <= "0%"), "100%", "200%", "300%", "400%", "500%", expression(phantom(x) >= "600%"))) +
  labs(x    = "Date",
       y    = "",
       fill = "Percent change \nfrom average") +
  scale_x_date(date_labels = "%b") +
  theme(panel.grid.major  = element_blank(),
        legend.position   = "right", 
        legend.background = element_rect(color = NA),
        legend.key.height = unit(2.5, "cm"),
        legend.direction  = "vertical",
        axis.text.x       = element_text(size = 10),
        axis.title        = element_text(size = 11),
        axis.text.y       = element_text(size = 8, hjust = 1))
##########################################################################################################
################## -- END FIGURE 3 ---------------------------------------------------------- ############
##########################################################################################################

##########################################################################################################
################## -- PERCENT CHANGE FROM AVERAGE OVER THE PANDEMIC ------------------------- ############
##########################################################################################################
# -- All the available data
viz_dat %>%
  filter(demo == "marginal") %>%
  mutate(smooth    = expected * (1 + fitted),
         smooth_se = sqrt(se^2 * expected_se^2 + expected^2 * se^2 + fitted^2 * expected_se^2)) %>%
  group_by(demo) %>%
  summarize(expected    = sum(expected),
            observed    = sum(observed),
            smooth      = sum(smooth),
            expected_se = sqrt(sum(expected_se^2)),
            smooth_se   = sqrt(sum(smooth_se^2))) %>%
  mutate(pc  = smooth /  expected - 1,
         pc2 = observed / expected - 1,
         se  = sqrt(smooth_se^2 / expected^2 + (smooth^2 * expected_se^2) / expected^4), 
         lwr = round(100 * (pc - 2 * se), 3),
         upr = round(100 * (pc + 2 * se), 3),
         pc  = round(100 * pc, 3),
         pc_lwr = pc-lwr,
         upr_pc = upr-pc) %>%
  select(demo, pc, lwr, upr, pc_lwr, upr_pc)

# -- First wave
viz_dat %>%
  filter(date >= "2020-04-29",
         date <= "2021-01-08") %>%
  mutate(smooth    = expected * (1 + fitted),
         smooth_se  = sqrt(se^2 * expected_se^2 + expected^2 * se^2 + fitted^2 * expected_se^2),
         smooth_se2 = sqrt(expected_se^2 + se^2 * expected_se^2 + expected^2 * se^2 + fitted^2 * expected_se^2)) %>%
  group_by(demo) %>%
  summarize(expected    = sum(expected),
            observed    = sum(observed),
            smooth      = sum(smooth),
            expected_se = sqrt(sum(expected_se^2)),
            smooth_se   = sqrt(sum(smooth_se^2))) %>%
  mutate(pc  = smooth /  expected - 1,
         pc2 = observed / expected - 1,
         se  = sqrt(smooth_se^2 / expected^2 + (smooth^2 * expected_se^2) / expected^4), 
         lwr = round(100 * (pc - 2 * se), 3),
         upr = round(100 * (pc + 2 * se), 3),
         pc  = round(100 * pc, 3),
         pc_lwr = pc-lwr,
         upr_pc = upr-pc) %>%
  select(demo, pc, lwr, upr, pc_lwr, upr_pc)

# -- Second wave
viz_dat %>%
  filter(date >= "2021-03-26") %>%
  mutate(smooth    = expected * (1 + fitted),
         smooth_se = sqrt(se^2 * expected_se^2 + expected^2 * se^2 + fitted^2 * expected_se^2)) %>%
  group_by(demo) %>%
  summarize(expected    = sum(expected),
            observed    = sum(observed),
            smooth      = sum(smooth),
            expected_se = sqrt(sum(expected_se^2)),
            smooth_se   = sqrt(sum(smooth_se^2))) %>%
  mutate(pc  = smooth /  expected - 1,
         pc2 = observed / expected - 1,
         se  = sqrt(smooth_se^2 / expected^2 + (smooth^2 * expected_se^2) / expected^4), 
         lwr = round(100 * (pc - 2 * se), 3),
         upr = round(100 * (pc + 2 * se), 3),
         pc  = round(100 * pc, 3),
         pc_lwr = pc-lwr,
         upr_pc = upr-pc) %>%
  select(demo, pc, lwr, upr, pc_lwr, upr_pc)
##########################################################################################################
################## -- END PERCENT CHANGE FROM AVERAGE OVER THE PANDEMIC --------------------- ############
##########################################################################################################

##########################################################################################################
################## -- TOTAL DEATH METRICS ----------------------------------------------- ################
##########################################################################################################
# -- Dates
dates <- tibble(date = seq(make_date(2019, 01, 01), make_date(2021, 12, 31), by = "days")) %>%
  mutate(week = ifelse(week(date) == 53, 52, week(date)),
         year = year(date)) %>%
  group_by(year, week) %>%
  filter(date == first(date)) %>%
  ungroup()

# -- Wrangling mortality data
counts <- dat %>%
  mutate(year = year(date),
         week = ifelse(week(date) == 53, 52, week(date))) %>%
  dplyr::select(-date) %>%
  left_join(dates, by = c("week", "year")) %>%
  group_by(date, village) %>%
  summarize(outcome = n()) %>%
  ungroup() %>%
  mutate(year = year(date)) %>%
  arrange(date) %>%
  filter(date >= "2019-01-01",
         date <= "2021-04-23")

# -- Getting rid of bad quality data
counts <- filter(counts, !village %in% c("Jetpur", "Modasa"))

# -- Total deaths recorded during the pandemic
counts %>%
  filter(date >= "2020-01-01") %>%
  summarize(total_deaths = sum(outcome))

# -- Total deaths per year (2019 to 2021)
counts %>%
  group_by(year) %>%
  summarize(total_deaths = sum(outcome))
##########################################################################################################
################## -- END TOTAL DEATH METRICS ------------------------------------------- ################
##########################################################################################################

##########################################################################################################
################## -- CUMULATIVE EXCESS DEATHS ------------------------------------------ ################
##########################################################################################################
# -- Putting cumulative excess deaths results together
mutate(ced_overall, demo = "marginal") %>%
  select_at(c("date", "demo", "ed", "ed_se",
              "intercept", "expected", "ced", "ced_se", "lwr",
              "upr", "ced_adj", "ced_adj_lwr", "ced_adj_upr")) %>%
  bind_rows(rename(ced_age, demo = agegroup)) %>%
  bind_rows(rename(ced_gender, demo = gender)) %>%
  mutate(demo = factor(demo, levels = c("(0,20]", "(20,40]", "(40,65]", "(65,Inf]", "female", "male", "marginal")),
         cexp = cumsum(expected)) %>%
  group_by(demo) %>%
  filter(date == max(date)) %>%
  ungroup() %>%
  mutate(n   = nchar(as.character(round(ced_se))) - 1,
         ced = round(ced, digits = -n),
         lwr = round(lwr, digits = -n),
         upr = round(upr, digits = -n),
         ced = prettyNum(ced, big.mark = ","),
         lwr = prettyNum(lwr, big.mark = ","),
         upr = prettyNum(upr, big.mark = ","),
         txt = paste0(ced, " [95% CI: ", lwr, " to ", upr, "]")) %>%
  arrange(demo) %>%
  select(demo, txt) %>%
  setNames(c("Group", "Cumulative excess deaths"))
##########################################################################################################
################## -- END CUMULATIVE EXCESS DEATHS -------------------------------------- ################
##########################################################################################################

##########################################################################################################
################## -- MUNICIPALITY METRICS -------------------------------------------------- ############
##########################################################################################################
# -- Putting everything together
res_municipality_all <- res_age %>%
  select(date, village, fitted, se, intercept, ced, ced_se, agegroup) %>%
  setNames(c("date", "village", "fitted", "se", "intercept", "ced", "ced_se", "demo")) %>%
  bind_rows(res_sex %>%
              select(date, village, fitted, se, intercept, ced, ced_se, gender) %>%
              setNames(c("date", "village", "fitted", "se", "intercept", "ced", "ced_se", "demo"))) %>%
  bind_rows(res_municipalities %>%
              select(date, village, fitted, se, intercept, ced, ced_se) %>%
              mutate(demo = "marginal"))

# -- Metrics used in section 3.3
res_municipality_all %>%
  group_by(village, demo) %>%
  filter(date == max(date)) %>%
  ungroup() %>%
  group_by(demo) %>%
  summarize(num_municipalities = n(),
            num_over_100       = sum(fitted >= 1),
            num_over_500       = sum(fitted >= 5),
            prop_over_100      = sum(fitted >= 1) / num_municipalities,
            prop_over_500      = sum(fitted >= 5) / num_municipalities)
##########################################################################################################
################## -- END MUNICIPALITY METRICS ---------------------------------------------- ############
##########################################################################################################

##########################################################################################################
################## -- SUPPLEMENTAL TABLE 1 ---------------------------------------------- ################
##########################################################################################################
municipality <- all_municipalities
marginal     <- municipality %in% marginal_municipalities
male         <- municipality %in% unique(res_males$village)
female       <- municipality %in% unique(res_females$village)
in20         <- municipality %in% unique(res_20$village)
in40         <- municipality %in% unique(res_40$village)
in65         <- municipality %in% unique(res_65$village)
inInf        <- municipality %in% unique(res_Inf$village)
tab <- data.frame(municipality, marginal, male, female, in20, in40, in65, inInf)
library(xtable)
xtable(tab)
##########################################################################################################
################## -- END SUPPLEMENTAL TABLE 1 ---------------------------------------------- ############
##########################################################################################################

##########################################################################################################
################## -- SUPP FIGURE 1 ---------------------------------------------------------- ###########
##########################################################################################################
# -- Dates 
dates <- tibble(date = seq(make_date(2019, 01, 01), make_date(2021, 12, 31), by = "days")) %>%
  mutate(week = ifelse(week(date) == 53, 52, week(date)),
         year = year(date)) %>%
  group_by(year, week) %>%
  filter(date == first(date)) %>%
  ungroup()

# -- Wrangling mortality data
counts <- dat %>%
  filter(!village %in% c("Jetpur", "Modasa")) %>%
  mutate(year = year(date),
         week = ifelse(week(date) == 53, 52, week(date))) %>%
  dplyr::select(-date) %>%
  left_join(dates, by = c("week", "year")) %>%
  group_by(date) %>%
  summarize(outcome = n()) %>%
  ungroup() %>%
  mutate(year = year(date)) %>%
  arrange(date) %>%
  filter(date <= "2021-05-31")

# -- Supp figure 1
counts %>%
  ggplot(aes(date, outcome, color = date <= "2021-04-23")) +
  geom_point(alpha = 0.70) +
  geom_vline(xintercept = ymd("2021-04-23"),
             linetype   = 2, 
             color      = "red3") +
  scale_color_manual(values = c("red3", "black"),
                     labels = c("Data not used", "Data used"),
                     name   = "") +
  scale_x_date(date_labels = "%b %Y") +
  scale_y_continuous(breaks = seq(0, 5000, by = 1000),
                     labels = scales::comma) +
  labs(x = "Date",
       y = "Number of deaths") +
  theme(legend.text     = element_text(face = "bold"),
        legend.position = c(0.50, 0.80))
##########################################################################################################
################## -- END SUPP FIGURE 1 ---------------------------------------------------------- #######
##########################################################################################################

##########################################################################################################
################## -- SUPP FIGURES 5 7 6 ---------------------------------------------------------- ######
##########################################################################################################
# -- Supp Figure 5
res_sex %>%
  mutate(fitted_new = ifelse(fitted <= 0, 0, fitted),
         fitted_new = ifelse(fitted >= 6, 6, fitted_new),
         gender = ifelse(gender == "male", "Male", "Female")) %>%
  ggplot(aes(date, reorder(village, intercept), fill = fitted_new)) +
  facet_wrap(~gender) +
  geom_tile(color = "black",
            size  = 0.10) +
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "Reds"),
                       breaks = 0:6,
                       labels = c(expression(phantom(x) <= "0%"), "100%", "200%", "300%", "400%", "500%", expression(phantom(x) >= "600%"))) +
  labs(x    = "Date",
       y    = "",
       fill = "Percent change \nfrom average") +
  scale_x_date(date_labels = "%b") +
  theme(panel.grid.major  = element_blank(),
        legend.position   = "right", 
        legend.background = element_rect(color = NA),
        legend.key.height = unit(2.5, "cm"),
        legend.direction  = "vertical",
        axis.text.x       = element_text(size = 10),
        axis.title        = element_text(size = 11),
        axis.text.y       = element_text(size = 8, hjust = 1),
        strip.background  = element_rect(fill = "#FBFCFC", color = NA),
        strip.text        = element_text(color = "black", hjust = 0))

# -- Supp Figure 6
res_age %>%
  mutate(fitted_new = ifelse(fitted <= 0, 0, fitted),
         fitted_new = ifelse(fitted >= 6, 6, fitted_new),
         agegroup   = case_when(agegroup == "(0,20]" ~ "Less than 20",
                                agegroup == "(20,40]" ~ "20 to 40",
                                agegroup == "(40,65]" ~ "40 to 65",
                                agegroup == "(65,Inf]" ~ "65 and over"),
         agegroup   = factor(agegroup, levels = c("Less than 20", "20 to 40", "40 to 65", "65 and over"))) %>% 
  ggplot(aes(date, reorder(village, intercept), fill = fitted_new)) +
  facet_wrap(~agegroup, ncol = 4) +
  geom_tile(color = "black",
            size  = 0.10) +
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "Reds"),
                       breaks = 0:6,
                       labels = c(expression(phantom(x) <= "0%"), "100%", "200%", "300%", "400%", "500%", expression(phantom(x) >= "600%"))) +
  labs(x    = "Date",
       y    = "",
       fill = "Percent change \nfrom average") +
  scale_x_date(date_labels = "%b") +
  theme(panel.grid.major  = element_blank(),
        legend.position   = "right", 
        legend.background = element_rect(color = NA),
        legend.key.height = unit(2.5, "cm"),
        legend.direction  = "vertical",
        axis.text.x       = element_text(size = 10),
        axis.title        = element_text(size = 11),
        axis.text.y       = element_text(size = 8, hjust = 1),
        strip.background  = element_rect(fill = "#FBFCFC", color = NA),
        strip.text        = element_text(color = "black", hjust = 0))
##########################################################################################################
################## -- END SUPP FIGURES 5 7 6 ---------------------------------------------------------- ##
##########################################################################################################