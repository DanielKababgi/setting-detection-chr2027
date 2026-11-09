# libaries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(nnet)
library(MASS)
library(mgcv)
library(ggthemes)

# Setup ----
options(scipen = 999)
set.seed(42)
my_colors <- c(
  "Other" = "#004D40",
  "Rural" = "#FFC107",
  "Urban" = "#1E88E5"
)

# Read in data ----
df_raw <- read_tsv("src/results-merged-meta.tsv")

## eyeball data ----
table(df_raw$prediction)
table(df_raw$gender)
table(df_raw$corpusSubset)

## subset for rows where 'pubYear' is not na ----
df <- df_raw[!is.na(df_raw$pubYear), ]
table(df$prediction)
ggplot(df, aes(x = pubYear)) +
  geom_histogram(binwidth = 10,
                 position = "dodge",) +
    labs(
    title="Diachronic distribution of texts",
    x = "Publication year",
    y = "Count"
  ) + theme_clean() 

# counts and proportions over time ----
## prep ----
df_binned <- df %>%
  filter(
    !is.na(pubYear),
    !is.na(prediction),
    pubYear >= 1750,
    pubYear <= 1950
  ) %>%
  mutate(
    period_start = 1750 + floor((pubYear - 1750) / 25) * 25
  )

proportions <- df_binned %>%
  count(period_start, prediction) %>%
  complete(
    period_start,
    prediction,
    fill = list(n = 0)
  ) %>%
  group_by(period_start) %>%
  mutate(
    total = sum(n),
    proportion = n / total
  ) %>%
  ungroup()

## plot: raw proportions in dataset ----
ggplot(
  proportions,
  aes(
    x = period_start,
    y = proportion,
    color = prediction,
    group = prediction
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(
   values = my_colors
  ) +
  scale_fill_viridis_d() + 
  scale_x_continuous(
    breaks = seq(1750, 1950, by = 25)
  ) +
  scale_y_continuous(
    labels = label_percent()
  ) +
  labs(
    x = "Publication year",
    y = "Proportion",
    color = "Prediction",
    title = "Proportion of setting over time from 1750–1925"
  ) +
  theme_clean()

# chi square test ----
# Is the distribution of categories independent over time period?
table_time <- table(df_binned$period_start, df_binned$prediction)
chisq.test(table_time)

# GAM analysis ----
## prep ---- 
df_gam <- df %>%
  filter(
    !is.na(pubYear),
    !is.na(prediction),
    pubYear >= 1750,
    pubYear <= 1950
  ) %>%
  mutate(
    prediction = factor(prediction)
  )

df_gam <- df_gam %>%
  mutate(
    prediction_num = as.integer(prediction) - 1
  )

table(df_gam$prediction, df_gam$prediction_num)

levels(df_gam$prediction)

K <- nlevels(df_gam$prediction) - 1

## GAM model & statistics ----
model_gam <- gam(
  list(prediction_num ~ s(pubYear, k = 10), ~ s(pubYear, k = 10)),
  family = multinom(K = K),
  data = df_gam,
  method = "REML"
)

summary(model_gam)

gam.check(model_gam)

### gam for referecen category "Rural"----
df_gam_rural <- df_gam %>%
  mutate(
    prediction = relevel(prediction, ref = "Rural"),
    prediction_num = as.integer(prediction) - 1
  )

table(
  df_gam_rural$prediction,
  df_gam_rural$prediction_num
)

model_gam_rural <- gam(
  list(
    prediction_num ~ s(pubYear, k = 10),
    ~ s(pubYear, k = 10)
  ),
  family = multinom(K = K),
  data = df_gam_rural,
  method = "REML"
)

summary(model_gam_rural)

## probability curves ----
### prep for probability curves ----
n_sims <- 2000

beta_hat <- coef(model_gam)

V <- vcov(
  model_gam,
  unconditional = TRUE
)

beta_draws <- mgcv::rmvn(
  n = n_sims,
  mu = beta_hat,
  V = V
)

categories <- levels(df_gam$prediction)

newdata <- data.frame(
  pubYear = seq(1750, 1950, by = 1)
)

prob_sims <- array(
  NA_real_,
  dim = c(
    nrow(newdata),
    length(categories),
    n_sims
  ),
  dimnames = list(
    year = newdata$pubYear,
    prediction = categories,
    simulation = NULL
  )
)

for (i in seq_len(n_sims)) {
  # Make a temporary copy of the fitted GAM
  model_sim <- model_gam
  
  # Replace fitted coefficients with posterior draw
  model_sim$coefficients <- beta_draws[i, ]
  
  # Predicted multinomial probabilities
  p <- predict(
    model_sim,
    newdata = newdata,
    type = "response"
  )
  
  prob_sims[, , i] <- p
}

fit_probs <- predict(
  model_gam,
  newdata = newdata,
  type = "response"
)

lower <- apply(
  prob_sims,
  c(1, 2),
  quantile,
  probs = 0.025
)

upper <- apply(
  prob_sims,
  c(1, 2),
  quantile,
  probs = 0.975
)

median_prob <- apply(
  prob_sims,
  c(1, 2),
  median
)

colnames(fit_probs) <- categories
colnames(lower) <- categories
colnames(upper) <- categories
colnames(median_prob) <- categories

# fitted probabilities
fit_long <- as.data.frame(fit_probs) %>%
  mutate(pubYear = newdata$pubYear) %>%
  pivot_longer(
    cols = -pubYear,
    names_to = "prediction",
    values_to = "fit"
  )

# lower interval
lower_long <- as.data.frame(lower) %>%
  mutate(pubYear = newdata$pubYear) %>%
  pivot_longer(
    cols = -pubYear,
    names_to = "prediction",
    values_to = "lower"
  )

# upper interval
upper_long <- as.data.frame(upper) %>%
  mutate(pubYear = newdata$pubYear) %>%
  pivot_longer(
    cols = -pubYear,
    names_to = "prediction",
    values_to = "upper"
  )

# posterior median
median_long <- as.data.frame(median_prob) %>%
  mutate(pubYear = newdata$pubYear) %>%
  pivot_longer(
    cols = -pubYear,
    names_to = "prediction",
    values_to = "median"
  )

plot_data <- fit_long %>%
  left_join(
    lower_long,
    by = c("pubYear", "prediction")
  ) %>%
  left_join(
    upper_long,
    by = c("pubYear", "prediction")
  ) %>%
  left_join(
    median_long,
    by = c("pubYear", "prediction")
  )

### plot: probability curves----
ggplot(
  plot_data,
  aes(
    x = pubYear,
    y = fit,
    color = prediction,
    fill = prediction
  )
) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.2,
    color = NA
  ) +
  geom_line(
    linewidth = 1
  ) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  scale_x_continuous(
    breaks = seq(1750, 1950, by = 25)
  ) +
  scale_y_continuous(
    labels = label_percent(),
    limits = c(0, 1)
  ) +
  labs(
    x = "Publication year",
    y = "Setting probability",
    color = "Prediction",
    fill = "Prediction",
    title = "Setting probability for setting over time, GAM"
  ) +
  theme_clean()

## difference Urban - Rural ----
### prep ----

# allowed values: Other, Urban, Rural
category1 <- "Rural" 
category2 <- "Other"

difference_sims <- prob_sims[, category1, ] -
  prob_sims[, category2, ]

difference_df <- data.frame(
  pubYear = newdata$pubYear,
  
  difference = apply(
    difference_sims,
    1,
    median
  ),
  
  lower = apply(
    difference_sims,
    1,
    quantile,
    probs = 0.025
  ),
  
  upper = apply(
    difference_sims,
    1,
    quantile,
    probs = 0.975
  )
)

difference_df <- difference_df %>%
  mutate(
    significant = case_when(
      lower > 0 ~ paste(category1, ">", category2),
      upper < 0 ~ paste(category2, ">", category1),
      TRUE ~ "No clear difference"
    )
  )

difference_df %>%
  filter(significant != "No clear difference")

dim(difference_sims)

fit_probs <- predict(
  model_gam,
  newdata = newdata,
  type = "response"
)

colnames(fit_probs) <- levels(df_gam$prediction)

fit_difference <- fit_probs[, category1] -
  fit_probs[, category2]

difference_sd <- apply(
  difference_sims,
  1,
  sd
)

# Avoid division by zero
difference_sd <- pmax(
  difference_sd,
  .Machine$double.eps
)

standardized_draws <- sweep(
  difference_sims,
  1,
  fit_difference,
  FUN = "-"
)

standardized_draws <- sweep(
  standardized_draws,
  1,
  difference_sd,
  FUN = "/"
)

max_deviation <- apply(
  abs(standardized_draws),
  2,
  max,
  na.rm = TRUE
)

critical_value <- quantile(
  max_deviation,
  probs = 0.95,
  na.rm = TRUE
)

critical_value

difference_df <- data.frame(
  pubYear = newdata$pubYear,
  difference = fit_difference,
  sd = difference_sd
) %>%
  mutate(
    lower = difference - critical_value * sd,
    upper = difference + critical_value * sd
  )

difference_df <- difference_df %>%
  mutate(
    significant = lower > 0 | upper < 0,
    
    direction = case_when(
      lower > 0 ~ paste(category1, " > ", category2),
      upper < 0 ~ paste(category2, " > ", category1),
      TRUE ~ "No clear difference"
    )
  )

difference_df %>%
  filter(significant)

### significant ranges ----
significant_ranges <- difference_df %>%
  mutate(
    sig_direction = if_else(
      significant,
      direction,
      NA_character_
    ),
    
    new_group =
      is.na(sig_direction) |
      sig_direction != dplyr::lag(
        sig_direction,
        default = first(sig_direction)
      ) |
      pubYear != dplyr::lag(
        pubYear,
        default = first(pubYear)
      ) + 1,
    
    group = cumsum(new_group)
  ) %>%
  
  filter(significant) %>%
  
  group_by(group, direction) %>%
  
  summarise(
    start_year = min(pubYear),
    end_year = max(pubYear),
    .groups = "drop"
  ) %>%
  
  dplyr::select(
    direction,
    start_year,
    end_year
  )

significant_ranges

### plot: significant ranges ----
ggplot(
  difference_df,
  aes(
    x = pubYear,
    y = difference
  )
) +
  
  # Highlight significant periods
  geom_rect(
    data = significant_ranges,
    aes(
      xmin = start_year,
      xmax = end_year,
      ymin = -Inf,
      ymax = Inf,
      fill = direction
    ),
    inherit.aes = FALSE,
    alpha = 0.15
  ) +
  
  # Simultaneous credible interval
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.2
  ) +
  
  # Estimated probability difference
  geom_line(
    linewidth = 1
  ) +
  
  # Zero = equal predicted probabilities
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_vline(
    data = significant_ranges,
    aes(xintercept = start_year),
    linetype = "dotted",
    alpha = 0.5
  ) +
  
  geom_vline(
    data = significant_ranges,
    aes(xintercept = end_year),
    linetype = "dotted",
    alpha = 0.5
  ) +
  
  scale_x_continuous(
    breaks = seq(1750, 1950, by = 25)
  ) +
  
  scale_y_continuous(
    labels = label_percent()
  ) +
  
  labs(
    x = "Publication year",
    y = "Difference in predicted probability",
    fill = "Significant range",
    title=paste(
      "Difference in predicted probability:",
      category1, " - ", category2
    )
  ) +
  
  theme_clean() #+
  #theme(
  #   legend.position = ""
  #)

## power analysis ----
simulate_gam_power <- function(
    model_gam,
    nsim = 1000,
    alpha = 0.05,
    effect_scale = 1
) {

  

  # 1. Get exactly the data used by the model

  dat <- model.frame(model_gam)
  
  if (!all(c("prediction_num", "pubYear") %in% names(dat))) {
    stop("The model data must contain prediction_num and pubYear.")
  }
  
  n <- nrow(dat)
  

  # 2. Predicted probabilities

  probs <- predict(
    model_gam,
    type = "response"
  )
  
  # Mean category probabilities:
  # represents a model with no publication-year effect
  mean_probs <- colMeans(probs)
  
  null_probs <- matrix(
    mean_probs,
    nrow = n,
    ncol = ncol(probs),
    byrow = TRUE
  )
  
  # 3. Scale the temporal effect

  # effect_scale = 0   -> no year effect
  # effect_scale = 0.5 -> half observed effect
  # effect_scale = 1   -> fitted/observed effect
  if (effect_scale < 0 || effect_scale > 1) {
    stop("effect_scale must be between 0 and 1.")
  }
  
  sim_probs <- null_probs +
    effect_scale * (probs - null_probs)
  
  # Storage for p-values
  p_values <- matrix(
    NA_real_,
    nrow = nsim,
    ncol = 2
  )
  
  colnames(p_values) <- c(
    "smooth_1",
    "smooth_2"
  )
  
  # 4. Simulation loop

  for (i in seq_len(nsim)) {
    
    # Generate one category (0, 1, 2)
    # for every observation
    y_sim <- apply(
      sim_probs,
      1,
      function(p) {
        sample(
          0:2,
          size = 1,
          prob = p
        )
      }
    )
    
    sim_dat <- dat
    sim_dat$prediction_num <- y_sim
    
    # Refit exactly the same model
    fit <- try(
      gam(
        list(
          prediction_num ~ s(pubYear, k = 10),
          ~ s(pubYear, k = 10)
        ),
        family = multinom(K = 2),
        data = sim_dat,
        method = "REML"
      ),
      silent = TRUE
    )
    
    # Only use successful fits
    if (!inherits(fit, "try-error")) {
      
      p <- summary(fit)$s.pv
      
      if (length(p) >= 2) {
        p_values[i, ] <- p[1:2]
      }
    }
  }
  
  # 5. Calculate power
  valid <- complete.cases(p_values)
  
  p_valid <- p_values[valid, , drop = FALSE]
  
  power_1 <- mean(p_valid[, 1] < alpha)
  power_2 <- mean(p_valid[, 2] < alpha)
  
  power_any <- mean(
    p_valid[, 1] < alpha |
      p_valid[, 2] < alpha
  )
  
  power_both <- mean(
    p_valid[, 1] < alpha &
      p_valid[, 2] < alpha
  )
  
  # Monte Carlo standard errors
  mcse_1 <- sqrt(
    power_1 * (1 - power_1) /
      nrow(p_valid)
  )
  
  mcse_2 <- sqrt(
    power_2 * (1 - power_2) /
      nrow(p_valid)
  )
  
  # 6. Return results
  list(
    n = n,
    simulations = nsim,
    successful_simulations = nrow(p_valid),
    effect_scale = effect_scale,
    
    power = data.frame(
      smooth = c(
        "s(pubYear)",
        "s.1(pubYear)"
      ),
      power = c(
        power_1,
        power_2
      ),
      mcse = c(
        mcse_1,
        mcse_2
      )
    ),
    
    power_any = power_any,
    power_both = power_both,
    
    p_values = p_values
  )
}

### run simulation ----
power_result <- simulate_gam_power(
  model_gam,
  nsim = 1000,
  alpha = 0.05,
  effect_scale = 1
)

power_result$power

### plot several simulations ----

power_25 <- simulate_gam_power(
  model_gam,
  nsim = 1000,
  effect_scale = 0.25
)

power_50 <- simulate_gam_power(
  model_gam,
  nsim = 1000,
  effect_scale = 0.50
)

power_75 <- simulate_gam_power(
  model_gam,
  nsim = 1000,
  effect_scale = 0.75
)

power_100 <- simulate_gam_power(
  model_gam,
  nsim = 1000,
  effect_scale = 1
)

power_curve <- data.frame(
  effect_scale = c(0.25, 0.50, 0.75, 1.00),
  
  smooth_1 = c(
    power_25$power$power[1],
    power_50$power$power[1],
    power_75$power$power[1],
    power_100$power$power[1]
  ),
  
  smooth_2 = c(
    power_25$power$power[2],
    power_50$power$power[2],
    power_75$power$power[2],
    power_100$power$power[2]
  )
)

power_curve

power_curve_long <- power_curve |>
  pivot_longer(
    cols = c(smooth_1, smooth_2),
    names_to = "smooth",
    values_to = "power"
  )

ggplot(
  power_curve_long,
  aes(
    x = effect_scale,
    y = power,
    color = smooth
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(
    yintercept = 0.80,
    linetype = "dashed"
  ) +
  scale_y_continuous(
    labels = percent,
    limits = c(0, 1)
  ) +
  labs(
    x = "Proportion of fitted temporal effect",
    y = "Statistical power",
    color = "Smooth"
  ) +
  theme_minimal()

### no temporal effect ----
type1_result <- simulate_gam_power(
  model_gam,
  nsim = 1000,
  alpha = 0.05,
  effect_scale = 0
)

type1_result$power
