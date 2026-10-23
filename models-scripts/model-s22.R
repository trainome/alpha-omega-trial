source("./R/packages.R")
source("./R/data-prep.R")
source("./R/helpers.R")

refit.models <- TRUE


dat <- data_load_s7s10()

d <- dat$d

# Moderator: baseline (T1) red blood cell omega-3 index. The participant index
# P in the Stan data is dense_rank(paste0("p", participant)) over criteria_dat
# (see data_load_s7s10()), the same ranking is used here to match the index to
# P. The index is standardised to mean 0 and SD 1 over participants with an
# observed baseline value, missing values (-999) are imputed in the model.
o3_base <- criteria_dat |>
  distinct(participant) |>
  mutate(P = dense_rank(paste0("p", participant))) |>
  left_join(
    AOData::omega3 |>
      filter(time == "T1") |>
      select(participant, omega3index),
    by = "participant"
  ) |>
  mutate(o3_z = (omega3index - mean(omega3index, na.rm = TRUE)) /
           sd(omega3index, na.rm = TRUE)) |>
  arrange(P)

stopifnot(nrow(o3_base) == d$N_P, all(o3_base$P == seq_len(d$N_P)))

d$o3_obs <- if_else(is.na(o3_base$o3_z), -999, o3_base$o3_z)
d$o3_miss <- if_else(is.na(o3_base$o3_z), 0L, 1L)

# Run model

# Compile model (with methods)
s22 <- cmdstan_model(
  stan_file = "models/s22.stan",
  exe_file = "models/s22.exe",
  compile_model_methods = TRUE,
  force_recompile = TRUE
)


# Fit and save the model with fit wrap. This will be possible to extract again
# in the qmd file.
fit_wrap(
  s22,
  DATA = d,
  fit = "fit.s22.post",
  init = make_inits(s22, DATA = d),
  refit = TRUE,
  ad = 0.95,
  iterwu = 2000
)
# Initialize methods for loo
fit.s22.post$init_model_methods()

# Use loo with moment matching
loo.s22 <- fit.s22.post$loo(
  variables = "log_lik",
  moment_match = TRUE,
  cores = 4
)

saveRDS(loo.s22, "fits/paper1/loo_s22.RDS")
