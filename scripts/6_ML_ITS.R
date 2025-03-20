#-------------------------------------------------------------------------------#
#---------Quasi-experimental methods for climate epidemiology-------------------# 
#-------------------------Drexel CCUH Workshop----------------------------------#
#-----------------------------Date:03/21/25-------------------------------------#
#--Tarik Benmarhnia (tbenmarhnia@health.ucsd.edu) & Yiqun Ma (yim022@ucsd.edu)--#
#-------------------------------------------------------------------------------#
#--------------------code written by Arnab K. Dey-------------------------------#
#------------this is a simplified script to only show key steps-----------------#
#------------see details in https://github.com/benmarhnia-lab/ML_ITS------------#

#---------------------------- Step 1: Model Tuning -----------------------------#
#------------------------------------ PHXGB ------------------------------------#
# load libraries ----------------------------------------------------------------
rm(list = ls())
set.seed(0112358)
pacman::p_load(here, tidymodels, tidyverse, modeltime, timetk, tictoc)

# load data ----------------
df_train_test <- read.csv(here("Data", "df-train-test-sf.csv")) |> mutate(date = as.Date(date))

# plit data into training and test sets -------------------------------------
set.seed(0112358)
splits_resp <- df_train_test |>
  time_series_split(
    assess = "24 months",
    cumulative = TRUE,
    date_var = date
  )

## resample data ----
set.seed(0112358)
resamples_kfold_resp <- training(splits_resp) |> 
  time_series_cv(
    assess = "12 months",     # Length of each assessment period
    initial = "5 years",     # Initial training period
    slice_limit = 10,        # Number of slices to create
    cumulative = TRUE       # Use expanding window
  )

# recipe for modeling ---------------------------------------------------
rec_obj_phxgb <- recipe(respiratory ~ ., training(splits_resp)) |>
  # Time series features 
  step_timeseries_signature(date) |>
  # Lags
  step_lag(pm25_diff, lag = 1:14) |>
  
  # Basic seasonal components
  step_fourier(date, period = 365, K = 3) |>   
  
  # cleaning steps
  step_rm(matches("(.iso$)|(.xts$)")) |>
  # step_rm(matches("county")) |>
  step_normalize(matches("(index.num$)|(_year$)")) |>
  step_dummy(all_nominal())

## review the recipe
rec_obj_phxgb |> prep() |> juice() |> colnames()

# specify models ---------------------------------------------------
model_phxgb_tune <- prophet_boost(
  mode = "regression",
  growth = tune(),
  changepoint_range = tune(),
  seasonality_yearly = tune(),
  prior_scale_changepoints = tune(),
  prior_scale_seasonality = tune(),
  #xgboost  
  mtry = tune(),
  min_n = tune(),
  tree_depth = tune(),
  learn_rate = tune(),
  loss_reduction = tune(),
  stop_iter = tune()
) |>
  set_engine("prophet_xgboost",
             set.seed = 0112358)

# generate grid for tuning ---------------------------------------------------
grid_phxgb_tune <- grid_space_filling(
  extract_parameter_set_dials(model_phxgb_tune) |>
    update(
      # Prophet parameters
      growth = growth(values = c("linear")),
      changepoint_range = changepoint_range(range = c(0.4, 0.8), trans = NULL), # Wider range
      seasonality_yearly = seasonality_yearly(values = c(TRUE)),
      prior_scale_changepoints = prior_scale_changepoints(
        range = c(0.01, 0.5),  
        trans = NULL
      ),
      prior_scale_seasonality = prior_scale_seasonality(
        range = c(0.1, 2.0),   
        trans = NULL
      ),
      
      # XGBoost parameters
      mtry = mtry(range = c(4, 25), trans = NULL),
      min_n = min_n(range = c(1L, 12L), trans = NULL),
      tree_depth = tree_depth(range = c(8, 20), trans = NULL),
      learn_rate = learn_rate(range = c(0.001, 0.1), trans = NULL),
      loss_reduction = loss_reduction(range = c(-12, 2), trans = log10_trans()),
      stop_iter = stop_iter(range = c(10L, 30L), trans = NULL)
    ),
  size = 100
)

# workflow for tuning ---------------------------------------------------
wflw_phxgb_tune <- workflow() |>
  add_model(model_phxgb_tune) |>
  add_recipe(rec_obj_phxgb)

# model tuning ---------------------------------------------------
tic(quite = FALSE)
set.seed(0112358)
tune_results_phxgb <- wflw_phxgb_tune |>
  tune_grid(
    resamples = resamples_kfold_resp,
    grid = grid_phxgb_tune,
    # control parameters
    control = control_grid(
      verbose = TRUE,
      allow_par = TRUE,
      save_pred = TRUE,
      save_workflow = TRUE,
      parallel_over = "resamples",
      event_level = "first",
      pkgs = c("tidymodels", "modeltime", "timetk")
    ),
    metrics = metric_set(rmse, rsq)
  )
toc() 

# save the results ---------------------------------------------------
rm(df_train_test)
save.image(here("Outputs", "model-tune-phxgb-final.RData"))

#---------------------------- Output Generation -----------------------------#
# Function to generate forecast intervals
generate_forecast_intervals <- function(
    model_spec,             
    training_data,    
    forecast_horizon_data,          
    n_iterations = 1000    
) {
  forecast_predictions <- list()
  
  pb <- txtProgressBar(min = 0, max = n_iterations, style = 3)
  
  forecast_dates <- unique(forecast_horizon_data$date)
  n_dates <- length(forecast_dates)
  
  forecast_matrix <- matrix(NA, nrow = n_dates, ncol = n_iterations)
  
  for(i in 1:n_iterations) {
    iter_seed <- sample.int(.Machine$integer.max, 1)
    set.seed(iter_seed)
    
    tryCatch({
      calibrated_models <- model_spec |>
        modeltime_calibrate(
          new_data = training_data,
          quiet = TRUE
        )
      
      final_models <- calibrated_models |>
        modeltime_refit(data = training_data)
      
      predictions <- final_models |>
        modeltime_forecast(
          new_data = forecast_horizon_data,
          actual_data = training_data,
          conf_interval = 0.95,
          conf_method = "conformal_default",
          set.seed = iter_seed,
          bootstrap_time = 100,
          allow_parallel = FALSE
        )
      
      pred_df <- predictions |>
        filter(.model_desc != "ACTUAL")
      forecast_matrix[, i] <- pred_df$.value
      
    }, error = function(e) {
      warning(sprintf("Iteration %d failed: %s", i, e$message))
      forecast_matrix[, i] <- NA
    })
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  forecast_matrix <- forecast_matrix[, colSums(!is.na(forecast_matrix)) > 0]
  
  results <- data.frame(
    .index = forecast_dates,
    .value = rowMeans(forecast_matrix, na.rm = TRUE),
    .pred_lo = apply(forecast_matrix, 1, quantile, probs = 0.025, na.rm = TRUE),
    .pred_hi = apply(forecast_matrix, 1, quantile, probs = 0.975, na.rm = TRUE)
  )
  
  attr(results, "n_successful") <- ncol(forecast_matrix)
  attr(results, "n_attempted") <- n_iterations
  
  return(results)
}

# ensure consistent numeric precision ----------------------------------------------
options(digits = 7)
options(scipen = 999)

# load data ---------------------------------------------------
df_preintervention <- read.csv(here("Data", "df-train-test-sf.csv")) |> mutate(date = as.Date(date))

df_all_cases <- read.csv(here("Data", "df-predict-sf.csv")) |> mutate(date = as.Date(date))

# forecast on all cases with bootstrapped CIs ------------------------------------------
forecast_cis <- generate_forecast_intervals(
  model_spec = model_tbl_best_phxgb,
  training_data = df_preintervention,
  forecast_horizon_data = df_all_cases,
  n_iterations = 1000
)
colnames(forecast_cis) <- c("date", "respiratory_pred", "conf_lo", "conf_hi")

# Merge with actuals ---------------------------------------------------
df_forecast <- df_all_cases |>
  rename(respiratory_actual = respiratory) |>
  select(date, respiratory_actual) |>
  left_join(forecast_cis, by = "date") 

# save final predictions ---------------------------------------------------
df_forecast |> saveRDS(here("Outputs", "final-preds-phxgb.rds"))

# save final predictions
df_forecast <- readRDS(here("Outputs", "final-preds-phxgb.rds"))

# subset
data.period <- df_forecast |>
  filter(date >= as.Date("2018-11-09") & date <= as.Date("2018-11-20")) |>
  dplyr::select(-date) |>
  mutate(period = "main event")

# calculate and format
data.period <- data.period|>
  group_by(period) |>
  summarise(observed = sum(respiratory_actual),
            expected = sum(respiratory_pred),
            expected_low = sum(conf_lo),
            expected_up = sum(conf_hi)) |>
  mutate(excess = observed - expected,
         excess_low = observed - expected_up,
         excess_up = observed - expected_low,
         excess_pct = excess / observed * 100,
         excess_low_pct = excess_low / observed * 100,
         excess_up_pct = excess_up / observed * 100) |>
  mutate(expected_CI = paste0(round(expected),
                              " (",
                              round(expected_low),
                              ", ",
                              round(expected_up),
                              ")"),
         excess_CI = paste0(round(excess),
                            " (",
                            round(excess_low),
                            ", ",
                            round(excess_up),
                            ")"),
         excess_pct_CI = paste0(round(excess_pct, 1),
                                " (",
                                round(excess_low_pct, 1),
                                ", ",
                                round(excess_up_pct, 1),
                                ")")) |>
  dplyr::select(period, observed, expected_CI, excess_CI, excess_pct_CI)

print(data.period)

