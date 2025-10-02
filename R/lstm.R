#' LSTM Forecasting for Financial Returns
#'
#' @param y Univariate time series of returns
#' @param h Forecast horizon
#' @param level Confidence levels for prediction intervals
#' @param lookback Number of previous periods to use for prediction
#' @param units Number of LSTM units
#' @param epochs Training epochs
#' @param batch_size Batch size for training
#'
#' @return Object of class 'forecast' with LSTM predictions
lstmf <- function(y, h = 10, level = c(80, 95), 
                         lookback = 20, units = 50, 
                         epochs = 100, batch_size = 32) {
  # Input validation
  if (!requireNamespace("keras", quietly = TRUE)) {
    stop("keras package required for LSTM forecasting")
  }
  # Ensure returns data (handle potential price input)
  if (any(y <= 0)) {
    warning("Negative values detected - ensure input is returns, not prices")
  }
  # Scale the returns data
  y_scaled <- scale(as.numeric(y))
  y_mean <- attr(y_scaled, "scaled:center")
  y_sd <- attr(y_scaled, "scaled:scale")
  # Create sequences for LSTM
  create_sequences <- function(data, lookback) {
    x <- array(dim = c(length(data) - lookback, lookback, 1))
    y <- array(dim = c(length(data) - lookback, 1))
    
    for (i in 1:(length(data) - lookback)) {
      x[i,,1] <- data[i:(i + lookback - 1)]
      y[i,1] <- data[i + lookback]
    }
    return(list(x = x, y = y))
  }
  
  sequences <- create_sequences(y_scaled, lookback)
  x_train <- sequences$x
  y_train <- sequences$y
  # Build LSTM model
  model <- keras::keras_model_sequential() %>%
    keras::layer_lstm(units = units, 
                      input_shape = c(lookback, 1),
                      return_sequences = FALSE) %>%
    keras::layer_dropout(rate = 0.2) %>%
    keras::layer_dense(units = 1)
  # Compile model
  model %>% keras::compile(
    optimizer = keras::optimizer_adam(learning_rate = 0.001),
    loss = "mse",
    metrics = c("mae")
  )
  # Train model
  history <- model %>% keras::fit(
    x_train, y_train,
    epochs = epochs,
    batch_size = batch_size,
    validation_split = 0.2,
    verbose = 0,
    callbacks = list(
      keras::callback_early_stopping(patience = 10, restore_best_weights = TRUE)
    )
  )
  # Generate forecasts recursively
  last_sequence <- tail(y_scaled, lookback)
  forecasts_scaled <- numeric(h)
  current_sequence <- last_sequence
  
  for (i in 1:h) {
    # Reshape for prediction
    input_seq <- array(current_sequence, dim = c(1, lookback, 1))
    # Predict next value
    pred <- predict(model, input_seq, verbose = 0)[1,1]
    forecasts_scaled[i] <- pred
    # Update sequence (remove first, add prediction)
    current_sequence <- c(current_sequence[-1], pred)
  }
  
  # Rescale forecasts back to original units
  forecasts <- forecasts_scaled * y_sd + y_mean
  # Create prediction intervals using empirical residuals
  # (Simplified approach - could be enhanced with conformal prediction)
  training_preds <- predict(model, x_train, verbose = 0)[,1]
  residuals <- as.numeric(y_train) - training_preds
  residual_sd <- sd(residuals, na.rm = TRUE)
  # Generate prediction intervals
  upper <- lower <- list()
  for (lvl in level) {
    z <- qnorm(1 - (1 - lvl/100)/2)
    error_margin <- z * residual_sd * sqrt(1:h)
    
    upper[[as.character(lvl)]] <- forecasts + error_margin
    lower[[as.character(lvl)]] <- forecasts - error_margin
  }
  # Return forecast object compatible with forecast package
  result <- list(
    method = "LSTM Returns",
    mean = ts(forecasts, frequency = frequency(y), 
              start = tsp(y)[2] + 1/tsp(y)[3]),
    x = y,
    fitted = ts(c(rep(NA, lookback), 
                  training_preds * y_sd + y_mean), 
                frequency = frequency(y), start = start(y)),
    residuals = ts(c(rep(NA, lookback), residuals * y_sd), 
                   frequency = frequency(y), start = start(y)),
    level = level,
    upper = do.call(cbind, upper),
    lower = do.call(cbind, lower)
  )
  
  class(result) <- "forecast"
  return(result)
}

# Example usage matching thetaf style:
# returns <- diff(log(eu_stock_prices))  # Convert prices to log returns
# lstm_forecast <- lstm_returns(returns, h = 10)
# plot(lstm_forecast)