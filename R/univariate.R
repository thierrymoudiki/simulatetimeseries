#' Simulate a univariate time series dataset 1
#'
#' @param n numerical, number of data points
#' @param trend string, "linear" or "quadratic"
#' @param seasonality string, "none" or "sinusoidal"
#' @param distribution string, "normal" and "student"
#' @param noise_sd numerical, standard deviation of noise
#'
#' @return a native time series object
#' @export
#'
#' @examples
#'
#' ts_data <-
#' simulate_time_series_1(
#'   n = 100L,
#'   trend = "quadratic",
#'   seasonality = "sinusoidal",
#'   noise_sd = 2500,
#'   distribution = "normal"
#' )
#' plot(ts_data, type = "l", main = "Simulated Time Series")
#'
simulate_time_series_1 <- function(n,
                                   trend = c("linear",
                                             "quadratic"),
                                   seasonality = c("none",
                                                   "sinusoidal"),
                                   distribution = c("normal",
                                                    "student"),
                                   noise_sd = 10) {
  trend <- match.arg(trend)
  seasonality <- match.arg(seasonality)
  distribution <- match.arg(distribution)
  # Generate time index
  times <- 1:n
  # Generate trend component
  if (!is.null(trend)) {
    if (trend == "linear") {
      trend_component <- seq(1, n, length.out = n)
    } else if (trend == "quadratic") {
      trend_component <- seq(1, n, length.out = n) ^ 2
    }
  } else {
    trend_component <- rep(0, n)
  }

  # Generate seasonal component
  if (!identical(seasonality, "none")) {
    if (identical(seasonality, "sinusoidal")) {
      if (trend == "linear")
      {
        seasonal_component <- sin(2 * pi * times / 12)
      } else {
        seasonal_component <- sin(2 * pi * (times ** 2) / 12)
      }
    } else {
      stop("Invalid seasonality type. Supported types: sinusoidal.")
    }
  } else {
    seasonal_component <- rep(0, n)
  }

  # Generate noise component
  if (distribution == "normal") {
    noise <- rnorm(n, mean = 0, sd = noise_sd)
  } else if (distribution == "student") {
    noise <- noise_sd * rt(n = 100, df = 3)
  }

  # Combine components to generate time series
  time_series <- trend_component + seasonal_component + noise

  return(ts(time_series))
}

#' Simulate a univariate time series dataset 2
#'
#' @param n numerical, number of data points
#' @param trend string, "linear" or "sinusoidal"
#' @param seasonality string, "none" or "sinusoidal"
#' @param noise_sd numerical, standard deviation of noise
#' @param ar autoregressive order
#' @param ma moving average order
#'
#' @return a native time series object
#' @export
#'
#' @examples
#'
#' ts_data <-
#' simulate_time_series_2(
#'   n = 100L,
#'   trend = "sinusoidal",
#'   seasonality = TRUE,
#'   noise_sd = runif(n = 1, min = 20, max=50)
#' )
#' plot(ts_data, type = "l", main = "Simulated Time Series")
#'
simulate_time_series_2 <- function(n,
                                   trend = c("linear",
                                             "sinusoidal"),
                                   seasonality = FALSE,
                                   noise_sd = 0.1,
                                   ar = 0,
                                   ma = 0) {

  trend <- match.arg(trend)
  # Generate base series
  series <- rnorm(n, mean = 0, sd = noise_sd)

  # Add trend
  times <- (1:n)
  if (trend == "linear") {
    series <- series + times * -2.9  # Adjust slope as desired
  } else if (trend == "sinusoidal") {
    series <- series + sin(2 * pi * times**2 / n)
  }

  # Add seasonality
  if (seasonality) {
    seasonal <- stats::rnorm(n, mean = 0, sd = noise_sd / 2)
    seasonal <-
      seasonal[seq(1, n, by = length(seasonal))]  # Repeat for length of series
    series <- series + seasonal
  }

  # Introduce autoregression (AR)
  if (ar > 0) {
    ar.coef <- stats::arima.sim(model = list(ar = ar),
                                n = n,
                                innov = noise_sd)
    series <- series + ar.coef * series[-1]  # Shift for causality
  }

  # Introduce moving average (MA)
  if (ma > 0) {
    ma.coef <- stats::arima.sim(model = list(ma = ma),
                                n = n,
                                innov = noise_sd)
    series <- series + ma.coef * rnorm(n, mean = 0, sd = noise_sd)
  }
  return(series)
}



#' Simulate a univariate time series dataset 3
#'
#' @param n numerical, number of data points
#'
#' @return a native time series object
#' @export
#'
#' @examples
#'
#' print(simulate_time_series_3(10))
#'
simulate_time_series_3 <- function(n=100)
{
  xi <- rnorm(n, mean = 0, sd = sqrt(0.01))
  eps <- rep(0, 100)
  for (i in 2:100) {
    eps[i] <- 0.99 * eps[i - 1] + xi[i]
  }
  trend <- seq_along(100)
  season_term <- 2 * pi * trend / 180
  return(ts(cos(season_term) + sin(season_term) + 0.01 * trend + eps))
}
