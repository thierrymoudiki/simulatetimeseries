% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/univariate.R
\name{simulate_time_series_2}
\alias{simulate_time_series_2}
\title{Simulate a univariate time series dataset 2}
\usage{
simulate_time_series_2(
  n,
  trend = c("linear", "sinusoidal"),
  seasonality = FALSE,
  noise_sd = 0.1,
  ar = 0,
  ma = 0,
  seed = 123
)
}
\arguments{
\item{n}{numerical, number of data points}

\item{trend}{string, "linear" or "sinusoidal"}

\item{seasonality}{string, "none" or "sinusoidal"}

\item{noise_sd}{numerical, standard deviation of noise}

\item{ar}{autoregressive order}

\item{ma}{moving average order}

\item{seed}{int, reproducibility seed}
}
\value{
a native time series object
}
\description{
Simulate a univariate time series dataset 2
}
\examples{

ts_data <-
simulate_time_series_2(
  n = 100L,
  trend = "sinusoidal",
  seasonality = TRUE,
  noise_sd = runif(n = 1, min = 20, max=50)
)
plot(ts_data, type = "l", main = "Simulated Time Series")

}
