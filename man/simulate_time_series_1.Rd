% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/univariate.R
\name{simulate_time_series_1}
\alias{simulate_time_series_1}
\title{Simulate a univariate time series dataset 1}
\usage{
simulate_time_series_1(
  n,
  trend = c("linear", "quadratic"),
  seasonality = c("none", "sinusoidal"),
  distribution = c("normal", "student"),
  noise_sd = 10,
  seed = 123
)
}
\arguments{
\item{n}{numerical, number of data points}

\item{trend}{string, "linear" or "quadratic"}

\item{seasonality}{string, "none" or "sinusoidal"}

\item{distribution}{string, "normal" and "student"}

\item{noise_sd}{numerical, standard deviation of noise}

\item{seed}{int, reproducibility seed}
}
\value{
a native time series object
}
\description{
Simulate a univariate time series dataset 1
}
\examples{

ts_data <-
simulate_time_series_1(
  n = 100L,
  trend = "quadratic",
  seasonality = "sinusoidal",
  noise_sd = 2500,
  distribution = "normal"
)
plot(ts_data, type = "l", main = "Simulated Time Series")

}
