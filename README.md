# simulatetimeseries

**Simulate complex synthetic time series for benchmarks**

Read vignettes too. 

```R
ts_data <- simulate_time_series_3(n = n_series)
ts_data2 <- simulate_time_series_3(n = n_series)
par(mfrow=c(1, 2))
plot(ts_data, type = "l", main = "Simulated Time Series")
plot(ts_data2, type = "l", main = "Simulated Time Series")
```
