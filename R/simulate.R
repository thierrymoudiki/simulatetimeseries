# Simulate using bootstrap -----


#' Simulate using bootstrap resampling
#'
#' Generates bootstrap samples from a numeric vector, optionally producing multiple replicates.
#'
#' @param x Numeric vector to resample.
#' @param n Integer. Number of samples per replicate (default: length of x).
#' @param p Integer. Number of replicates (default: 1).
#' @param seed Integer. Random seed for reproducibility (default: 123).
#'
#' @return A vector or matrix of bootstrap samples.
#'
#' @examples
#' x <- rnorm(10)
#' rbootstrap(x, n = 10, p = 3)
#'
#' @export
rbootstrap <- function(x,
                       n = length(x),
                       p = 1,
                       seed = 123) {
  if (p <= 1)
  {
    set.seed(seed)
    return(sample(x, size = n, replace = TRUE))
  } else {
    return(sapply(1:p,
                  function(i) {
                    set.seed(seed + i - 1)
                    sample(x, size = n, replace = TRUE)
                  }))
  }
}

# Simulate Gaussian kernel density -----

#' Simulate Gaussian Kernel Density
#'
#' Generates samples from a Gaussian kernel density estimate of a numeric vector.
#'
#' @param x Numeric vector to estimate density from.
#' @param n Integer. Number of samples per replicate (default: length of x).
#' @param p Integer. Number of replicates (default: 1).
#' @param seed Integer. Random seed for reproducibility (default: 123).
#' @param method Character. Sampling method: "antithetic" or "traditional".
#'
#' @return A vector or matrix of samples from the estimated density.
#'
#' @examples
#' x <- rnorm(10)
#' rgaussiandens(x, n = 10, p = 3)
#'
#' @export
rgaussiandens <- function(x,
                          n = length(x),
                          p = 1,
                          seed = 123,
                          method = c("antithetic",
                                     "traditional")) {
  z <- try(stats::density(x, bw = "sj", kernel = "gaussian"),
           silent = TRUE)

  if (inherits(z, "try-error"))
    z <- stats::density(x, kernel = "gaussian")

  width <- z$bw # Kernel width

  method <- match.arg(method)

  rkernel <- function(n, seed) {
    set.seed(seed)
    if (!identical(method, "antithetic"))
    {
      return(stats::rnorm(n, sd = width))
    } else {
      half_n <- n %/% 2
      eps <- stats::rnorm(half_n, sd = width)
      if (2 * length(eps) < n)
      {
        return(c(eps,-eps, stats::rnorm(1, sd = width)))
      }
      return(sample(c(eps,-eps),
                    replace = FALSE))
    }
  }  # Kernel sampler

  if (p <= 1)
  {
    set.seed(seed)
    return(sample(x, n, replace = TRUE) + rkernel(n, seed))    # Here's the entire algorithm
  } else {
    return(sapply(1:p,
                  function(i) {
                    set.seed(seed + i - 1)
                    sample(x, n, replace = TRUE) + rkernel(n, seed + i - 1)
                  }))
  }
}

# Simulate using surrogate -----

#' Simulate using surrogate data
#'
#' @export
#'
rsurrogate <- function(x,
                       n = length(x),
                       p = 1,
                       seed = 123) {
  if (n > length(x))
  {
    stop("For surrogates, must have number of predictions < number of training observations")
  }
  if (p <= 1)
  {
    set.seed(seed)
    res <- tseries::surrogate(x, ns = p,
                              fft = TRUE)
    return(res[seq_len(n), ])
  } else {
    res <- sapply(1:p,
                  function(i) {
                    set.seed(seed + i - 1)
                    tseries::surrogate(x, ns = p,
                                       fft = TRUE)
                  })
    return(res[seq_len(n), ])
  }
}

# fitdistr ------

#' Simulate from parametric distribution
#'
#' @export
#'
rfitdistr <- function(x, n=length(x), p=1)
{
  mean_x <- mean(x)
  sd_x <- sd(x)
  scaled_x <- (x - mean_x)/sd_x
  simulate_function <- misc::fit_param_dist(as.numeric(scaled_x), 
                                            verbose = FALSE)
  res <- simulate_function(n*p)
  if (p <= 1)
  {
    return(res)
  } else {
    return(matrix(res, nrow = n, ncol = p))
  }
}

# GAN ------

#' Generate Synthetic Data using GANs
#'
#' Creates synthetic data using Generative Adversarial Networks with predefined
#' architectures optimized for different data types.
#'
#' @param x Matrix or data.frame. Input data to learn distribution from.
#' @param n Integer. Number of synthetic samples (rows) to generate.
#' @param p Integer. Number of features (columns) to generate. Currently must match input dimension.
#' @param type_input Character. Type of data distribution: "unimodal", "mixture", or "financial".
#'
#' @return A matrix of synthetic data with n rows and p columns.
#'
#' @examples
#' \dontrun{
#' # Unimodal normal distribution
#' real_data <- matrix(rnorm(1000, mean = 2, sd = 4))
#' synthetic <- rgan(real_data, n = 500, p = 1, type_input = "unimodal")
#' 
#' # Mixture distribution
#' mixture <- rbinom(1000, 1, 0.5)
#' real_mixture <- matrix(mixture * rnorm(1000, 2, 1) + (1-mixture) * rnorm(1000, 8, 2))
#' synthetic <- rgan(real_mixture, n = 500, p = 1, type_input = "mixture")
#' 
#' # Financial data
#' synthetic <- rgan(EuStockMarkets[1:500,1], n = 500, p = 1, type_input = "financial")
#' }
#'
#' @export
rgan <- function(x, n, p = 1, type_input = c("unimodal", "mixture", "financial")) {
  # Validate inputs
  type_input <- match.arg(type_input)
  
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  
  if (ncol(x) != p) {
    stop("Input data dimension (ncol(x)) must match p parameter. Currently p must equal ncol(x).")
  }
  
  if (n <= 0) {
    stop("n must be positive")
  }
  # Define architectures based on type
  architectures <- switch(type_input,
                          "unimodal" = list(
                            generator = function(latent_dim = 1) {
                              keras_model_sequential(input_shape = latent_dim) |> 
                                layer_dense(units = 1, activation = "linear")
                            },
                            discriminator = function(dim = 1) {
                              keras_model_sequential(input_shape = dim) |> 
                                layer_dense(units = 8, activation = "relu") |>
                                layer_dense(units = 1, activation = "sigmoid")
                            },
                            n_iter = 5,
                            epochs_per_iter = 30
                          ),
                          "mixture" = list(
                            generator = function(latent_dim = 1) {
                              keras_model_sequential(input_shape = latent_dim) |> 
                                layer_dense(units = 8, activation = "tanh") |>
                                layer_dense(units = 1, activation = "linear")
                            },
                            discriminator = function(dim = 1) {
                              keras_model_sequential(input_shape = dim) |> 
                                layer_dense(units = 16, activation = "relu") |>
                                layer_dense(units = 8, activation = "relu") |>
                                layer_dense(units = 1, activation = "sigmoid")
                            },
                            n_iter = 8,
                            epochs_per_iter = 50
                          ),
                          "financial" = list(
                            generator = function(latent_dim = 1) {
                              keras_model_sequential(input_shape = latent_dim) |> 
                                layer_dense(units = 16, activation = "tanh") |>
                                layer_dense(units = 8, activation = "tanh") |>
                                layer_dense(units = 1, activation = "linear")
                            },
                            discriminator = function(dim = 1) {
                              keras_model_sequential(input_shape = dim) |> 
                                layer_dense(units = 32, activation = "relu") |>
                                layer_dense(units = 16, activation = "relu") |>
                                layer_dense(units = 8, activation = "relu") |>
                                layer_dense(units = 1, activation = "sigmoid")
                            },
                            n_iter = 10,
                            epochs_per_iter = 80
                          )
  )
  # Train the GAN
  result <- train_gan(
    train_dat = x,
    generator_fn = architectures$generator,
    discriminator_fn = architectures$discriminator,
    n_iter = architectures$n_iter,
    epochs_per_iter = architectures$epochs_per_iter,
    num_resamples = n,
    seed = 123L
  )
  # Return only the synthetic data
  return(result$resamples)
}

# Optional: Add a helper function to automatically detect the best type
#' @export
rgan_auto <- function(x, n, p = 1) {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  # Simple automatic type detection
  kurt <- moments::kurtosis(x)
  skew <- moments::skewness(x)
  
  if (kurt < 4 && abs(skew) < 1) {
    type <- "unimodal"
    message("Auto-detected: unimodal distribution")
  } else if (kurt > 6) {
    type <- "financial" 
    message("Auto-detected: financial-like distribution (high kurtosis)")
  } else {
    type <- "mixture"
    message("Auto-detected: mixture distribution")
  }

  rgan(x, n, p, type)
}