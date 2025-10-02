# Enable JIT compilation for performance
compiler::enableJIT(3)

# Define GAN class once (reused for both examples)
gan <- keras3::new_model_class(
  classname = "GAN",
  initialize = function(discriminator, generator, latent_dim) {
    super$initialize()
    self$discriminator <- discriminator
    self$generator <- generator
    self$latent_dim <- latent_dim
    self$d_loss_metric <- metric_mean(name = "d_loss")
    self$g_loss_metric <- metric_mean(name = "g_loss")
  },
  compile = function(d_optimizer, g_optimizer, loss_fn) {
    super$compile()
    self$d_optimizer <- d_optimizer
    self$g_optimizer <- g_optimizer
    self$loss_fn <- loss_fn
  },
  metrics = keras3::mark_active(function() {
    list(self$d_loss_metric, self$g_loss_metric)
  }),
  train_step = function(real_data) {
    batch_size <- tf$shape(real_data)[1]
    random_latent_vectors <- tf$random$normal(shape = c(batch_size, self$latent_dim))
    generated_data <- self$generator(random_latent_vectors)
    combined_data <- tf$concat(list(generated_data, real_data), axis = 0L)
    labels <- tf$concat(list(tf$ones(tuple(batch_size, 1L)),
                             tf$zeros(tuple(batch_size, 1L))), axis = 0L)
    labels <- labels + tf$random$uniform(tf$shape(labels), maxval = 0.01)
    
    with(tf$GradientTape() %as% tape, {
      predictions <- self$discriminator(combined_data)
      d_loss <- self$loss_fn(labels, predictions)
    })
    
    grads <- tape$gradient(d_loss, self$discriminator$trainable_weights)
    self$d_optimizer$apply_gradients(
      zip_lists(grads, self$discriminator$trainable_weights)
    )
    
    random_latent_vectors <- tf$random$normal(shape = c(batch_size, self$latent_dim))
    misleading_labels <- tf$zeros(tuple(batch_size, 1L))
    
    with(tf$GradientTape() %as% tape, {
      predictions <- random_latent_vectors |>
        self$generator() |>
        self$discriminator()
      g_loss <- self$loss_fn(misleading_labels, predictions)
    })
    
    grads <- tape$gradient(g_loss, self$generator$trainable_weights)
    self$g_optimizer$apply_gradients(
      zip_lists(grads, self$generator$trainable_weights)
    )
    
    self$d_loss_metric$update_state(d_loss)
    self$g_loss_metric$update_state(g_loss)
    list(d_loss = self$d_loss_metric$result(),
         g_loss = self$g_loss_metric$result())
  }
)

# Reusable training function

#' Train a Simple GAN Model
#'
#' Trains a simple Generative Adversarial Network (GAN) using user-supplied generator and discriminator functions.
#'
#' @param train_dat Matrix or data.frame. Training data for the GAN.
#' @param generator_fn Function. Returns a generator model given latent dimension.
#' @param discriminator_fn Function. Returns a discriminator model given input dimension.
#' @param n_iter Integer. Number of training iterations (default: 5).
#' @param epochs_per_iter Integer. Number of epochs per iteration (default: 30).
#' @param num_resamples Integer. Number of synthetic samples to generate after training (default: 500).
#' @param seed Integer. Random seed for reproducibility (default: 123).
#'
#' @return A list containing the trained GAN model (`model`), a matrix of synthetic resamples (`resamples`), and elapsed training time in seconds (`time`).
#'
#' @details
#' Requires both 'tensorflow' and 'keras3' packages to be installed. The generator and discriminator functions must return valid Keras models.
#'
#' @examples
#' # Example usage (requires keras3 and tensorflow)
#' # train_gan(train_dat = matrix(rnorm(100)), generator_fn = my_gen, discriminator_fn = my_disc)
#'
#' @export
train_gan <- function(train_dat, generator_fn, discriminator_fn, 
                      n_iter = 5, epochs_per_iter = 30, 
                      num_resamples = 500, seed=123L) {
  
  is_tf_available <- misc::is_package_available("tensorflow")
  is_keras3_available <- misc::is_package_available("keras3")
  if (!is_tf_available | !is_keras3_available)
    stop("both 'tensorflow' and 'keras3' must be installed")
  set.seed(seed)
  set_random_seed(seed)
  tf$random$set_seed(seed)
  
  latent_dim <- as.integer(1)
  nsyn <- nrow(train_dat)
  # Build and compile model
  mod <- gan(
    discriminator = discriminator_fn(dim = ncol(train_dat)), 
    generator = generator_fn(latent_dim = latent_dim), 
    latent_dim = latent_dim
  )
  mod$compile(
    d_optimizer = optimizer_adam(beta_1 = 0.5),
    g_optimizer = optimizer_adam(beta_1 = 0.5),
    loss_fn = loss_binary_crossentropy()
  )
  # Training loop
  start <- proc.time()[3]
  pb <- utils::txtProgressBar(max = n_iter, style = 3L)
  
  for (i in seq_len(n_iter)) {
    mod |> fit(train_dat, epochs = epochs_per_iter, batch_size = 32, verbose = 0)
    newdat <- mod$generator(tf$random$normal(shape = c(nsyn, latent_dim)))
    utils::setTxtProgressBar(pb, i)
  }
  elapsed <- proc.time()[3] - start
  close(pb)
  cat(sprintf("\nElapsed %.2f seconds\n", elapsed))
  # Generate resamples
  latent_vectors <- matrix(rnorm(n = num_resamples * latent_dim), 
                           nrow = num_resamples, ncol = latent_dim)
  resamples_matrix <- as.matrix(mod$generator(latent_vectors))
  
  return(list(model = mod, resamples = resamples_matrix, time = elapsed))
}