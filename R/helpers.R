# Helper functions ############################################################

# Wrapper function for sampling from Stan models.
#
# When models are run and fits are stored the function retrieves the fit
# from storage. This will make it easier to work with the qmd document.
#
#
# model, a compiled Stan model object
# DATA, the data set used to fit the model
# fit, name of the fitted model
# iter, number of iterations
# iterwu, number of iterations in war-up
# chains, number of chains
# pchains, number of parrallel chains,
# refresh, progress refresh
fit_wrap <- function(
  model,
  DATA = d,
  fit = "fit.m1",
  iter = 2000,
  iterwu = 1000,
  chains = 4,
  pchains = 4,
  refresh = 1000,
  refit = refit.models
) {
  fit.filepath <- paste0("fits/paper1/", fit, ".RDS")

  # If refit is set to true in the document, all model fits will be refitted
  if (refit) file.remove(fit.filepath)

  # Sample from the model if a saved version does not exixts.
  if (!file.exists(fit.filepath)) {
    fittedmodel <- model$sample(
      data = DATA, # The list of variables
      iter_sampling = iter, # Number of iterations for sampling
      iter_warmup = iterwu, # Number of warm-up iterations
      chains = chains, # Number of Markov chains
      parallel_chains = pchains, # Number of parallel chains
      refresh = refresh
    ) # Number of iterations before each print
    # Save the model for easy access
    fittedmodel$save_object(file = fit.filepath) # <-- use this
  } else {
    fittedmodel <- readRDS(fit.filepath)
  }
  # Assign the model object to the environment.
  assign(fit, fittedmodel, envir = parent.frame())
}

#
