# functions for running stochastic character mapping.

# =============================================================================

#' Prune tree for running `fitMk()`
#'
#' @param .phy tree
#' @param .dat tipdata as a named vector where names are tip labels in tree
#'
#' @return
#' @export
#'
#' @examples
prune_tree_for_fitmk <- function(.phy, .dat) {
  stopifnot(all(names(.dat) %in% .phy$tip.label))

  .phy_pruned <- drop.tip(.phy, which(!(.phy$tip.label %in% names(.dat))))
  .phy_pruned <- ladderize(.phy_pruned)

  return(.phy_pruned)
}

# =============================================================================
#' Multistart fitMk
#'
#' Start fitMk from multiple starting points to avoid local optimum. Returns the
#'  best of .nstart fits from randomized starting values, which is more likely
#'  to be the global optimum than a single-start fit."
#'
#' @param .phy tree
#' @param .dat tipdata
#' @param .model model to fit
#' @param .nstart number of random starts
#'
#' @returns
#' @export
#'
#' @examples
fit_mk_multistart <- function(.phy, .dat, .model, .nstart = 12) {
  fits <- map(seq_len(.nstart), function(i) {
    fitMk(.phy, .dat, model = .model, pi = "estimated", rand_start = TRUE)
  })
  ll <- map_dbl(fits, function(f) f$logLik)

  best <- fits[[which.max(ll)]]
  attr(best, "multistart_logLik") <- ll
  best
}


# =============================================================================

#' Fit models
#'
#' Currently only does ER, SYM, and ARD models. Need to modify in future to take
#'    any set of provided models.
#'
#' @param .phy tree
#' @param .dat tip data as a named vector where names are tip labels in tree
#'
#' @return results of `anova()`
#' @export
#'
#' @examples
#'
fit_models <- function(.phy, .dat) {
  .dat_sorted <- .dat[.phy$tip.label]
  stopifnot(identical(names(.dat_sorted), .phy$tip.label))

  # fit models
  message("    fitting models...")
  fit_er <- fit_mk_multistart(.phy, .dat_sorted, .model = "ER")
  fit_sym <- fit_mk_multistart(.phy, .dat_sorted, .model = "SYM")
  fit_ard <- fit_mk_multistart(.phy, .dat_sorted, .model = "ARD")

  # model comparison
  message("    running model comparison...")
  fit_aov <- anova(fit_er, fit_sym, fit_ard)

  return(fit_aov)
}

# =============================================================================

#' Fit single model
#'
#' @param .phy tree
#' @param .dat tip data as a named vector where names are tip labels in tree
#' @param .model either one of ER, SYM, ARD, or a custom model as a matrix
#'
#' @return results of `fitMK()`
#' @export
#'
#' @examples
#'
fit_single_model <- function(.phy, .dat, .model) {
  .dat_sorted <- .dat[.phy$tip.label]
  stopifnot(identical(names(.dat_sorted), .phy$tip.label))

  message("    fitting model...")
  fit <- fit_mk_multistart(.phy, .dat_sorted, .model = .model)

  return(fit)
}

# =============================================================================

#' Run simmap
#'
#' @param .models_aov outout of `fit_models()`
#' @param .nsim number of simulations
#'
#' @return simmap object
#' @export
#'
#' @examples
run_simmap <- function(.models_aov, .nsim = 1000) {
  simmap(.models_aov, nsim = .nsim)
}

# =============================================================================

#' Wrapper to run simmap and get ancestral character estimates
#'
#' @param .models_aov
#' @param .nsim
#'
#' @return
#' @export
#'
#' @examples
run_simmap_and_get_ace <- function(.models_aov, .nsim = 1000) {
  .simmap <- run_simmap(.models_aov = .models_aov, .nsim = .nsim)
  .simmap_summary <- summary(.simmap)
  .ace <- .simmap_summary$ace
  return(.ace)
}

###############################################################################
