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
  fit_er  <- fitMk(.phy, .dat_sorted, model = "ER",  pi = "estimated")
  fit_sym <- fitMk(.phy, .dat_sorted, model = "SYM", pi = "estimated")
  fit_ard <- fitMk(.phy, .dat_sorted, model = "ARD", pi = "estimated")
  
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
  fit <- fitMk(.phy, .dat_sorted, model = .model,  pi = "estimated")
  
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
