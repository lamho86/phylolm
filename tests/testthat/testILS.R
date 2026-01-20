test_that("GC fit", {
  #' @title Variance matrix for a BM with ILS on a tree
  #'
  #' @description
  #' Computes the un-scaled variance matrix of a trait evolving as a BM with ILS on
  #' a phylogenetic tree.
  #' The ancestral population is assumed to have a given distribution
  #' with mean `m_0` and variance `sig2_0`.
  #' If $N$ is the population size, for a BM process with variance `2*N*sig2_BM`,
  #' the variance matrix of the trait at the tips of the tree is given by
  #' `2*N*sig2_BM * vcv_ils(tree, lambda)`, where
  #' `lambda = sig2_0 / (2*N*sig2_BM)`.
  #' Take `lambda=0` for a fixed ancestral population distribution (`sig2_0=0`),
  #' and `lambda=1` for an ancestral distribution at equilibrium (`sig2_0=2*N*sig2_BM`).
  #'
  #' @param tree a phylogenetic tree of class ape::phylo.
  #' Branch lengths must be in coalescent units.
  #' @param lambda the ratio between the ancestral population variance and the
  #' (scaled) BM variance.
  #'
  #' @return The un-scaled variance matrix of the tree.
  #'
  vcv_ils <- function(tree, lambda) {
    s_all <- vcv(tree)
    lambda_all <- 1 + (lambda - 1) * (1 - shared_time_ratio(s_all))
    q_s_diag <- proba_coal(diag(s_all))
    sig_diag <- q_s_diag + lambda * (1 - q_s_diag)
    return(lambda_all * s_all + diag(sig_diag))
  }

  #' @title Coalescence probability (q)
  #'
  #' @description
  #' Probability that two individual coalesce after a time `t` in coalescent unit.
  #'
  #' @param t time in coalescent unit
  #'
  #' @return probability of coalescence
  #'
  proba_coal <- function(t) {
    return(1 - exp(-t))
  }

  #' @title Expected shared time (r)
  #'
  #' @description
  #' Two individuals from a population of size $N$ share a common ancestor
  #' for an average of `l*shared_time_ratio(l/(2*N))` generations within the last `l`
  #' generations.
  #'
  #' @param t time in coalescent unit
  #'
  #' @return The proportion of time shared over the last `t` coalescent units.
  #'
  shared_time_ratio <- function(t) {
    res <- 1 - proba_coal(t) / t
    res[is.na(res)] <- 0 # case when t = 0
    return(res)
  }

  set.seed(1289)
  tree <- rphylo(15, 1, 0.1, fossils = TRUE)
  # direct computation
  varils <- vcv_ils(tree, 0.2)
  # using branch length transform
  treeils <- transf.branch.lengths(tree, model = "GC", parameters = list("lambda_GC" = 0.2))
  varils2 <- vcv(treeils$tree)
  # equal ?
  expect_equal(varils, varils2)
})

test_that("GC fit", {
  set.seed(1289)
  ## simulate tree
  tree <- rphylo(15, 1, 0.1, fossils = TRUE)
  ## simulate trait
  b0 <- 0; b1 <- 1;
  x <- rTrait(n = 1, phy = tree, model = "BM",
              parameters = list(ancestral.state = 0, sigma2 = 10))
  y <- b0 + b1 * x +
    rTrait(n = 1, phy = tree, model = "lambda",
           parameters = list(ancestral.state = 0, sigma2 = 1, lambda = 0.5))
  # adding measurement errors and bootstrap
  z <- y + rnorm(20, 0, 1)
  dat <- data.frame(trait = z, pred = x)

  ## Fit - Fixed values
  fit <- phylolm(trait ~ pred, data = dat, phy = tree, model = "GC")
  expect_equal(fit$sigma2, 1.141, tolerance = 1e-3)
  expect_equal(fit$optpar, 5.546e-08, tolerance = 1e-3)

  fit <- phylolm(trait ~ pred, data = dat, phy = tree, model = "GC", measurement_error = TRUE)
  expect_equal(fit$sigma2, 0.362, tolerance = 1e-3)
  expect_equal(fit$optpar, 1.920e-06, tolerance = 1e-3)
})
