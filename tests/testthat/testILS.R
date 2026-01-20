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

test_that("GC fit with replicates", {
  ## Tree
  set.seed(1289)
  ntips <- 5
  tree <- ape::rphylo(ntips, 0.1, 0)

  ## Replicates
  reps <- sample(1:3, ntips, replace = TRUE)
  rep_ids <- make.unique(rep(tree$tip.label, times = reps), sep = "_")
  ## match species to samples
  data <- data.frame(species = rep(tree$tip.label, times = reps),
                     sample = rep_ids)
  ## add replicates on the tree
  tree_rep <- add.replicates(tree, data)

  ## Simulate data
  traits <- rTrait(2, tree, model = "OU", parameters = list(alpha = 0.1))
  traits_rep <- matrix(rep(traits, times = rep(reps, 2)), nrow = sum(reps))
  traits_rep <- matrix(rnorm(prod(dim(traits_rep)), traits_rep, 0.5), nrow = sum(reps))
  rownames(traits_rep) <- rep_ids
  colnames(traits_rep) <- paste0("trait", 1:2)
  traits_rep <- as.data.frame(traits_rep)
  traits_rep$trait2 <- 2 * traits_rep$trait1 + traits_rep$trait2

  ## Fit - BM errors
  fit <- phylolm(trait1 ~ 1, data = traits_rep, phy = tree_rep, model = "BM", measurement_error = TRUE, REML = TRUE)
  ## comparison with julia
  expect_equal(fit$logLik, -13.7119218037, tolerance = 1e-5)
  expect_equal(fit$sigma2, 0.447819, tolerance = 1e-5)
  expect_equal(fit$sigma2_error, 0.204118, tolerance = 1e-5)

  ## Fit - GC - lambda fixed
  fit <- phylolm(trait1 ~ 1, data = traits_rep, phy = tree_rep, model = "GC", REML = TRUE,
                 upper.bound = list(lambda_GC = 1), lower.bound = list(lambda_GC = 1), starting.value = list(lambda_GC = 1))
  ## comparison with julia
  expect_equal(fit$logLik, -14.0489140317, tolerance = 1e-5)
  expect_equal(fit$sigma2, 0.304069, tolerance = 1e-5)

  fit <- phylolm(trait1 ~ 1, data = traits_rep, phy = tree_rep, model = "GC", REML = TRUE,
                 upper.bound = list(lambda_GC = 0), lower.bound = list(lambda_GC = 0), starting.value = list(lambda_GC = 0))
  ## comparison with julia
  expect_equal(fit$logLik, -14.0316023914, tolerance = 1e-5)
  expect_equal(fit$sigma2, 0.305006, tolerance = 1e-5)

  fit <- phylolm(trait1 ~ trait2, data = traits_rep, phy = tree_rep, model = "GC", REML = TRUE,
                 upper.bound = list(lambda_GC = 1), lower.bound = list(lambda_GC = 1), starting.value = list(lambda_GC = 1))
  ## comparison with julia
  expect_equal(fit$logLik, -9.3567699602, tolerance = 1e-5)
  expect_equal(fit$sigma2, 0.0855647, tolerance = 1e-5)

  ## julia code
  # tree = "((t4:6.964056834,(t5:2.862271346,t2:2.862271346):4.101785488):12.14330626,(t3:3.331896183,t1:3.331896183):15.77546691);"
  # tree = readnewick(tree)
  #
  # df = DataFrame(
  #   tipnames  = Vector{String}(["t1","t1","t2","t2","t3","t3","t4","t4","t4","t5"]),
  #   trait1    = Vector{Float64}([0.3606089,1.4910356,1.2622511,1.7206378,1.8314243,1.3526612,-3.5467733,-3.1195793,-2.9613277,0.6563839]),
  #   trait2    = Vector{Float64}([3.859196,6.197472,2.847903,3.205740,3.903874,3.928627,-6.149997,-6.252434,-5.375226,4.034418]))
  #
  #
  # phylolm(@formula(trait1 ~ 1), df, tree; model = "BM", withinspecies_var = true)
  # phylolm(@formula(trait1 ~ 1), df, tree; model = "gaussiancoalescent", paramlist=Dict(:lambda => (start=1, fixed=true)))
  # phylolm(@formula(trait1 ~ trait2), df, tree; model = "gaussiancoalescent", paramlist=Dict(:lambda => (start=1, fixed=true)))

})

test_that("GC fit with replicates and measurement error", {
  ## Tree
  set.seed(1289)
  ntips <- 100
  tree <- ape::rphylo(ntips, 0.1, 0)

  ## Replicates
  reps <- sample(3:5, ntips, replace = TRUE)
  rep_ids <- make.unique(rep(tree$tip.label, times = reps), sep = "_")
  ## match species to samples
  data <- data.frame(species = rep(tree$tip.label, times = reps),
                     sample = rep_ids)
  ## add replicates on the tree
  tree_rep <- add.replicates(tree, data)

  ## Simulate data
  traits <- rTrait(2, tree, model = "OU", parameters = list(alpha = 0.1))
  traits_rep <- matrix(rep(traits, times = rep(reps, 2)), nrow = sum(reps))
  traits_rep <- matrix(rnorm(prod(dim(traits_rep)), traits_rep, 2), nrow = sum(reps))
  rownames(traits_rep) <- rep_ids
  colnames(traits_rep) <- paste0("trait", 1:2)
  traits_rep <- as.data.frame(traits_rep)

  ## Fit - BM errors
  fit <- phylolm(trait1 ~ 1, data = traits_rep, phy = tree_rep, model = "GC", measurement_error = TRUE, REML = TRUE)
  expect_equal(fit$sigma2_error, 3.592375, tolerance = 1e-5)

})

