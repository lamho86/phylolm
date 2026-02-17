context("addReplicates")

test_that("addReplicates", {
  skip_if_not_installed("phytools")
  ## Tree
  set.seed(1289)
  ntips <- 10
  tree <- ape::rphylo(ntips, 0.1, 0)

  ## Simulate data
  reps <- sample(0:3, ntips, replace = TRUE)
  traits <- data.frame(g1 = rnorm(sum(reps), 1, 0.5),
                       species = unname(do.call(c, mapply(function(x, y) rep(x, each = y), tree$tip.label, reps)))
  )
  traits$id <- mapply(paste0, traits$species, do.call(c, sapply(reps, function(x) if (x > 0) paste0("_", 1:x))))

  ## Replicates
  expect_warning(tree_rep <- add.replicates(tree, traits, species_id = "species", sample_id = "id"),
                 "Species 't4', 't9' are in the tree but not in the data. They will be droped from the final tree.")
  expect_equal(length(tree_rep$tip.label), sum(reps))
  expect_equal(ape::is.ultrametric(tree_rep), TRUE)

  ## Extra data
  traits_extra <- rbind(traits, data.frame(g1 = rnorm(3, 3, 3),
                                           species = rep("r1", 3),
                                           id = paste0("r1_", 1:3)))
  expect_error(add.replicates(tree, traits_extra, species_id = "species", sample_id = "id"),
               "are in the data but not in the tree.")

  ## Wrong names - species
  traits <- data.frame(g1 = rnorm(sum(reps), 1, 0.5),
                       labels = unname(do.call(c, mapply(function(x, y) rep(x, each = y), tree$tip.label, reps)))
  )
  traits$id <- mapply(paste0, traits$labels, do.call(c, sapply(reps, function(x) if (x > 0) paste0("_", 1:x))))

  expect_error(add.replicates(tree, traits, species_id = "species", sample_id = "id"),
               "data frame should contain a column named species with species names for each sample")
  expect_warning(add.replicates(tree, traits, species_id = "labels", sample_id = "id"),
                 "Species 't4', 't9' are in the tree but not in the data. They will be droped from the final tree.")

  ## Wrong names - ids
  traits <- data.frame(g1 = rnorm(sum(reps), 1, 0.5),
                       species = unname(do.call(c, mapply(function(x, y) rep(x, each = y), tree$tip.label, reps)))
  )
  traits$unique_ids <- mapply(paste0, traits$species, do.call(c, sapply(reps, function(x) if (x > 0) paste0("_", 1:x))))

  expect_error(add.replicates(tree, traits, species_id = "species", sample_id = "id"),
               "data frame should contain a column named id with sample ids")
  expect_warning(add.replicates(tree, traits, species_id = "species", sample_id = "unique_ids"),
                 "Species 't4', 't9' are in the tree but not in the data. They will be droped from the final tree.")

  ## from a vector
  traits$unique_ids <- make.unique(traits$species, sep = "_")
  expect_warning(tree_1 <- add.replicates(tree, traits, species_id = "species", sample_id = "unique_ids"),
                 "Species 't4', 't9' are in the tree but not in the data. They will be droped from the final tree.")
  expect_warning(tree_2 <- add.replicates(tree, traits$species),
                 "Species 't4', 't9' are in the tree but not in the data. They will be droped from the final tree.")
  expect_true(isTRUE(all.equal(tree_1, tree_2)))

  ## Replicated ids
  traits$unique_ids[2] <- traits$unique_ids[3]

  expect_error(add.replicates(tree, traits, species_id = "species", sample_id = "unique_ids"),
               "should be unique identifiers of the samples")


})

test_that("phylolm - rep", {
  skip_if_not_installed("phytools")

  ## Tree
  set.seed(1289)
  ntips <- 20
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
  traits <- rTrait(3, tree, model = "OU", parameters = list(alpha = 0.1))
  traits_rep <- matrix(rep(traits, times = rep(reps, 3)), nrow = sum(reps))
  traits_rep <- matrix(rnorm(prod(dim(traits_rep)), traits_rep, 0.5), nrow = sum(reps))
  rownames(traits_rep) <- rep_ids
  colnames(traits_rep) <- paste0("trait", 1:3)
  traits_rep <- as.data.frame(traits_rep)

  ## fit trait 1
  fit1 <- phylolm(trait1 ~ 1, traits_rep, tree_rep, model = "OUfixedRoot", measurement_error = TRUE)

  ## Rphylopars
  # traits_rphylopars <- cbind(data$species, traits_rep)
  # colnames(traits_rphylopars)[1] <- "species"
  # fit_phylopars <- Rphylopars::phylopars(traits_rphylopars[, 1:2, drop = F], tree, model = "OU",
  #                                        REML = FALSE, phylo_correlated = FALSE, pheno_correlated = FALSE)
  # fit_phylopars$pars$phenocov
  # fit_phylopars$model$alpha

  expect_equal(fit1$sigma2_error, 0.2117206, tolerance = 1e-3) # value given by Rphylopars
  expect_equal(fit1$optpar, 0.09693093, tolerance = 1e-3) # value given by Rphylopars

  ## fit trait 1 vs 2 and 3
  fit123 <- phylolm(trait1 ~ trait2 + trait3, traits_rep, tree_rep, model = "OUfixedRoot", measurement_error = TRUE)
  expect_equal(fit123$sigma2_error, 0.18, tolerance = 1e-2)

})
