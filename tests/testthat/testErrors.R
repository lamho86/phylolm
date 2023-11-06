test_that("Multivariate Response throws an error", {
  set.seed(12891026)
  ## Tree
  ntips <- 10
  tree <- ape::rphylo(ntips, 0.1, 0)
  ## multivariate data
  ntraits <- 2
  y_data <- phylolm::rTrait(ntraits, tree, model = "BM", parameters = list(sigma2 = 1))
  ## Condition
  cond <- sample(c(0, 1), ntips, replace = TRUE)
  names(cond) <- tree$tip.label
  y_data[cond == 1, ] <- y_data[cond == 1, ] + 5 ## strong effect
  ## Data frame
  all_dat <- data.frame(y = y_data, cond = cond)

  ## Univariate fits (various specifications, all equal)
  fit1 <- phylolm(y_data[, 1] ~ cond, phy = tree)
  fit2 <- phylolm(y_data[, 1, drop = FALSE] ~ cond, phy = tree)
  fit3 <- phylolm(y.1 ~ cond, phy = tree, data = all_dat)

  expect_equal(fit1[!(names(fit1) %in% c("call", "formula", "model.frame"))],
               fit2[!(names(fit1) %in% c("call", "formula", "model.frame"))])
  expect_equal(fit1[!(names(fit1) %in% c("call", "formula", "model.frame"))],
               fit3[!(names(fit1) %in% c("call", "formula", "model.frame"))])

  ## Multivariate fits throw errors
  expect_error(phylolm(y_data ~ cond, phy = tree),
               "'phylolm' can only handle a simple \\(univariate\\) response vector y.")

  expect_error(phylolm(cbind(y.1, y.2) ~ cond, phy = tree, data = all_dat),
               "'phylolm' can only handle a simple \\(univariate\\) response vector y.")
})
