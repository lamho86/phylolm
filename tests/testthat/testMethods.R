test_that("Test confint", {
  set.seed(12891026)
  ## Tree
  ntips <- 10
  tree <- ape::rphylo(ntips, 0.1, 0)
  ## data
  y_data <- phylolm::rTrait(1, tree, model = "BM", parameters = list(sigma2 = 1))
  ## Condition
  cond <- sample(c(0, 1), ntips, replace = TRUE)
  names(cond) <- tree$tip.label
  y_data[cond == 1] <- y_data[cond == 1] + 5 ## strong effect

  ## Fit
  fit <- phylolm(y_data ~ cond, phy = tree)

  ## Confidence Interval
  level <- 0.95
  ci <- confint(fit, level = level)

  ## Confidence Interval - Manual
  quants <- c((1 - level) / 2, 1 - (1 - level) / 2)
  ci_manual <- coef(fit) + summary(fit)$coefficients[, 2] %o% qt(quants, fit$n - fit$d)

  ## Equal
  expect_equivalent(ci, ci_manual)

  ## One parameter at a time, changing level
  level <- 0.90
  ci <- confint(fit, parm = "cond", level = level)
  quants <- c((1 - level) / 2, 1 - (1 - level) / 2)
  ci_manual <- coef(fit)[2] + summary(fit)$coefficients[2, 2] * qt(quants, fit$n - fit$d)
  expect_equivalent(ci, ci_manual)
})
