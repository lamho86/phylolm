test_that("REML estimates", {
  skip_if_not_installed("nlme")

  ## Taken from example
  set.seed(123456)
  tre = rphylo(60, 0.1, 0.01)
  tre_h <- max(vcv(tre))
  taxa = sort(tre$tip.label)
  b0 = 0; b1 = 1; b2 = -2;
  x <- rTrait(n = 1, phy = tre, model = "BM",
              parameters = list(ancestral.state = 0, sigma2 = 10))
  f <- sample(c(0, 1), length(taxa), replace = TRUE)
  names(f) <- names(x)
  y <- b0 + b1*x + b2 * f +
    rTrait(n = 1, phy = tre, model = "lambda",
           parameters = list(
             ancestral.state = 0, sigma2 = 1, lambda = 0.5))
  # adding measurement errors
  z <- y + rnorm(length(taxa), 0, 1)
  dat = data.frame(trait = z[taxa], pred = x[taxa], fac = as.factor(f[taxa]), species = taxa)

  ## Fit - BM - ML
  fit = phylolm(trait ~ pred + fac, data = dat, phy = tre, model = "BM")
  fit_gls <- nlme::gls(trait ~ pred + fac, data = dat, correlation = corBrownian(1, tre, form = ~species), method = "ML")

  expect_equal(as.numeric(logLik(fit)), fit_gls$logLik)
  expect_equal(fit$coefficients, fit_gls$coefficients)
  expect_equal(fit$vcov, fit_gls$varBeta)
  expect_equal(fit$n, fit_gls$dims$N)
  expect_equal(fit$d, fit_gls$dims$p)
  expect_equal(summary(fit)$coefficients[, "t.value"], summary(fit_gls)$tTable[, "t-value"])
  expect_equal(fit$fitted[taxa], fit_gls$fitted[taxa])
  expect_equal(fit$residuals[taxa], fit_gls$residuals[taxa])
  expect_equal(AIC(fit), AIC(fit_gls))
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(fit_gls)))
  expect_equal(fit$sigma2, fit_gls$sigma^2 / tre_h)

  ## Fit - BM - REML
  fit = phylolm(trait ~ pred, data = dat, phy = tre, model = "BM", REML = TRUE)
  fit_gls <- nlme::gls(trait ~ pred, data = dat, correlation = corBrownian(1, tre, form = ~species), method = "REML")

  expect_equal(as.numeric(logLik(fit)), fit_gls$logLik)
  expect_equal(fit$coefficients, fit_gls$coefficients)
  expect_equal(fit$vcov, fit_gls$varBeta)
  expect_equal(fit$n, fit_gls$dims$N)
  expect_equal(fit$d, fit_gls$dims$p)
  expect_equal(summary(fit)$coefficients[, "t.value"], summary(fit_gls)$tTable[, "t-value"])
  expect_equal(fit$fitted[taxa], fit_gls$fitted[taxa])
  expect_equal(fit$residuals[taxa], fit_gls$residuals[taxa])
  expect_equal(AIC(fit), AIC(fit_gls))
  expect_equal(as.numeric(logLik(fit)), logLik(fit_gls)[1])
  expect_equal(fit$sigma2, fit_gls$sigma^2 / tre_h)

  ## Fit - lambda - ML
  fit = phylolm(trait ~ pred + fac, data = dat, phy = tre, model = "lambda")
  fit_gls <- nlme::gls(trait ~ pred + fac, data = dat, correlation = corPagel(1, tre, form = ~species), method = "ML")

  expect_equal(as.numeric(logLik(fit)), fit_gls$logLik)
  expect_equal(fit$coefficients, fit_gls$coefficients, tolerance = 1e-5)
  expect_equal(fit$vcov, fit_gls$varBeta, tolerance = 1e-5)
  expect_equal(fit$n, fit_gls$dims$N)
  expect_equal(fit$d, fit_gls$dims$p)
  expect_equal(summary(fit)$coefficients[, "t.value"], summary(fit_gls)$tTable[, "t-value"], tolerance = 1e-5)
  expect_equal(fit$fitted[taxa], fit_gls$fitted[taxa], tolerance = 1e-5)
  expect_equal(fit$residuals[taxa], fit_gls$residuals[taxa], tolerance = 1e-5)
  expect_equal(AIC(fit), AIC(fit_gls))
  expect_equal(as.numeric(logLik(fit)), logLik(fit_gls)[1])
  expect_equal(fit$optpar, fit_gls$modelStruct$corStruct[1], tolerance = 1e-5)
  expect_equal(fit$sigma2, fit_gls$sigma^2 / tre_h, tolerance = 1e-5)

  ## Fit - lambda - REML
  fit = phylolm(trait ~ pred, data = dat, phy = tre, model = "lambda", REML = TRUE)
  fit_gls <- nlme::gls(trait ~ pred, data = dat, correlation = corPagel(1, tre, form = ~species), method = "REML")

  expect_equal(fit$logLik, fit_gls$logLik)
  expect_equal(fit$coefficients, fit_gls$coefficients, tolerance = 1e-6)
  expect_equal(fit$vcov, fit_gls$varBeta, tolerance = 1e-6)
  expect_equal(fit$n, fit_gls$dims$N)
  expect_equal(fit$d, fit_gls$dims$p)
  expect_equal(summary(fit)$coefficients[, "t.value"], summary(fit_gls)$tTable[, "t-value"], tolerance = 1e-6)
  expect_equal(fit$fitted[taxa], fit_gls$fitted[taxa], tolerance = 1e-6)
  expect_equal(fit$residuals[taxa], fit_gls$residuals[taxa], tolerance = 1e-6)
  expect_equal(AIC(fit), AIC(fit_gls))
  expect_equal(as.numeric(logLik(fit)), logLik(fit_gls)[1])
  expect_equal(fit$optpar, fit_gls$modelStruct$corStruct[1], tolerance = 1e-6)
  expect_equal(fit$sigma2, fit_gls$sigma^2 / tre_h, tolerance = 1e-6)

  ## Fit - OU - ML
  fit = phylolm(trait ~ pred + fac, data = dat, phy = tre, model = "OUrandomRoot")
  fit_gls <- nlme::gls(trait ~ pred + fac, data = dat, correlation = corMartins(0.1, tre, form = ~species), method = "ML")

  expect_equal(fit$logLik, fit_gls$logLik)
  expect_equal(fit$coefficients, fit_gls$coefficients, tolerance = 1e-6)
  expect_equal(fit$vcov, fit_gls$varBeta, tolerance = 1e-6)
  expect_equal(fit$n, fit_gls$dims$N)
  expect_equal(fit$d, fit_gls$dims$p)
  expect_equal(summary(fit)$coefficients[, "t.value"], summary(fit_gls)$tTable[, "t-value"], tolerance = 1e-5)
  expect_equal(fit$fitted[taxa], fit_gls$fitted[taxa], tolerance = 1e-6)
  expect_equal(fit$residuals[taxa], fit_gls$residuals[taxa], tolerance = 1e-6)
  expect_equal(AIC(fit), AIC(fit_gls))
  expect_equal(as.numeric(logLik(fit)), logLik(fit_gls)[1])
  expect_equal(fit$optpar, fit_gls$modelStruct$corStruct[1], tolerance = 1e-6)
  expect_equal(fit$sigma2 / 2 / fit$optpar, fit_gls$sigma^2, tolerance = 1e-6)

  ## Fit - OU - REML
  fit = phylolm(trait ~ pred, data = dat, phy = tre, model = "OUrandomRoot", REML = TRUE)
  fit_gls <- nlme::gls(trait ~ pred, data = dat, correlation = corMartins(0.1, tre, form = ~species), method = "REML")

  expect_equal(fit$logLik, fit_gls$logLik)
  expect_equal(fit$coefficients, fit_gls$coefficients, tolerance = 1e-6)
  expect_equal(fit$vcov, fit_gls$varBeta, tolerance = 1e-6)
  expect_equal(fit$n, fit_gls$dims$N)
  expect_equal(fit$d, fit_gls$dims$p)
  expect_equal(summary(fit)$coefficients[, "t.value"], summary(fit_gls)$tTable[, "t-value"], tolerance = 1e-5)
  expect_equal(fit$fitted[taxa], fit_gls$fitted[taxa], tolerance = 1e-6)
  expect_equal(fit$residuals[taxa], fit_gls$residuals[taxa], tolerance = 1e-6)
  expect_equal(AIC(fit), AIC(fit_gls))
  expect_equal(as.numeric(logLik(fit)), logLik(fit_gls)[1])
  expect_equal(fit$optpar, fit_gls$modelStruct$corStruct[1], tolerance = 1e-6)
  expect_equal(fit$sigma2 / 2 / fit$optpar, fit_gls$sigma^2, tolerance = 1e-6)

})

test_that("REML estimates - non ultrametric tree", {
  # skip_if_not_installed("mvMORPH")

  ## Taken from example
  set.seed(123456)
  tre = rphylo(100, 0.1, 0.05, fossils = TRUE)
  tre_h <- max(vcv(tre))
  taxa = sort(tre$tip.label)
  b0 = 0; b1 = 1; b2 = -2;
  x <- rTrait(n = 1, phy = tre, model = "BM",
              parameters = list(ancestral.state = 0, sigma2 = 10))
  f <- sample(c(0, 1), length(taxa), replace = TRUE)
  names(f) <- names(x)
  y <- b0 + b1*x + b2 * f +
    rTrait(n = 1, phy = tre, model = "lambda",
           parameters = list(
             ancestral.state = 0, sigma2 = 1, lambda = 0.5))
  # adding measurement errors
  z <- y + rnorm(length(taxa), 0, 1)
  dat = data.frame(trait = z[taxa], pred = x[taxa], fac = as.factor(f[taxa]), species = taxa)
  # Dummy variable for mvMORPH
  datmv <- list(trait = cbind(z[taxa], rnorm(length(taxa))),
                pred = x[taxa], fac = as.factor(f[taxa]), species = taxa)

  ## Fit - BM - ML
  fit = phylolm(trait ~ pred + fac, data = dat, phy = tre, model = "BM")
  # fit_mvMORPH <- mvMORPH::mvgls(trait ~ pred + fac, data = datmv,
  #                               tree = tre, model = "BM",
  #                               method = "LL", REML = FALSE,
  #                               param = list(constraint = "diagonal"))

  ## check against values obtained with mvMORPH
  # fit_mvMORPH$coefficients[, 1]
  expect_equal(fit$coefficients, c(-0.6197818, 0.8807205, -0.1857472), tolerance = 1e-3, check.attributes = FALSE)
  # fit_mvMORPH$fitted[taxa, 1][c(1, 4, 39)]
  expect_equal(fit$fitted[taxa][c(1, 4, 39)], c(18.348420, -3.333741, 21.302602), tolerance = 1e-3, check.attributes = FALSE)
  # fit_mvMORPH$residuals[taxa, 1][c(24, 32, 49)]
  expect_equal(fit$residuals[taxa][c(24, 32, 49)], c(-8.862541, -18.115390, 4.652933), tolerance = 1e-3, check.attributes = FALSE)
  # fit_mvMORPH$sigma$S[1,1]
  expect_equal(fit$sigma2, 11.60188, tolerance = 1e-3, check.attributes = FALSE)

  ## Fit - BM - REML
  fit = phylolm(trait ~ pred, data = dat, phy = tre, model = "BM", REML = TRUE)
  # fit_mvMORPH <- mvMORPH::mvgls(trait ~ pred, data = datmv,
  #                               tree = tre, model = "BM",
  #                               method = "LL", REML = TRUE,
  #                               param = list(constraint = "diagonal"))

  ## check against values obtained with mvMORPH
  # fit_mvMORPH$coefficients[, 1]
  expect_equal(fit$coefficients, c(-0.6753924, 0.8805762), tolerance = 1e-5, check.attributes = FALSE)
  # fit_mvMORPH$fitted[taxa, 1][c(1, 4, 39)]
  expect_equal(fit$fitted[taxa][c(1, 4, 39)], c(18.47542, -3.20319, 21.24340 ), tolerance = 1e-5, check.attributes = FALSE)
  # fit_mvMORPH$residuals[taxa, 1][c(24, 32, 49)]
  expect_equal(fit$residuals[taxa][c(24, 32, 49)], c(-8.807384, -18.248363, 4.709006 ), tolerance = 1e-5, check.attributes = FALSE)
  # fit_mvMORPH$sigma$S[1,1]
  expect_equal(fit$sigma2, 11.74626, tolerance = 1e-5, check.attributes = FALSE)

  ## Fit - OU - ML
  # fit_mvMORPH <- mvMORPH::mvgls(trait ~ pred, data = datmv,
  #                               tree = tre, model = "OU",
  #                               method = "LL", REML = FALSE,
  #                               param = list(constraint = "diagonal"))
  # fit_mvMORPH$param
  fit <- phylolm(trait ~ pred, data = dat, phy = tre, model = "OUrandomRoot",
                 lower.bound = list(alpha = 0.5134518),
                 upper.bound = list(alpha = 0.5134518),
                 starting.value = list(alpha = 0.5134518))

  ## check against values obtained with mvMORPH
  # fit_mvMORPH$param
  expect_equal(fit$optpar, 0.5134518, tolerance = 1e-5, check.attributes = FALSE)
  # fit_mvMORPH$coefficients[, 1]
  expect_equal(fit$coefficients, c(-4.261479, 1.036382), tolerance = 1e-4, check.attributes = FALSE)
  # fit_mvMORPH$fitted[taxa, 1][c(1, 4, 39)]
  expect_equal(fit$fitted[taxa][c(1, 4, 39)], c(18.277810, -7.236537, 21.535548), tolerance = 1e-4, check.attributes = FALSE)
  # fit_mvMORPH$residuals[taxa, 1][c(24, 32, 49)]
  expect_equal(fit$residuals[taxa][c(24, 32, 49)], c(-4.731292, -11.599618, 7.795520 ), tolerance = 1e-4, check.attributes = FALSE)
  # fit_mvMORPH$sigma$S[1,1]
  expect_equal(fit$sigma2, 50.40184, tol = 1e-4)

  ## Fit - OU - REML
  # fit_mvMORPH <- mvMORPH::mvgls(trait ~ pred, data = datmv,
  #                               tree = tre, model = "OU",
  #                               method = "LL", REML = TRUE,
  #                               param = list(constraint = "diagonal"))
  # fit_mvMORPH$param
  fit <- phylolm(trait ~ pred, data = dat, phy = tre, model = "OUrandomRoot",
                 REML = TRUE,
                 lower.bound = list(alpha = 0.5134518),
                 upper.bound = list(alpha = 0.5134518),
                 starting.value = list(alpha = 0.5134518))

  ## check against values obtained with mvMORPH
  # fit_mvMORPH$param
  expect_equal(fit$optpar, 0.5134518, tolerance = 1e-5, check.attributes = FALSE)
  # fit_mvMORPH$coefficients[, 1]
  expect_equal(fit$coefficients, c(-4.261479, 1.036382), tolerance = 1e-4, check.attributes = FALSE)
  # fit_mvMORPH$fitted[taxa, 1][c(1, 4, 39)]
  expect_equal(fit$fitted[taxa][c(1, 4, 39)], c(18.277810, -7.236537, 21.535548), tolerance = 1e-4, check.attributes = FALSE)
  # fit_mvMORPH$residuals[taxa, 1][c(24, 32, 49)]
  expect_equal(fit$residuals[taxa][c(24, 32, 49)], c(-4.731292, -11.599618, 7.795520 ), tolerance = 1e-4, check.attributes = FALSE)
  # fit_mvMORPH$sigma$S[1,1]
  expect_equal(fit$sigma2, 51.01649, tol = 1e-4)

})
