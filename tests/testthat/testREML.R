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
