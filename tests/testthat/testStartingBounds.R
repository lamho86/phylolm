test_that("Setting upper, lower and starting values", {
  ## Taken from example
  set.seed(123456)
  tre = rcoal(60)
  taxa = sort(tre$tip.label)
  b0=0; b1=1;
  x <- rTrait(n=1, phy=tre,model="BM",
              parameters=list(ancestral.state=0,sigma2=10))
  y <- b0 + b1*x +
    rTrait(n=1,phy=tre,model="lambda",parameters=list(
      ancestral.state=0,sigma2=1,lambda=0.5))
  # adding measurement errors and bootstrap
  z <- y + rnorm(60,0,1)
  dat = data.frame(trait=z[taxa],pred=x[taxa])

  ## Fit - Fixed values
  fit = phylolm(trait~pred,data=dat,phy=tre,model="lambda",
                starting.value = list(lambda = 0.5),
                upper.bound = list(lambda = 0.5),
                lower.bound = list(lambda = 0.5))

  expect_equal(fit$optpar, 0.5)

  ## Fit - Fixed values - unnamed (backward compatibility)
  fit = phylolm(trait~pred,data=dat,phy=tre,model="lambda",
                starting.value = 0.5,
                upper.bound = 0.5,
                lower.bound = 0.5)

  expect_equal(fit$optpar, 0.5)

  ## Fit - bounds
  fit <- expect_warning(phylolm(trait~pred,data=dat,phy=tre,model="kappa",
                                starting.value = list(kappa = 0.5),
                                upper.bound = list(kappa = 2),
                                lower.bound = list(kappa = 0.2)),
                        "the estimation of kappa matches the upper/lower bound for this parameter")
  
  expect_true(fit$optpar <= 2)
  expect_true(fit$optpar >= 0.2)

  ## Fit - bounds default
  fit <- expect_warning(phylolm(trait~pred,data=dat,phy=tre,model="kappa",
                                starting.value = list(kappa = 0.5),
                                upper.bound = list(kappa = 2)),
                        "the estimation of kappa matches the upper/lower bound for this parameter")
  
  expect_true(fit$optpar <= 2)
  expect_true(fit$optpar >= 10^{-6}) ## default lower bound

  ## Fit - bounds - measurement error
  fit <- expect_warning(phylolm(trait~pred,data=dat,phy=tre,model="BM", measurement_error = TRUE,
                                starting.value = list(sigma2_error = 0.5),
                                upper.bound = list(sigma2_error = 0.51),
                                lower.bound = list(sigma2_error = 0.49)),
                        "the estimation of sigma2_error matches the upper/lower bound for this parameter")

  expect_true(fit$sigma2_error <= 0.51 * fit$sigma2)
  expect_true(fit$sigma2_error >= 0.49 * fit$sigma2)

  ## Fit - bounds - measurement error
  expect_warning(phylolm(trait~pred,data=dat,phy=tre,model="OUfixedRoot", measurement_error = TRUE,
                         starting.value = list(alpha = 1, sigma2_error = 0.5),
                         upper.bound = list(alpha = 1, sigma2_error = 0.51),
                         lower.bound = list(alpha = 1, sigma2_error = 0.49)),
                 "the estimation of alpha matches the upper/lower bound for this parameter")
  fit <- expect_warning(phylolm(trait~pred,data=dat,phy=tre,model="OUfixedRoot", measurement_error = TRUE,
                                starting.value = list(alpha = 1, sigma2_error = 0.5),
                                upper.bound = list(alpha = 1, sigma2_error = 0.51),
                                lower.bound = list(alpha = 1, sigma2_error = 0.49)),
                        "the estimation of sigma2_error matches the upper/lower bound for this parameter")

  expect_true(fit$sigma2_error <= 0.51 * fit$sigma2 / (2 * 1))
  expect_true(fit$sigma2_error >= 0.49 * fit$sigma2 / (2 * 1))
  expect_equal(fit$optpar, 1)

  ## Fit - bounds - measurement error
  expect_warning(phylolm(trait~pred,data=dat,phy=tre,model="OUfixedRoot", measurement_error = TRUE,
                         starting.value = list(alpha = 0.5, sigma2_error = 0.5),
                         upper.bound = list(alpha = 0.6, sigma2_error = 0.55),
                         lower.bound = list(alpha = 0.4, sigma2_error = 0.49)),
                 "the estimation of alpha matches the upper/lower bound for this parameter")
  fit <- expect_warning(phylolm(trait~pred,data=dat,phy=tre,model="OUfixedRoot", measurement_error = TRUE,
                                starting.value = list(alpha = 0.5, sigma2_error = 0.5),
                                upper.bound = list(alpha = 0.6, sigma2_error = 0.55),
                                lower.bound = list(alpha = 0.4, sigma2_error = 0.49)),
                        "the estimation of sigma2_error matches the upper/lower bound for this parameter")

  expect_true(fit$sigma2_error <= 0.55 * fit$sigma2 / (2 * fit$optpar) + 1e-8)
  expect_true(fit$sigma2_error >= 0.49 * fit$sigma2 / (2 * fit$optpar) + 1e-8)
  expect_true(fit$optpar <= 0.6)
  expect_true(fit$optpar >= 0.4)

  ## Fit - default bounds - measurement error
  fit = phylolm(trait~pred,data=dat,phy=tre,model="delta", measurement_error = TRUE,
                starting.value = list(delta = 0.5, sigma2_error = 5 * 10^{15}),
                upper.bound = list(delta = 0.6),
                lower.bound = list(delta = 0.4, sigma2_error = 10^{15}))

  expect_true(fit$sigma2_error <= 10^{16} * fit$sigma2 / (2 * fit$optpar))
  expect_true(fit$sigma2_error >= 10^{15} * fit$sigma2 / (2 * fit$optpar))
  expect_true(fit$optpar <= 0.6)
  expect_true(fit$optpar >= 0.4)

  ## Fit - default bounds - measurement error - unnamed length one (backward compatibility)
  fit <- expect_warning(phylolm(trait~pred,data=dat,phy=tre,model="OUrandomRoot", measurement_error = TRUE,
                                starting.value = 0.5,
                                upper.bound = 0.5,
                                lower.bound = 0.5),
                        "the estimation of alpha matches the upper/lower bound for this parameter")

  expect_true(fit$sigma2_error <= 10^{16} * fit$sigma2 / (2 * fit$optpar))
  expect_true(fit$sigma2_error >= 10^{-16} * fit$sigma2 / (2 * fit$optpar))
  expect_equal(fit$optpar, 0.5)

  ## Fit - BM - unnamed length one
  fit = phylolm(trait~pred,data=dat,phy=tre,model="BM", measurement_error = FALSE,
                starting.value = 0.5,
                upper.bound = 0.5,
                lower.bound = 0.5)

  expect_equal(fit$sigma2_error, 0.0)
  expect_equal(fit$optpar, NULL)

  ## Fit - BM
  fit = phylolm(trait~pred,data=dat,phy=tre,model="BM", measurement_error = FALSE,
                starting.value = list(sigma2_error = 0.5),
                upper.bound = list(sigma2_error = 0.5),
                lower.bound = list(sigma2_error = 0.5))

  expect_equal(fit$sigma2_error, 0.0)
  expect_equal(fit$optpar, NULL)

})
