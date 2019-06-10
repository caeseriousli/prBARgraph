library("testthat")

test_that("Poisson regression w/o intercept is the same (small sample).", {
  set.seed(10)
  nobs  <- 20
  beta1 <- c(1, 0, 0.5, 0, 0, 0.7, 0.8, 0, 0.5)
  ncovs <- length(beta1)
  X <- matrix(rnorm(nobs * ncovs), nrow = nobs)
  y <- rpois(nobs, lambda = exp(X %*% beta1))
  a <- glm.fit(X, y, family = poisson())

  b <- poisBAR(y, X, lambda = 0, xi = 0, delta = 0,
          eps = 1E-6, tol = 1E-6,
          max.iter = 1000)
  expect_equal(as.vector(coef(a)), as.vector(b$coef), tolerance = 1E-6)
})

test_that("Poisson regression w/o intercept is the same (large sample).", {
  set.seed(10)
  nobs  <- 500
  beta1 <- c(1, 0, 0.5, 0, 0, 0.7, 0.8, 0, 0.5)
  ncovs <- length(beta1)
  X <- matrix(rnorm(nobs * ncovs), nrow = nobs)
  y <- rpois(nobs, lambda = exp(X %*% beta1))
  a <- glm.fit(X, y, family = poisson())

  b <- poisBAR(y, X, lambda = 0, xi = 0, delta = 0,
               eps = 1E-6, tol = 1E-6,
               max.iter = 1000)
  expect_equal(as.vector(coef(a)), as.vector(b$coef), tolerance = 1E-6)
})

