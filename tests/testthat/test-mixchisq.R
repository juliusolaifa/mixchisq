test_that("weights are normalized", {
  expect_equal(
    pmixchisq(1, df = c(2, 3), w = c(1, 1)),
    pmixchisq(1, df = c(2, 3), w = c(0.5, 0.5))
  )
})

test_that("pmixchisq agrees with explicit formula without atom", {
  x <- 3
  expected <- 0.25 * pchisq(x, df = 2) + 0.75 * pchisq(x, df = 5)
  expect_equal(pmixchisq(x, df = c(2, 5), w = c(0.25, 0.75)), expected)
})

test_that("pmixchisq handles atom at zero", {
  expect_equal(pmixchisq(0, df = c(0, 1), w = c(0.5, 0.5)), 0.5)
  expect_equal(pmixchisq(-1, df = c(0, 1), w = c(0.5, 0.5)), 0)
})

test_that("qmixchisq returns zero below the atom", {
  q <- qmixchisq(c(0.1, 0.5), df = c(0, 1), w = c(0.5, 0.5))
  expect_equal(q, c(0, 0))
})

test_that("rmixchisq returns zeros when df=0 component exists", {
  set.seed(42)
  x <- rmixchisq(200, df = c(0, 1), w = c(0.5, 0.5))
  expect_true(any(x == 0))
  expect_true(all(x >= 0))
})
