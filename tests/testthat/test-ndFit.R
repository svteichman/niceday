test_that("ndFit works as expected", {

  # fix this to test something more meaningful!
  W <- data.frame(matrix(rpois(20 * 30, 100), nrow = 20, ncol = 30))
  A <- rep(0:1, each = 10)
  X <- data.frame(rnorm(20, 0, 1), sample(c(1:3), size = 20, replace = TRUE))

  # fit niceday model
  niceday_res <- ndFit(W = W, A = A, X = X)

  expect_type(niceday_res, "list")

  # also add tests that things fail when expected (if you give the wrong inputs), and that test different
  # argument values
})
