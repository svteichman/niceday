# implement Poisson HAL
SL.hal9001.pois <- function (Y, X, newX, family, obsWeights, id, max_degree = 2,
                              smoothness_orders = 1, num_knots = 5, ...) {
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (!is.null(newX) & !is.matrix(newX)) {
    newX <- as.matrix(newX)
  }
  hal_fit <- hal9001::fit_hal(Y = Y, X = X, family = "poisson",
                              weights = obsWeights, id = id, max_degree = max_degree,
                              smoothness_orders = smoothness_orders, num_knots = num_knots,
                              ...)
  if (!is.null(newX)) {
    pred <- stats::predict(hal_fit, new_data = newX)
  }
  else {
    pred <- stats::predict(hal_fit, new_data = X)
  }
  fit <- list(object = hal_fit)
  class(fit) <- "SL.hal9001"
  out <- list(pred = pred, fit = fit)
  return(out)
}
