#' @export

# implement Tweedie (essentially Poisson) xgboost
SL.xgboost.pois <- function (Y, X, newX, family, obsWeights, id, nrounds = 50, params = list(),
                         max_depth = 5, min_child_weight = 6, eta = 0.1,
                         nthread = 1, verbose = 0, save_period = NULL, ...) {
  require("xgboost")
  if (packageVersion("xgboost") < "0.6") {
    stop("SL.xgboost requires xgboost version >= 0.6, try help('SL.xgboost') for details")
  }
  if (!is.matrix(X)) {
    X = model.matrix(~. - 1, X)
  }
  #max_thread <- parallel::detectCores()
  xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = rep(1, NROW(X)))
  model = xgboost::xgboost(data = xgmat,
                           objective = "reg:tweedie",
                           eval_metric = "rmse",
                           nrounds = nrounds,
                           max_depth = max_depth,
                           min_child_weight = min_child_weight,
                           eta = eta,
                           verbose = verbose,
                           nthread = nthread,
                           params = list(tweedie_variance_power = 1.01),
                           save_period = save_period,
                           early_stopping_rounds = 10)
  if (!is.matrix(newX)) {
    newX = model.matrix(~. - 1, newX)
  }
  pred = predict(model, newdata = newX)
  fit = list(object = model)
  class(fit) = c("SL.xgboost.pois")
  out = list(pred = pred, fit = fit)
  return(out)
}
