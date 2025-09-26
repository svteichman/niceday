# implement binomial xgboost
SL.xgboost.binom <- function (Y, X, newX, family, obsWeights, id, nrounds = 50, params = list(),
                             max_depth = 5, min_child_weight = 6, eta = 0.1,
                             nthread = 1, verbose = 0, save_period = NULL, ...) {

  if (utils::packageVersion("xgboost") < "0.6") {
    stop("SL.xgboost requires xgboost version >= 0.6, try help('SL.xgboost') for details")
  }
  if (!is.matrix(X)) {
    X = stats::model.matrix(~. - 1, X)
  }
  #max_thread <- parallel::detectCores()
  xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = rep(1, NROW(X)))
  model = xgboost::xgboost(data = xgmat,
                           objective = "binary:logistic",
                           eval_metric = "logloss",
                           nrounds = nrounds,
                           max_depth = max_depth,
                           min_child_weight = min_child_weight,
                           eta = eta,
                           verbose = verbose,
                           nthread = nthread,
                           params = params,
                           save_period = save_period,
                           early_stopping_rounds = 10)
  if (!is.matrix(newX)) {
    newX = stats::model.matrix(~. - 1, newX)
  }
  pred = stats::predict(model, newdata = newX)
  fit = list(object = model)
  class(fit) = c("SL.xgboost.binom")
  out = list(pred = pred, fit = fit)
  return(out)
}
