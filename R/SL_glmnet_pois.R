# implement poisson GLM with lasso regularization
SL.glmnet.pois <- function (Y, X, newX, family, obsWeights, id, alpha = 1, nfolds = 10,
                          nlambda = 100, useMin = TRUE, loss = "deviance", ...)
{
  if (!is.matrix(X)) {
    X <- stats::model.matrix(~-1 + ., X)
    newX <- stats::model.matrix(~-1 + ., newX)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights,
                             lambda = NULL, type.measure = loss, nfolds = nfolds,
                             family = "poisson", alpha = alpha, nlambda = nlambda,
                             ...)
  pred <- stats::predict(fitCV, newx = newX, type = "response",
                  s = ifelse(useMin, "lambda.min", "lambda.1se"))
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.glmnet"
  out <- list(pred = pred, fit = fit)
  return(out)
}
