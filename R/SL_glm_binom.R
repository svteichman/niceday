#' @export

# re-implement binomial GLM with more iterations and larger epsilon
SL.glm.binom <- function (Y, X, newX, family, obsWeights, model = TRUE, ...) {
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
  fit.glm <- glm(Y ~ ., data = X, family = binomial(link = "logit"),
                 weights = obsWeights, model = model,
                 control = list(maxit = 1e3, epsilon = 1e-4))
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- "SL.glm.binom"
  out <- list(pred = pred, fit = fit)
  return(out)
}
