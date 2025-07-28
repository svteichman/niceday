#' @export

# implement ranger for binomial family
SL.ranger.binom <- function (Y, X, newX, obsWeights,
                             num.trees = 100,
                             min.node.size = 1,
                             probability = TRUE,
                             mtry = floor(sqrt(ncol(X))), write.forest = TRUE,
                             sample.fraction = ifelse(replace, 1, 0.632), 
                             replace = TRUE, num.threads = 1, verbose = FALSE, ...) {
  require("ranger")
  Y = factor(Y, levels = c(1, 0))
  if (is.matrix(X)) {
    X = data.frame(X)
  }
  fit <- ranger::ranger(`_Y` ~ ., data = cbind(`_Y` = Y, X), 
                        num.trees = num.trees, mtry = mtry, min.node.size = min.node.size, 
                        replace = replace, sample.fraction = sample.fraction, 
                        case.weights = obsWeights, write.forest = write.forest, 
                        probability = probability, num.threads = num.threads, 
                        verbose = verbose)
  pred <- predict(fit, data = newX)$predictions[, "1"]
  fit <- list(object = fit, verbose = verbose)
  class(fit) <- c("SL.ranger")
  out <- list(pred = pred, fit = fit)
  return(out)
}
