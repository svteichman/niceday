#' @export

# implement ranger for Poisson family (i.e., not probabilities)
SL.ranger.pois <- function(Y, X, newX, family, obsWeights, num.trees = 500, mtry = floor(sqrt(ncol(X))), 
                           write.forest = TRUE, replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632), 
                           num.threads = 1, verbose = T, ...) {
  require("ranger")
  if (is.matrix(X)) {
    X = data.frame(X)
  }
  fit <- ranger::ranger(`_Y` ~ ., data = cbind(`_Y` = Y, X), 
                        num.trees = num.trees, mtry = mtry, min.node.size = 5, 
                        replace = replace, sample.fraction = sample.fraction, 
                        case.weights = obsWeights, write.forest = write.forest, 
                        probability = FALSE, num.threads = num.threads, 
                        verbose = verbose)
  pred <- predict(fit, data = newX)$predictions
  fit <- list(object = fit, verbose = verbose)
  class(fit) <- c("SL.ranger")
  out <- list(pred = pred, fit = fit)
  return(out)
}
