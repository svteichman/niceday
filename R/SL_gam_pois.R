#' @export

# implement Generalized Additive Model (GAM) under Poisson family
SL.gam.pois <- function (Y, X, newX, family, obsWeights, deg.gam = 2, cts.num = 4,
                          ...) {
  if (!requireNamespace("gam")) {
    stop("SL.gam.qpois requires the gam package, but it isn't available")
  }
  if (!"package:gam" %in% search()) {attachNamespace("gam")}
  cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))
  if (sum(!cts.x) > 0) {
    gam.model <- as.formula(paste("Y~", paste(paste("s(",
                            colnames(X[, cts.x, drop = FALSE]), ",", deg.gam,
                            ")", sep = ""), collapse = "+"), "+",
                            paste(colnames(X[, !cts.x, drop = FALSE]), collapse = "+")))
  } else {
    gam.model <- as.formula(paste("Y~", paste(paste("s(",
                            colnames(X[, cts.x, drop = FALSE]), ",", deg.gam,
                            ")", sep = ""), collapse = "+")))
  }
  if (sum(!cts.x) == length(cts.x)) {
    gam.model <- as.formula(paste("Y~", paste(colnames(X),
                                              collapse = "+"), sep = ""))
  }
  fit.gam <- gam::gam(gam.model, data = data.frame(Y, X), family = quasipoisson(link = "log"),
                      control = gam::gam.control(maxit = 50, bf.maxit = 50),
                      weights = obsWeights)
  if (packageVersion("gam") >= "1.15") {
    pred <- gam::predict.Gam(fit.gam, newdata = data.frame(newX), type = "response")
  } else {
    stop("This SL.gam.qpois wrapper requires gam version >= 1.15, please update the gam package with 'update.packages('gam')'")
  }
  fit <- list(object = fit.gam)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.gam.pois")
  return(out)
}


