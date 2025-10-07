#' Fit niceday model
#'
#' @param W Matrix of multivariate outcome variables.
#' @param A Formula indicating binary covariate of interest.
#' @param X Formula indicating covariates to adjust for. Default is no adjustment, `~ 1`
#' @param data Data frame containing any covariates in formulas for `A` and `X`.
#' @param alpha Desired asymptotic type I error rate. Default is `0.05`.
#' @param uniform_CI Generate simultaneous confidence intervals across al categories? Default is `TRUE`.
#' @param d Parameter for smoothed median centering. Default is `0.1`.
#' @param gtrunc Truncation parameter passed to `est_nuis`, bounding the estimated propensity scores away from `0` and `1` by `gtrunc`. Default is `0.01`.
#' @param bs_rep Number of replicates for bootstrap sampling of uniform confidence intervals. Default is `1e5`.
#' @param num_crossval_folds Number of folds for cross-validation. Default is `10`.
#' @param num_crossfit_folds Number of folds for cross-fitting. Default is `10`.
#' @param enforce_pos_reg Should estimates of \eqn{E[W_j|A=a,X]} to forced to be strictly positive? Default is `FALSE`.
#' @param sl.lib.pi Libraries used to estimate the propensity score nuisance function. Default is `c("SL.mean", "SL.lm", "SL.glm.binom", "SL.xgboost.binom")`.
#' @param sl.lib.m Libraries used to estimate conditional mean mu_j nuisance functions. Default is `c("SL.mean", "SL.lm", "SL.glm.qpois", "SL.xgboost.pois")`.
#' @param sl.lib.q Libraries used to estimate conditional probabilities on nonzero observations. Default is the input to `sl.lib.pi`.
#' @param nuis Optional list of estimated nuisance functions (from function `est_nuis()`). Default is `NULL`.
#' @param verbose Do you want to receive updates as this function runs? Default is `FALSE`.
#' @param cross_fit Should cross-fitting be run? Default is `TRUE`.
#'
#' @return A list containing elements `coef`, `nuisances`, `crossfit_nuisances`, `variance`, and `call`. `coef` contains estimates,
#' standard errors, and confidence intervals for all categories. `nuisances` and `crossfit_nuisances` contain estimated nuisance
#' parameters and nuisance parameters from the crossfitting procedure (included when adjustment covariates are included).
#' `variance` contains the estimated covariance matrix. `call` contains the call to `ndFit()`.
#'
#' @examples
#' data(EcoZUR_meta)
#' data(EcoZUR_count)
#' ndFit(W = EcoZUR_count[, 1:50], # consider only the first 50 taxa to run quickly
#'       data = EcoZUR_meta,
#'       A = ~ Diarrhea,
#'       X = ~ sex + age_months,
#'       num_crossval_folds = 2, # use more folds in practice
#'       num_crossfit_folds = 2, # for cross validation and cross fitting
#'       sl.lib.pi = c("SL.mean"), # choosing single learner for the example to run quickly,
#'       sl.lib.m = c("SL.mean"))  # in practice would use other options as well
#'
#' @export
ndFit <- function(W,
                  A,
                  X = ~ 1,
                  data,
                  alpha = 0.05,
                  uniform_CI = TRUE,
                  d = 0.1,
                  gtrunc = 0.01,
                  bs_rep = 1e5,
                  num_crossval_folds = 10,
                  num_crossfit_folds = 10,
                  enforce_pos_reg = FALSE,
                  sl.lib.pi = c("SL.mean",
                                "SL.lm",
                                "SL.glm.binom",
                                "SL.xgboost.binom"),
                  sl.lib.m = c("SL.mean",
                               "SL.lm",
                               "SL.glm.qpois",
                               "SL.xgboost.pois"),
                  sl.lib.q = sl.lib.pi,
                  nuis = NULL,
                  verbose = FALSE,
                  cross_fit = TRUE) {

  call <- match.call(expand.dots = FALSE)

  # perform checks to ensure valid data input
  if (!is.data.frame(W)) {
    stop("W must be an nxJ data.frame")
  }

  if (!inherits(A, "formula")) {
    stop("A must be a formula")
  }

  # transform A from formula into a vector
  A <- stats::model.matrix(A, data)
  if (ncol(A) > 2) {
    stop("The `A` formula must only include a single covariate")
  }
  A <- A[, 2]
  if (any(!(A %in% c(0, 1)))) {
    stop("The covariate A must be binary")
  }

  if (nrow(data) != nrow(W)) {
    stop("The number of rows in `data` does not match the number of rows of `W`")
  }

  if (any(colMeans(W[A == 1, ]) == 0) | any(colMeans(W[A == 0, ]) == 0)) {
   stop(paste("One or more columns in W has a subgroup mean count of 0.",
              "You must remove it prior to input."))
  }

  if (mean(A == 1) == 0 | mean(A == 0) == 0) {
    stop("One of the subgroups of interest (0 or 1) is entirely absent.")
  }

  # transform X formula into design matrix
  X <- stats::model.matrix(X, data)
  if (ncol(X) == 1) {
    adjust_covariates = FALSE
  } else {
    adjust_covariates = TRUE
  }

  if (!is.null(nuis)) {
    # confirm that nuisance is provided in required format
  }

  # construct components of estimator
  n <- nrow(W)
  J <- ncol(W)

  #########################
  # Case 1: No adjustment #
  #########################

  # first, construct estimator in the simple case
  cov_W_A1 <- stats::cov(W[A == 1, ])
  E_W_A1 <-  colMeans(W[A == 1, ])
  P_A1 <- mean(A == 1)

  cov_W_A0 <- stats::cov(W[A == 0, ])
  E_W_A0 <-  colMeans(W[A == 0, ])
  P_A0 <- mean(A == 0)

  # calculate log-fold differences in observed data
  est <- log(E_W_A1) - log(E_W_A0)
  psi_hat_ABC_simp <- est

  # estimate variance of estimator
  Sigmahat <- (diag(1 / E_W_A1) %*% cov_W_A1 %*% diag(1 / E_W_A1) / P_A1 +
                   diag(1 / E_W_A0) %*% cov_W_A0 %*% diag(1 / E_W_A0) / P_A0)
  corrhat_est <- stats::cov2cor(Sigmahat)

  # standard error
  se_hat_psi_hat_ABC_simp <- sqrt(diag(Sigmahat) / n)

  lower_log_noadj_marg <- psi_hat_ABC_simp - stats::qnorm(1 - alpha / 2) * se_hat_psi_hat_ABC_simp
  upper_log_noadj_marg <- psi_hat_ABC_simp + stats::qnorm(1 - alpha / 2) * se_hat_psi_hat_ABC_simp

  # OPTIONAL: alternative code for estimating upper-alpha quantile for simultaneous coverage
  #
  # root <- uniroot(f =
  #                   function(t) {
  #                     mvtnorm::pmvnorm(lower = rep(-t, J),
  #                                      upper = rep(t, J),
  #                                      mean = rep(0, J),
  #                                      sigma = corrhat_est,
  #                                      algorithm = GenzBretz(maxpts = 100, abseps = 0.001))[1] - 0.951
  #                   },
  #                 interval = c(0.95 * stats::qnorm(1 - alpha / 2),
  #                              1.05 * stats::qnorm(1 - alpha / 2 / J)),
  #                 tol = 0.01)
  # (critical_value <- min(max(root$root + root$estim.prec,
  #                            stats::qnorm(1 - alpha / 2)),
  #                        stats::qnorm(1 - alpha / 2 / J)))

  if (uniform_CI) {
    q <- stats::quantile(apply(MASS::mvrnorm(n = 1e5, mu = rep(0, J),
                                      Sigma = corrhat_est),
                        1,
                        function(r){max(abs(r))}), 1 - alpha)
    cimult <- pmin(pmax(q, stats::qnorm(1 - alpha / 2)), stats::qnorm(1 - alpha / 2 / J))

    lower_log_noadj_sim <- psi_hat_ABC_simp - cimult * se_hat_psi_hat_ABC_simp
    upper_log_noadj_sim <- psi_hat_ABC_simp + cimult * se_hat_psi_hat_ABC_simp
  } else {
    lower_log_noadj_sim <- NA
    upper_log_noadj_sim <- NA
  }

  # apply smooth transformation of estimates
  g_plugin <- pseudohuber_center(x = est, d = d)
  grad_g_plugin <- dpseudohuber_center_dx(x = est, d = d)

  # construct estimator
  est_g <- est - g_plugin
  psi_hat_ABC_g_simp <- est_g

  # estimate covariance of estimates
  grad_h_mat <- t(diag(1, length(grad_g_plugin)) - grad_g_plugin)
  Sigmahat_g <- grad_h_mat %*% Sigmahat %*% t(grad_h_mat)
  corrhat_g <- stats::cov2cor(Sigmahat_g)

  # standard error
  se_hat_psi_hat_ABC_g_simp <- sqrt(diag(Sigmahat_g) / n)

  lower_g_log_noadj_marg <- psi_hat_ABC_g_simp - stats::qnorm(1 - alpha / 2) * se_hat_psi_hat_ABC_g_simp
  upper_g_log_noadj_marg <- psi_hat_ABC_g_simp + stats::qnorm(1 - alpha / 2) * se_hat_psi_hat_ABC_g_simp

  if (uniform_CI) {
    q <- stats::quantile(apply(MASS::mvrnorm(n = 1e5, mu = rep(0, J),
                                      Sigma = corrhat_g),
                        1,
                        function(r){max(abs(r))}), 1 - alpha)
    cimult <- pmin(pmax(q, stats::qnorm(1 - alpha / 2)), stats::qnorm(1 - alpha / 2 / J))

    lower_g_log_noadj_sim <- psi_hat_ABC_g_simp - cimult * se_hat_psi_hat_ABC_g_simp
    upper_g_log_noadj_sim <- psi_hat_ABC_g_simp + cimult * se_hat_psi_hat_ABC_g_simp
  } else {
    lower_g_log_noadj_sim <- NA
    upper_g_log_noadj_sim <- NA
  }

  results_simple <- list(Psi1 = data.frame(category = 1:J,
                                          est = psi_hat_ABC_simp,
                                          se = se_hat_psi_hat_ABC_simp,
                                          lower_marg = lower_log_noadj_marg,
                                          upper_marg = upper_log_noadj_marg,
                                          lower_sim = lower_log_noadj_sim,
                                          upper_sim = upper_log_noadj_sim),

                         Psi1g = data.frame(category = 1:J,
                                            est = psi_hat_ABC_g_simp,
                                            se = se_hat_psi_hat_ABC_g_simp,
                                            lower_marg = lower_g_log_noadj_marg,
                                            upper_marg = upper_g_log_noadj_marg,
                                            lower_sim = lower_g_log_noadj_sim,
                                            upper_sim = upper_g_log_noadj_sim))

  ###########################
  # Case 2: With adjustment #
  ###########################

  if (adjust_covariates) {
    # estimate nuisances E[W_j|A=1,X], E[W_j|A=0,X], and P(A=1|X), if not given
    if (is.null(nuis)) {
      nuis <- est_nuis(W = W,
                       A = A,
                       X = X,
                       num_crossfit_folds = num_crossfit_folds,
                       num_crossval_folds = num_crossval_folds,
                       sl.lib.pi = sl.lib.pi,
                       sl.lib.m = sl.lib.m,
                       sl.lib.q = sl.lib.q,
                       enforce_pos_reg = enforce_pos_reg,
                       gtrunc = gtrunc,
                       verbose = verbose)
    }

    # estimate log-fold ratio log(E[E[W_j|A=1,X]] / E[E[W_j|A=0,X]]) + c
    psi_ABCD <- est_psi_ABCD(W = W,
                             A = A,
                             X = X,
                             nuis = nuis,
                             alpha = alpha,
                             bs_rep = bs_rep,
                             d = d,
                             uniform_CI = uniform_CI,
                             verbose = verbose)

    results_adjust <- list(Psi2 = data.frame(category = 1:J,
                                             est = psi_ABCD$res$est,
                                             se = psi_ABCD$res$se,
                                             lower_marg = psi_ABCD$res$lower_marg,
                                             upper_marg = psi_ABCD$res$upper_marg,
                                             lower_sim = psi_ABCD$res$lower_sim,
                                             upper_sim = psi_ABCD$res$upper_sim,
                                             type = psi_ABCD$res$type),

                           Psi2g = data.frame(category = 1:J,
                                               est = psi_ABCD$res_g$est,
                                               se = psi_ABCD$res_g$se,
                                               lower_marg = psi_ABCD$res_g$lower_marg,
                                               upper_marg = psi_ABCD$res_g$upper_marg,
                                               lower_sim = psi_ABCD$res_g$lower_sim,
                                               upper_sim = psi_ABCD$res_g$upper_sim,
                                               type = psi_ABCD$res_g$type))
    cf_nuis <- psi_ABCD$cf_nuis
  } else {
    results_adjust <- NULL
    nuis <- NULL
    cf_nuis <- NULL
  }

  # now return results
  if (adjust_covariates) {
    results <- list(coef = results_adjust,
                    nuisances = nuis,
                    crossfit_nuisances = cf_nuis,
                    variance = list("Sigmahat" = psi_ABCD$Sigmahat$Sigmahat_log_AIPW,
                                    "Sigmahat_g" = psi_ABCD$Sigmahat$Sigmahat_g_log_AIPW))
  } else {
    results <- list(coef = results_simple,
                    variance = list("Sigmahat" = Sigmahat, "Sigmahat_g" = Sigmahat_g))
  }
  results$call <- call

  return(structure(results, class = "ndFit"))
}
