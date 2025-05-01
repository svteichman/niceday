#' @export

ndFit <- function(W, # matrix of responses
                  A, # binary case/control covariate
                  X, # matrix of covariates for adjustment

                  alpha = 0.05, # asymptotic type-I error control

                  uniform_CI = TRUE,

                  d = 0.1, # parameter for 'smooth median' centering

                  gtrunc = 0.01,

                  bs_rep = 1e5, # number of replicates for bootstrap sampling of uniform confidence intervals

                  num_crossval_folds = 10,
                  num_crossfit_folds = 10,

                  enforce_pos_reg = FALSE, # force estimates of E[W_j|A=a,X] to be strictly positive

                  sl.lib.pi = c("SL.mean",
                                "SL.lm",
                                "SL.glm.binom",
                                "SL.xgboost.binom"),
                  sl.lib.m = c("SL.mean",
                               "SL.lm",
                               "SL.glm.qpois",
                               "SL.xgboost.pois"),

                  sl.lib.q = sl.lib.pi,

                  adjust_covariates = TRUE,

                  nuis = NULL,

                  verbose = FALSE,

                  cross_fit = TRUE) {

  # perform checks to ensure valid data input
  if (!is.data.frame(W)) {
    stop("W must be an nxJ data.frame")
  }

  if (!is.vector(A)) {
    stop("A must be a length n vector")
  }

  if (length(A) != nrow(W)) {
    stop("The length of A does not match the number of rows of W")
  }

  if (any(colMeans(W[A == 1, ]) == 0) | any(colMeans(W[A == 0, ]) == 0)) {
   stop(paste("One or more columns in W has a subgroup mean count of 0.",
              "You must remove it prior to input."))
  }

  if (mean(A == 1) == 0 | mean(A == 0) == 0) {
    stop("One of the subgroups of interest (0 or 1) is entirely absent.")
  }

  if (adjust_covariates) {
    if (!(is.data.frame(X) | is.matrix(X))) {
      stop("User requested covariate adjustment. The input X must be ",
           "provided as a matrix or a data.frame object.")
    }
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
  covhat_est <- (cov_W_A1 / (E_W_A1 %*% t(E_W_A1)) / P_A1 +
                   cov_W_A0 / (E_W_A0 %*% t(E_W_A0)) / P_A0)
  corrhat_est <- covhat_est /
    (sqrt(diag(covhat_est)) %*% t(sqrt(diag(covhat_est))))

  # standard error
  se_hat_psi_hat_ABC_simp <- sqrt(diag(covhat_est) / n)

  lower_log_noadj_marg <- psi_hat_ABC_simp - qnorm(1 - alpha / 2) * se_hat_psi_hat_ABC_simp
  upper_log_noadj_marg <- psi_hat_ABC_simp + qnorm(1 - alpha / 2) * se_hat_psi_hat_ABC_simp

  if (uniform_CI) {
    q <- quantile(apply(MASS::mvrnorm(n = 1e5, mu = rep(0, J),
                                      Sigma = corrhat_est),
                        1,
                        function(r){max(abs(r))}), 1 - alpha)
    cimult <- pmin(pmax(q, qnorm(1 - alpha / 2)), qnorm(1 - alpha / 2 / J))

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
  Sigmahat_g <- t(grad_h_mat) %*% covhat_est %*% grad_h_mat
  corrhat_g <- Sigmahat_g /
    (sqrt(diag(Sigmahat_g)) %*% t(sqrt(diag(Sigmahat_g))))

  # standard error
  se_hat_psi_hat_ABC_g_simp <- sqrt(diag(Sigmahat_g) / n)

  lower_g_log_noadj_marg <- psi_hat_ABC_g_simp - qnorm(1 - alpha / 2) * se_hat_psi_hat_ABC_g_simp
  upper_g_log_noadj_marg <- psi_hat_ABC_g_simp + qnorm(1 - alpha / 2) * se_hat_psi_hat_ABC_g_simp

  if (uniform_CI) {
    q <- quantile(apply(MASS::mvrnorm(n = 1e5, mu = rep(0, J),
                                      Sigma = corrhat_g),
                        1,
                        function(r){max(abs(r))}), 1 - alpha)
    cimult <- pmin(pmax(q, qnorm(1 - alpha / 2)), qnorm(1 - alpha / 2 / J))

    lower_g_log_noadj_sim <- psi_hat_ABC_g_simp - cimult * se_hat_psi_hat_ABC_g_simp
    upper_g_log_noadj_sim <- psi_hat_ABC_g_simp + cimult * se_hat_psi_hat_ABC_g_simp
  } else {
    lower_g_log_noadj_sim <- NA
    upper_g_log_noadj_sim <- NA
  }

  results_simple <- list(ABC = data.frame(taxon = 1:J,
                                          est = psi_hat_ABC_simp,
                                          se = se_hat_psi_hat_ABC_simp,
                                          lower_marg = lower_log_noadj_marg,
                                          upper_marg = upper_log_noadj_marg,
                                          lower_sim = lower_log_noadj_sim,
                                          upper_sim = upper_log_noadj_sim),

                         ABC_g = data.frame(taxon = 1:J,
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

    results_adjust <- list(ABCD = data.frame(taxon = 1:J,
                                             est = psi_ABCD$res$est,
                                             se = psi_ABCD$res$se,
                                             lower_marg = psi_ABCD$res$lower_marg,
                                             upper_marg = psi_ABCD$res$upper_marg,
                                             lower_sim = psi_ABCD$res$lower_sim,
                                             upper_sim = psi_ABCD$res$upper_sim,
                                             type = psi_ABCD$res$type),

                           ABCD_g = data.frame(taxon = 1:J,
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
  results <- list(noadjust = results_simple,
                  adjust = results_adjust,
                  nuis = nuis,
                  cf_nuis = cf_nuis)

  return(results)
}
