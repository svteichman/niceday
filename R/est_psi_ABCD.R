#' @export

est_psi_ABCD <- function(W,
                         A,
                         X,
                         nuis,
                         d = 0.1,
                         gtrunc = 0.01,
                         alpha = 0.05,
                         bs_rep = 1e5,
                         uniform_CI = TRUE,
                         use_TMLE = TRUE,
                         TMLE_fluctuation = "poisson",
                         quiet = FALSE) {

  n <- nrow(W)
  J <- ncol(W)

  #################################################
  # Step 1:                                       #
  # Construct one-step estimator of E[E[W|A=a,X]] #
  #################################################

  fold_list <- nuis$fold_list
  nfold <- length(fold_list)

  # create arrays to store estimated influence function applied to data across folds
  if_AIPW1_array <- array(NA,
                          dim = list(n, J, nfold),
                          dimnames = c(list(paste0("samp", 1:n)),
                                       list(paste0("tax", 1:J)),
                                       list(paste0("fold", 1:nfold))))
  if_AIPW0_array <- if_AIPW1_array

  # create arrays to store estimated AIPW estimator of E[E[W|A=a,X]] across folds
  est_AIPW1_mat <- array(NA,
                         dim = list(J, nfold),
                         dimnames = c(list(paste0("tax", 1:J)),
                                      list(paste0("fold", 1:nfold))))
  est_AIPW0_mat <- est_AIPW1_mat

  # create matrices to store cross-fitted nuisances for P(A=1|X), E[W_j|W_j>0,A,X], and P(W_j>0|A,X)
  cf_est_m1 <- matrix(NA, nrow = n, ncol = J,
                       dimnames = c(list(paste0("samp", 1:n)),
                                    list(paste0("tax", 1:J))))
  cf_est_q0 <- cf_est_q1 <- cf_est_m0 <- cf_est_m1
  cf_esp_pi <- rep(NA, n)

  # now iterate across folds
  for (k in 1:nfold) {
    # identify which samples were used to learn the nuisances
    samp_subset_k <- unname(fold_list[[k]])
    samp_subset_comp_k <- sort(unname(unlist(fold_list[-k])))

    # map nuisance to variables
    est_m1_k <- nuis$mat_m1[,,k]
    est_m0_k <- nuis$mat_m0[,,k]

    est_q1_k <- nuis$mat_q1[,,k]
    est_q0_k <- nuis$mat_q0[,,k]

    est_pi_k <- nuis$mat_pi[,,k]

    # store the cross-fitted nuisances
    cf_est_m1[samp_subset_k, ] <- est_m1_k[samp_subset_k,]
    cf_est_m0[samp_subset_k, ] <- est_m0_k[samp_subset_k,]

    cf_est_q1[samp_subset_k, ] <- est_q1_k[samp_subset_k,]
    cf_est_q0[samp_subset_k, ] <- est_q0_k[samp_subset_k,]

    cf_esp_pi[samp_subset_k] <- est_pi_k[samp_subset_k]

    # calculate plug-in estimator
    plugin_AIPW1_k <- colMeans(est_m1_k[samp_subset_comp_k, , drop = FALSE] *
                                 est_q1_k[samp_subset_comp_k, , drop = FALSE])

    plugin_AIPW0_k <- colMeans(est_m0_k[samp_subset_comp_k, , drop = FALSE] *
                                 est_q0_k[samp_subset_comp_k, , drop = FALSE])

    # calculate the influence functions for E[E[W|A=a,X]]
    if_AIPW1_k <- t(t((A == 1) * (W - est_m1_k * est_q1_k) / est_pi_k + est_m1_k * est_q1_k) - plugin_AIPW1_k)
    if_AIPW0_k <- t(t((A == 0) * (W - est_m0_k * est_q0_k) / (1 - est_pi_k) + est_m0_k * est_q0_k) - plugin_AIPW0_k)

    # save results
    if_AIPW1_array[,,k] <- as.matrix(if_AIPW1_k)
    if_AIPW0_array[,,k] <- as.matrix(if_AIPW0_k)
  }

  # calculate average value of the estimated influence function applied to each data point across cross-fitted slices
  avg_if_AIPW1 <- apply(if_AIPW1_array, c(1, 2), mean)
  avg_if_AIPW0 <- apply(if_AIPW0_array, c(1, 2), mean)

  # confirm that the infleunce function contains only real-valued numbers
  if (!all(is.numeric(c(avg_if_AIPW0, avg_if_AIPW1)) &
           is.finite(c(avg_if_AIPW0, avg_if_AIPW1)) &
           !is.nan(c(avg_if_AIPW0, avg_if_AIPW1)) &
           !is.na(c(avg_if_AIPW0, avg_if_AIPW1)))) {
    stop("The estimated influence function contains non-real-valued numbers.")
  }

  #####################################################################
  # Step 2:                                                           #
  # Use TMLE to update the conditional mean regressions such that the #
  # empirical mean of the estimated influence functions are 0         #
  #####################################################################

  est_AIPW0_tmle <- est_AIPW1_tmle <- rep(NA, J)

  if (!quiet) {
    cat("Beginning TMLE\n")
    pb <- txtProgressBar(min = 0, max = J, style = 3)
  }

  for (j in 1:J) {
    ###############################################################################
    # Part A: apply ad-hoc perturbation to conditional mean regression if W_j > 0 #
    ###############################################################################
    mhat_1j <- cf_est_m1[, j]
    mhat_0j <- cf_est_m0[, j]
    mhat_Aj <- ifelse(A == 1, mhat_1j, mhat_0j)

    # we cannot estimate E[W_j|W_j>0,A,X] by zero, so replace by smallest prediction
    which_shift_m <- W[, j] > 0 & mhat_Aj == 0
    shift_m <- min(c(c(W[W[, j] > 0, j]),
                   c(nuis$mat_m0[ , j, ])[c(nuis$mat_m0[ , j, ]) > 0],
                   c(nuis$mat_m1[ , j, ])[c(nuis$mat_m1[ , j, ]) > 0]))
    mhat_1j[which_shift_m] <- shift_m
    mhat_0j[which_shift_m] <- shift_m
    mhat_Aj[which_shift_m] <- shift_m

    m_epsilon <- coef(glm(formula = clever.resp ~ -1 + offset(clever.offset) + clever.covar1 + clever.covar2,
                          weights = clever.weight,
                          subset = clever.subset,
                          data =
                            data.frame(clever.resp   = W[, j],
                                       clever.weight = (1 / ifelse(A == 1, cf_esp_pi, 1 - cf_esp_pi)) /
                                         mean(1 / ifelse(A == 1, cf_esp_pi, 1 - cf_esp_pi)),
                                       clever.offset = log(mhat_Aj),
                                       clever.covar1 = A,
                                       clever.covar2 = 1 - A,
                                       clever.subset = W[, j] > 0),
                          family = quasipoisson(link = "log"),
                          control = glm.control(epsilon = 1e-10,
                                                maxit = 1e3)))

    mhat_1j_update <- (mhat_1j) * exp(m_epsilon[1] * A)
    mhat_0j_update <- (mhat_0j) * exp(m_epsilon[2] * (1 - A))
    mhat_Aj_update <- ifelse(A == 1, mhat_1j_update, mhat_0j_update)

    ####################################################################
    # Part B: apply perturbation to full iterated expectation estimate #
    ####################################################################

    if (all(W[, j] > 0)) {
      qhat_1j_update <- rep(1, n)
      qhat_0j_update <- rep(1, n)
      qhat_Aj_update <- rep(1, n)
    } else {
      qhat_1j <- cf_est_q1[, j]
      qhat_0j <- cf_est_q0[, j]
      qhat_Aj <- ifelse(A == 1, qhat_1j, qhat_0j)

      # we cannot estimate P(W_j>0|A,X) by zero or one, so truncate by a small amount
      which_trunc_q <- (W[, j] > 0 & qhat_Aj == 0) | (W[, j] == 0 & qhat_Aj == 1)
      qpreds <- c(c(nuis$mat_q0[ , j, ]), c(nuis$mat_q1[ , j, ]))
      trunc_q <- min(c(c(qpreds[qpreds > 0], 1 - qpreds[qpreds < 1]),
                       1 - mean(W[, j] > 0), mean(W[, j] > 0),
                       1e-4))
      qhat_1j[which_trunc_q] <- pmax(pmin(qhat_1j[which_trunc_q], 1 - trunc_q), trunc_q)
      qhat_0j[which_trunc_q] <- pmax(pmin(qhat_0j[which_trunc_q], 1 - trunc_q), trunc_q)
      qhat_Aj <- ifelse(A == 1, qhat_1j, qhat_0j)

      q_epsilon <- coef(glm(formula = clever.resp ~ -1 + offset(clever.offset) + clever.covar1 + clever.covar2,
                            weights = clever.weight,
                            subset = clever.subset,
                            data =
                              data.frame(clever.resp   = ifelse(W[, j] > 0, 1, 0),
                                         clever.weight = (mhat_Aj_update / ifelse(A == 1, cf_esp_pi, 1 - cf_esp_pi)) /
                                           mean(mhat_Aj_update / ifelse(A == 1, cf_esp_pi, 1 - cf_esp_pi)),
                                         clever.offset = qlogis(qhat_Aj),
                                         clever.covar1 = A,
                                         clever.covar2 = 1 - A,
                                         clever.subset = (W[, j] > 0 & qhat_Aj < 1) |
                                                         (W[, j] == 0 & qhat_Aj > 0)),
                            family = quasibinomial(link = "logit"),
                            control = glm.control(epsilon = 1e-10,
                                                  maxit = 1e3)))

      qhat_1j_update <- plogis(qlogis(qhat_1j) + q_epsilon[1] * A)
      qhat_0j_update <- plogis(qlogis(qhat_0j) + q_epsilon[2] * (1 - A))
      qhat_Aj_update <- ifelse(A == 1, qhat_1j_update, qhat_0j_update)

      # ggplot(data = data.frame(pred = qhat_1j_update,
      #                          actual = factor(ifelse(W[, j] > 0, 1, 0)),
      #                          treatment = factor(A))) +
      #   geom_jitter(aes(x = treatment, y = pred, color = actual), width = 0.2, alpha = 0.5) +
      #   labs(title = "Logistic Regression Fit",
      #        x = "Treatment (A)",
      #        y = "Probability / Outcome") +
      #   theme_minimal()
    }

    if (any(abs(c(mean((A == 0) * (W[, j] - mhat_0j_update * qhat_0j_update) / (1 - cf_esp_pi)),
                  mean((A == 1) * (W[, j] - mhat_1j_update * qhat_1j_update) / cf_esp_pi))) > 1e-3)) {
      warning(paste("The TMLE fluctuation in taxon", j ,"seems to have poor fit.\n",
                    "A = 0:", mean((A == 0) * (W[, j] - mhat_0j_update * qhat_0j_update) / (1 - cf_esp_pi)),"\n",
                    "A = 1:", mean((A == 1) * (W[, j] - mhat_1j_update * qhat_1j_update) / cf_esp_pi),"\n"))
    }

    est_AIPW1_tmle[j] <- mean(mhat_1j_update * qhat_1j_update)
    est_AIPW0_tmle[j] <- mean(mhat_0j_update * qhat_0j_update)

    # update progress bar
    if (!quiet) {
      setTxtProgressBar(pb, j)
    }
  }

  ##################################################################
  # Step 2:                                                        #
  # Apply Delta Method to learn log(E[E[W|A=1,X]] / E[E[W|A=0,X]]) #
  ##################################################################

  # construct the TMLE estimator
  est_log_AIPW <- log(est_AIPW1_tmle) - log(est_AIPW0_tmle)
  est_AIPW1 <- est_AIPW1_tmle
  est_AIPW0 <- est_AIPW0_tmle

  avg_if_log_AIPW <- t(t(avg_if_AIPW1) / est_AIPW1) - t(t(avg_if_AIPW0) / est_AIPW0)

  Sigmahat_log_AIPW <- cov(avg_if_log_AIPW)

  ###########################################################################
  # Step 3:                                                                 #
  # Apply Delta Method to learn log(E[E[W|A=1,X]] / E[E[W|A=0,X]]) - g(...) #
  ###########################################################################

  est_g_log_AIPW <- est_log_AIPW - niceday::pseudohuber_center(x = est_log_AIPW, d = d)

  est_grad_g_log_AIPW <- niceday::dpseudohuber_center_dx(x = est_log_AIPW, d = d)

  avg_if_g_log_AIPW <- avg_if_log_AIPW - (avg_if_log_AIPW %*% est_grad_g_log_AIPW) %*% rep(1, J)

  Sigmahat_g_log_AIPW <- cov(avg_if_g_log_AIPW)

  ##########################################################################
  # Step 4:                                                                #
  # Construct marginal and uniform confidence intervals for each parameter #
  ##########################################################################

  # produce confidence intervals and p-values
  se_est_log_AIPW <- sqrt(diag(Sigmahat_log_AIPW)) / sqrt(n)

  lower_log_AIPW_marg <- est_log_AIPW - qnorm(1 - alpha / 2) * se_est_log_AIPW
  upper_log_AIPW_marg <- est_log_AIPW + qnorm(1 - alpha / 2) * se_est_log_AIPW

  if (uniform_CI) {
    q <- quantile(apply(MASS::mvrnorm(n = bs_rep, mu = rep(0, J),
                                      Sigma = cor(avg_if_log_AIPW)),
                        1,
                        function(r){max(abs(r))}), 1 - alpha)
    cimult <- pmin(pmax(q, qnorm(1 - alpha / 2)), qnorm(1 - alpha / 2 / J))

    lower_log_AIPW_sim <- est_log_AIPW - cimult * se_est_log_AIPW
    upper_log_AIPW_sim <- est_log_AIPW + cimult * se_est_log_AIPW
  } else {
    lower_log_AIPW_sim <- NA
    upper_log_AIPW_sim <- NA
  }

  pval_log_AIPW <- 2 * (1 - pnorm(abs(est_log_AIPW / se_est_log_AIPW)))

  # produce confidence intervals and p-values
  se_est_g_log_AIPW <- sqrt(diag(Sigmahat_g_log_AIPW)) / sqrt(n)

  lower_g_log_AIPW_marg <- est_g_log_AIPW - qnorm(1 - alpha / 2) * se_est_g_log_AIPW
  upper_g_log_AIPW_marg <- est_g_log_AIPW + qnorm(1 - alpha / 2) * se_est_g_log_AIPW

  if (uniform_CI) {
    q <- quantile(apply(MASS::mvrnorm(n = bs_rep, mu = rep(0, J),
                                      Sigma = cor(avg_if_g_log_AIPW)),
                        1,
                        function(r){max(abs(r))}), 1 - alpha)
    cimult_g <- pmin(pmax(q, qnorm(1 - alpha / 2)), qnorm(1 - alpha / 2 / J))

    lower_g_log_AIPW_sim <- est_g_log_AIPW - cimult_g * se_est_g_log_AIPW
    upper_g_log_AIPW_sim <- est_g_log_AIPW + cimult_g * se_est_g_log_AIPW
  } else {
    lower_g_log_AIPW_sim <- NA
    upper_g_log_AIPW_sim <- NA
  }

  pval_g_log_AIPW <- 2 * (1 - pnorm(abs(est_g_log_AIPW / se_est_g_log_AIPW)))

  pnames <- paste0("log(E[E[V", 1:ncol(W), "|A=1,X]] / E[E[V", 1:ncol(W),"|A=0,X]])")

  res <- data.frame(taxon = 1:ncol(W),
                    parameter = pnames,
                    est = est_log_AIPW,
                    se = se_est_log_AIPW,
                    lower_marg = lower_log_AIPW_marg,
                    upper_marg = upper_log_AIPW_marg,
                    lower_sim = lower_log_AIPW_sim,
                    upper_sim = upper_log_AIPW_sim,
                    pval = pval_log_AIPW,
                    type = ifelse(use_TMLE, "TMLE", "OSE"),
                    row.names = NULL)

  res_g <- data.frame(taxon = 1:ncol(W),
                      parameter = paste0(pnames, " - g(...)"),
                      est = est_g_log_AIPW,
                      se = se_est_g_log_AIPW,
                      lower_marg = lower_g_log_AIPW_marg,
                      upper_marg = upper_g_log_AIPW_marg,
                      lower_sim = lower_g_log_AIPW_sim,
                      upper_sim = upper_g_log_AIPW_sim,
                      pval = pval_g_log_AIPW,
                      type = ifelse(use_TMLE, "TMLE", "OSE"),
                      row.names = NULL)

  return(list(res = res,
              res_g = res_g,
              nuis = list(cf_esp_pi = cf_esp_pi,
                          cf_est_m0 = cf_est_m0,
                          cf_est_m1 = cf_est_m1,
                          cf_est_q0 = cf_est_q0,
                          cf_est_q1 = cf_est_q1),
              crossfit_if = list(avg_if_AIPW1 = avg_if_AIPW1,
                                 avg_if_AIPW0 = avg_if_AIPW0,
                                 avg_if_log_AIPW = avg_if_log_AIPW,
                                 avg_if_g_log_AIPW = avg_if_g_log_AIPW)))
}
