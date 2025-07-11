#' @export

est_psi_ABCD <- function(W,
                         A,
                         X,
                         nuis,
                         d = 0.1,
                         alpha = 0.05,
                         bs_rep = 1e5,
                         uniform_CI = TRUE,
                         verbose = FALSE,
                         tmle_fluctuation = "poisson") {

  n <- nrow(W)
  J <- ncol(W)

  ####################################################
  # Step 1:                                          #
  # Extract cross-fitted nuisances for E[E[W|A=a,X]] #
  ####################################################

  fold_list <- nuis$fold_list
  nfold <- length(fold_list)

  # create matrices to store cross-fitted nuisances
  cf_est_m1 <- matrix(NA, nrow = n, ncol = J,
                       dimnames = c(list(paste0("samp", 1:n)),
                                    list(paste0("tax", 1:J))))
  cf_est_q0 <- cf_est_q1 <- cf_est_m0 <- cf_est_m1
  cf_est_pi <- rep(NA, n)

  # now iterate across folds
  for (k in 1:nfold) {
    # identify which samples were used to learn the nuisances
    samp_subset_k <- unname(fold_list[[k]])
    samp_subset_comp_k <- sort(unname(unlist(fold_list[-k])))

    # store the cross-fitted nuisances
    cf_est_m1[samp_subset_k, ] <- nuis$mat_m1[samp_subset_k, , k, drop = FALSE]
    cf_est_m0[samp_subset_k, ] <- nuis$mat_m0[samp_subset_k, , k, drop = FALSE]

    cf_est_q1[samp_subset_k, ] <- nuis$mat_q1[samp_subset_k, , k, drop = FALSE]
    cf_est_q0[samp_subset_k, ] <- nuis$mat_q0[samp_subset_k, , k, drop = FALSE]

    cf_est_pi[samp_subset_k] <- nuis$mat_pi[samp_subset_k, , k, drop = TRUE]
  }

  #####################################################################
  # Step 2:                                                           #
  # Use TMLE to update the conditional mean regressions such that the #
  # empirical mean of the estimated influence functions are 0         #
  #####################################################################

  est_AIPW0_tmle <- est_AIPW1_tmle <- rep(NA, J)

  mhat_1_update <- matrix(NA, nrow = n, ncol = J,
                          dimnames = c(list(paste0("samp", 1:n)),
                                       list(paste0("tax", 1:J))))
  qhat_0_update <- qhat_1_update <- mhat_0_update <- mhat_1_update

  if (verbose %in% c(TRUE, "development")) {
    cat("Beginning TMLE\n")
    pb <- txtProgressBar(min = 0, max = J, style = 3)
  }

  tmle_perturbation <- matrix(NA, nrow = J, ncol = 6,
                              dimnames = c(list(paste0("tax", 1:J)),
                                           list(c("m1", "m0", "q1", "q0", "zstat1", "zstat0"))))

  for (j in 1:J) {
    ###############################################################################
    # Part A: apply ad-hoc perturbation to conditional mean regression if W_j > 0 #
    ###############################################################################
    mhat_1j <- cf_est_m1[, j]
    mhat_0j <- cf_est_m0[, j]
    mhat_Aj <- ifelse(A == 1, mhat_1j, mhat_0j)

    # shift falsifiable estimates of E[W_j|W_j>0,A,X]
    which_shift_m <- W[, j] > 0 & mhat_Aj == 0
    if(any(which_shift_m)){message("shift required")}
    shift_m <- min(c(c(W[W[, j] > 0, j]),
                   c(nuis$mat_m0[ , j, ])[c(nuis$mat_m0[ , j, ]) > 0],
                   c(nuis$mat_m1[ , j, ])[c(nuis$mat_m1[ , j, ]) > 0]))
    mhat_1j[which_shift_m] <- mhat_1j[which_shift_m] + shift_m
    mhat_0j[which_shift_m] <- mhat_0j[which_shift_m] + shift_m
    mhat_Aj[which_shift_m] <- shift_m

    if (tmle_fluctuation == "poisson") {

      m_epsilon <- coef(glm(formula = clever.resp ~ -1 + offset(clever.offset) + clever.covar1 + clever.covar2,
                            weights = clever.weight,
                            subset = clever.subset,
                            data =
                              data.frame(clever.resp   = W[, j],
                                         clever.weight = (1 / ifelse(A == 1, cf_est_pi, 1 - cf_est_pi)) /
                                           mean(1 / ifelse(A == 1, cf_est_pi, 1 - cf_est_pi)),
                                         clever.offset = log(mhat_Aj),
                                         clever.covar1 = A,
                                         clever.covar2 = 1 - A,
                                         clever.subset = W[, j] > 0),
                            family = quasipoisson(link = "log"),
                            control = glm.control(epsilon = 1e-10,
                                                  maxit = 1e3)))

      mhat_1j_update <- mhat_1j * exp(m_epsilon[1] * A)
      mhat_0j_update <- mhat_0j * exp(m_epsilon[2] * (1 - A))
      mhat_Aj_update <- ifelse(A == 1, mhat_1j_update, mhat_0j_update)

    } else if (tmle_fluctuation == "logistic") {

      ub <- max(c(W[, j], mhat_0j, mhat_1j)) * 1.01
      m_epsilon <- coef(glm(formula = clever.resp ~ -1 + offset(clever.offset) + clever.covar1 + clever.covar2,
                            weights = clever.weight,
                            subset = clever.subset,
                            data =
                              data.frame(clever.resp   = W[, j] / ub,
                                         clever.weight = (1 / ifelse(A == 1, cf_est_pi, 1 - cf_est_pi)) /
                                           mean(1 / ifelse(A == 1, cf_est_pi, 1 - cf_est_pi)),
                                         clever.offset = qlogis(mhat_Aj / ub),
                                         clever.covar1 = A,
                                         clever.covar2 = 1 - A,
                                         clever.subset = W[, j] > 0),
                            family = quasibinomial(link = "logit"),
                            control = glm.control(epsilon = 1e-10,
                                                  maxit = 1e3)))

      m_epsilon <- ifelse(is.na(m_epsilon), 0, m_epsilon)
      mhat_1j_update <- plogis(qlogis(mhat_1j / ub) + m_epsilon[1] * A) * ub
      mhat_0j_update <- plogis(qlogis(mhat_0j / ub) + m_epsilon[2] * (1 - A)) * ub
      mhat_Aj_update <- ifelse(A == 1, mhat_1j_update, mhat_0j_update)

    } else {
      stop('Invalid choice for TMLE fluctuation, it must be either "poisson" or "logistic".')
    }

    ####################################################################
    # Part B: apply perturbation to full iterated expectation estimate #
    ####################################################################

    if (all(W[, j] > 0)) {
      qhat_1j_update <- rep(1, n)
      qhat_0j_update <- rep(1, n)
      qhat_Aj_update <- rep(1, n)
      q_epsilon <- c(0, 0)
    } else {
      qhat_1j <- cf_est_q1[, j]
      qhat_0j <- cf_est_q0[, j]
      qhat_Aj <- ifelse(A == 1, qhat_1j, qhat_0j)

      # truncate falsifiable estimates of P(W_j>0|A,X)
      which_trunc_q <- (W[, j] > 0 & qhat_Aj == 0) | (W[, j] == 0 & qhat_Aj == 1)
      if(any(which_trunc_q)){message("trunc required")}
      qpreds <- c(c(nuis$mat_q0[ , j, ]), c(nuis$mat_q1[ , j, ]))
      trunc_q <- min(c(c(qpreds[qpreds > 0], 1 - qpreds[qpreds < 1]),
                       1 - mean(W[, j] > 0), mean(W[, j] > 0),
                       1 / sum(which_trunc_q) / 2))
      qhat_1j[which_trunc_q] <- pmax(pmin(qhat_1j[which_trunc_q], 1 - trunc_q), trunc_q)
      qhat_0j[which_trunc_q] <- pmax(pmin(qhat_0j[which_trunc_q], 1 - trunc_q), trunc_q)
      qhat_Aj <- ifelse(A == 1, qhat_1j, qhat_0j)

      q_epsilon <- coef(glm(formula = clever.resp ~ -1 + offset(clever.offset) + clever.covar1 + clever.covar2,
                            weights = clever.weight,
                            subset = clever.subset,
                            data =
                              data.frame(clever.resp   = ifelse(W[, j] > 0, 1, 0),
                                         clever.weight = (mhat_Aj_update / ifelse(A == 1, cf_est_pi, 1 - cf_est_pi)) /
                                           mean(mhat_Aj_update / ifelse(A == 1, cf_est_pi, 1 - cf_est_pi)),
                                         clever.offset = qlogis(qhat_Aj),
                                         clever.covar1 = A,
                                         clever.covar2 = 1 - A,
                                         clever.subset = !(qhat_Aj == 1 | qhat_Aj == 0)),
                            family = quasibinomial(link = "logit"),
                            control = glm.control(epsilon = 1e-10,
                                                  maxit = 1e3)))

      # if perturbation coefficient is NA, make no change
      q_epsilon <- ifelse(is.na(q_epsilon), 0, q_epsilon)
      qhat_1j_update <- plogis(qlogis(qhat_1j) + q_epsilon[1] * A)
      qhat_0j_update <- plogis(qlogis(qhat_0j) + q_epsilon[2] * (1 - A))
    }

    if (any(abs(c(mean((A == 0) * (W[, j] - mhat_0j_update * qhat_0j_update) / (1 - cf_est_pi)),
                  mean((A == 1) * (W[, j] - mhat_1j_update * qhat_1j_update) / cf_est_pi))) > 1e-2)) {
      warning(paste("The TMLE fluctuation in taxon", j ,"seems to have poor fit.\n",
                    "A = 0:", mean(((1 - A) / (1 - cf_est_pi)) * (W[, j] - mhat_0j_update * qhat_0j_update)),"\n",
                    "A = 1:", mean((A / cf_est_pi) * (W[, j] - mhat_1j_update * qhat_1j_update)),"\n"))
    }

    # mhat_1j_update <- cf_est_m1[, j]
    # mhat_0j_update <- cf_est_m0[, j]
    # qhat_1j_update <- cf_est_q1[, j]
    # qhat_0j_update <- cf_est_q0[, j]
    est_AIPW1_tmle[j] <- mean(mhat_1j_update * qhat_1j_update) # + mean((A / cf_est_pi) * (W[, j] - mhat_1j_update * qhat_1j_update))
    est_AIPW0_tmle[j] <- mean(mhat_0j_update * qhat_0j_update) # + mean(((1 - A) / (1 - cf_est_pi)) * (W[, j] - mhat_0j_update * qhat_0j_update))

    z1 <- sqrt(n) * mean((A / cf_est_pi) * (W[, j] - cf_est_m1[, j] * cf_est_q1[, j])) /
      sd((A / cf_est_pi) * (W[, j] - cf_est_m1[, j] * cf_est_q1[, j]))
    z0 <- sqrt(n) * mean(((1 - A) / (1 - cf_est_pi)) * (W[, j] - cf_est_m0[, j] * cf_est_q0[, j])) /
      sd(((1 - A) / (1 - cf_est_pi)) * (W[, j] - cf_est_m0[, j] * cf_est_q0[, j]))
    tmle_perturbation[j, ] <- c(m_epsilon, q_epsilon, z1, z0)

    mhat_1_update[, j] <- mhat_1j_update
    qhat_1_update[, j] <- qhat_1j_update

    mhat_0_update[, j] <- mhat_0j_update
    qhat_0_update[, j] <- qhat_0j_update

    # update progress bar
    if (verbose %in% c(TRUE, "development")) {
      setTxtProgressBar(pb, j)
    }
  }
  cat("\n")

  # plot(x = tmle_perturbation[, "m0"], y = tmle_perturbation[, "m1"],
  #      xlab = "epsilon_m0", ylab = "epsilon_m1",
  #      xlim = range(c(tmle_perturbation[, c("m0", "m1")])),
  #      ylim = range(c(tmle_perturbation[, c("m0", "m1")]))); abline(h = 0); abline(v = 0)

  # plot(x = tmle_perturbation[, "q0"], y = tmle_perturbation[, "q1"],
  #      xlab = "epsilon_q0", ylab = "epsilon_q1",
  #      xlim = c(-1, 1) * max(abs(range(c(tmle_perturbation[, c("q0", "q1")])))),
  #      ylim = c(-1, 1) * max(abs(range(c(tmle_perturbation[, c("q0", "q1")]))))); abline(h = 0); abline(v = 0)

  # hist(tmle_perturbation[, "zstat1"], breaks = J/2, xlim = c(-1, 1) * max(abs(tmle_perturbation[, "zstat1"])))
  # hist(tmle_perturbation[, "zstat0"], breaks = J/2, xlim = c(-1, 1) * max(abs(tmle_perturbation[, "zstat0"])))

  ##################################################################
  # Step 2:                                                        #
  # Apply Delta Method to learn log(E[E[W|A=1,X]] / E[E[W|A=0,X]]) #
  ##################################################################

  # construct the TMLE estimator
  est_log_AIPW <- log(est_AIPW1_tmle) - log(est_AIPW0_tmle)

  # calculate the influence functions evaluated at each point using TMLE estimator
  if_AIPW1_tmle <- t(t((A / cf_est_pi) * (W - mhat_1_update * qhat_1_update) +
                         mhat_1_update * qhat_1_update) - est_AIPW1_tmle)

  if_AIPW0_tmle <- t(t((1 - A) / (1 - cf_est_pi) * (W - mhat_0_update * qhat_0_update) +
                         mhat_0_update * qhat_0_update) - est_AIPW0_tmle)

  if_log_AIPW_tmle <- t(t(if_AIPW1_tmle) / est_AIPW1_tmle) - t(t(if_AIPW0_tmle) / est_AIPW0_tmle)

  Sigmahat_log_AIPW <- cov(if_log_AIPW_tmle)

  ###########################################################################
  # Step 3:                                                                 #
  # Apply Delta Method to learn log(E[E[W|A=1,X]] / E[E[W|A=0,X]]) - g(...) #
  ###########################################################################

  est_g_log_AIPW <- est_log_AIPW - niceday::pseudohuber_center(x = est_log_AIPW, d = d)

  est_grad_g_log_AIPW <- niceday::dpseudohuber_center_dx(x = est_log_AIPW, d = d)

  if_g_log_AIPW_tmle <- if_log_AIPW_tmle - (if_log_AIPW_tmle %*% est_grad_g_log_AIPW) %*% rep(1, J)

  Sigmahat_g_log_AIPW <- cov(if_g_log_AIPW_tmle)

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
                                      Sigma = cor(if_log_AIPW_tmle)), 1,
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
                                      Sigma = cor(if_g_log_AIPW_tmle)), 1,
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
                    type = "TMLE",
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
                      type = "TMLE",
                      row.names = NULL)

  return(list(res = res,
              res_g = res_g,
              cf_nuis = list(cf_est_pi = cf_est_pi,
                             cf_est_m0 = cf_est_m0,
                             cf_est_m1 = cf_est_m1,
                             cf_est_q0 = cf_est_q0,
                             cf_est_q1 = cf_est_q1),
              cf_if = list(if_AIPW1_tmle = if_AIPW1_tmle,
                           if_AIPW0_tmle = if_AIPW0_tmle,
                           if_log_AIPW_tmle = if_log_AIPW_tmle,
                           if_g_log_AIPW_tmle = if_g_log_AIPW_tmle),
              Sigmahat = list(Sigmahat_log_AIPW = Sigmahat_log_AIPW,
                              Sigmahat_g_log_AIPW = Sigmahat_g_log_AIPW)))
}
