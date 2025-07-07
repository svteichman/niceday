#' @export

est_nuis <- function(W,
                     A,
                     X,
                     num_crossfit_folds = 10,
                     num_crossval_folds = 10,
                     gtrunc = 0.01,
                     sl.lib.pi = c("SL.mean",
                                   "SL.lm",
                                   "SL.glm.binom",
                                   "SL.xgboost.binom"),
                     sl.lib.m = c("SL.mean",
                                  "SL.lm",
                                  "SL.glm.qpois",
                                  "SL.xgboost.pois"),
                     sl.lib.q = sl.lib.pi,
                     allow_warnings = TRUE,
                     enforce_pos_reg = FALSE,
                     verbose = FALSE) {
  require("SuperLearner")

  n <- nrow(W)
  J <- ncol(W)

  # number of folds used for cross-fitting
  nfold <- num_crossfit_folds
  # if one fold, skip cross-fitting
  if (nfold == 1) {
    fold_list <- list(1:n)
  # otherwise, use n-fold cross-fitting
  } else if (nfold > 1) {
    fold_list <- lapply(split(sample(1:n, size = n, replace = FALSE),
                              cut(1:n, breaks = nfold, labels = FALSE)), sort)
  } else {
    stop("num_crossfit_folds must be a natural number")
  }

  mat_pi <- array(NA,
                 dim = list(n, 1, nfold),
                 dimnames = c(list(paste0("samp", 1:n)),
                              list("notax"),
                              list(paste0("fold", 1:nfold))))
  check_pi <- rep(NA, nfold)

  mat_m1 <- array(NA,
                   dim = list(n, J, nfold),
                   dimnames = c(list(paste0("samp", 1:n)),
                                list(paste0("tax", 1:J)),
                                list(paste0("fold", 1:nfold))))
  mat_q0 <- mat_q1 <- mat_m0 <- mat_m1
  check_m1 <- array(NA,
                    dim = list(J, nfold),
                    dimnames = c(list(paste0("tax", 1:J)),
                                 list(paste0("fold", 1:nfold))))
  check_q0 <- check_q1 <- check_m0 <- check_m1

  if(verbose %in% c(TRUE, "development")){cat("Beginning propensity score estimation\n")}

  ##################################
  ### Estimate propensity score: ###
  ###  P(A=1|X)                  ###
  ##################################

  for (k in 1:nfold) {
    # identify which samples are used to estimate the propensity score
    samp_subset <- 1:n %in% unname(fold_list[[k]])
    if (nfold == 1) {
      samp_subset_comp <- rep(TRUE, n)
    } else {
      samp_subset_comp <- 1:n %in% sort(unname(unlist(fold_list[-k])))
    }

    fit_pi <- function(SL.library, subset_id, method) {
      #suppressMessages(
        withCallingHandlers(expr = {
          out <- c(SuperLearner(Y = A[subset_id],
                         X = setNames(object = data.frame(X[subset_id, , drop = FALSE]),
                                      nm = paste0("X", 1:ncol(X))),
                         newX = setNames(object = data.frame(X),
                                         nm = paste0("X", 1:ncol(X))),
                         method = method,
                         family = binomial(link = "logit"),
                         SL.library = SL.library,
                         cvControl = list(V = min(num_crossval_folds,
                                                  sum(A[subset_id] == 0),
                                                  sum(A[subset_id] == 1)),
                                          stratifyCV = TRUE))$SL.predict)
        }, warning = function(w) {
          if (allow_warnings) {invokeRestart("muffleWarning")} else {stop(conditionMessage(w))}
        })
      #)
      return(out)
    }

    # use fall-back learners if SuperLearner fails
    est_pi <- tryCatch(
      # Attempt 1:
      # try to fit cross-fitted, user-specified SuperLearner models
      {flag_pi <<- 1 # ; cat("Attempt 1\n")
       fit_pi(SL.library = sl.lib.pi,
              subset_id = samp_subset_comp,
              method = "method.NNloglik")},
      error = function(e1) {
        tryCatch(
          # Attempt 2:
          # try to fit cross-fitted intercept and GLM SuperLearner models
          {flag_pi <<- 2 # ; cat("Attempt 2\n")
           fit_pi(SL.library = c("SL.mean", "SL.glm.binom"),
                  subset_id = samp_subset_comp,
                  method = "method.NNloglik")},
          error = function(e2) {
            tryCatch(
              # Attempt 3:
              # try to fit full-data intercept and GLM SuperLearner models
              {flag_pi <<- 3 # ; cat("Attempt 3\n")
               fit_pi(SL.library = c("SL.mean", "SL.glm.binom"),
                      subset_id = rep(TRUE, n),
                      method = "method.NNloglik")},
              error = function(e3) {
                # Attempt 4:
                # if all else fails, take the full-data empirical mean
                {flag_pi <<- 4 # ; cat("Attempt 4\n")
                 rep(mean(A == 1), n)}
              }
            )
          }
        )
      }
    )

    # save the results
    mat_pi[1:n, 1, k] <- pmin(pmax(est_pi, gtrunc), 1 - gtrunc)
    check_pi[k] <- flag_pi
  }

  ######################################
  ### Estimate conditional mean:     ###
  ###  E[W_j|A=1,X] and E[W_j|A=0,X] ###
  ######################################

  if (verbose %in% c(TRUE, "development")) {
    cat("Beginning conditional mean estimation\n")
    pb <- txtProgressBar(min = 0, max = J, style = 3)
  }

  for (j in 1:J) {
    for (k in 1:nfold) {
      # identify which samples are used to estimate the conditional mean
      samp_subset <- 1:n %in% unname(fold_list[[k]])
      if (nfold == 1) {
        samp_subset_comp <- rep(TRUE, n)
      } else {
        samp_subset_comp <- 1:n %in% sort(unname(unlist(fold_list[-k])))
      }

      fit_m <- function(SL.library, subset_id, method) {
        #suppressMessages(
          withCallingHandlers(expr = {
            out <- c(SuperLearner(Y = W[subset_id, j, drop = TRUE],
                           X = setNames(object = data.frame(X[subset_id, , drop = FALSE]),
                                        nm = paste0("X", 1:ncol(X))),
                           newX = setNames(object = data.frame(X),
                                           nm = paste0("X", 1:ncol(X))),
                           method = method,
                           family = gaussian(link = "identity"),
                           SL.library = SL.library,
                           cvControl = list(V = min(num_crossval_folds,
                                                    sum(subset_id))))$SL.predict)
            if (all(out) == 0) {stop("Estimates cannot all be zero.")}
          }, warning = function(w) {
            if (allow_warnings) {invokeRestart("muffleWarning")} else {stop(conditionMessage(w))}
          })
        #)
        return(out)
      }

      # for each taxon j, estimate P(W_j>0|A=1,X)
      fit_q <- function(SL.library, subset_id, method) {
        #suppressMessages(
          withCallingHandlers(expr = {
            out <- c(SuperLearner(Y = ifelse(W[subset_id, j] > 0, 1, 0),
                           X = setNames(object = data.frame(X[subset_id, , drop = FALSE]),
                                        nm = paste0("X", 1:ncol(X))),
                           newX = setNames(object = data.frame(X),
                                           nm = paste0("X", 1:ncol(X))),
                           method = method,
                           family = binomial(link = "logit"),
                           SL.library = SL.library,
                           cvControl = list(V = min(num_crossval_folds,
                                                    sum(W[subset_id, j] > 0),
                                                    sum(W[subset_id, j] == 0)),
                                            stratifyCV = TRUE))$SL.predict)
            if (all(out) == 0) {stop("Estimates cannot all be zero.")}
          }, warning = function(w) {
            if (allow_warnings) {invokeRestart("muffleWarning")} else {stop(conditionMessage(w))}
          })
        #)
        return(out)
      }

      # for each taxon j, estimate E[W_j|W_j>0,A=1,X]
      est_m1 <- tryCatch(
        # Attempt 1:
        # try to fit cross-fitted, user-specified SuperLearner models
        {flag_m1 <<- 1 # ; cat("Attempt 1\n")
         if (length(unique(W[samp_subset_comp & A == 1 & W[, j] > 0, j])) <= 1) {stop("Skip")}
         fit_m(SL.library = sl.lib.m,
               subset_id = samp_subset_comp & A == 1 & W[, j] > 0,
               method = "method.NNLS")},
        error = function(e1) {
          tryCatch(
            # Attempt 2:
            # try to fit cross-fitted intercept and GLM SuperLearner models
            {flag_m1 <<- 2 # ; cat("Attempt 2\n")
             if (length(unique(W[samp_subset_comp & A == 1 & W[, j] > 0, j])) <= 1) {stop("Skip")}
             fit_m(SL.library = c("SL.mean", "SL.glm.qpois"),
                   subset_id = samp_subset_comp & A == 1 & W[, j] > 0,
                   method = "method.NNLS")},
            error = function(e2) {
              tryCatch(
                # Attempt 3:
                # try to fit full-data intercept and GLM SuperLearner models
                {flag_m1 <<- 3 # ; cat("Attempt 3\n")
                 if (length(unique(W[A == 1 & W[, j] > 0, j])) <= 1) {stop("Skip")}
                 fit_m(SL.library = c("SL.mean", "SL.glm.qpois"),
                       subset_id = c(A == 1 & W[, j] > 0),
                       method = "method.NNLS")},
                 error = function(e3) {
                  # Attempt 4:
                  # if all else fails, take the full-data empirical mean
                  {flag_m1 <<- 4 # ; cat("Attempt 4\n")
                   rep(mean(W[A == 1 & W[, j] > 0, j]), n)}
                }
              )
            }
          )
        }
      )

      # for each taxon j, estimate E[W_j|W_j>0,A=1,X]
      est_m0 <- tryCatch(
        # Attempt 1:
        # try to fit cross-fitted, user-specified SuperLearner models
        {flag_m0 <<- 1 # ; cat("Attempt 1\n")
         if (length(unique(W[samp_subset_comp & A == 0 & W[, j] > 0, j])) <= 1) {stop("Skip")}
         fit_m(SL.library = sl.lib.m,
               subset_id = samp_subset_comp & A == 0 & W[, j] > 0,
               method = "method.NNLS")},
        error = function(e1) {
          tryCatch(
            # Attempt 2:
            # try to fit cross-fitted intercept and GLM SuperLearner models
            {flag_m0 <<- 2 # ; cat("Attempt 2\n")
             if (length(unique(W[samp_subset_comp & A == 0 & W[, j] > 0, j])) <= 1) {stop("Skip")}
             fit_m(SL.library = c("SL.mean", "SL.glm.qpois"),
                   subset_id = samp_subset_comp & A == 0 & W[, j] > 0,
                   method = "method.NNLS")},
            error = function(e2) {
              tryCatch(
                # Attempt 3:
                # try to fit full-data intercept and GLM SuperLearner models
                {flag_m0 <<- 3 # ; cat("Attempt 3\n")
                 if (length(unique(W[A == 0 & W[, j] > 0, j])) <= 1) {stop("Skip")}
                 fit_m(SL.library = c("SL.mean", "SL.glm.qpois"),
                       subset_id = c(A == 0 & W[, j] > 0),
                       method = "method.NNLS")},
                error = function(e3) {
                  # Attempt 4:
                  # if all else fails, take the full-data empirical mean
                  {flag_m0 <<- 4 # ; cat("Attempt 4\n")
                   rep(mean(W[A == 0 & W[, j] > 0, j]), n)}
                }
              )
            }
          )
        }
      )

      # for each taxon j, estimate P(W_j>0|A=1,X)
      est_q1 <- tryCatch(
        # Attempt 1:
        # try to fit cross-fitted, user-specified SuperLearner models
        {flag_q1 <<- 1 # ; cat("Attempt 1\n")
         if (sum(W[samp_subset_comp & A == 1, j]  > 0) <= 1 |
             sum(W[samp_subset_comp & A == 1, j] == 0) <= 1) {stop("Skip")}
         fit_q(SL.library = sl.lib.q,
               subset_id = samp_subset_comp & A == 1,
               method = "method.NNloglik")},
        error = function(e1) {
          tryCatch(
            # Attempt 2:
            # try to fit cross-fitted intercept and GLM SuperLearner models
            {flag_q1 <<- 2 # ; cat("Attempt 2\n")
             if (sum(W[samp_subset_comp & A == 1, j]  > 0) <= 1 |
                 sum(W[samp_subset_comp & A == 1, j] == 0) <= 1) {stop("Skip")}
             fit_q(SL.library = c("SL.mean", "SL.glm.binom"),
                   subset_id = samp_subset_comp & A == 1,
                   method = "method.NNloglik")},
            error = function(e2) {
              tryCatch(
                # Attempt 3:
                # try to fit full-data intercept and GLM SuperLearner models
                {flag_q1 <<- 3 # ; cat("Attempt 3\n")
                 if (sum(W[A == 1, j]  > 0) <= 1 |
                     sum(W[A == 1, j] == 0) <= 1) {stop("Skip")}
                 fit_q(SL.library = c("SL.mean", "SL.glm.binom"),
                       subset_id = c(A == 1),
                       method = "method.NNLS")},
                error = function(e3) {
                  # Attempt 4:
                  # if all else fails, take the full-data empirical mean
                  {flag_q1 <<- 4 # ; cat("Attempt 4\n")
                   rep(mean(W[A == 1, j] > 0), n)}
                }
              )
            }
          )
        }
      )

      # for each taxon j, estimate P(W_j>0|A=0,X)
      est_q0 <- tryCatch(
        # Attempt 1:
        # try to fit cross-fitted, user-specified SuperLearner models
        {flag_q0 <<- 1 # ; cat("Attempt 1\n")
         if (sum(W[samp_subset_comp & A == 0, j]  > 0) <= 1 |
             sum(W[samp_subset_comp & A == 0, j] == 0) <= 1) {stop("Skip")}
         fit_q(SL.library = sl.lib.q,
               subset_id = samp_subset_comp & A == 0,
               method = "method.NNloglik")},
        error = function(e1) {
          tryCatch(
            # Attempt 2:
            # try to fit cross-fitted intercept and GLM SuperLearner models
            {flag_q0 <<- 2 # ; cat("Attempt 2\n")
             if (sum(W[samp_subset_comp & A == 0, j]  > 0) <= 1 |
                 sum(W[samp_subset_comp & A == 0, j] == 0) <= 1) {stop("Skip")}
             fit_q(SL.library = c("SL.mean", "SL.glm.binom"),
                   subset_id = samp_subset_comp & A == 0,
                   method = "method.NNloglik")},
            error = function(e2) {
              tryCatch(
                # Attempt 3:
                # try to fit full-data intercept and GLM SuperLearner models
                {flag_q0 <<- 3 # ; cat("Attempt 3\n")
                 if (sum(W[A == 0, j]  > 0) <= 1 |
                     sum(W[A == 0, j] == 0) <= 1) {stop("Skip")}
                 fit_q(SL.library = c("SL.mean", "SL.glm.binom"),
                       subset_id = c(A == 0),
                       method = "method.NNLS")},
                error = function(e3) {
                  # Attempt 4:
                  # if all else fails, take the full-data empirical mean
                  {flag_q0 <<- 4 # ; cat("Attempt 4\n")
                   rep(mean(W[A == 0, j] > 0), n)}
                }
              )
            }
          )
        }
      )

      # perform truncation to avoid numerical instability
      est_m1 <- pmin(pmax(est_m1, 0), 1.5 * max(W[A == 1, j]))
      est_m0 <- pmin(pmax(est_m0, 0), 1.5 * max(W[A == 0, j]))
      est_q1 <- pmin(pmax(est_q1, 0), 1)
      est_q0 <- pmin(pmax(est_q0, 0), 1)

      # save estimates
      mat_m1[1:n, j, k] <- as.numeric(est_m1)
      check_m1[j, k] <- flag_m1

      mat_q1[1:n, j, k] <- as.numeric(est_q1)
      check_q1[j, k] <- flag_q1

      mat_m0[1:n, j, k] <- as.numeric(est_m0)
      check_m0[j, k] <- flag_m0

      mat_q0[1:n, j, k] <- as.numeric(est_q0)
      check_q0[j, k] <- flag_q0
    }

    # update progress bar
    if (verbose %in% c(TRUE, "development")) {
      setTxtProgressBar(pb, j)
    }
  }

  cat("\n")

  nuis <- list(mat_pi = mat_pi,
               mat_m0 = mat_m0,
               mat_m1 = mat_m1,
               mat_q0 = mat_q0,
               mat_q1 = mat_q1,
               fold_list = fold_list,
               checks = list(check_pi = check_pi,
                             check_m0 = check_m0,
                             check_m1 = check_m1,
                             check_q0 = check_q0,
                             check_q1 = check_q1))

  # warn user about important issues that may arise with nuisance estimates
  if (any(is.na(mat_pi)) > 0) {
    warning("The propensity score estimates contain missing values.")
  }

  if (any(mat_pi == 0 & !is.na(mat_pi)) |
      any(mat_pi == 1 & !is.na(mat_pi))) {
    warning("The propensity score estimates contain predictions of 0 or 1.")
  }

  if (any(mat_pi < 0 & !is.na(mat_pi)) |
      any(mat_pi > 1 & !is.na(mat_pi))) {
    warning("The propensity score estimates contain predictions outside of (0,1).")
  }

  if (any(is.na(mat_m0)) |
      any(is.na(mat_m1)) |
      any(is.na(mat_q0)) |
      any(is.na(mat_q1))) {
    warning("The conditional mean regression estimates contain missing values.")
  }

  if (any(mat_m0 < 0 & !is.na(mat_m0)) |
      any(mat_m1 < 0 & !is.na(mat_m1)) |
      any(mat_q0 < 0 & !is.na(mat_q0)) |
      any(mat_q1 < 0 & !is.na(mat_q1))) {
    warning("The conditional mean regression estimates contain negative values.")
  }

  if (any(mat_m0 == 0 & !is.na(mat_m0)) |
      any(mat_m1 == 0 & !is.na(mat_m1)) |
      any(mat_q0 == 0 & !is.na(mat_q0)) |
      any(mat_q1 == 0 & !is.na(mat_q1))) {
    warning("The conditional mean regression estimates contain zero values.")
  }

  if (any(mat_m0 == Inf & !is.na(mat_m0)) |
      any(mat_m1 == Inf & !is.na(mat_m1)) |
      any(mat_q0 == Inf & !is.na(mat_q0)) |
      any(mat_q1 == Inf & !is.na(mat_q1))) {
    warning("The conditional mean regression estimates contain infinite values.")
  }

  return(nuis)
}
