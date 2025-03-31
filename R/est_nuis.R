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
                                  "SL.xgboost.qpois"),
                     sl.lib.q = sl.lib.pi,
                     enforce_pos_reg = FALSE,
                     quiet = FALSE) {
  require("SuperLearner")

  n <- nrow(W)
  J <- ncol(W)
  nfold <- num_crossfit_folds
  if (nfold < 2) {stop("The nfold argument must be an integer greater than or equal to 2.")}
  fold_list <- lapply(split(sample(1:n, size = n, replace = FALSE),
                            cut(1:n, breaks = nfold, labels = FALSE)), sort)

  mat_pi <- array(NA,
                 dim = list(n, 1, nfold),
                 dimnames = c(list(paste0("samp", 1:n)),
                              list("notax"),
                              list(paste0("fold", 1:nfold))))

  mat_m1 <- array(NA,
                   dim = list(n, J, nfold),
                   dimnames = c(list(paste0("samp", 1:n)),
                                list(paste0("tax", 1:J)),
                                list(paste0("fold", 1:nfold))))
  mat_q0 <- mat_q1 <- mat_m0 <- mat_m1

  if(!quiet){cat("Beginning propensity score estimation\n")}

  ##################################
  ### Estimate propensity score: ###
  ###  P(A=1|X)                  ###
  ##################################

  for (k in 1:nfold) {
    # identify which samples are used to estimate the propensity score
    samp_subset <- 1:n %in% unname(fold_list[[k]])
    samp_subset_comp <- 1:n %in% sort(unname(unlist(fold_list[-k])))

    est_pi <- tryCatch(expr =
      c(SuperLearner(Y = A[samp_subset_comp],
                     X = setNames(object = data.frame(X[samp_subset_comp, , drop = FALSE]),
                                  nm = paste0("X", 1:ncol(X))),
                     newX = setNames(object = data.frame(X),
                                     nm = paste0("X", 1:ncol(X))),
                     method = "method.NNLS",
                     family = binomial(link = "logit"),
                     SL.library = sl.lib.pi,
                     cvControl = list(V = min(num_crossval_folds,
                                              sum(A[samp_subset_comp] == 0),
                                              sum(A[samp_subset_comp] == 1)),
                                      stratifyCV = TRUE))$SL.predict),
    error = function(e){
      #cat("An error occurred:\n", e$message, "\n")
      rep(mean(A[samp_subset_comp] == 1), n)
    },
    warning = function(w){
      #cat("A warning occurred:\n", w$message, "\n")
      rep(mean(A[samp_subset_comp] == 1), n)
    })

    # save the results
    mat_pi[1:n, 1, k] <- pmin(pmax(est_pi, gtrunc), 1 - gtrunc)
  }

  ######################################
  ### Estimate conditional mean:     ###
  ###  E[W_j|A=1,X] and E[W_j|A=0,X] ###
  ######################################

  if (!quiet) {
    cat("Beginning conditional mean estimation\n")
    pb <- txtProgressBar(min = 0, max = J, style = 3)
  }

  for (j in 1:J) {
    for (k in 1:nfold) {
      # identify which samples are used to estimate the conditional mean
      samp_subset <- 1:n %in% unname(fold_list[[k]])
      samp_subset_comp <- 1:n %in% sort(unname(unlist(fold_list[-k])))

      # identify subset that satisfies A=a
      samp_subset_comp_q1 <- samp_subset_comp & A == 1
      samp_subset_comp_q0 <- samp_subset_comp & A == 0

      # identify subset that also satisfies W_j>0 and A=a
      samp_subset_comp_m1 <- samp_subset_comp & A == 1 & W[, j] > 0
      samp_subset_comp_m0 <- samp_subset_comp & A == 0 & W[, j] > 0

      # for each taxon j, estimate E[W_j|W_j>0,A=1,X]
      if (sum(samp_subset_comp_m1) > 0) {
        est_m1 <- tryCatch(expr =
                             c(SuperLearner(Y = W[samp_subset_comp_m1, j, drop = TRUE],
                                            X = setNames(object = data.frame(X[samp_subset_comp_m1, , drop = FALSE]),
                                                         nm = paste0("X", 1:ncol(X))),
                                            newX = setNames(object = data.frame(X),
                                                            nm = paste0("X", 1:ncol(X))),
                                            method = "method.NNLS",
                                            family = gaussian(link = "identity"),
                                            SL.library = sl.lib.m,
                                            cvControl = list(V = min(num_crossval_folds,
                                                                     sum(samp_subset_comp_m0),
                                                                     sum(samp_subset_comp_m1))))$SL.predict),
                           error = function(e){
                             #cat("An error occurred:\n", e$message, "\n")
                             rep(mean(W[samp_subset_comp_m1, j]), n)
                           },
                           warning = function(w){
                             #cat("A warning occurred:\n", w$message, "\n")
                             rep(mean(W[samp_subset_comp_m1, j]), n)
                           })
      } else {
        warning(paste0("Cross-fitting for E[W_", j,"|W_", j,">0,A=1,X] is impossible.\n",
                       "  Reverted to overall mean(W[A==1 & W[, j] > 0, j]).\n"))
        est_m1 <- rep(mean(W[A == 1 & W[, j] > 0, j]), n)
      }

      # for each taxon j, estimate P(W_j>0|A=1,X)
      if (sum(samp_subset_comp_q1) > 0) {
        est_q1 <- tryCatch(expr =
                             c(SuperLearner(Y = ifelse(W[samp_subset_comp_q1, j] > 0, 1, 0),
                                            X = setNames(object = data.frame(X[samp_subset_comp_q1, , drop = FALSE]),
                                                         nm = paste0("X", 1:ncol(X))),
                                            newX = setNames(object = data.frame(X),
                                                            nm = paste0("X", 1:ncol(X))),
                                            method = "method.NNLS",
                                            family = binomial(link = "logit"),
                                            SL.library = sl.lib.q,
                                            cvControl = list(V = min(num_crossval_folds,
                                                                     sum(samp_subset_comp_q0),
                                                                     sum(samp_subset_comp_q1))))$SL.predict),
                           error = function(e){
                             #cat("An error occurred:\n", e$message, "\n")
                             rep(mean(W[samp_subset_comp_q1, j] > 0), n)
                           },
                           warning = function(w){
                             #cat("A warning occurred:\n", w$message, "\n")
                             rep(mean(W[samp_subset_comp_q1, j] > 0), n)
                           })
      } else {
        warning(paste0("Cross-fitting for P(W_", j, ">0|A=1,X) is impossible.\n",
                       "  Reverted to overall mean(W[A == 1, j] > 0).\n"))
        est_q1 <- rep(mean(W[A == 1, j] > 0), n)
      }


      # for each taxon j, estimate E[W_j|W_j>0,A=0,X]
      if (sum(samp_subset_comp_m0) > 0) {
        est_m0 <- tryCatch(expr =
           c(SuperLearner(Y = W[samp_subset_comp_m0, j, drop = TRUE],
                          X = setNames(object = data.frame(X[samp_subset_comp_m0, , drop = FALSE]),
                                       nm = paste0("X", 1:ncol(X))),
                          newX = setNames(object = data.frame(X),
                                          nm = paste0("X", 1:ncol(X))),
                          method = "method.NNLS",
                          family = gaussian(link = "identity"),
                          SL.library = sl.lib.m,
                          cvControl = list(V = min(num_crossval_folds,
                                                   sum(samp_subset_comp_m0),
                                                   sum(samp_subset_comp_m1))))$SL.predict),
         error = function(e){
           #cat("An error occurred:\n", e$message, "\n")
           rep(mean(W[samp_subset_comp_m0, j]), n)
         },
         warning = function(w){
           #cat("A warning occurred:\n", w$message, "\n")
           rep(mean(W[samp_subset_comp_m0, j]), n)
        })
      } else {
        warning(paste0("Cross-fitting for E[W_", j,"|W_", j,">0,A=0,X] is impossible.\n",
                       "  Reverted to overall mean(W[A==0 & W[, j] > 0, j]).\n"))
        est_m0 <- rep(mean(W[A == 0 & W[, j] > 0, j]), n)
      }

      # for each taxon j, estimate P(W_j>0|A=0,X)
      if (sum(samp_subset_comp_q0) > 0) {
        est_q0 <- tryCatch(expr =
           c(SuperLearner(Y = ifelse(W[samp_subset_comp_q0, j] > 0, 1, 0),
                          X = setNames(object = data.frame(X[samp_subset_comp_q0, , drop = FALSE]),
                                       nm = paste0("X", 1:ncol(X))),
                          newX = setNames(object = data.frame(X),
                                          nm = paste0("X", 1:ncol(X))),
                          method = "method.NNLS",
                          family = binomial(link = "logit"),
                          SL.library = sl.lib.q,
                          cvControl = list(V = min(num_crossval_folds,
                                                   sum(samp_subset_comp_q0),
                                                   sum(samp_subset_comp_q1))))$SL.predict),
         error = function(e){
           #cat("An error occurred:\n", e$message, "\n")
           rep(mean(W[samp_subset_comp_q0, j] > 0), n)
         },
         warning = function(w){
           #cat("A warning occurred:\n", w$message, "\n")
           rep(mean(W[samp_subset_comp_q0, j] > 0), n)
        })
      } else {
        warning(paste0("Cross-fitting for P(W_", j, ">0|A=0,X) is impossible.\n",
                       "  Reverted to overall mean(W[A == 0, j] > 0).\n"))
        est_q0 <- rep(mean(W[A == 0, j] > 0), n)
      }

      # perform truncation to avoid numerical instability
      est_m1 <- pmin(pmax(est_m1, 0), max(W[A == 1 & W[, j] > 0, j]))
      est_m0 <- pmin(pmax(est_m0, 0), max(W[A == 0 & W[, j] > 0, j]))
      est_q1 <- pmin(pmax(est_q1, 0), 1)
      est_q0 <- pmin(pmax(est_q0, 0), 1)

      # save estimates
      mat_m1[1:n, j, k] <- as.numeric(est_m1)
      mat_q1[1:n, j, k] <- as.numeric(est_q1)

      mat_m0[1:n, j, k] <- as.numeric(est_m0)
      mat_q0[1:n, j, k] <- as.numeric(est_q0)
    }

    # update progress bar
    if (!quiet) {
      setTxtProgressBar(pb, j)
    }
  }

  cat("\n")

  nuis <- list(mat_pi = mat_pi,
               mat_m0 = mat_m0,
               mat_m1 = mat_m1,
               mat_q0 = mat_q0,
               mat_q1 = mat_q1,
               fold_list = fold_list)

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
