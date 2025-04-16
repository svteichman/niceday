#' @export

true_nuis <- function(W,
                      A,
                      X,
                      muV_AX,
                      muW_AX,
                      piA_X,
                      drawV,
                      num_crossfit_folds = 10) {

  n <- nrow(W)
  J <- ncol(W)
  nfold <- num_crossfit_folds
  if (nfold < 2) {stop("The nfold argument must be an integer greater than or equal to 2.")}

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

  ###################################
  ### Calculate propensity score: ###
  ###  P(A=1|X)                   ###
  ###################################

  for (i in 1:n) {
    for (k in 1:nfold) {
      mat_pi[i, 1, k] <- piA_X(x = X[i, ])
    }
  }

  pzerovec <- withr::with_seed(9, sample(seq(0.1, 0.9, length = J)))
  if (drawV == "zinb") {
    q_oracle <- function(a, x, j) {
      return(
        1 - (pzerovec[j] + (2 / (2 + muV_AX(a = a, x = x, j) / (1 - pzerovec[j])))^2 * (1 - pzerovec[j]))
      )
    }
    m_oracle <- function(a, x, j) {
      return(
        muW_AX(a, x, j) / q_oracle(a, x, j)
      )
    }
  } else if (drawV == "poisson") {
    q_oracle <- function(a, x, j) {
      return(
        1 - exp(-muV_AX(a, x, j))
      )
    }
    m_oracle <- function(a, x, j) {
      return(
        muW_AX(a, x, j) / q_oracle(a, x, j)
      )
    }
  } else {
    stop("Invalid input for drawV")
  }

  for (j in 1:J) {
    for (i in 1:n) {
      for (k in 1:nfold) {
        mat_m0[i, j, k] <- m_oracle(a = 0, x = X[i, ], j = j)
        mat_m1[i, j, k] <- m_oracle(a = 1, x = X[i, ], j = j)
        mat_q0[i, j, k] <- q_oracle(a = 0, x = X[i, ], j = j)
        mat_q1[i, j, k] <- q_oracle(a = 1, x = X[i, ], j = j)
      }
    }
  }

  nuis <- list(mat_pi = mat_pi,
               mat_m0 = mat_m0,
               mat_m1 = mat_m1,
               mat_q0 = mat_q0,
               mat_q1 = mat_q1,
               fold_list = rep(list(1:n), num_crossfit_folds))

  return(nuis)
}
