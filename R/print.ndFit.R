#' Print function
#'
#' @param x Object of class \code{ndFit}
#' @param n The number of coefficient estimates to be printed (ordered by largest absolute value to smallest)
#' @param ... No optional arguments are accepted at this time.
#'
#' @return \code{NULL}. Displays printed model summary.
#'
#' @method print ndFit
#'
#' @export
print.ndFit <- function(x,
                        n = 20,
                        ...) {

  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("\nCoefficient estimates with the largest magnitudes (Psi1/Psi2):\n")
  n_or_nrow <- min(n, nrow(x$coef[[1]]))
  sorted_ind <- order(abs(x$coef[[1]]$est), decreasing = TRUE)[1:n_or_nrow]
  coefs_tab <- x$coef[[1]][sorted_ind,, drop = FALSE]
  cols_NA <- which(colMeans(is.na(x$coef[[1]])) == 1)
  if (length(cols_NA) > 0) {
    coefs_tab <- coefs_tab[, -(cols_NA)]
  }
  print(data.frame(coefs_tab))

  cat("\n")

  cat("\nCoefficient estimates with the largest magnitudes (Psi1g/Psi2g):\n")
  n_or_nrow <- min(n, nrow(x$coef[[2]]))
  sorted_ind <- order(abs(x$coef[[2]]$est), decreasing = TRUE)[1:n_or_nrow]
  coefs_tab <- x$coef[[2]][sorted_ind,, drop = FALSE]
  cols_NA <- which(colMeans(is.na(x$coef[[2]])) == 1)
  if (length(cols_NA) > 0) {
    coefs_tab <- coefs_tab[, -(cols_NA)]
  }
  print(data.frame(coefs_tab))

  cat("\n")

  message("To obtain the entire coefficient table, use the command `ndFit_object[[1]]$coef` for the Psi1 or Psi2 parameter and `ndFit_object[[2]]$coef` for the Psi1g or Psi2g parameter.")
}
