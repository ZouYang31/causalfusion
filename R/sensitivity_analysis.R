# ------------------------------------------------------------------------------#
#' Sensitivity Analysis for Algorithm Hyperparameters (Constraint Slack)
#'
#' This function performs a sensitivity analysis over the algorithm's hyperparameters \code{eta_X} and \code{eta_Z},
#' which control the tolerance for the constraints in the optimization procedure. Specifically, these
#' parameters bind the error for covariates in the target domain (\code{X}) and
#' the reference domain (\code{Z}).
#'
#' @param eta_values A numeric vector of candidate values for both \code{eta_X} and \code{eta_Z}.
#'        The function evaluates all combinations via a grid search.
#' @param B A list of combinations of a budgeting vector.
#' @param F A numeric matrix of outcomes in the reference domain. The first row is the treated unit.
#' @param Y A numeric matrix of outcomes in the target domain.
#' @param X A numeric matrix of covariates in the target domain.
#' @param Z A numeric matrix of covariates in the reference domain.
#' @param t_max Integer. Number of time periods in reference domain.
#' @param dr Integer. Number of covariates in \code{Z}.
#' @param dt Integer. Number of covariates in \code{X}.
#' @param i_max Integer. Total number of units (treated + control).
#' @param NSE_Z_baseline Numeric. Baseline normalized squared error for \code{Z}.
#' @param NSE_X_baseline Numeric. Baseline normalized squared error for \code{X}.
#'
#' @return A data frame with one row per \code{(eta_X, eta_Z)} combination, including:
#' \describe{
#'   \item{eta_X}{Hyperparameter controlling constraint slack in the target domain (\code{X}).}
#'   \item{eta_Z}{Hyperparameter controlling constraint slack in the reference domain (\code{Z}).}
#'   \item{diff_outcome_ref}{Absolute difference between observed and synthetic means in the reference domain.}
#'   \item{diff_outcome_tar}{Absolute difference between observed and synthetic means in the target domain.}
#'   \item{ett_synthetic}{Estimated treatment effect on the treated (ETT) from the synthetic control.}
#' }
#'
#' @details
#' The hyperparameters \code{eta_X} and \code{eta_Z} are introduced as algorithmic constants that define
#' the allowable slack in the covariate reconstruction constraints. By varying these parameters,
#' the function quantifies how sensitive the synthetic control results are to the tightness of the
#' optimization constraints.
#'
#' @seealso \code{\link{find_best_B}}, \code{\link{sc_ett}}
#'
#' @examples
#' \dontrun{
#' eta_vals <- c(0.01, 0.05, 0.1, 0.2)
#' result <- sensitivity_param(
#'   eta_values = eta_vals,
#'   F = F_matrix,
#'   Y = Y_matrix,
#'   X = X_matrix,
#'   Z = Z_matrix,
#'   t_max = 10,
#'   dr = 4,
#'   dt = 5,
#'   i_max = 20,
#'   NSE_Z_baseline = 0.05,
#'   NSE_X_baseline = 0.07
#' )
#' }
#'
#' @export

sensitivity_param <- function(eta_values, B,
                              F, Y, X, Z,
                              t_max, dr, dt, i_max,
                              NSE_Z_baseline,
                              NSE_X_baseline){

  # Check eta_values
  if (!is.numeric(eta_values) || length(eta_values) == 0) {
    stop("'eta_values' must be a non-empty numeric vector.")
  }

  # Check matrix dimensions
  if (!is.matrix(F) || !is.matrix(Y) || !is.matrix(X) || !is.matrix(Z)) {
    stop("F, Y, X, and Z must all be numeric matrices.")
  }

  message("Starting sensitivity analysis over eta_X and eta_Z combinations...")

  eta_combinations <- expand.grid(eta_X = eta_values, eta_Z = eta_values)

  # Initialize a dataframe to store results
  results <- eta_combinations %>%
    mutate(
      diff_outcome_ref = NA,
      diff_outcome_tar = NA,
      ett_synthetic = NA
    )

  J <- i_max - 1
  F_treated <- F[1, ]
  F_control <- F[2:(J+1), ]
  Z_treated <- Z[1, ]
  Z_control <- Z[2:(J+1), ]
  X_treated <- X[1, ]
  X_control <- X[2:(J+1), ]
  Y_treated <- Y[1, ]
  Y_control <- Y[2:(J+1), ]


  # Loop over all 16 combinations of eta_X and eta_Z
  for (i in seq_len(nrow(results))) {
    eta_X <- results$eta_X[i]
    eta_Z <- results$eta_Z[i]

    # Call the function with different eta_X and eta_Z values
    b_list <- find_best_B(
      B = B,
      F_treated = F_treated,
      F_control = F_control,
      X_treated = X_treated,
      X_control = X_control,
      Z_treated = Z_treated,
      Z_control = Z_control,
      t_max = t_max,
      dr = dr,
      dt = dt,
      i_max = i_max,
      NSE_Z_baseline = NSE_Z_baseline,
      NSE_X_baseline = NSE_X_baseline,
      eta_Z = eta_Z,
      eta_X = eta_X
    )

    # Solve for weight vector (w) using the optimal B values
    w <- b_list$best_w_star

    # Sensitivity matrix
    diff_outcome_ref <- abs(mean(F_treated) - mean(w %*% F_control))

    diff_outcome_tar <- abs(mean(Y_treated) - mean(w %*% Y_control))

    # Causal effect
    ett_synthetic <- sc_ett(Y_treated, Y = Y, w = w, J = J)

    # Store results
    results$diff_outcome_ref[i] <- round(diff_outcome_ref, 3)
    results$diff_outcome_tar[i] <- round(diff_outcome_tar, 3)

    results$ett_synthetic[i] <- round(ett_synthetic, 3)
  }


  return(results)
}
