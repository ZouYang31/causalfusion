##########################################
# Causal Data Fusion Algorithm
#
#########################################


#' Find Optimal Weighting Parameters B for Synthetic Control
#'
#' This function searches over a list of candidate B parameter sets (weightings for domains F, X, and Z)
#' to find the one that yields the lowest loss using the Gurobi-based optimizer.
#'
#' @param B A list of candidate parameter sets. Each element should be a list with named elements \code{b_F}, \code{b_X}, and \code{b_Z}.
#' @param F_treated A numeric vector for the treated unit’s outcome.
#' @param F_control A matrix for the control units' outcomes (rows = units, cols = time points).
#' @param X_treated A numeric vector of covariates for the treated unit in the target domain.
#' @param X_control A matrix of covariates for the control units in the target domain.
#' @param Z_treated A numeric vector of covariates for the treated unit in the reference domain.
#' @param Z_control A matrix of covariates for the control units in the reference domain.
#' @param t_max An integer representing the number of time points for the \code{F} matrix.
#' @param dr An integer representing the number of covariates in \code{Z}.
#' @param dt An integer representing the number of covariates in \code{X}.
#' @param i_max An integer indicating the total number of units (treated + controls).
#' @param NSE_Z_baseline Best possibly normalized squared error for \code{Z}, used to define the quadratic constraint.
#' @param NSE_X_baseline Best possibly normalized squared error for \code{X}, used to define the quadratic constraint.
#' @param eta_Z A numeric parameter for the \code{Z} constraint. Default is \code{0.1}.
#' @param eta_X A numeric parameter for the \code{X} constraint. Default is \code{0.1}.
#'
#' @return A list containing:
#' \describe{
#'   \item{best_B_list}{The parameter set B that achieved the smallest loss.}
#'   \item{least_loss_B}{The minimum loss value achieved across all candidates.}
#'   \item{best_w_star}{The synthetic control weights corresponding to the best \code{B}.}
#' }
#'
#' @details
#' For each candidate \code{B} in the input list, this function calls \code{optimize_w_b_gurobi_B_list()}
#' to compute the optimal weights and loss. The candidate \code{B} that achieves the lowest loss on
#' the outcome dimension \code{F} is selected as the best. This can be useful in hyperparameter
#' tuning or model selection for synthetic control algorithms with domain-weighted objectives.
#'
#' @seealso \code{\link{optimize_w_b_gurobi_B_list}}
#'
#' @examples
#' \dontrun{
#' B_candidates <- list(
#'   list(b_F = 0.4, b_X = 0.3, b_Z = 0.3),
#'   list(b_F = 0.6, b_X = 0.2, b_Z = 0.2)
#' )
#' result <- find_best_B(
#'   B = B_candidates,
#'   F_treated = rnorm(10),
#'   F_control = matrix(rnorm(30), nrow = 3),
#'   X_treated = rnorm(5),
#'   X_control = matrix(rnorm(15), nrow = 3),
#'   Z_treated = rnorm(4),
#'   Z_control = matrix(rnorm(12), nrow = 3),
#'   t_max = 10, dr = 4, dt = 5, i_max = 4,
#'   NSE_Z_baseline = 0.1, NSE_X_baseline = 0.2,
#'   eta_Z = 0.1, eta_X = 0.1
#' )
#' }
#'
#' @export

find_best_B <- function(B, F_treated, F_control, X_treated, X_control, Z_treated, Z_control,
                        t_max, dr, dt, i_max,
                        NSE_Z_baseline, NSE_X_baseline, eta_Z, eta_X) {

  # Initialize variables to store the least loss and corresponding B
  least_loss_B <- Inf
  best_B_list <- NULL
  best_w_star <- NULL

  # Iterate over all combinations in result_list
  for (B_list in B) {
    # Calculate the loss for the current B_list
    result <- optimize_w_b_gurobi_B_list(
      par = B_list,
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
    loss_B <- result$loss.B

    w_star <- result$solution.w

    # Update least loss and best B if the current loss is smaller
    if (!is.na(loss_B) && loss_B < least_loss_B) {
      least_loss_B <- loss_B
      best_B_list <- B_list
      best_w_star <- w_star
    }
  }

  # Return the best B_list and the corresponding least loss as a list
  return(list(best_B_list = best_B_list, least_loss_B = least_loss_B, best_w_star = best_w_star))
}




#' Optimize Synthetic Weights with Quadratic Constraints (Gurobi)
#'
#' This function computes optimal synthetic control weights \code{w} by solving a constrained
#' quadratic optimization problem using the \code{gurobi} solver. The objective minimizes the sum of weighted
#' NSEs (\code{F}, \code{X}, \code{Z}), subject to
#' constraints on the allowable reconstruction error in \code{X} and \code{Z}.
#'
#' @param par A named list of weight parameters: \code{b_F}, \code{b_X}, and \code{b_Z}. These control the contribution of each data domain in the objective.
#' @param F_treated A numeric vector for the treated unit’s outcome.
#' @param F_control A matrix for the control units' outcomes (rows = units, cols = time points).
#' @param X_treated A numeric vector of covariates for the treated unit in the target domain.
#' @param X_control A matrix of covariates for the control units in the target domain.
#' @param Z_treated A numeric vector of covariates for the treated unit in the reference domain.
#' @param Z_control A matrix of covariates for the control units in the reference domain.
#' @param t_max An integer representing the number of time points for the \code{F} matrix.
#' @param dr An integer representing the number of covariates in \code{Z}.
#' @param dt An integer representing the number of covariates in \code{X}.
#' @param i_max An integer indicating the total number of units (treated + controls).
#' @param NSE_Z_baseline Best possibly normalized squared error for \code{Z}, used to define the quadratic constraint.
#' @param NSE_X_baseline Best possibly normalized squared error for \code{X}, used to define the quadratic constraint.
#' @param eta_Z A numeric parameter for the \code{Z} constraint. Default is \code{0.1}.
#' @param eta_X A numeric parameter for the \code{X} constraint. Default is \code{0.1}.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{loss.B}{The loss value from the \code{F} domain, based on the optimized weights.}
#'   \item{solution.w}{A numeric vector of optimized weights \code{w}. If optimization fails, returns \code{NA}.}
#' }
#'
#' @details
#' The function solves a constrained optimization problem:
#' \deqn{w^* = \arg\min_w \left[ b_F \cdot \text{NSE}_F(F_1, w) + b_X \cdot \text{NSE}_X(X_1, w) + b_Z \cdot \text{NSE}_Z(Z_1, w) \right]}
#' subject to quadratic constraints on \code{X} and \code{Z} reconstruction errors.
#' Constraints ensure that the error does not exceed a slack-adjusted baseline error.
#'
#' The optimization uses Gurobi's `quadcon` feature for enforcing quadratic constraints and
#' assumes convexity of the problem. The output is a set of weights \code{w} that minimize
#' the total loss while satisfying the constraints.
#'
#' @import gurobi
#' @examples
#' \dontrun{
#' par <- list(b_F = 0.3, b_X = 0.4, b_Z = 0.3)
#' F_treated <- rnorm(10)
#' F_control <- matrix(rnorm(30), nrow = 3)
#' X_treated <- rnorm(5)
#' X_control <- matrix(rnorm(15), nrow = 3)
#' Z_treated <- rnorm(4)
#' Z_control <- matrix(rnorm(12), nrow = 3)
#' result <- optimize_w_b_gurobi_B_list(
#'   par, F_treated, F_control,
#'   X_treated, X_control,
#'   Z_treated, Z_control,
#'   t_max = 10, dr = 4, dt = 5, i_max = 4,
#'   NSE_Z_baseline = 0.1, NSE_X_baseline = 0.2
#' )
#' }
#'
#' @export
optimize_w_b_gurobi_B_list <- function(par = stop("B missing"),
                                            F_treated = stop("F_treated missing"), F_control = stop("F_control missing"),
                                       X_treated = stop("X missing"), X_control = stop("X missing"),
                                       Z_treated = stop("Z missing"), Z_control = stop("Z missing"),
                                       t_max = stop("t_max missing"),
                                       dr = stop("dr missing"), dt = stop("dt missing"), i_max = stop("i_max missing"),
                                       NSE_Z_baseline = stop("NSE_Z_baseline missing"),
                                       NSE_X_baseline = stop("NSE_X_baseline missing"),
                                       eta_Z = 0.1, eta_X = 0.1) {

  J <- i_max - 1

  # Extract treated and control matrices
  F_treated <- F_treated[1:t_max]
  F_control <- F_control[, 1:t_max]

  # Normalize B so that sum(B) = 1
  b_F <- par$b_F
  b_Z <- par$b_Z
  b_X <- par$b_X

  # Quadratic matrix for the objective function
  H <- b_F/t_max * (F_control %*% t(F_control)) +
    b_Z/dr * (Z_control %*% t(Z_control)) +
    b_X/dt * (X_control %*% t(X_control))

  # Linear term for the objective function
  c <- -2 * (b_F/t_max * t(F_treated) %*% t(F_control) +
               b_Z/dr * t(Z_treated) %*% t(Z_control) +
               b_X/dt * t(X_treated) %*% t(X_control))

  ## -------- constraints -------- ##
  scale_n = 1

  # Quadratic constraints for Z and X
  Q_Z <- Z_control %*% t(Z_control) * scale_n
  L_Z <- -2 * t(matrix(Z_treated)) %*% t(Z_control) * scale_n
  #rhs_Z <- dr *  eta_Z * NSE_Z_baseline - sum(Z_treated^2)
  rhs_Z <- dr * ((1 + eta_Z) * (1 + NSE_Z_baseline * scale_n) - 1) - sum(Z_treated^2) * scale_n

  Q_X <- X_control %*% t(X_control) * scale_n
  L_X <- -2 * t(matrix(X_treated)) %*% t(X_control) * scale_n
  #rhs_X <- dt * eta_X * NSE_X_baseline - sum(X_treated^2)
  rhs_X <- dt * ((1 + eta_X) * (1 + NSE_X_baseline * scale_n) - 1) - sum(X_treated^2) * scale_n

  ## -------- Gurobi Model -------- ##
  model <- list()

  # Objective function (quadratic)
  model$Q <- H
  model$obj <- as.vector(c)
  model$modelsense <- "min"

  # Quadratic constraints (Z and X)
  model$quadcon <- list(
    list(Qc = Q_Z, q = as.vector(L_Z), rhs = rhs_Z),
    list(Qc = Q_X, q = as.vector(L_X), rhs = rhs_X)
  )

  # Linear constraint: sum(w) = 1
  model$A <- matrix(1, nrow = 1, ncol = J)  # Row of ones for the sum constraint
  model$rhs <- c(1)
  model$sense <- c("=")


  # Variable bounds
  model$lb <- rep(0, J)  # Non-negativity bounds
  model$ub <- rep(1, J)  # No upper bound

  # Solve the model using Gurobi
  params <- list(OutputFlag = 0, TimeLimit = 100) # debug model false

  # Extract weights
  result <- tryCatch({
    gurobi(model, params = params)
  }, error = function(e) {
    return(NA)
  })

  # Check if result is NA due to timeout or other issues
  if (is.null(result) || !is.list(result) || is.null(result$status) || result$status != "OPTIMAL") {
    return(list(loss.B = NA, solution.w = NA))
  }

  # Extract weights
  solution.w <- as.numeric(result$x)
  if (is.null(solution.w) || any(is.na(solution.w))) {
    return(list(loss.B = NA, solution.w = NA))
  }

  solution.w[solution.w < 0.005] <- 0

  # Calculate the loss for w
  if (any(is.na(solution.w))) {
    loss.w <- NA
    loss.B <- NA
  } else {
    loss.w <- as.numeric(
      b_F * (sum((matrix(F_treated) - t(F_control) %*% solution.w)^2) / t_max) +
        b_Z * (sum((matrix(Z_treated) - t(Z_control) %*% solution.w)^2) / dr) +
        b_X * (sum((matrix(X_treated) - t(X_control) %*% solution.w)^2) / dt)
    )

    # Return the loss for B
    loss.B <- as.numeric(sum((matrix(F_treated) - t(F_control) %*% solution.w)^2) / t_max)
  }

  return(list(loss.B = loss.B, solution.w = solution.w))
}




