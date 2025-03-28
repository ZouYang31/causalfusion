##########################################
# Functions for Optimization
#
#########################################
#' Compute Normalized Squared Error (NSE)
#'
#' This function calculates the normalized squared error (NSE) between
#' a target unit and control units.
#' It can be applied to different data matrices: \code{F}, \code{X}, or \code{Z},
#' depending on the value of the \code{target} parameter.
#'
#' @param w A numeric weight vector of length \code{J}, where \code{J = i_max - 1}.
#' @param F A numeric matrix representing outcome variables in the reference domain.
#' @param Z A numeric matrix of covariates in reference domain.
#' @param X A numeric matrix for covariates in target domain.
#' @param t_max An integer. The number of time points used for outcome in reference domain (used for \code{F}).
#' @param dr An integer. The number of covariates in the reference domain \code{Z}.
#' @param dt An integer. The number of covariates in the target domain \code{X}.
#' @param i_max An integer. Total number of units, including treated and donor units.
#' @param target A character string. Indicates which matrix to compute the error for: \code{"F"}, \code{"X"}, or \code{"Z"}. Default is \code{"X"}.
#'
#' @return A numeric scalar representing the normalized squared error.
#'
#' @examples
#' w <- runif(3)
#' F <- matrix(rnorm(40), nrow = 4)
#' X <- matrix(rnorm(40), nrow = 4)
#' Z <- matrix(rnorm(40), nrow = 4)
#' NSE_x(w, F, X, Z, t_max = 10, dr = 10, dt = 10, i_max = 4, target = "X")
#'
#' @export


NSE_x <- function(w, F, X, Z, t_max, dr, dt, i_max, target = "X") {
  J <- i_max - 1
  if (target == "F") {
    return(sum((F[1, 1:t_max] - t(w) %*% F[2:(J+1), 1:t_max])^2) / t_max)
  } else if (target == "X") {
    return(sum((X[1, ] - t(w) %*% X[2:(J+1), ])^2) / dt)
  } else if (target == "Z") {
    return(sum((Z[1, ] - t(w) %*% Z[2:(J+1), ])^2) / dr)
  } else {
    stop("Invalid target. Choose 'F', 'X', or 'Z'.")
  }
}


#' Optimize Synthetic Weights via Quadratic Programming (ipop)
#'
#' This function estimates synthetic control weights using the
#' \code{ipop()} function from the \code{kernlab} package to solve a constrained
#' quadratic programming problem. The weights are optimized to match the treated unit
#' to a convex combination of control units based on covariates \code{X} or \code{Z}.
#'
#' @param X A numeric matrix of covariates in the target domain or reference domain. First row is the treated unit,
#'        remaining rows are control units.
#' @param n_X Integer. Number of covariates in the target domain or reference domain, used for scaling.
#' @param i_max Integer. Total number of units (treated + control).
#' @param margin.ipop Numeric. Margin for constraint violation tolerance in \code{ipop}. Default is \code{1e-4}.
#' @param sigf.ipop Integer. Number of significant figures in \code{ipop} solution. Default is \code{5}.
#' @param bound.ipop Numeric. Box constraint bound in \code{ipop}. Default is \code{10}.
#'
#' @return A list containing:
#' \describe{
#'   \item{weights}{A numeric vector of optimized synthetic control weights.}
#'   \item{loss}{The squared loss between treated and weighted control units.}
#' }
#'
#' @details
#' This function uses quadratic programming to find weights \eqn{w} such that
#' \eqn{X_{treated} \approx X_{control} \cdot w}, subject to constraints:
#' \itemize{
#'   \item \eqn{\sum w = 1}
#'   \item \eqn{w \geq 0}
#'   \item \eqn{w \leq 1}
#' }
#' The optimization problem is solved using \code{ipop()} from the \code{kernlab} package.
#'
#' @import kernlab
#' @examples
#' X <- matrix(rnorm(40), nrow = 4)
#' Z <- matrix(rnorm(40), nrow = 4)
#' optimize_w_ipop(X = X, n_X = 4, i_max = 4)
#'
#' @export
optimize_w_ipop <- function(X, n_X, i_max,
                            margin.ipop = 1e-4, sigf.ipop = 5, bound.ipop = 10) {

  J <- i_max - 1

  # Extract treated and control data
  X_treated <- matrix(X[1, ])
  X_control <- t(X[2:(J+1), ])

  V <- diag(rep(1 / n_X, n_X))  # Scaling factor for X

  # Define quadratic programming components
  H <- t(X_control) %*% V %*% X_control  # Quadratic term

  c <- -1 * t(X_treated) %*% V %*% X_control  # Linear term
  A <- t(rep(1, J))  # Equality constraint (sum w = 1)
  b <- 1
  l <- rep(0, J)  # Lower bound (w >= 0)
  u <- rep(1, J)  # Upper bound (w < 1)
  r <- 0

  # Solve the quadratic programming problem
  res <- ipop(c = c, H = H, A = A, b = b, l = l, u = u, r = r,
              bound = bound.ipop, margin = margin.ipop, sigf = sigf.ipop)

  # Extract weights and compute losses
  solution.w <- as.numeric(primal(res))

  # Set the threshold for w,
  solution.w[solution.w < 0.005] <- 0
  loss.w <- sum((X_treated - X_control %*% solution.w)^2)

  # Return the optimized weights and loss
  return(list(weights = solution.w, loss = loss.w))
}


