##########################################
# Effect of Treatment on the Treated
#########################################


#' Calculate Synthetic Control Data Fusion Estimated Treatment Effect (ETT)
#'
#' This function computes the estimated causal effect using synthetic control data fusion method.
#' The ETT is calculated as the difference between the outcomes for treated unit and weighted sum of the outcomes of the control units.
#' @details
#' \deqn{
#' \hat{\psi}^{sc} =  Y_{1} - \sum_{i=2}^{J+1} w_i Y_{i}.
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{\hat\psi^{eq0}} is the causal estimator using a synthetic control data fusion method.
#'   \item \eqn{Y_{1}} is the observed outcome for the treated unit
#'   \item \eqn{Y_{i}} is the outcome for control unit \eqn{i}
#'   \item \eqn{w_i} is the synthetic control weight for control unit \eqn{i}
#' }
#' @param Y_treated A numeric vector of outcomes for the treated unit.
#' @param Y A matrix of outcomes for control units including treated unit in the first row.
#' @param w A numeric weight vector for constructing the synthetic control.
#' @param J The number of control units used for the synthetic control.
#'
#' @return Estimated treatment effect (ETT) as a single numeric value.
#' @export
# ------------------------------------------------------------------------------#
# Function: Calculate Synthetic Estimated Treatment Effect (ETT)
# ------------------------------------------------------------------------------#
sc_ett <- function(Y_treated, Y, w, J) {
  # Calculate the synthetic treated unit
  synthetic_treated_Y <- t(w) %*% Y[2:(J + 1), ]

  # Compute the Estimated Treatment Effect (ETT)
  ett_synthetic <- Y_treated - synthetic_treated_Y

  # Return the ETT result
  return(mean(ett_synthetic))
}

#' Estimated Treatment Effect (ETT) using Linear Equi-Confounding
#'
#' This function computes the estimated treatment effect using a linear equi-confounding method.
#' The ETT is calculated to adjust the naive difference in outcomes between the treated and control units in the target domain by subtracting the corresponding difference in the reference domain.
#'
#' @details
#' \deqn{
#' \hat\psi^{eq1}=\left\{Y_1-F_1\right\}-\frac{1}{J}\sum_{i=2}^{J+1}\left\{Y_i- F_i\right\}.
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{\hat\psi^{eq1}} is the causal estimator using a linear equi-confounding method.
#'   \item \eqn{Y_{1}} is the observed outcome for the treated unit of target domain
#'   \item \eqn{Y_{i}} is the outcome for the control unit\eqn{i} of target domain.
#'   \item \eqn{F_{1}} is the observed outcome for the treated unit of reference domain
#'   \item \eqn{F_{i}} is the outcome for control for the control unit\eqn{i} of reference domain.
#'   \item \eqn{J} is the number of control units.
#'
#' }
#' @param Y_treated A numeric vector of outcomes for the treated group of target domain.
#' @param Y_control A numeric vector of outcomes for the control group of target domain.
#' @param F_treated A numeric vector of outcomes measurements for the treated group of reference domain.
#' @param F_control A numeric vector of outcomes measurements for the control group of reference domain.
#'
#' @return Estimated treatment effect (ETT) as a single numeric value.
#' @export
# ------------------------------------------------------------------------------#
# Function: Compute Estimated Treatment Effect (ETT) using Equi-Confounding Data Fusion
# ------------------------------------------------------------------------------#
linear_equi_confounding_ett <- function(Y_treated, Y_control, F_treated, F_control) {

  # Compute means for observed and control units
  Y_1 <- mean(Y_treated)  # Observed target unit
  Y_0 <- mean(Y_control)
  F_1 <- mean(F_treated)
  F_0 <- mean(F_control)

  # Compute the Estimated Treatment Effect (ETT) for Equi-Confounding
  ett <- (Y_1 - Y_0) - (F_1 - F_0)

  return(ett)
}

#' Estimated Treatment Effect (ETT) using Logarithmic Equi-Confounding
#'
#' This function computes the estimated treatment effect using a logarithmic equi-confounding method.
#' The ETT is calculated to adjust the target domain's treated outcome using outcome ratios from a reference domain.
#' @details
#' \deqn{
#' \hat\psi^{eq2}=Y_1 - \frac{F_1}{ \sum_{i=2}^{J+1} F_i} \sum_{i=2}^{J+1} Y_i.
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{\hat\psi^{eq2}} is the causal estimator using a logarithmic equi-confounding method.
#'   \item \eqn{Y_{1}} is the observed outcome for the treated unit of target domain
#'   \item \eqn{Y_{i}} is the outcome for the control unit \eqn{i} of target domain.
#'   \item \eqn{F_{1}} is the observed outcome for the treated unit of reference domain
#'   \item \eqn{F_{i}} is the outcome for control for the control unit \eqn{i} of reference domain.
#'   \item \eqn{J} is the number of control units.
#'
#' }
#' @param Y_treated A numeric vector of outcomes for the treated group of target domain.
#' @param Y_control A numeric vector of outcomes for the control group of target domain.
#' @param F_treated A numeric vector of outcomes measurements for the treated group of reference domain.
#' @param F_control A numeric vector of outcomes measurements for the control group of reference domain.
#'
#' @return Estimated treatment effect (ETT) as a single numeric value.
#' @export
log_equi_confounding_ett <- function(Y_treated, Y_control, F_treated, F_control) {

  # Normalize using log of the mean of control group
  nor_a <- log(mean(Y_control))

  # Compute log-transformed values
  log_Y_1 <- log(mean(Y_treated)) / nor_a  # Observed target unit
  log_Y_0 <- log(mean(Y_control)) / nor_a
  log_B_1 <- log(mean(F_treated)) / nor_a
  log_B_0 <- log(mean(F_control)) / nor_a

  # Compute the Estimated Treatment Effect (ETT) for Logarithm Equi-Confounding
  estimate <- mean(Y_treated) - (mean(F_treated) / mean(F_control)) * mean(Y_control)

  return(estimate)
}




