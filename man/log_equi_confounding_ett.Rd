% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/etts.R
\name{log_equi_confounding_ett}
\alias{log_equi_confounding_ett}
\title{Estimated Treatment Effect (ETT) using Logarithmic Equi-Confounding}
\usage{
log_equi_confounding_ett(Y_treated, Y_control, F_treated, F_control)
}
\arguments{
\item{Y_treated}{A numeric vector of outcomes for the treated group of target domain.}

\item{Y_control}{A numeric vector of outcomes for the control group of target domain.}

\item{F_treated}{A numeric vector of outcomes measurements for the treated group of reference domain.}

\item{F_control}{A numeric vector of outcomes measurements for the control group of reference domain.}
}
\value{
Estimated treatment effect (ETT) as a single numeric value.
}
\description{
This function computes the estimated treatment effect using a logarithmic equi-confounding method.
The ETT is calculated to adjust the target domain's treated outcome using outcome ratios from a reference domain.
}
\details{
\deqn{
\hat\psi^{eq2}=Y_1 - \frac{F_1}{ \sum_{i=2}^{J+1} F_i} \sum_{i=2}^{J+1} Y_i.
}

where:
\itemize{
\item \eqn{\hat\psi^{eq2}} is the causal estimator using a logarithmic equi-confounding method.
\item \eqn{Y_{1}} is the observed outcome for the treated unit of target domain
\item \eqn{Y_{i}} is the outcome for the control unit \eqn{i} of target domain.
\item \eqn{F_{1}} is the observed outcome for the treated unit of reference domain
\item \eqn{F_{i}} is the outcome for control for the control unit \eqn{i} of reference domain.
\item \eqn{J} is the number of control units.

}
}
