% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/etts.R
\name{sc_ett}
\alias{sc_ett}
\title{Calculate Synthetic Control Data Fusion Estimated Treatment Effect (ETT)}
\usage{
sc_ett(Y_treated, Y, w, J)
}
\arguments{
\item{Y_treated}{A numeric vector of outcomes for the treated unit.}

\item{Y}{A matrix of outcomes for control units including treated unit in the first row.}

\item{w}{A numeric weight vector for constructing the synthetic control.}

\item{J}{The number of control units used for the synthetic control.}
}
\value{
Estimated treatment effect (ETT) as a single numeric value.
}
\description{
This function computes the estimated causal effect using synthetic control data fusion method.
The ETT is calculated as the difference between the outcomes for treated unit and weighted sum of the outcomes of the control units.
}
\details{
\deqn{
\hat{\psi}^{sc} =  Y_{1} - \sum_{i=2}^{J+1} w_i Y_{i}.
}

where:
\itemize{
\item \eqn{\hat\psi^{eq0}} is the causal estimator using a synthetic control data fusion method.
\item \eqn{Y_{1}} is the observed outcome for the treated unit
\item \eqn{Y_{i}} is the outcome for control unit \eqn{i}
\item \eqn{w_i} is the synthetic control weight for control unit \eqn{i}
}
}
