

#' Placebo Tests for Synthetic Control Data Fusion Method
#'
#' This function runs placebo tests by treating each control unit as if it were the treated unit
#' and estimating the synthetic control weights. It compares the actual outcome against the
#' synthetic counterfactual outcome for both the reference outcome \code{F} and the target outcome \code{Y}.
#'
#' @param F A numeric matrix representing the reference domain outcome for treated and control units.
#'          The first row is the treated unit, and the rest are control units.
#' @param Y A numeric matrix representing the target domain outcome.
#' @param t_max Integer. Number of pre-treatment time periods for \code{F}.
#' @param s_max Integer. Number of post-treatment time periods for \code{Y}.
#' @param i_max Integer. Total number of units.
#' @param w A numeric weight vector representing the synthetic control weights for the treated unit.
#' @param X A numeric matrix of covariates in the target domain.
#' @param Z A numeric matrix of covariates in the reference domain.
#' @param dr Integer. Number of covariates in \code{Z}.
#' @param dt Integer. Number of covariates in \code{X}.
#' @param eta_Z Numeric. Default is \code{0.1}.
#' @param eta_X Numeric. Default is \code{0.1}.
#'
#' @return A list containing:
#' \describe{
#'   \item{plot_F}{A \code{ggplot} object showing placebo effects in the reference domain (\code{F}).}
#'   \item{plot_Y}{A \code{ggplot} object showing placebo effects in the target domain (\code{Y}).}
#' }
#' The plots visualize the difference between observed and synthetic outcomes across placebo-treated control units.
#'
#' @details
#' For each control unit, the function constructs a synthetic version of that unit using the other control units.
#' It calculates the synthetic control effect as the difference between observed and synthetic outcomes.
#' These placebo effects are plotted alongside the actual treated unit’s effect to assess significance.
#' Synthetic control weights are optimized using `optimize_w_ipop()` and `find_best_B()` with constraints on covariate fit.
#'
#' @import ggplot2
#' @importFrom scales alpha
#' @importFrom utils head
#' @importFrom stats reshape
#' @export

# ## ----         placebo test -------------  ##
placebo_test<- function(F, Y, t_max, s_max, i_max,
                         w, X, Z, dr, dt,
                         eta_Z = 0.1, eta_X = 0.1){

  # Select all control units as placebo treated units
  placebo_units <- 2:i_max
  J <- i_max - 1

  synthetic_treated_F <- t(w) %*% F[2:(J+1), ]
  synthetic_treated_Y <- t(w) %*% Y[2:(J+1), ]

  # Initialize lists to store synthetic placebo outcomes
  synthetic_placebo_list_F <- list()
  synthetic_placebo_list_Y <- list()

  F_treated <- F[1, ]
  F_control <- F[2:(J+1), ]
  # Calculate the baseline X and Z
  result_X <- optimize_w_ipop(F_treated = F_treated, F_control = F_control,
                              X, Z, t_max, dr, dt, i_max, target = "X")
  wX <- result_X$weights
  NSE_X_baseline <- NSE_x(wX, F, X, Z, t_max, dr, dt, i_max, target = "X")


  result_Z <- optimize_w_ipop(F_treated = F_treated, F_control = F_control,
                              X, Z, t_max, dr, dt, i_max, target = "Z")
  wZ <- result_Z$weights
  NSE_Z_baseline <- NSE_x(wZ, F, X, Z, t_max, dr, dt, i_max, target = "Z")


  for (placebo_unit in placebo_units) {
    print(placebo_unit)
    F_placebo_treated <- F[placebo_unit, ]
    F_placebo_control <- F[-placebo_unit, ]

    Y_placebo_treated <- Y[placebo_unit, ]
    Y_placebo_control <- Y[-placebo_unit, ]

    X_placebo_treated <- X[placebo_unit, ]
    X_placebo_control <- X[-placebo_unit, ]

    Z_placebo_treated <- Z[placebo_unit, ]
    Z_placebo_control <- Z[-placebo_unit, ]

    b_list <- find_best_B(B, F_treated = F_placebo_treated, F_control = F_placebo_control,
                          X_treated = X_placebo_treated, X_control = X_placebo_control,
                          Z_treated = Z_placebo_treated, Z_control = Z_placebo_control,
                          t_max, dr, dt, i_max,
                          NSE_Z_baseline = NSE_Z_baseline, NSE_X_baseline = NSE_X_baseline, eta_Z= eta_Z, eta_X= eta_X)
    # Solve for W
    w_placebo <- b_list$best_w_star

    if (is.null(w_placebo)) {
      # Define default behavior when w_placebo is NULL
      warning(paste("w_placebo is NULL for unit", placebo_unit, ". Skipping this unit."))
      next # Skip to the next iteration of the loop
    }

    synthetic_placebo_treated_F <- t(w_placebo) %*% F[-placebo_unit, ]
    effect_F <- F[placebo_unit, ] - synthetic_placebo_treated_F

    synthetic_placebo_list_F[[paste0("Unit_", placebo_unit)]] <- effect_F

    synthetic_placebo_treated_Y <- t(w_placebo) %*% Y[-placebo_unit, ]

    effect_Y <- Y[placebo_unit, ] - synthetic_placebo_treated_Y
    synthetic_placebo_list_Y[[paste0("Unit_", placebo_unit)]] <- effect_Y
  }

  # Combine the synthetic placebo outcomes into data frames
  synthetic_placebo_df_F <- do.call(rbind, lapply(names(synthetic_placebo_list_F), function(unit) {
    data.frame(
      year = 1:t_max,
      Outcome = as.numeric(synthetic_placebo_list_F[[unit]]),
      Unit = unit
    )
  }))

  synthetic_placebo_df_Y <- do.call(rbind, lapply(names(synthetic_placebo_list_Y), function(unit) {
    data.frame(
      year = 1:s_max,
      Outcome = as.numeric(synthetic_placebo_list_Y[[unit]]),
      Unit = unit
    )
  }))

  # Highlight the treated unit (Unit_1)
  synthetic_placebo_df_F <- rbind(synthetic_placebo_df_F, data.frame(
    year = 1:t_max,
    Outcome = as.numeric(F[1,] - synthetic_treated_F),
    Unit = "Synthetic_Treated_F"
  ))

  synthetic_placebo_df_Y <- rbind(synthetic_placebo_df_Y, data.frame(
    year = 1:s_max,
    Outcome = as.numeric(Y[1,] - synthetic_treated_Y),
    Unit = "Synthetic_Treated_Y"
  ))

  color_map_F <- c("Synthetic_Treated_F" = "black")
  line_type_map_F <- c("Synthetic_Treated_F" = "solid")

  for (unit in unique(synthetic_placebo_df_F$Unit)) {
    if (!unit %in% names(color_map_F)) {
      color_map_F[unit] <- scales::alpha("grey", 0.6)
      line_type_map_F[unit] <- "solid"
    }
  }

  color_map_Y <- c("Synthetic_Treated_Y" = "black")
  line_type_map_Y <- c("Synthetic_Treated_Y" = "solid")

  for (unit in unique(synthetic_placebo_df_Y$Unit)) {
    if (!unit %in% names(color_map_Y)) {
      color_map_Y[unit] <- scales::alpha("grey", 0.6)
      line_type_map_Y[unit] <- "solid"
    }
  }

  # Custom legend breaks and labels
  legend_breaks <- c("Synthetic_Treated_F", "Unit_2")  # Add "Unit_2" to represent control units
  legend_labels <- c("Chelsea", "19 Cities")

  # Create the plots
  plot_F <- ggplot(synthetic_placebo_df_F, aes(x = year, y = Outcome, color = Unit, linetype = Unit)) +
    geom_line(data = subset(synthetic_placebo_df_F, !Unit %in% c("Synthetic_Treated_F")),
              aes(color = Unit, linetype = Unit), size = 0.8, alpha = 0.6) +  # Grey control lines first
    geom_line(data = subset(synthetic_placebo_df_F, Unit == "Synthetic_Treated_F"),
              aes(color = Unit, linetype = Unit), size = 0.8) +  # Solid black line on top
    geom_hline(yintercept = 0, linetype = "dashed", color = "blue", size = 0.5) + # Line at y = 0
    scale_color_manual(values = color_map_F,
                       breaks = legend_breaks,
                       labels = legend_labels) +
    scale_linetype_manual(values = line_type_map_F,
                          breaks = legend_breaks,
                          labels = legend_labels) +
    scale_y_continuous(limits = c(-0.1, 0.1)) +
    labs(title = NULL, x = "Time", y = "Causal Effect", color = "Legend", linetype = "Legend") +
    theme_minimal() +
    theme(text = element_text(family = "sans"),
          legend.position = c(0.95, 0.1),
          legend.justification = c(1, 1),
          legend.text = element_text(size = 14),
          legend.key.size = unit(1, "lines"),
          legend.key.width = unit(2, "lines"),
          legend.margin = margin(0, 0, 0, 0),
          legend.background = element_rect(fill = scales::alpha("grey", 0.5)),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.ticks.length = unit(0.25, "cm")) +
    guides(color = guide_legend(title = NULL, override.aes = list(linetype = c("solid", "solid"),
                                                                  color = c("black", "grey")),
                                order = 1),
           linetype = guide_legend(title = NULL, override.aes = list(color = c("black", "grey")),
                                   order = 1))

  legend_breaks <- c("Synthetic_Treated_Y", "Unit_2")  # Add "Unit_2" to represent control units
  legend_labels <- c("Chelsea", "19 Cities")

  plot_Y <- ggplot(synthetic_placebo_df_Y, aes(x = year, y = Outcome, color = Unit, linetype = Unit)) +
    geom_line(data = subset(synthetic_placebo_df_Y, !Unit %in% c("Synthetic_Treated_Y")),
              aes(color = Unit, linetype = Unit), size = 0.8, alpha = 0.6) +  # Grey control lines first
    geom_line(data = subset(synthetic_placebo_df_Y, Unit == "Synthetic_Treated_Y"),
              aes(color = Unit, linetype = Unit), size = 0.8) +  # Solid black line on top
    geom_hline(yintercept = 0, linetype = "dashed", color = "blue", size = 0.5) + # Line at y = 0
    scale_color_manual(values = color_map_Y,
                       breaks = legend_breaks,
                       labels = legend_labels) +
    scale_linetype_manual(values = line_type_map_Y,
                          breaks = legend_breaks,
                          labels = legend_labels) +
    labs(title = NULL, x = "Time", y = "Causal Effect", color = "Legend", linetype = "Legend") +
    theme_minimal() +
    theme(text = element_text(family = "sans"),
          legend.position = c(0.95, 0.1),
          legend.justification = c(1, 1),
          legend.text = element_text(size = 14),
          legend.key.size = unit(1, "lines"),
          legend.key.width = unit(2, "lines"),
          legend.margin = margin(0, 0, 0, 0),
          legend.background = element_rect(fill = scales::alpha('grey', 0.5)),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.ticks.length = unit(0.25, "cm")) +
    guides(color = guide_legend(title = NULL, override.aes = list(linetype = c("solid", "solid"),
                                                                  color = c("black", "grey")),
                                order = 1),
           linetype = guide_legend(title = NULL, override.aes = list(color = c("black", "grey")),
                                   order = 1))

  return(list(plot_F = plot_F, plot_Y = plot_Y))
}
