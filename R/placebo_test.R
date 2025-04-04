
utils::globalVariables(c("year", "Outcome", "Unit"))
#' Placebo Tests for Synthetic Control Data Fusion Method
#'
#' This function performs placebo tests by treating each control unit as a pseudo-treated unit and estimating
#' its synthetic control using all other units. It then compares the observed and synthetic outcomes to
#' compute placebo causal effects in both the reference domain (\code{F}) and the target domain (\code{Y}).
#' Results are returned as two customizable ggplot objects.
#'
#' @param B A list of combinations of a budgeting vector.
#' @param F A numeric matrix for the reference domain outcome. Rows = units; columns = time periods.
#'          The first row is the treated unit.
#' @param Y A numeric matrix for the target domain outcome (same structure as \code{F}).
#' @param t_max Integer. Number of time periods in the reference domain (\code{F}).
#' @param s_max Integer. Number of time periods in the target domain (\code{Y}).
#' @param i_max Integer. Total number of units (treated + controls).
#' @param w A numeric weight vector representing the synthetic control weights for the treated unit.
#' @param X A matrix of covariates in the target domain.
#' @param Z A matrix of covariates in the reference domain.
#' @param dr Integer. Number of covariates in \code{Z}.
#' @param dt Integer. Number of covariates in \code{X}.
#' @param eta_Z Numeric. Default is \code{0.1}.
#' @param eta_X Numeric. Default is \code{0.1}.
#' @param Ylab Character string. Y-axis label for both plots. Default is \code{"Causal Effect"}.
#' @param Xlab Character string. X-axis label for both plots. Default is \code{"Time"}.
#' @param ref_lim Numeric vector of length 2. Y-axis limits for the reference domain plot (\code{plot_F}).
#'                Default is \code{c(-0.1, 0.1)}.
#' @param tar_lim Numeric vector of length 2. Y-axis limits for the target domain plot (\code{plot_Y}).
#'                Default is \code{c(-0.5, 0.5)}.
#' @param Tunit_lab Character string. Legend label for the treated unit. Default is \code{"Target Unit"}.
#' @param Cunit_lab Character string. Legend label for the control units. Default is \code{"Control Units"}.
#'
#' @return A list containing:
#' \describe{
#'   \item{plot_F}{A \code{ggplot} object showing placebo effects in the reference domain (\code{F}).}
#'   \item{plot_Y}{A \code{ggplot} object showing placebo effects in the target domain (\code{Y}).}
#' }
#' The plots visualize the difference between observed and synthetic outcomes across placebo-treated control units.
#'
#' @details
#' This function uses synthetic control methods to estimate placebo effects by rotating treatment assignment
#' across control units. It calculates the baseline normalized squared error (NSE) using covariate fit in both
#' X and Z domains, then uses those baselines with slack parameters (\code{eta_X}, \code{eta_Z}) to optimize
#' weights and compare effects. The resulting placebo plots visualize the effect trajectories across units.
#'
#' The treated unit's effect is plotted in black, and control units are shown in grey.
#' Labels and y-axis ranges are fully customizable via function arguments.
#'
#' @seealso \code{\link{optimize_w_ipop}}, \code{\link{find_best_B}}, \code{\link{NSE_x}}
#'
#' @examples
#' \dontrun{
#' placebo_result <- placebo_test(F = F_matrix, Y = Y_matrix,
#'                                t_max = 10, s_max = 10, i_max = 20,
#'                                w = w_star, X = X_matrix, Z = Z_matrix,
#'                                dr = 5, dt = 5,
#'                                eta_Z = 0.1, eta_X = 0.1,
#'                                Ylab = "Effect", Xlab = "Time",
#'                                Tunit_lab = "Treated", Cunit_lab = "Controls")
#' print(placebo_result$plot_F)
#' print(placebo_result$plot_Y)
#' }
#'
#' @import ggplot2
#' @importFrom scales alpha
#' @export

# ## ----         placebo test -------------  ##
placebo_test<- function(B, F, Y, t_max, s_max, i_max,
                        w, X, Z, dr, dt,
                        eta_Z = 0.1, eta_X = 0.1,
                        Ylab = c("Causal Effect"),
                        Xlab = c("Time"),
                        ref_lim = c(-0.1, 0.1),
                        tar_lim = c(-0.5, 0.5),
                        Tunit_lab = c("Target Unit"),
                        Cunit_lab = c("Control Units") ){

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
  result_X <- optimize_w_ipop(X = X, n_X = dt, i_max = i_max)
  wX <- result_X$weights
  NSE_X_baseline <- NSE_x(wX, F, X, Z, t_max, dr, dt, i_max, target = "X")


  result_Z <- optimize_w_ipop(X = Z, n_X = dr, i_max = i_max)
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

    b_list <- find_best_B(B = B, F_treated = F_placebo_treated, F_control = F_placebo_control,
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
      color_map_F[unit] <-  scales::alpha("grey", 0.6)
      line_type_map_F[unit] <- "solid"
    }
  }

  color_map_Y <- c("Synthetic_Treated_Y" = "black")
  line_type_map_Y <- c("Synthetic_Treated_Y" = "solid")

  for (unit in unique(synthetic_placebo_df_Y$Unit)) {
    if (!unit %in% names(color_map_Y)) {
      color_map_Y[unit] <-  scales::alpha("grey", 0.6)
      line_type_map_Y[unit] <- "solid"
    }
  }

  # Custom legend breaks and labels
  legend_breaks <- c("Synthetic_Treated_F", "Unit_2")  # Add "Unit_2" to represent control units
  legend_labels <- c(Tunit_lab, Cunit_lab)

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
    scale_y_continuous(limits = ref_lim) +
    labs(title = NULL, x = Xlab, y = Ylab, color = "Legend", linetype = "Legend") +
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
           linetype = guide_legend(title = NULL,
                                   order = 1))

  legend_breaks <- c("Synthetic_Treated_Y", "Unit_2")  # Add "Unit_2" to represent control units
  legend_labels <- c(Tunit_lab, Cunit_lab)

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
    scale_y_continuous(limits = tar_lim) +
    labs(title = NULL, x = Xlab, y = Ylab, color = "Legend", linetype = "Legend") +
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
           linetype = guide_legend(title = NULL,
                                   order = 1))


  return(list(plot_F = plot_F, plot_Y = plot_Y))
}
