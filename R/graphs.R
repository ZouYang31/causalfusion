#################################################################
# Generate Reference and Target Domain Plot for Synthetic Control
#################################################################

#' Plot Synthetic Reference Outcome Over Time
#'
#' This function creates a line plot comparing the treated target unit and its synthetic target unit
#' with other control units over time of the reference domain.
#'
#' @param F A matrix of outcomes of reference domain (rows = units, columns = time).
#' @param w A numeric vector of weights for synthetic control.
#' @param t_max Maximum number of time periods.
#' @param i_max Total number of units (including treated and controls).
#' @param J Number of control units.
#' @param title Optional plot title.
#' @param x_label Label for the x-axis. Default is "Time".
#' @param y_label Label for the y-axis. Default is "Outcomes in Target Domain".
#' @param legend_labels Character vector for custom legend labels. Default is "Target "
#' @return A graph representing the comparison of outcomes for the observed target unit and synthetic target unit in the reference domain
#' @export
synthetic_reference_plot_real_data <- function(F, w, t_max, i_max, J,
                                               title = NULL,
                                               x_label = "Time",
                                               y_label = "Outcome",
                                               legend_labels = c("Target Unit", "Synthetic Target Unit", "Control Units")) {

  F_treated <- F[1, ]
  synthetic_treated_F <- t(w) %*% F[2:(J+1), ]

  data_all <- data.frame(
    year = rep(1:t_max, i_max),
    Outcome = as.vector(t(F)),
    Unit = rep(paste0("Unit_", 1:i_max), each = t_max)
  )

  data_synthetic <- data.frame(
    year = 1:t_max,
    Outcome = as.vector(synthetic_treated_F),
    Unit = rep("Synthetic_control", t_max)
  )

  data_all <- rbind(data_all, data_synthetic)

  color_map <- c("Unit_1" = "black", "Synthetic_control" = "black")
  line_type_map <- c("Unit_1" = "solid", "Synthetic_control" = "dashed")

  for (unit in unique(data_all$Unit)) {
    if (!unit %in% names(color_map)) {
      color_map[unit] <- scales::alpha("grey", 0.5)
      line_type_map[unit] <- "solid"
    }
  }
  legend_breaks = c("Unit_1", "Synthetic_control", "Unit_2")
  # Match override aesthetics to number of legend items
  override_colors <- c("black", "black", "grey")[1:length(legend_breaks)]
  override_linetypes <- c("solid", "dashed", "solid")[1:length(legend_breaks)]

  p <- ggplot(data_all, aes(x = year, y = Outcome, group = Unit, color = Unit, linetype = Unit)) +
    geom_line(data = subset(data_all, !Unit %in% c("Unit_1", "Synthetic_control")),
              size = 0.8, alpha = 0.5) +
    geom_line(data = subset(data_all, Unit == "Synthetic_control"), size = 0.8) +
    geom_line(data = subset(data_all, Unit == "Unit_1"), size = 0.8) +
    scale_color_manual(values = color_map,
                       breaks = legend_breaks,
                       labels = legend_labels) +
    scale_linetype_manual(values = line_type_map,
                          breaks = legend_breaks,
                          labels = legend_labels) +
    labs(title = title,
         x = x_label,
         y = y_label,
         color = "Legend",
         linetype = "Legend") +
    theme_minimal() +
    theme(
      text = element_text(family = "sans"),
      legend.position = c(0.95, 0.15),
      legend.justification = c(1, 1),
      legend.text = element_text(size = 14),
      legend.key.size = unit(1, "lines"),
      legend.key.width = unit(2, "lines"),
      legend.margin = margin(0, 0, 0, 0),
      legend.background = element_rect(fill = scales::alpha('grey', 0.5)),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.ticks.length = unit(0.25, "cm")
    ) +
    guides(
      color = guide_legend(
        title = NULL,
        override.aes = list(
          linetype = override_linetypes,
          color = override_colors
        ),
        order = 1
      ),
      linetype = "none"
    )

  return(p)
}



#' Plot Synthetic Target Outcome Over Time
#'
#' This function creates a line plot comparing the treated target unit and its synthetic target unit
#' with other control units over time of the target domain.
#'
#' @param Y A matrix of outcomes for control units including treated unit in the first row.
#' @param w A numeric weight vector for constructing the synthetic control.
#' @param s_max Maximum number of time periods.
#' @param i_max Total number of units (including treated and controls).
#' @param J The number of control units.
#' @param title Optional plot title.
#' @param x_label Label for the x-axis. Default is "Time".
#' @param y_label Label for the y-axis. Default is "Outcomes in Target Domain".
#' @param legend_labels Character vector for custom legend labels. Default is "Target Unit", "Synthetic Target Unit", "Control Units"
#' @return A graph representing the comparison of outcomes for the observed target unit and synthetic target unit in the target domain
#' @export
synthetic_target_plot_real_data <- function(Y, w, s_max, i_max, J,
                                            title = NULL,
                                            x_label = "Time",
                                            y_label = "Outcomes in Target Domain",
                                            legend_labels = c("Target Unit", "Synthetic Target Unit", "Control Units")) {
  # Compute synthetic treated outcome
  synthetic_treated_Y <- t(w) %*% Y[2:(J + 1), ]

  # Create full data set
  data_all <- data.frame(
    year = rep(1:s_max, i_max),
    Outcome = as.vector(t(Y)),
    Unit = rep(paste0("Unit_", 1:i_max), each = s_max)
  )

  # Add synthetic treated
  data_synthetic <- data.frame(
    year = 1:s_max,
    Outcome = as.vector(synthetic_treated_Y),
    Unit = rep("Synthetic_control", s_max)
  )

  # Add observed treated unit (row 1)
  data_observed <- data.frame(
    year = 1:s_max,
    Outcome = as.vector(Y[1, ]),
    Unit = rep("Observed_target", s_max)
  )

  # Combine all data
  data_all <- rbind(data_all, data_synthetic, data_observed)

  # Define colors and line types
  color_map <- c("Observed_target" = "black", "Synthetic_control" = "black", "Unit_1" = "black")
  line_type_map <- c("Observed_target" = "solid", "Synthetic_control" = "dashed", "Unit_1" = "solid")

  # Set others to grey
  for (unit in unique(data_all$Unit)) {
    if (!unit %in% names(color_map)) {
      color_map[unit] <- scales::alpha("grey", 0.5)
      line_type_map[unit] <- "solid"
    }
  }

  legend_breaks <- c("Observed_target", "Synthetic_control", "Unit_2")
  override_colors <- c("black", "black", "grey")[1:length(legend_breaks)]
  override_linetypes <- c("solid", "dashed", "solid")[1:length(legend_breaks)]

  # Plot
  p <- ggplot(data_all, aes(x = year, y = Outcome, group = Unit, color = Unit, linetype = Unit)) +
    geom_line(data = subset(data_all, !Unit %in% c("Observed_target", "Synthetic_control")),
              size = 0.8, alpha = 0.5) +
    geom_line(data = subset(data_all, Unit == "Synthetic_control"), size = 0.8) +
    geom_line(data = subset(data_all, Unit == "Observed_target"), size = 0.8) +
    scale_color_manual(values = color_map,
                       breaks = legend_breaks,
                       labels = legend_labels) +
    scale_linetype_manual(values = line_type_map,
                          breaks = legend_breaks,
                          labels = legend_labels) +
    labs(title = title,
         x = x_label,
         y = y_label,
         color = "Legend",
         linetype = "Legend") +
    theme_minimal() +
    theme(
      text = element_text(family = "sans"),
      legend.position = c(0.95, 0.2),
      legend.justification = c(1, 1),
      legend.text = element_text(size = 14),
      legend.key.size = unit(1, "lines"),
      legend.key.width = unit(2, "lines"),
      legend.margin = margin(0, 0, 0, 0),
      legend.background = element_rect(fill = scales::alpha('grey', 0.5)),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.ticks.length = unit(0.25, "cm")
    ) +
    guides(
      color = guide_legend(title = NULL,
                           override.aes = list(linetype = override_linetypes,
                                               color = override_colors),
                           order = 1),
      linetype = "none"
    )

  return(p)
}


