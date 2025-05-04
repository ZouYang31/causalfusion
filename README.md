# causalfusion <img src="https://img.shields.io/badge/R-package-blue.svg" align="right" />

**CausalFusion** is an R package for estimating causal effects when pre-intervention data in the target domain is missing or incomplete. It introduces **three data fusion methods** that leverage auxiliary panel data from related reference domains to estimate treatment effects in the target domain.

These methods overcome the limitations of conventional synthetic control by recovering counterfactual outcomes even in the absence of pre-treatment information.

### 🔍 Included Methods

- **Linear Equi-Confounding**  
- **Logarithmic Equi-Confounding**  
- **Synthetic Control Data Fusion**

Each method encodes structural assumptions across domains and solves for synthetic control weights using constrained optimization with interpretable hyperparameters.

📄 Read the full theory paper on arXiv "Causal Data Fusion for Panel Data without Pre-Intervention Period": [arXiv:2410.16391](https://arxiv.org/abs/2410.16391)

---

## 🚀 Installation

You can install the development version directly from GitHub:

```{r,eval=F}
# install.packages("devtools")
devtools::install_github("ZouYang31/causalfusion")
```
You also need to download the gurobi package and license from gurobi website. Academic has free version license to use. 

```{r,eval=F}
# Setup the gurobi function
Sys.setenv(GRB_LICENSE_FILE = "your_path_to_gurobi_license/gurobi.lic")
install.packages("/Library/gurobi1200/macos_universal2/R/gurobi_12.0-0_R_4.4.1.tgz", repos = NULL, type = "source")
```

---

## Usage 

The package contains functions for data cleaning and generating the synthetic and negative control. 
| Function  | Description
|:---------|:----------|
| `NSE_x()` | This function calculates the normalized squared error (NSE) between a target unit and control units. |
| `optimize_w_ipop()` | Calculate the best possible NSEs for covariates. |
| `optimize_w_b_gurobi_B_list()` | This function computes optimal synthetic control weights \code{w} by solving a constrained quadratic optimization problem using the \code{gurobi} solver. |
| `find_best_B()` | This function searches over a list of candidate B parameter sets (weightings for domains F, X, and Z) to find the one that yields the lowest loss using the Gurobi-based optimizer. |


## Usage Example
Suppose you have one data frames, which includes your data for reference and target domain.
```{r,message=F}
library(gurobi)
library(causalfusion)
```
Data preprocessing
```{r,cache=T}
# 1. Separate the reference and target domain. 
df_tar <- panel_data[, c("unit", "month", "X1", "X2", "Y")]
df_ref <- panel_data[, c("unit", "month", "Z1", "Z2", "Z3", "F")]

# 2. Normalize the covariates in both domains.
df_tar <- normalize_columns(
  df_tar,
  c("X1","X2")
)

df_ref <- normalize_columns(
  df_ref,
  c("Z1","Z2", "Z3")
)

# 3. Reorder to target unit to the top.
df_tar <- reorder_unit_first(df_tar, "unit_VV4Q")
df_ref <- reorder_unit_first(df_ref, "unit_VV4Q")
```
Next, we prepare for the matrix F, Y, Z, X. 
```{r,cache=T}
# 4. Get the covariates and outcome matrix.
F <- outcome_matrix(df_ref, colname_outcome_var = "F",
                    colname_unit = "unit", colname_time = "month")

Y <- outcome_matrix(df_tar, colname_outcome_var = "Y",
                    colname_unit = "unit", colname_time = "month")

# Aggregate covariates 
Z_mean <- aggregate(cbind(Z1, Z2, Z3) ~ unit, data = panel_data, FUN = mean)

Z <- covariates_matrix(Z_mean, c("Z1", "Z2", "Z3"), "unit")

X_mean <- aggregate(cbind(X1, X2) ~ unit, data = panel_data, FUN = mean)

X <- covariates_matrix(X_mean, c("X1", "X2"), "unit")
```

Calculate the treated and control units
```{r,cache=T}
i_max <- 20   # Total number of cities
J <- i_max - 1 # Control units
w <- rep(0, J)    # Initialize weight vector

Y_treated <- Y[1, ]
Y_control <- Y[2:(J + 1), ]

F_treated <- F[1, ]
F_control <- F[2:(J + 1), ]

Z_treated <- Z[1, ]
Z_control <- Z[2:(J + 1), ]

X_treated <- X[1, ]
X_control <- X[2:(J + 1), ]

# Number of unique time points in vaccination data
t_max <- length(Y_treated)
s_max <- length(F_treated)

# Number of covariates in reference and target domains
dr <- length(Z_treated)   # Reference domain covariates
dt <- length(X_treated)  # Target domain covariates

B <- generate_b_list(step = 0.01, min_value = 0.01)

# Calculate the baseline for NSE_Z and NSE_X
result_X <- optimize_w_ipop(X = X, n_X = dt, i_max = num_city,
                            margin.ipop = 1e-4, sigf.ipop = 5, bound.ipop = 10)

wX <- result_X$weights
NSE_X_baseline <- NSE_x(w = wX, F= F, X= X, Z= Z, t_max, dr, dt, i_max=num_city, target = "X")

result_Z <- optimize_w_ipop(X = Z, n_X = dr, i_max = num_city,
                            margin.ipop = 1e-4, sigf.ipop = 5, bound.ipop = 10)
wZ <- result_Z$weights
NSE_Z_baseline <- NSE_x(w = wZ,  F = F, X = X, Z = Z, t_max, dr, dt, i_max=num_city, target = "Z")
```

Find the best possible B from optimization process.
```{r,cache=T}
# Optimize Best B
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
  i_max = num_city,
  NSE_Z_baseline = NSE_Z_baseline,
  NSE_X_baseline = NSE_X_baseline,
  eta_Z = 0.1,
  eta_X = 0.1
)


# Solve for weight vector (w) using the optimal B values
w <- b_list$best_w_star
```

Calculate the causal effect using 3 methods. 

```{r,cache=T}
# Calculate the Estimated Treatment Effect (ETT)
ett_synthetic <- sc_ett(Y_treated, Y = Y, w = w, J = J)
ett_linear <- linear_equi_confounding_ett(Y_treated, Y_control, F_treated, F_control)
ett_log <- log_equi_confounding_ett(Y_treated, Y_control, F_treated, F_control)
```


## Debugging

Spot an issue? Please let me know by posting an issue.