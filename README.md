# causalfusion <img src="https://img.shields.io/badge/R-package-blue.svg" align="right" />

**CausalFusion** is an R package for estimating causal effects when pre-intervention data in the target domain is missing or incomplete. It introduces **three data fusion methods** that leverage auxiliary panel data from related reference domains to estimate treatment effects in the target domain.

These methods overcome the limitations of conventional synthetic control by recovering counterfactual outcomes even in the absence of pre-treatment information.

### 🔍 Included Methods

- **Linear Equi-Confounding**  
- **Logarithmic Equi-Confounding**  
- **Synthetic Control Data Fusion**

Each method encodes structural assumptions across domains and solves for synthetic control weights using constrained optimization with interpretable hyperparameters.

📄 Read the full theory paper on arXiv: [arXiv:2410.16391](https://arxiv.org/abs/2410.16391)

---

## 🚀 Installation

You can install the development version directly from GitHub:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install causalfusion from GitHub
devtools::install_github("ZouYang31/causalfusion")

## 📚 Load the Package

```r
library(causalfusion)
