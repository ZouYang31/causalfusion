# CausalPanelFusion <img src="https://img.shields.io/badge/R-package-blue.svg" align="right" />

**CausalPanelFusion** is an R package for estimating causal effects when pre-intervention data in the target domain is missing or incomplete. It introduces **two data fusion methods** that leverage auxiliary panel data from related reference domains to estimate treatment effects in the target domain.

These methods overcome the limitations of conventional synthetic control by recovering counterfactual outcomes even in the absence of pre-treatment information.

### 🔍 Included Methods

- **Equi-Confounding Data Fusion**   
- **Synthetic Control Data Fusion**

📄 Read the full theory paper on arXiv: [Causal Data Fusion for Panel Data without Pre-intervention Period](https://arxiv.org/abs/2410.16391)

---

## 🚀 Installation

You can install the development version directly from GitHub:

```r
# install.packages("devtools")
devtools::install_github("ZouYang31/causalfusion")

