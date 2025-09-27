# CausalPanelFusion <img src="https://img.shields.io/badge/R-package-blue.svg" align="right" />

**CausalPanelFusion** is an R package for estimating causal effects when pre-intervention data in the target domain is not available. It contains **two data fusion methods** that leverage auxiliary panel data from related reference domains to estimate treatment effects in the target domain.

Unlike the conventional panel data methods, our proposed methods can recover counterfactual outcomes even in the absence of pre-intervention data in the target domain.

### 🔍 Included Methods

- **Equi-Confounding Data Fusion**   
- **Synthetic Control Data Fusion**

📄 Read the full paper on arXiv: [Causal Data Fusion for Panel Data without a Pre-intervention Period](https://arxiv.org/abs/2410.16391)

---

## 🚀 Installation

You can install the development version directly from GitHub:

```r
# install.packages("devtools")
devtools::install_github("ZouYang31/causalfusion")

