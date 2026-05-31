# Feature Inventory: JASP Regression

Upstream inspected: `/private/tmp/jaspRegression` cloned from https://github.com/jasp-stats/jaspRegression on 2026-05-31.

License found in upstream DESCRIPTION: GPL (>= 2).

No upstream source code has been copied into this scaffold yet.

| JASP Feature | Included | Notes |
|--------------|----------|-------|
| Correlation | Shell only | Upstream files: `R/correlation.R`, `R/correlationWrapper.R`, `inst/qml/Correlation.qml`; implementation pending. |
| Linear Regression | Shell only | Upstream files: `R/regressionlinear.R`, wrapper, QML; implementation pending. |
| Logistic Regression | Shell only | Upstream files: `R/regressionlogistic.R`, wrapper, QML; implementation pending. |
| Generalized Linear Model | Shell only | Upstream files: `R/generalizedlinearmodel.R`, wrapper, QML; implementation pending. |
| Bayesian Correlation | Shell only | Upstream files: `R/correlationbayesian.R`; implementation pending. |
| Bayesian Linear Regression | Shell only | Upstream files: `R/regressionlinearbayesian.R`; implementation pending. |
| Bayesian Logistic Regression | Shell only | Upstream files: `R/regressionlogisticbayesian.R`; implementation pending. |
| Model comparison | Not yet | Must preserve JASP model-building behavior and outputs. |
| Diagnostics and residual checks | Not yet | Preserve JASP calculations; add jamovi-native guidance. |
| Plots | Not yet | Preserve JASP plot content while rendering in jamovi style. |
| APA reporting | Shell only | Placeholder section exists; statistical narrative pending. |
| About and citation panel | Included | Implemented as static HTML per analysis for now. |

