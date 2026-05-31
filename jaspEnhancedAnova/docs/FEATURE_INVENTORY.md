# Feature Inventory: JASP Enhanced ANOVA

Upstream inspected: `/private/tmp/jaspAnova` cloned from https://github.com/jasp-stats/jaspAnova on 2026-05-31.

License found in upstream DESCRIPTION: GPL (>= 2).

No upstream source code has been copied into this scaffold yet. The current implementation is a feature-parity scaffold: menus, options, and output placeholders are represented so the later statistical port has an explicit contract.

| JASP Feature | Included | Notes |
|--------------|----------|-------|
| ANOVA analysis entry | Shell only | Upstream files: `R/anova.R`, `R/anovaWrapper.R`, `inst/qml/Anova.qml`; implementation pending. |
| Repeated Measures ANOVA analysis entry | Shell only | Upstream files: `R/anovarepeatedmeasures.R`, wrapper, QML, help; implementation pending. |
| ANCOVA analysis entry | Shell only | Upstream files: `R/ancova.R`, wrapper, QML; implementation pending. |
| MANOVA analysis entry | Shell only | Upstream files: `R/manova.R`, wrapper, QML; implementation pending. |
| Bayesian ANOVA analysis entry | Shell only | Upstream files: `R/anovabayesian.R`, common Bayesian helpers; implementation pending. |
| Bayesian Repeated Measures ANOVA analysis entry | Shell only | Upstream files: `R/anovarepeatedmeasuresbayesian.R`; implementation pending. |
| Bayesian ANCOVA analysis entry | Shell only | Upstream files: `R/ancovabayesian.R`; implementation pending. |
| Dependent variable | UI scaffold | Large jamovi assignment box present. |
| Fixed factors | UI scaffold | Large jamovi assignment box present. |
| Covariates | UI scaffold | Large jamovi assignment box present. |
| Random factors | UI scaffold | Present as requested; must verify upstream support and document limitations if unsupported. |
| WLS weights | UI scaffold | Present as requested; calculation support pending. |
| Descriptive statistics | UI scaffold | Option and output placeholder present. |
| omega-squared | UI scaffold | Option present; JASP-equivalent calculation pending. |
| partial omega-squared | UI scaffold | Option present; JASP-equivalent calculation pending. |
| eta-squared | UI scaffold | Option present; JASP-equivalent calculation pending. |
| partial eta-squared | UI scaffold | Option present; JASP-equivalent calculation pending. |
| Confidence intervals for effect sizes | UI scaffold | Option and confidence level present. |
| User-defined confidence level | UI scaffold | Numeric confidence level option present. |
| Vovk-Sellke maximum p-ratio | UI scaffold | Option present. |
| Homogeneity correction: none | UI scaffold | List option present. |
| Homogeneity correction: Brown-Forsythe | UI scaffold | List option present. |
| Homogeneity correction: Welch | UI scaffold | List option present. |
| Q-Q plot of residuals | UI scaffold | Option present. |
| Residual plots | UI scaffold | Added enhancement option. |
| Leverage plots | UI scaffold | Added enhancement option. |
| Influence diagnostics | UI scaffold | Added enhancement option. |
| Cook's Distance | UI scaffold | Added enhancement option. |
| Mahalanobis distance | UI scaffold | Added enhancement option. |
| Custom contrasts | UI scaffold | Option and syntax field present. |
| Planned contrasts | UI scaffold | Option present. |
| Polynomial contrasts | UI scaffold | Option present. |
| Repeated contrasts | UI scaffold | Option present. |
| Helmert contrasts | UI scaffold | Option present. |
| Difference contrasts | UI scaffold | Option present. |
| Deviation contrasts | UI scaffold | Option present. |
| Contrast confidence intervals | UI scaffold | Option present. |
| Cohen's d for contrasts | UI scaffold | Option present. |
| Order restricted hypotheses syntax editor | UI scaffold | Text field present; syntax validation pending. |
| Multiple hypothesis models | UI scaffold | Option present. |
| Model comparison | UI scaffold | Option present. |
| Reference model selection | UI scaffold | Unconstrained, complement, and null options present. |
| Weight ratios | UI scaffold | Option present. |
| Relative weights matrix | UI scaffold | Option present. |
| Compare model coefficients | UI scaffold | Option present. |
| Model summary output | Output placeholder | Result section present. |
| Marginal means output | UI scaffold | Options and output placeholder present. |
| Informed hypothesis tests | UI scaffold | Option present. |
| Heterogeneity correction | UI scaffold | Option present. |
| Bootstrap confidence intervals | UI scaffold | Option present. |
| Bootstrap samples | UI scaffold | Numeric sample count present. |
| Tukey post hoc | UI scaffold | Option present. |
| Scheffe post hoc | UI scaffold | Option present. |
| Bonferroni post hoc | UI scaffold | Option present. |
| Holm post hoc | UI scaffold | Option present. |
| Sidak post hoc | UI scaffold | Option present. |
| Games-Howell post hoc | UI scaffold | Option present. |
| Dunnett post hoc | UI scaffold | Option present. |
| Post hoc confidence intervals | UI scaffold | Option present. |
| Significant comparison flags | UI scaffold | Option present. |
| Letter-based grouping tables | UI scaffold | Option present. |
| Post hoc effect sizes | UI scaffold | Option present. |
| Post hoc bootstrap estimates | UI scaffold | Option present. |
| Conditional comparisons for interactions | UI scaffold | Option present. |
| Descriptive plots | UI scaffold | Factor assignment, horizontal axis, separate lines, separate plots, and error bars present. |
| Bar plots | UI scaffold | Factor assignment, horizontal axis, separate plots, error bars, and zero-axis option present. |
| Raincloud plots | UI scaffold | Factor assignment, horizontal axis, separate plots, and horizontal display option present. |
| Estimated marginal means table | UI scaffold | Option present. |
| Marginal means pairwise comparisons | UI scaffold | Option present. |
| Marginal means confidence intervals | UI scaffold | Option present. |
| Compare marginal means to zero | UI scaffold | Option present. |
| CI adjustment: none | UI scaffold | Option present. |
| CI adjustment: Bonferroni | UI scaffold | Option present. |
| CI adjustment: Holm | UI scaffold | Option present. |
| CI adjustment: Sidak | UI scaffold | Option present. |
| CI adjustment: Tukey | UI scaffold | Option present. |
| Simple effects analysis | UI scaffold | Simple effect factor and two moderators present. |
| Simple effects tables | UI scaffold | Option present. |
| Simple effects post hoc comparisons | UI scaffold | Option present. |
| Simple effects confidence intervals | UI scaffold | Option present. |
| Simple effects effect sizes | UI scaffold | Option present. |
| Kruskal-Wallis test | UI scaffold | Option present. |
| Kruskal-Wallis epsilon-squared | UI scaffold | Option present. |
| Kruskal-Wallis eta-squared | UI scaffold | Option present. |
| Kruskal-Wallis confidence intervals | UI scaffold | Option present. |
| Dunn's post hoc tests | UI scaffold | Option present. |
| Append residuals | UI scaffold | Option present. |
| Raw residuals | UI scaffold | Option present. |
| Studentized residuals | UI scaffold | Option present. |
| Standardized residuals | UI scaffold | Option present. |
| Custom residual column name | UI scaffold | Text field present. |
| Append predictions | UI scaffold | Option present. |
| Custom prediction column name | UI scaffold | Text field present. |
| APA 7 narrative text | UI scaffold | Option and output placeholder present. |
| APA 7 tables | UI scaffold | Option present. |
| Copy-to-clipboard report | UI scaffold | Option present; actual clipboard behavior must be evaluated against jamovi capabilities. |
| Teaching mode | UI scaffold | Options and output placeholder present. |
| Explanation of each statistic | UI scaffold | Option present. |
| Effect-size interpretation | UI scaffold | Option present. |
| Assumption guidance | UI scaffold | Option present. |
| Recommended remedies | UI scaffold | Option present. |
| Publication-ready tables | UI scaffold | Option present. |
| Publication-ready figures | UI scaffold | Option present. |
| Export to Word | UI scaffold | Option present; export implementation pending. |
| Export to PDF | UI scaffold | Option present; export implementation pending. |
| Analysis syntax export | UI scaffold | Option present. |
| R syntax export | UI scaffold | Option present. |
| Save analysis configuration | UI scaffold | Option present. |
| About and citation panel | Included | Implemented as static HTML per analysis for now. |

## Repeated Measures ANOVA Feature Contract

| JASP Feature | Included | Notes |
|--------------|----------|-------|
| Repeated measures factor construction | UI scaffold | Factor names, levels, and free-form factor specification fields present; true dynamic cell grid pending. |
| Multiple repeated measures factors | UI scaffold | Supported through repeated measures factor specification field and repeated measures cell assignment. |
| Custom factor names | UI scaffold | Text field present. |
| Unlimited factor levels | UI scaffold | Text field present; practical UI validation pending. |
| Repeated measures cells | UI scaffold | Large assignment box present. |
| Multi-factor repeated designs | UI scaffold | Cell assignment and factor specification present. |
| Nested repeated designs | UI scaffold | Specification field present; implementation pending. |
| Between-subject factors | UI scaffold | Large assignment box present. |
| Covariates | UI scaffold | Large assignment box present. |
| Optional grouping factor | UI scaffold | Assignment box present for nonparametric follow-up workflows. |
| Descriptive statistics | UI scaffold | Option and output placeholder present. |
| omega-squared | UI scaffold | Option present. |
| partial omega-squared | UI scaffold | Option present. |
| eta-squared | UI scaffold | Option present. |
| partial eta-squared | UI scaffold | Option present. |
| generalized eta-squared | UI scaffold | Option present. |
| Effect-size confidence intervals | UI scaffold | Option and confidence level present. |
| Vovk-Sellke maximum p-ratio | UI scaffold | Option present. |
| Repeated measures model builder | UI scaffold | Automatic/manual model, RM effects, interactions, higher-order interactions, between effects, and mixed interactions present. |
| Sum of Squares Type I | UI scaffold | Option present. |
| Sum of Squares Type II | UI scaffold | Option present. |
| Sum of Squares Type III | UI scaffold | Option present. |
| Pool error terms | UI scaffold | Global and follow-up-specific options present. |
| Mauchly's Test of Sphericity | UI scaffold | Option present. |
| Sphericity correction: none | UI scaffold | Option present. |
| Sphericity correction: Greenhouse-Geisser | UI scaffold | Option present. |
| Sphericity correction: Huynh-Feldt | UI scaffold | Option present. |
| Levene's Test | UI scaffold | Option present. |
| Homogeneity tests | UI scaffold | Option present. |
| Q-Q plot of residuals | UI scaffold | Option present. |
| Residual diagnostics | UI scaffold | Option present. |
| Repeated measures contrasts: none/simple/deviation/difference/Helmert/polynomial/repeated | UI scaffold | List option present. |
| Contrast confidence intervals | UI scaffold | Option present. |
| Cohen's d for contrasts | UI scaffold | Option present. |
| Contrast effect sizes | UI scaffold | Option present. |
| Order restricted hypotheses syntax editor | UI scaffold | Text field present. |
| Multiple user-defined hypothesis models | UI scaffold | Toggle and model specification field present. |
| Model comparison and reference model | UI scaffold | Reference model, complement model, weight ratios, relative weights matrix, and coefficient comparison present. |
| Bootstrapping | UI scaffold | Toggle, sample size, and CI options present. |
| Post hoc repeated measures factors | UI scaffold | Option present. |
| Post hoc between-subject factors | UI scaffold | Option present. |
| Post hoc interactions | UI scaffold | Option present. |
| Holm/Bonferroni/Tukey/Scheffe/Sidak corrections | UI scaffold | Options present. |
| Post hoc confidence intervals | UI scaffold | Option present. |
| Significant comparison flags | UI scaffold | Option present. |
| Letter-based grouping tables | UI scaffold | Option present. |
| Post hoc effect sizes | UI scaffold | Option present. |
| Conditional interaction comparisons | UI scaffold | Option present. |
| Descriptive plots | UI scaffold | Horizontal axis, separate lines, separate plots, error bars, y-axis label, and average-unused-factors option present. |
| Bar plots | UI scaffold | Horizontal axis, separate plots, error bars, zero-axis, y-axis label, and average-unused-factors option present. |
| Raincloud plots | UI scaffold | Horizontal axis, separate plots, and y-axis label present. |
| Marginal means | UI scaffold | Terms, bootstrap, compare-to-zero, CI, adjustment, and pool-error options present. |
| Simple main effects | UI scaffold | Simple factor, two moderators, tests, pairwise comparisons, effect sizes, CIs, and pool error terms present. |
| Friedman-style repeated measures analysis | UI scaffold | Option present. |
| Conover post hoc tests | UI scaffold | Option present. |
| Raw/studentized/standardized residual export | UI scaffold | Options and custom column name present. |
| Prediction export | UI scaffold | Option and custom column name present. |
| APA reporting | UI scaffold | Narrative, tables, and assumption summaries present. |
| Teaching mode | UI scaffold | Sphericity, corrections, effect size, and interaction explanations present. |
| Publication mode | UI scaffold | Publication-ready tables/figures and Word/PDF export options present. |
| Added diagnostics | UI scaffold | Residual plots, influence diagnostics, Cook's distance, and Mahalanobis distance present. |
| Reproducibility | UI scaffold | Analysis syntax, R syntax, and configuration-save options present. |

## Mandatory Porting Rule

Feature parity is mandatory. Any missing JASP feature must be documented, justified, and accompanied by a recommended alternative or implementation plan.
