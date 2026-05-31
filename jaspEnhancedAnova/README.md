# JASP Enhanced ANOVA for jamovi

This is an independent jamovi module for adapting selected open-source JASP ANOVA functionality into a jamovi-native workflow.

Design principle: JASP functionality. jamovi usability.

Current status: focused development build. Only the two active analyses below are exposed in jamovi:

- JASP Enhanced ANOVA
- JASP Enhanced Repeated Measures ANOVA

Other analysis shells have been disabled for now and moved to `disabled/` so they can be restored later without cluttering the jamovi menu.

Upstream reference: https://github.com/jasp-stats/jaspAnova

License: GPL (>= 2), matching the upstream JASP module license declared in its DESCRIPTION file.

## Active Analyses

### JASP Enhanced ANOVA

Implemented for the current build:

- classical fixed-effects ANOVA
- WLS weights
- descriptives
- eta/partial eta/omega/partial omega effect sizes
- residual diagnostics
- custom one-factor contrast inference from numeric contrast rows
- post hoc comparisons
- marginal means and simple effects
- Kruskal-Wallis fallback
- Q-Q, residual, and raincloud-style plots
- generated residual and prediction columns
- APA/reproducibility text

### JASP Enhanced Repeated Measures ANOVA

Implemented for the current build:

- wide-format repeated-measures analysis using selected cells as one within-subject factor
- repeated-measures descriptives
- within-factor contrast inference
- post hoc paired comparisons
- marginal means and simple effects by between-subject factor
- Friedman test
- Q-Q, residual, and raincloud-style plots
- generated subject-level mean residual and prediction columns
- APA/reproducibility text

## Disabled Analyses

The following analysis shells are retained in `disabled/` but are not currently built into the jamovi module:

- ANCOVA
- MANOVA
- Bayesian ANOVA
- Bayesian Repeated Measures ANOVA
- Bayesian ANCOVA

To re-enable one later, move its `.a.yaml`, `.u.yaml`, and `.r.yaml` files from `disabled/jamovi/` back into `jamovi/`, and move its `.b.R` implementation from `disabled/R/` back into `R/`. Then run the build command in `DEPLOYMENT.md`.

## Known Limitations

This is not yet full JASP parity. Remaining work includes exact JASP multi-factor custom contrast grammar, full multi-factor/nested repeated-measures construction, Bayesian/order-restricted Bayes factors, complete sphericity epsilon reporting, and independent verification against JASP fixtures.

## Attribution

This module adapts functionality from open-source JASP analyses while providing a jamovi-native user experience. Please cite both jamovi and JASP when using this module.

The jamovi project. jamovi [Computer Software]. https://www.jamovi.org

JASP Team. JASP [Computer Software]. https://jasp-stats.org
