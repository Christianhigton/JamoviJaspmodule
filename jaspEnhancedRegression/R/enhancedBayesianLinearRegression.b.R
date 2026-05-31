#' @importFrom jmvcore .
enhancedBayesianLinearRegressionClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
    "enhancedBayesianLinearRegressionClass",
    inherit = enhancedBayesianLinearRegressionBase,
    private = list(.run = function() {
        self$results$status$setContent(.je_pending_html("JASP Enhanced Bayesian Linear Regression", "jaspRegression/R/regressionlinearbayesian.R"))
        self$results$about$setContent(.je_about_html("https://github.com/jasp-stats/jaspRegression"))
    })
)
