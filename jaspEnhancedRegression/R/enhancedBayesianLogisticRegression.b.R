#' @importFrom jmvcore .
enhancedBayesianLogisticRegressionClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
    "enhancedBayesianLogisticRegressionClass",
    inherit = enhancedBayesianLogisticRegressionBase,
    private = list(.run = function() {
        self$results$status$setContent(.je_pending_html("JASP Enhanced Bayesian Logistic Regression", "jaspRegression/R/regressionlogisticbayesian.R"))
        self$results$about$setContent(.je_about_html("https://github.com/jasp-stats/jaspRegression"))
    })
)
