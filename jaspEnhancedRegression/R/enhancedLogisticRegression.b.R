#' @importFrom jmvcore .
enhancedLogisticRegressionClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
    "enhancedLogisticRegressionClass",
    inherit = enhancedLogisticRegressionBase,
    private = list(.run = function() {
        self$results$status$setContent(.je_pending_html("JASP Enhanced Logistic Regression", "jaspRegression/R/regressionlogistic.R"))
        self$results$about$setContent(.je_about_html("https://github.com/jasp-stats/jaspRegression"))
    })
)
