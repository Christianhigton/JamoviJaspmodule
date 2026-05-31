#' @importFrom jmvcore .
enhancedLinearRegressionClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
    "enhancedLinearRegressionClass",
    inherit = enhancedLinearRegressionBase,
    private = list(.run = function() {
        self$results$status$setContent(.je_pending_html("JASP Enhanced Linear Regression", "jaspRegression/R/regressionlinear.R"))
        self$results$about$setContent(.je_about_html("https://github.com/jasp-stats/jaspRegression"))
    })
)
