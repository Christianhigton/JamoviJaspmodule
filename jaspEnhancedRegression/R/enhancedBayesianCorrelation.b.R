#' @importFrom jmvcore .
enhancedBayesianCorrelationClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
    "enhancedBayesianCorrelationClass",
    inherit = enhancedBayesianCorrelationBase,
    private = list(.run = function() {
        self$results$status$setContent(.je_pending_html("JASP Enhanced Bayesian Correlation", "jaspRegression/R/correlationbayesian.R"))
        self$results$about$setContent(.je_about_html("https://github.com/jasp-stats/jaspRegression"))
    })
)
