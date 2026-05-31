#' @importFrom jmvcore .
enhancedCorrelationClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
    "enhancedCorrelationClass",
    inherit = enhancedCorrelationBase,
    private = list(.run = function() {
        self$results$status$setContent(.je_pending_html("JASP Enhanced Correlation", "jaspRegression/R/correlation.R"))
        self$results$about$setContent(.je_about_html("https://github.com/jasp-stats/jaspRegression"))
    })
)
