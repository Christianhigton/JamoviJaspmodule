#' @importFrom jmvcore .
enhancedBayesianAnovaClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
    "enhancedBayesianAnovaClass",
    inherit = enhancedBayesianAnovaBase,
    private = list(.run = function() {
        self$results$status$setContent(.je_pending_html("JASP Enhanced Bayesian ANOVA", "jaspAnova/R/anovabayesian.R"))
        self$results$about$setContent(.je_about_html("https://github.com/jasp-stats/jaspAnova"))
    })
)
