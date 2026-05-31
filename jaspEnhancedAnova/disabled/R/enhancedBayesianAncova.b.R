#' @importFrom jmvcore .
enhancedBayesianAncovaClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
    "enhancedBayesianAncovaClass",
    inherit = enhancedBayesianAncovaBase,
    private = list(.run = function() {
        self$results$status$setContent(.je_pending_html("JASP Enhanced Bayesian ANCOVA", "jaspAnova/R/ancovabayesian.R"))
        self$results$about$setContent(.je_about_html("https://github.com/jasp-stats/jaspAnova"))
    })
)
