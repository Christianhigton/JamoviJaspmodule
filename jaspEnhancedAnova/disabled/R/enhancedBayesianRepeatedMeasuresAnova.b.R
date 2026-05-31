#' @importFrom jmvcore .
enhancedBayesianRepeatedMeasuresAnovaClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
    "enhancedBayesianRepeatedMeasuresAnovaClass",
    inherit = enhancedBayesianRepeatedMeasuresAnovaBase,
    private = list(.run = function() {
        self$results$status$setContent(.je_pending_html("JASP Enhanced Bayesian Repeated Measures ANOVA", "jaspAnova/R/anovarepeatedmeasuresbayesian.R"))
        self$results$about$setContent(.je_about_html("https://github.com/jasp-stats/jaspAnova"))
    })
)
