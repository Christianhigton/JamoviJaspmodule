#' @importFrom jmvcore .
enhancedGeneralizedLinearModelClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
    "enhancedGeneralizedLinearModelClass",
    inherit = enhancedGeneralizedLinearModelBase,
    private = list(.run = function() {
        self$results$status$setContent(.je_pending_html("JASP Enhanced Generalized Linear Model", "jaspRegression/R/generalizedlinearmodel.R"))
        self$results$about$setContent(.je_about_html("https://github.com/jasp-stats/jaspRegression"))
    })
)
