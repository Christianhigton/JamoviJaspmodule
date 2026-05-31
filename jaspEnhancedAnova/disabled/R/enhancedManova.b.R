#' @importFrom jmvcore .
enhancedManovaClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
    "enhancedManovaClass",
    inherit = enhancedManovaBase,
    private = list(.run = function() {
        self$results$status$setContent(.je_pending_html("JASP Enhanced MANOVA", "jaspAnova/R/manova.R"))
        self$results$about$setContent(.je_about_html("https://github.com/jasp-stats/jaspAnova"))
    })
)
