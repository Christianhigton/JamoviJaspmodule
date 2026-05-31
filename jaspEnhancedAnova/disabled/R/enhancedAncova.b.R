#' @importFrom jmvcore .
enhancedAncovaClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
    "enhancedAncovaClass",
    inherit = enhancedAncovaBase,
    private = list(.run = function() {
        self$results$status$setContent(.je_pending_html("JASP Enhanced ANCOVA", "jaspAnova/R/ancova.R"))
        self$results$about$setContent(.je_about_html("https://github.com/jasp-stats/jaspAnova"))
    })
)
