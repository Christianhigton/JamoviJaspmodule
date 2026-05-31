#' @importFrom jmvcore .
.je_about_html <- function(source_repo) {
    paste0(
        "<h3>JASP Enhanced for jamovi</h3>",
        "<p>This module adapts functionality from open-source JASP analyses ",
        "while providing a jamovi-native user experience.</p>",
        "<p>Designed for users who prefer the workflow and usability of jamovi ",
        "while retaining the extensive analytical capabilities available in JASP.</p>",
        "<p>Please cite both jamovi and JASP when using this module.</p>",
        "<h4>References</h4>",
        "<p>The jamovi project. jamovi [Computer Software]. ",
        "https://www.jamovi.org</p>",
        "<p>JASP Team. JASP [Computer Software]. https://jasp-stats.org</p>",
        "<p>Adapted from: ", source_repo, "</p>"
    )
}

.je_pending_html <- function(analysis, upstream) {
    paste0(
        "<p><strong>", analysis, " is scaffolded but not implemented yet.</strong></p>",
        "<p>The upstream JASP implementation has been identified at <code>",
        upstream, "</code>. The next step is to port options, calculations, ",
        "tables, plots, diagnostics, Bayesian behavior where applicable, and ",
        "verification tests without dropping functionality.</p>"
    )
}
