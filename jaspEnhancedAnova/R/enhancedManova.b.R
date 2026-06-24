#' @importFrom jmvcore .
enhancedManovaClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
    "enhancedManovaClass",
    inherit = enhancedManovaBase,
    private = list(
        .run = function() {
            self$results$about$setContent(.je_about_html("https://github.com/jasp-stats/jaspAnova"))

            deps    <- .je_chr_vec(self$options$dependent)
            factors <- .je_chr_vec(self$options$fixedFactors)
            random  <- .je_chr_vec(self$options$randomFactors)
            covs    <- .je_chr_vec(self$options$covariates)

            if (length(deps) < 2 || (length(factors) == 0 && length(random) == 0 && length(covs) == 0)) {
                self$results$status$setContent(
                    "<p>Select at least two dependent variables and at least one fixed factor, random factor, or covariate to run MANCOVA.</p>"
                )
                private$.setPending("Waiting for required variables.")
                return()
            }

            vars <- unique(c(deps, factors, random, covs))
            vars <- vars[nzchar(vars)]
            dat <- self$data[, vars, drop = FALSE]
            dat <- dat[stats::complete.cases(dat), , drop = FALSE]

            if (nrow(dat) < 4) {
                self$results$status$setContent("<p>Not enough complete cases to run MANCOVA (need at least 4).</p>")
                private$.setPending("Insufficient complete cases.")
                return()
            }

            all_factors <- unique(c(factors, random))
            for (f in all_factors) dat[[f]] <- as.factor(dat[[f]])

            rhs_terms <- .je_model_terms_from_options(self$options$modelTerms, all_factors, covs)
            rhs <- if (length(rhs_terms) == 0) {
                if (isTRUE(self$options$includeIntercept)) "1" else "0"
            } else {
                paste(rhs_terms, collapse = " + ")
            }
            if (!isTRUE(self$options$includeIntercept) && rhs != "0")
                rhs <- paste("0", rhs, sep = " + ")

            response <- paste0("cbind(", paste(deps, collapse = ", "), ")")
            formula <- stats::as.formula(paste(response, "~", rhs))
            fit <- tryCatch(stats::manova(formula, data = dat), error = function(e) e)

            if (inherits(fit, "error")) {
                self$results$status$setContent(paste0(
                    "<p>The MANCOVA model could not be fitted: ", .je_escape(fit$message), "</p>"
                ))
                private$.setPending("Model fit failed.")
                return()
            }

            tests <- .je_manova_tests(fit, self$options)
            univ <- .je_manova_univariate(fit, self$options)
            desc <- .je_manova_descriptives(dat, deps, all_factors, self$options)
            assumptions <- .je_manova_assumptions(dat, deps, all_factors, fit, self$options)
            apa <- .je_manova_apa(tests, deps, all_factors, covs)

            self$results$status$setContent(paste0(
                "<p><strong>Engine:</strong> True MANCOVA via <code>stats::manova()</code>, using ",
                "multivariate omnibus tests and per-dependent-variable univariate follow-up tables.</p>",
                if (length(random) > 0)
                    "<p><em>Random factors are included as model factors in this jamovi-native port; mixed-effects random variance components are not estimated in this pass.</em></p>"
                else "",
                "<p><strong>Dependencies for full JASP parity:</strong> ",
                "<code>biotools</code> (Box's M) - ", if (requireNamespace("biotools", quietly = TRUE)) "installed" else "<strong>missing; internal fallback used when possible</strong>",
                "; <code>mvnormtest</code> (multivariate normality) - ", if (requireNamespace("mvnormtest", quietly = TRUE)) "installed" else "<strong>missing; residual Shapiro fallback used</strong>",
                ".</p>"
            ))
            self$results$modelSummary$setContent(.je_manova_model_html(formula, nrow(dat), deps, all_factors, covs))
            self$results$multivariateTests$setContent(.je_manova_tests_html(tests))
            self$results$univariateTables$setContent(univ)
            self$results$descriptivesSection$setContent(desc)
            self$results$assumptionsSection$setContent(assumptions)
            self$results$apa$setContent(apa)
            self$results$reproducibility$setContent(.je_repro_html(formula, self$options))
        },

        .setPending = function(reason) {
            self$results$modelSummary$setContent(paste0("<p>", .je_escape(reason), "</p>"))
            self$results$multivariateTests$setContent(paste0("<p>", .je_escape(reason), "</p>"))
            self$results$univariateTables$setContent(paste0("<p>", .je_escape(reason), "</p>"))
            self$results$descriptivesSection$setContent(paste0("<p>", .je_escape(reason), "</p>"))
            self$results$assumptionsSection$setContent(paste0("<p>", .je_escape(reason), "</p>"))
            self$results$apa$setContent(paste0("<p>", .je_escape(reason), "</p>"))
            self$results$reproducibility$setContent(paste0("<p>", .je_escape(reason), "</p>"))
        }
    )
)

.je_manova_tests <- function(fit, options) {
    selected <- c(
        Pillai = isTRUE(options$testPillai),
        Wilks = isTRUE(options$testWilks),
        `Hotelling-Lawley` = isTRUE(options$testHotellingLawley),
        Roy = isTRUE(options$testRoy)
    )
    selected <- names(selected)[selected]
    if (length(selected) == 0) selected <- "Pillai"

    rows <- list()
    for (test in selected) {
        tab <- tryCatch(summary(fit, test = test)$stats, error = function(e) NULL)
        if (is.null(tab)) next
        df <- as.data.frame(tab, check.names = FALSE)
        df$Term <- rownames(df)
        rownames(df) <- NULL
        df <- df[df$Term != "Residuals", , drop = FALSE]
        if (nrow(df) == 0) next
        value_col <- setdiff(names(df), c("Term", "Df", "approx F", "num Df", "den Df", "Pr(>F)"))[1]
        out <- data.frame(
            Test = test,
            Term = df$Term,
            Df = df$Df,
            Statistic = if (!is.na(value_col)) .je_fmt(df[[value_col]]) else NA,
            `approx F` = .je_fmt(df$`approx F`),
            `num Df` = .je_fmt(df$`num Df`),
            `den Df` = .je_fmt(df$`den Df`),
            p = .je_p(df$`Pr(>F)`),
            check.names = FALSE
        )
        if (isTRUE(options$vovkSellke))
            out$`VS-MPR` <- .je_fmt(.je_vovk_sellke(df$`Pr(>F)`))
        rows[[length(rows) + 1]] <- out
    }

    if (length(rows) == 0) data.frame() else do.call(rbind, rows)
}

.je_manova_tests_html <- function(tests) {
    if (nrow(tests) == 0)
        return("<p>No multivariate tests could be calculated.</p>")
    .je_table_html(tests)
}

.je_manova_univariate <- function(fit, options) {
    if (!isTRUE(options$anovaTables))
        return("<p>Enable univariate ANOVA tables to show per-dependent-variable follow-up tests.</p>")

    tabs <- tryCatch(stats::summary.aov(fit), error = function(e) NULL)
    if (is.null(tabs))
        return("<p>Univariate follow-up tables could not be calculated.</p>")

    chunks <- character()
    for (dv in names(tabs)) {
        tab <- if (is.data.frame(tabs[[dv]]))
            as.data.frame(tabs[[dv]], check.names = FALSE)
        else
            as.data.frame(tabs[[dv]][[1]], check.names = FALSE)
        tab$Term <- rownames(tab)
        rownames(tab) <- NULL
        tab <- tab[, c("Term", setdiff(names(tab), "Term")), drop = FALSE]
        if ("Pr(>F)" %in% names(tab))
            tab$`Pr(>F)` <- .je_p(tab$`Pr(>F)`)
        for (nm in intersect(c("Sum Sq", "Mean Sq", "F value"), names(tab)))
            tab[[nm]] <- .je_fmt(tab[[nm]])
        chunks <- c(chunks, paste0("<h4>", .je_escape(dv), "</h4>", .je_table_html(tab)))
    }

    paste(chunks, collapse = "")
}

.je_manova_descriptives <- function(dat, deps, factors, options) {
    if (!isTRUE(options$descriptives))
        return("<p>Enable descriptive statistics to show dependent-variable summaries.</p>")

    rows <- list()
    group <- if (length(factors) > 0) interaction(dat[, factors, drop = FALSE], drop = TRUE, sep = " x ") else factor("Overall")
    for (dv in deps) {
        split_y <- split(dat[[dv]], group)
        for (g in names(split_y)) {
            x <- split_y[[g]]
            rows[[length(rows) + 1]] <- data.frame(
                Dependent = dv,
                Group = g,
                N = length(x),
                Mean = .je_fmt(mean(x)),
                SD = .je_fmt(stats::sd(x)),
                check.names = FALSE
            )
        }
    }

    .je_table_html(do.call(rbind, rows))
}

.je_manova_assumptions <- function(dat, deps, factors, fit, options) {
    chunks <- character()

    if (isTRUE(options$boxMTest)) {
        chunks <- c(chunks, .je_box_m_html(dat, deps, factors))
    } else {
        chunks <- c(chunks, "<p>Enable Box's M test to assess covariance homogeneity.</p>")
    }

    if (isTRUE(options$shapiroTest)) {
        chunks <- c(chunks, .je_multivariate_normality_html(fit, deps))
    } else {
        chunks <- c(chunks, "<p>Enable multivariate Shapiro-Wilk test to assess residual normality.</p>")
    }

    paste(chunks, collapse = "")
}

.je_box_m_html <- function(dat, deps, factors) {
    if (length(factors) == 0)
        return("<p>Box's M test requires at least one fixed factor.</p>")

    group <- interaction(dat[, factors, drop = FALSE], drop = TRUE)
    if (nlevels(group) < 2)
        return("<p>Box's M test requires at least two groups.</p>")

    if (requireNamespace("biotools", quietly = TRUE)) {
        bt <- tryCatch(biotools::boxM(dat[, deps, drop = FALSE], group), error = function(e) NULL)
        if (!is.null(bt)) {
            out <- data.frame(
                Test = "Box's M",
                Statistic = .je_fmt(unname(bt$statistic)),
                df = .je_fmt(unname(bt$parameter)),
                p = .je_p(bt$p.value),
                check.names = FALSE
            )
            return(paste0("<h4>Covariance homogeneity</h4>", .je_table_html(out)))
        }
    }

    res <- .je_box_m_fallback(dat[, deps, drop = FALSE], group)
    if (is.null(res))
        return("<p>Box's M test could not be calculated. Install <code>biotools</code> for the JASP-equivalent implementation.</p>")

    out <- data.frame(
        Test = "Box's M (fallback)",
        Statistic = .je_fmt(res$statistic),
        df = .je_fmt(res$df),
        p = .je_p(res$p),
        check.names = FALSE
    )
    paste0("<h4>Covariance homogeneity</h4>", .je_table_html(out))
}

.je_box_m_fallback <- function(y, group) {
    groups <- levels(group)
    p <- ncol(y)
    covs <- list()
    ns <- numeric()
    for (g in groups) {
        yg <- y[group == g, , drop = FALSE]
        if (nrow(yg) <= p) return(NULL)
        covs[[g]] <- stats::cov(yg)
        ns[g] <- nrow(yg)
    }

    pooled <- Reduce("+", Map(function(s, n) (n - 1) * s, covs, ns)) / (sum(ns) - length(groups))
    det_pooled <- determinant(pooled, logarithm = TRUE)$modulus[1]
    if (!is.finite(det_pooled)) return(NULL)
    m <- (sum(ns) - length(groups)) * det_pooled -
        sum(mapply(function(s, n) (n - 1) * determinant(s, logarithm = TRUE)$modulus[1], covs, ns))
    cfac <- ((2 * p^2 + 3 * p - 1) / (6 * (p + 1) * (length(groups) - 1))) *
        (sum(1 / (ns - 1)) - 1 / (sum(ns) - length(groups)))
    chi <- (1 - cfac) * m
    df <- (length(groups) - 1) * p * (p + 1) / 2
    list(statistic = chi, df = df, p = stats::pchisq(chi, df = df, lower.tail = FALSE))
}

.je_multivariate_normality_html <- function(fit, deps) {
    res <- tryCatch(stats::residuals(fit), error = function(e) NULL)
    if (is.null(res))
        return("<p>Residuals were not available for normality checks.</p>")
    res <- as.matrix(res)

    if (requireNamespace("mvnormtest", quietly = TRUE) && nrow(res) >= 3 && nrow(res) <= 5000) {
        mt <- tryCatch(mvnormtest::mshapiro.test(t(res)), error = function(e) NULL)
        if (!is.null(mt)) {
            out <- data.frame(
                Test = "Multivariate Shapiro-Wilk",
                W = .je_fmt(unname(mt$statistic)),
                p = .je_p(mt$p.value),
                check.names = FALSE
            )
            return(paste0("<h4>Residual normality</h4>", .je_table_html(out)))
        }
    }

    rows <- list()
    for (i in seq_along(deps)) {
        x <- res[, i]
        if (length(x) < 3 || length(x) > 5000) next
        sh <- stats::shapiro.test(x)
        rows[[length(rows) + 1]] <- data.frame(
            Dependent = deps[i],
            W = .je_fmt(unname(sh$statistic)),
            p = .je_p(sh$p.value),
            check.names = FALSE
        )
    }
    if (length(rows) == 0)
        return("<p>Residual Shapiro-Wilk checks skipped (n outside 3-5000).</p>")
    paste0(
        "<h4>Residual normality</h4>",
        "<p><em>Install <code>mvnormtest</code> for a multivariate Shapiro-Wilk test. Per-dependent-variable residual checks are shown below.</em></p>",
        .je_table_html(do.call(rbind, rows))
    )
}

.je_manova_model_html <- function(formula, n, deps, factors, covs) {
    data.frame(
        Item = c("Formula", "Complete cases", "Dependent variables", "Factors", "Covariates"),
        Value = c(
            paste(deparse(formula), collapse = " "),
            n,
            paste(deps, collapse = ", "),
            if (length(factors) > 0) paste(factors, collapse = ", ") else "None",
            if (length(covs) > 0) paste(covs, collapse = ", ") else "None"
        ),
        check.names = FALSE
    ) |>
        .je_table_html()
}

.je_manova_apa <- function(tests, deps, factors, covs) {
    if (nrow(tests) == 0)
        return("<p>No APA statement is available because no multivariate tests were calculated.</p>")

    first <- tests[1, , drop = FALSE]
    paste0(
        "<p>A MANCOVA was conducted with dependent variables ",
        .je_escape(paste(deps, collapse = ", ")),
        if (length(factors) > 0) paste0(", factor(s) ", .je_escape(paste(factors, collapse = ", "))) else "",
        if (length(covs) > 0) paste0(", and covariate(s) ", .je_escape(paste(covs, collapse = ", "))) else "",
        ". The first reported multivariate test was ", .je_escape(first$Test),
        " for ", .je_escape(first$Term), ", F(", .je_escape(first$`num Df`),
        ", ", .je_escape(first$`den Df`), ") = ", .je_escape(first$`approx F`),
        ", p = ", .je_escape(first$p), ".</p>"
    )
}
