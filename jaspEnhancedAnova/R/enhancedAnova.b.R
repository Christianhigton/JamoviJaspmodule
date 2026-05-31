#' @importFrom jmvcore .
enhancedAnovaClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
    "enhancedAnovaClass",
    inherit = enhancedAnovaBase,
    private = list(
        .plotState = NULL,
        .run = function() {
            self$results$about$setContent(.je_about_html(
                "https://github.com/jasp-stats/jaspAnova"
            ))

            dep <- .je_chr(self$options$dependent)
            factors <- .je_chr_vec(self$options$fixedFactors)
            covariates <- .je_chr_vec(self$options$covariates)
            weights <- .je_chr(self$options$wlsWeights)

            if (is.null(dep) || dep == "" || (length(factors) == 0 && length(covariates) == 0)) {
                self$results$status$setContent(
                    "<p>Select a dependent variable and at least one fixed factor or covariate to run the implemented classical ANOVA engine.</p>"
                )
                .je_set_pending_outputs(self, "Waiting for required variables.")
                return()
            }

            vars <- unique(c(dep, factors, covariates, weights))
            vars <- vars[nzchar(vars)]
            dat <- self$data[, vars, drop = FALSE]
            names(dat) <- vars
            dat <- dat[stats::complete.cases(dat), , drop = FALSE]
            row_nums <- as.integer(rownames(dat))
            if (anyNA(row_nums))
                row_nums <- seq_len(nrow(dat))

            if (nrow(dat) < 3) {
                self$results$status$setContent("<p>There are not enough complete cases to run ANOVA.</p>")
                .je_set_pending_outputs(self, "Insufficient complete cases.")
                return()
            }

            for (f in factors)
                dat[[f]] <- as.factor(dat[[f]])

            rhs <- c(factors, covariates)
            formula <- stats::as.formula(paste(dep, "~", paste(rhs, collapse = " * ")))

            fit <- tryCatch({
                if (!is.null(weights) && weights %in% names(dat)) {
                    stats::lm(formula, data = dat, weights = dat[[weights]])
                } else {
                    stats::lm(formula, data = dat)
                }
            }, error = function(e) e)

            if (inherits(fit, "error")) {
                self$results$status$setContent(paste0(
                    "<p>The ANOVA model could not be fitted: ",
                    .je_escape(fit$message),
                    "</p>"
                ))
                .je_set_pending_outputs(self, "Model fit failed.")
                return()
            }

            aov_tab <- as.data.frame(stats::anova(fit))
            aov_tab$Term <- rownames(aov_tab)
            rownames(aov_tab) <- NULL

            effect_tab <- .je_effect_sizes(aov_tab)
            desc_tab <- .je_descriptives(dat, dep, factors)
            assumptions <- .je_assumptions(dat, dep, factors, fit)
            posthoc <- .je_posthoc(dat, dep, factors, self$options)
            kruskal <- .je_kruskal(dat, dep, factors)
            marginal <- .je_marginal_means(dat, dep, factors, covariates, self$options)
            simple <- .je_simple_effects(dat, dep, factors, self$options)
            contrast <- .je_contrasts(dat, dep, factors, fit, self$options)
            bootstrap <- .je_bootstrap_effects(dat, dep, factors, covariates, self$options)
            apa <- .je_apa(aov_tab, effect_tab, dep, factors)

            self$results$status$setContent(
                paste0(
                    "<p><strong>Implemented:</strong> classical fixed-effects ANOVA via R <code>lm()</code>/<code>anova()</code>, ",
                    "descriptives, core effect sizes, residual diagnostics, one-way Tukey/post hoc tests, ",
                    "custom contrast inference for one-factor contrasts, marginal means, simple effects, bootstrap effect-size intervals, ",
                    "Kruskal-Wallis fallback, generated residual/prediction columns, raincloud-style plots, and APA draft text.</p>",
                    "<p><strong>Still pending for full JASP parity:</strong> Bayesian/order-restricted Bayes factors and exact JASP multi-factor custom contrast grammar.</p>"
                )
            )

            self$results$modelSummary$setContent(.je_model_html(aov_tab, formula, nrow(dat), weights))
            self$results$descriptivesSection$setContent(.je_table_html(desc_tab))
            self$results$effectSizesSection$setContent(.je_table_html(effect_tab))
            self$results$assumptionsSection$setContent(assumptions)
            self$results$contrasts$setContent(contrast)
            self$results$orderRestrictions$setContent(.je_order_restricted_status(self$options))
            self$results$postHocSection$setContent(posthoc)
            self$results$marginalMeansSection$setContent(marginal)
            self$results$simpleEffectsSection$setContent(simple)
            self$results$nonparametric$setContent(kruskal)
            self$results$plots$setContent(.je_plot_status(dat, dep, factors, self$options))
            self$results$savedColumns$setContent(.je_saved_columns_html(self$options))
            self$results$apa$setContent(paste0(apa, bootstrap))
            self$results$teaching$setContent(.je_teaching_html(effect_tab, assumptions))
            self$results$publication$setContent(.je_publication_html(self$options))
            self$results$reproducibility$setContent(.je_repro_html(formula, self$options))

            private$.plotState <- list(data = dat, dep = dep, factors = factors, fit = fit)
            if (isTRUE(self$options$qqPlotResiduals))
                self$results$qqPlot$setState(private$.plotState)
            if (isTRUE(self$options$residualPlots))
                self$results$residualPlot$setState(private$.plotState)
            if (isTRUE(self$options$raincloudPlots))
                self$results$raincloudPlot$setState(private$.plotState)

            private$.populateOutputs(fit, row_nums)
        },
        .populateOutputs = function(fit, row_nums) {
            if (isTRUE(self$options$saveResiduals) && isTRUE(self$options$saveRawResiduals) &&
                self$options$residsOV && self$results$residsOV$isNotFilled()) {
                self$results$residsOV$setRowNums(row_nums)
                self$results$residsOV$setValues(stats::residuals(fit))
            }
            if (isTRUE(self$options$saveResiduals) && isTRUE(self$options$saveStudentizedResiduals) &&
                self$options$studentizedResidsOV && self$results$studentizedResidsOV$isNotFilled()) {
                self$results$studentizedResidsOV$setRowNums(row_nums)
                self$results$studentizedResidsOV$setValues(stats::rstudent(fit))
            }
            if (isTRUE(self$options$saveResiduals) && isTRUE(self$options$saveStandardizedResiduals) &&
                self$options$standardizedResidsOV && self$results$standardizedResidsOV$isNotFilled()) {
                self$results$standardizedResidsOV$setRowNums(row_nums)
                self$results$standardizedResidsOV$setValues(stats::rstandard(fit))
            }
            if (isTRUE(self$options$savePredictions) && self$options$predictOV &&
                self$results$predictOV$isNotFilled()) {
                self$results$predictOV$setRowNums(row_nums)
                self$results$predictOV$setValues(stats::fitted(fit))
            }
        },
        .qqPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state) || !requireNamespace("ggplot2", quietly = TRUE))
                return(FALSE)
            res <- stats::rstandard(image$state$fit)
            df <- as.data.frame(stats::qqnorm(res, plot.it = FALSE))
            ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
                ggplot2::geom_abline(slope = 1, intercept = 0, colour = theme$color[1]) +
                ggplot2::geom_point(colour = theme$color[1]) +
                ggplot2::labs(x = "Theoretical Quantiles", y = "Standardized Residuals") +
                ggtheme
        },
        .residualPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state) || !requireNamespace("ggplot2", quietly = TRUE))
                return(FALSE)
            df <- data.frame(Fitted = stats::fitted(image$state$fit), Residuals = stats::rstandard(image$state$fit))
            ggplot2::ggplot(df, ggplot2::aes(x = Fitted, y = Residuals)) +
                ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = theme$color[1]) +
                ggplot2::geom_point(colour = theme$color[1], alpha = 0.75) +
                ggplot2::labs(x = "Fitted Values", y = "Standardized Residuals") +
                ggtheme
        },
        .raincloudPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state) || !requireNamespace("ggplot2", quietly = TRUE))
                return(FALSE)
            dat <- image$state$data
            dep <- image$state$dep
            factors <- image$state$factors
            if (length(factors) == 0)
                return(FALSE)
            x <- factors[1]
            dat[[x]] <- as.factor(dat[[x]])
            ggplot2::ggplot(dat, ggplot2::aes_string(x = x, y = dep, fill = x)) +
                ggplot2::geom_violin(trim = FALSE, alpha = 0.35, colour = NA) +
                ggplot2::geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.65) +
                ggplot2::geom_jitter(width = 0.08, alpha = 0.55, size = 1.6) +
                ggplot2::labs(x = x, y = dep) +
                ggplot2::guides(fill = "none") +
                ggtheme
        }
    )
)

.je_chr <- function(x) {
    if (is.null(x) || length(x) == 0)
        return(NULL)
    as.character(x[[1]])
}

.je_chr_vec <- function(x) {
    if (is.null(x) || length(x) == 0)
        return(character())
    x <- as.character(unlist(x, use.names = FALSE))
    x[nzchar(x)]
}

.je_escape <- function(x) {
    x <- as.character(x)
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x
}

.je_fmt <- function(x, digits = 3) {
    ifelse(is.na(x), "", formatC(x, digits = digits, format = "f"))
}

.je_p <- function(p) {
    ifelse(is.na(p), "", ifelse(p < .001, "&lt; .001", sub("^0", "", formatC(p, digits = 3, format = "f"))))
}

.je_table_html <- function(dat) {
    if (is.null(dat) || nrow(dat) == 0)
        return("<p>No results available for the current options.</p>")
    headers <- paste0("<th>", .je_escape(names(dat)), "</th>", collapse = "")
    rows <- apply(dat, 1, function(row) {
        paste0("<tr>", paste0("<td>", .je_escape(row), "</td>", collapse = ""), "</tr>")
    })
    paste0("<table><thead><tr>", headers, "</tr></thead><tbody>", paste(rows, collapse = ""), "</tbody></table>")
}

.je_model_html <- function(aov_tab, formula, n, weights) {
    dat <- data.frame(
        Term = aov_tab$Term,
        df = .je_fmt(aov_tab$Df),
        `Sum Sq` = .je_fmt(aov_tab$`Sum Sq`),
        `Mean Sq` = .je_fmt(aov_tab$`Mean Sq`),
        F = .je_fmt(aov_tab$`F value`),
        p = .je_p(aov_tab$`Pr(>F)`),
        check.names = FALSE
    )
    weight_note <- if (!is.null(weights) && nzchar(weights)) paste0("<p>WLS weights: <code>", .je_escape(weights), "</code></p>") else ""
    paste0(
        "<p>Model: <code>", .je_escape(deparse(formula)), "</code></p>",
        "<p>Complete cases analysed: ", n, "</p>",
        weight_note,
        .je_table_html(dat)
    )
}

.je_descriptives <- function(dat, dep, factors) {
    if (length(factors) == 0) {
        x <- dat[[dep]]
        return(data.frame(Group = "Overall", N = length(x), Mean = .je_fmt(mean(x)), SD = .je_fmt(stats::sd(x)), check.names = FALSE))
    }

    group <- interaction(dat[, factors, drop = FALSE], drop = TRUE, sep = " : ")
    split_y <- split(dat[[dep]], group)
    data.frame(
        Group = names(split_y),
        N = vapply(split_y, length, integer(1)),
        Mean = .je_fmt(vapply(split_y, mean, numeric(1))),
        SD = .je_fmt(vapply(split_y, stats::sd, numeric(1))),
        check.names = FALSE
    )
}

.je_effect_sizes <- function(aov_tab) {
    resid_row <- which(aov_tab$Term == "Residuals")
    if (length(resid_row) == 0)
        return(data.frame())
    sse <- aov_tab$`Sum Sq`[resid_row]
    mse <- aov_tab$`Mean Sq`[resid_row]
    ss_total <- sum(aov_tab$`Sum Sq`, na.rm = TRUE)
    terms <- aov_tab[-resid_row, , drop = FALSE]
    if (nrow(terms) == 0)
        return(data.frame())

    eta <- terms$`Sum Sq` / ss_total
    partial_eta <- terms$`Sum Sq` / (terms$`Sum Sq` + sse)
    omega <- pmax(0, (terms$`Sum Sq` - terms$Df * mse) / (ss_total + mse))
    partial_omega <- pmax(0, (terms$`Sum Sq` - terms$Df * mse) / (terms$`Sum Sq` + sse + mse))

    data.frame(
        Term = terms$Term,
        `eta squared` = .je_fmt(eta),
        `partial eta squared` = .je_fmt(partial_eta),
        `omega squared` = .je_fmt(omega),
        `partial omega squared` = .je_fmt(partial_omega),
        check.names = FALSE
    )
}

.je_assumptions <- function(dat, dep, factors, fit) {
    out <- character()
    res <- stats::residuals(fit)
    if (length(res) >= 3 && length(res) <= 5000) {
        sh <- stats::shapiro.test(res)
        out <- c(out, paste0("<p>Residual normality: Shapiro-Wilk W = ", .je_fmt(unname(sh$statistic)), ", p = ", .je_p(sh$p.value), ".</p>"))
    } else {
        out <- c(out, "<p>Residual normality: Shapiro-Wilk test skipped because the sample size is outside the supported range.</p>")
    }

    if (length(factors) > 0) {
        group <- interaction(dat[, factors, drop = FALSE], drop = TRUE, sep = " : ")
        med <- stats::ave(dat[[dep]], group, FUN = stats::median)
        lev_fit <- stats::lm(abs(dat[[dep]] - med) ~ group)
        lev <- as.data.frame(stats::anova(lev_fit))
        p <- lev$`Pr(>F)`[1]
        out <- c(out, paste0("<p>Homogeneity: Brown-Forsythe/Levene-style median test F = ", .je_fmt(lev$`F value`[1]), ", p = ", .je_p(p), ".</p>"))
    } else {
        out <- c(out, "<p>Homogeneity: no fixed factor was supplied, so groupwise homogeneity was not tested.</p>")
    }

    paste(out, collapse = "")
}

.je_posthoc <- function(dat, dep, factors, options) {
    if (!isTRUE(options$postHoc))
        return("<p>Post hoc tests are disabled.</p>")
    if (length(factors) == 0)
        return("<p>Post hoc tests require at least one fixed factor.</p>")
    if (length(factors) > 2)
        return("<p>Post hoc testing currently supports one or two fixed factors. Higher-order interaction post hoc tests are still pending.</p>")

    chunks <- character()

    for (factor in factors) {
        formula <- stats::as.formula(paste(dep, "~", factor))
        fit <- stats::aov(formula, data = dat)
        tk <- tryCatch(stats::TukeyHSD(fit), error = function(e) e)
        if (!inherits(tk, "error")) {
            tab <- as.data.frame(tk[[1]])
            tab$Comparison <- rownames(tab)
            rownames(tab) <- NULL
            dat_out <- data.frame(
                Factor = factor,
                Comparison = tab$Comparison,
                Difference = .je_fmt(tab$diff),
                Lower = .je_fmt(tab$lwr),
                Upper = .je_fmt(tab$upr),
                p = .je_p(tab$`p adj`),
                check.names = FALSE
            )
            chunks <- c(chunks, paste0("<h4>Tukey HSD: ", .je_escape(factor), "</h4>", .je_table_html(dat_out)))
        }

        pw <- tryCatch(stats::pairwise.t.test(dat[[dep]], dat[[factor]], p.adjust.method = .je_posthoc_adjust(options)), error = function(e) e)
        if (!inherits(pw, "error")) {
            chunks <- c(chunks, paste0("<h4>Pairwise t tests: ", .je_escape(factor), " (", .je_escape(.je_posthoc_adjust(options)), ")</h4>", .je_pairwise_html(pw$p.value)))
        }
    }

    if (length(factors) == 2 && isTRUE(options$conditionalComparisons)) {
        chunks <- c(chunks, .je_conditional_posthoc(dat, dep, factors, options))
    }

    if (length(chunks) == 0)
        return("<p>No post hoc output could be calculated for the current model.</p>")

    paste0(
        "<p>Implemented post hoc output includes Tukey HSD and adjusted pairwise t tests for main-effect factors. ",
        "Scheffe, Dunnett, Games-Howell, and letter groupings remain pending for exact JASP parity.</p>",
        paste(chunks, collapse = "")
    )
}

.je_kruskal <- function(dat, dep, factors) {
    if (length(factors) != 1)
        return("<p>Kruskal-Wallis is currently implemented for one grouping factor only.</p>")
    kt <- stats::kruskal.test(stats::as.formula(paste(dep, "~", factors[1])), data = dat)
    paste0("<p>Kruskal-Wallis chi-square = ", .je_fmt(unname(kt$statistic)), ", df = ", kt$parameter, ", p = ", .je_p(kt$p.value), ".</p>")
}

.je_apa <- function(aov_tab, effect_tab, dep, factors) {
    test_rows <- aov_tab[aov_tab$Term != "Residuals", , drop = FALSE]
    if (nrow(test_rows) == 0)
        return("<p>No APA text is available for this model.</p>")
    row <- test_rows[1, ]
    eta <- if (nrow(effect_tab) > 0) effect_tab$`partial eta squared`[1] else ""
    fac <- if (length(factors) > 0) paste(factors, collapse = ", ") else "the model"
    paste0(
        "<p>A classical ANOVA tested the effect of ", .je_escape(fac), " on ", .je_escape(dep), ". ",
        "The first model term was F(", .je_fmt(row$Df), ", ", .je_fmt(aov_tab$Df[aov_tab$Term == "Residuals"][1]), ") = ",
        .je_fmt(row$`F value`), ", p = ", .je_p(row$`Pr(>F)`), ", partial eta squared = ", eta, ".</p>"
    )
}

.je_teaching_html <- function(effect_tab, assumptions) {
    paste0(
        "<p>The model summary tests whether mean differences are larger than expected from residual variation. ",
        "Partial eta squared describes the proportion of explainable variance associated with a term after accounting for residual error.</p>",
        assumptions
    )
}

.je_marginal_means <- function(dat, dep, factors, covariates, options) {
    if (!isTRUE(options$marginalMeans))
        return("<p>Estimated marginal means are disabled.</p>")
    terms <- .je_chr_vec(options$marginalMeanTerms)
    terms <- intersect(terms, factors)
    if (length(terms) == 0)
        terms <- factors
    if (length(terms) == 0)
        return("<p>Estimated marginal means require at least one fixed factor.</p>")

    group <- interaction(dat[, terms, drop = FALSE], drop = TRUE, sep = " : ")
    split_y <- split(dat[[dep]], group)
    alpha <- 1 - options$ciWidth / 100
    out <- lapply(names(split_y), function(g) {
        x <- split_y[[g]]
        se <- stats::sd(x) / sqrt(length(x))
        crit <- stats::qt(1 - alpha / 2, df = max(1, length(x) - 1))
        data.frame(
            Term = paste(terms, collapse = " : "),
            Level = g,
            N = length(x),
            Mean = .je_fmt(mean(x)),
            SE = .je_fmt(se),
            Lower = .je_fmt(mean(x) - crit * se),
            Upper = .je_fmt(mean(x) + crit * se),
            check.names = FALSE
        )
    })
    html <- paste0("<p>Marginal means are computed as observed cell means in this build. Covariate-adjusted emmeans-style estimates require an emmeans-backed implementation.</p>", .je_table_html(do.call(rbind, out)))

    if (isTRUE(options$marginalMeansPairwise) && length(terms) == 1) {
        adj <- .je_adjust_method(options$marginalMeansCiAdjustment)
        pw <- tryCatch(stats::pairwise.t.test(dat[[dep]], dat[[terms]], p.adjust.method = adj), error = function(e) e)
        if (!inherits(pw, "error"))
            html <- paste0(html, "<h4>Pairwise marginal mean comparisons</h4><p>Adjustment used: ", .je_escape(adj), ".</p>", .je_pairwise_html(pw$p.value))
    }
    html
}

.je_simple_effects <- function(dat, dep, factors, options) {
    if (!isTRUE(options$simpleEffects))
        return("<p>Simple effects analysis is disabled.</p>")
    simple <- .je_chr(options$simpleEffectFactor)
    mod1 <- .je_chr(options$moderatorFactor1)
    mod2 <- .je_chr(options$moderatorFactor2)
    if (is.null(simple) || !(simple %in% factors))
        return("<p>Select a simple effect factor from the fixed factors.</p>")
    moderators <- intersect(c(mod1, mod2), factors)
    moderators <- setdiff(moderators[nzchar(moderators)], simple)
    if (length(moderators) == 0)
        return("<p>Select at least one moderator factor for simple effects.</p>")

    strata <- interaction(dat[, moderators, drop = FALSE], drop = TRUE, sep = " : ")
    pieces <- lapply(levels(strata), function(level) {
        sub <- dat[strata == level, , drop = FALSE]
        if (length(unique(sub[[simple]])) < 2 || nrow(sub) < 3)
            return(paste0("<h4>", .je_escape(level), "</h4><p>Not enough data for this stratum.</p>"))
        fit <- stats::aov(stats::as.formula(paste(dep, "~", simple)), data = sub)
        tab <- as.data.frame(stats::anova(fit))
        tab$Term <- rownames(tab)
        rownames(tab) <- NULL
        html <- paste0("<h4>", .je_escape(level), "</h4>", .je_model_html(tab, stats::formula(fit), nrow(sub), NULL))
        if (isTRUE(options$simpleEffectsPostHoc) || isTRUE(options$simpleEffectsEffectSizes)) {
            html <- paste0(html, .je_posthoc(sub, dep, simple, options))
        }
        html
    })
    paste0("<p>Simple effects are implemented as one-way ANOVAs within moderator strata.</p>", paste(pieces, collapse = ""))
}

.je_contrasts <- function(dat, dep, factors, fit, options) {
    requested <- isTRUE(options$plannedContrasts) || isTRUE(options$polynomialContrasts) ||
        isTRUE(options$helmertContrasts) || isTRUE(options$differenceContrasts) ||
        isTRUE(options$deviationContrasts) || isTRUE(options$repeatedContrasts) ||
        isTRUE(options$customContrasts)
    if (!requested)
        return("<p>No contrast options are enabled.</p>")
    if (length(factors) == 0)
        return("<p>Contrasts require at least one fixed factor.</p>")

    chunks <- character()
    for (factor in factors) {
        levs <- levels(dat[[factor]])
        mats <- list()
        if (isTRUE(options$helmertContrasts))
            mats$Helmert <- stats::contr.helmert(length(levs))
        if (isTRUE(options$polynomialContrasts))
            mats$Polynomial <- stats::contr.poly(length(levs))
        if (isTRUE(options$deviationContrasts))
            mats$Deviation <- stats::contr.sum(length(levs))
        if (isTRUE(options$differenceContrasts) || isTRUE(options$repeatedContrasts))
            mats$Difference <- .je_contr_diff(length(levs))

        if (isTRUE(options$customContrasts)) {
            custom <- .je_parse_custom_contrasts(options$contrastSyntax, length(levs))
            if (nrow(custom) > 0) {
                mats$Custom <- t(custom)
                chunks <- c(chunks, .je_custom_contrast_tests(dat, dep, factor, custom, fit))
            } else {
                chunks <- c(chunks, paste0(
                    "<h4>Custom contrasts: ", .je_escape(factor), "</h4>",
                    "<p>Enter one contrast per line with ", length(levs), " numeric weights matching levels: ",
                    .je_escape(paste(levs, collapse = ", ")), ".</p>"
                ))
            }
        }
        for (nm in names(mats)) {
            mat <- as.data.frame(mats[[nm]])
            names(mat) <- paste0("C", seq_len(ncol(mat)))
            mat$Level <- levs
            mat <- mat[, c("Level", setdiff(names(mat), "Level")), drop = FALSE]
            chunks <- c(chunks, paste0("<h4>", .je_escape(nm), " contrasts: ", .je_escape(factor), "</h4>", .je_table_html(mat)))
        }
    }
    paste0("<p>Contrast matrices are generated for built-in contrast types. Custom one-factor contrast inference is implemented using the ANOVA residual MSE.</p>", paste(chunks, collapse = ""))
}

.je_parse_custom_contrasts <- function(text, n_levels) {
    if (is.null(text) || !nzchar(text))
        return(matrix(numeric(), nrow = 0, ncol = n_levels))
    lines <- unlist(strsplit(text, "\n", fixed = TRUE), use.names = FALSE)
    rows <- lapply(lines, function(line) {
        line <- trimws(sub(".*:", "", line))
        if (!nzchar(line))
            return(NULL)
        vals <- suppressWarnings(as.numeric(unlist(strsplit(gsub(",", " ", line), "\\s+"))))
        vals <- vals[!is.na(vals)]
        if (length(vals) != n_levels)
            return(NULL)
        vals
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) == 0)
        return(matrix(numeric(), nrow = 0, ncol = n_levels))
    do.call(rbind, rows)
}

.je_custom_contrast_tests <- function(dat, dep, factor, contrasts, fit) {
    levs <- levels(dat[[factor]])
    means <- tapply(dat[[dep]], dat[[factor]], mean)
    ns <- tapply(dat[[dep]], dat[[factor]], length)
    aov_tab <- as.data.frame(stats::anova(fit))
    resid_row <- which(rownames(aov_tab) == "Residuals")
    if (length(resid_row) == 0)
        return("<p>Custom contrasts could not be tested because residual error was unavailable.</p>")
    mse <- aov_tab$`Mean Sq`[resid_row]
    df <- aov_tab$Df[resid_row]
    rows <- lapply(seq_len(nrow(contrasts)), function(i) {
        w <- contrasts[i, ]
        est <- sum(w * means[levs])
        se <- sqrt(mse * sum((w ^ 2) / ns[levs]))
        t <- est / se
        p <- 2 * stats::pt(abs(t), df = df, lower.tail = FALSE)
        data.frame(
            Contrast = paste0("C", i),
            Weights = paste(w, collapse = ", "),
            Estimate = .je_fmt(est),
            SE = .je_fmt(se),
            t = .je_fmt(t),
            df = .je_fmt(df),
            p = .je_p(p),
            check.names = FALSE
        )
    })
    paste0("<h4>Custom contrast tests: ", .je_escape(factor), "</h4>", .je_table_html(do.call(rbind, rows)))
}

.je_order_restricted_status <- function(options) {
    if (!isTRUE(options$orderRestricted) && !isTRUE(options$modelComparison) && !isTRUE(options$informedHypothesisTests))
        return("<p>Bayesian/order-restricted hypothesis testing is disabled.</p>")
    paste0(
        "<p><strong>Captured hypothesis syntax:</strong></p><pre>",
        .je_escape(options$orderRestrictedSyntax),
        "</pre><p>Order-restricted and Bayesian model comparison still require the JASP informed-hypothesis engine. ",
        "The UI and syntax capture are implemented, but Bayes factors/weights are not yet calculated.</p>"
    )
}

.je_bootstrap_effects <- function(dat, dep, factors, covariates, options) {
    if (!isTRUE(options$bootstrapCi) && !isTRUE(options$postHocBootstrap))
        return("")
    rhs <- c(factors, covariates)
    if (length(rhs) == 0)
        return("<p>Bootstrap CIs require model terms.</p>")
    reps <- min(as.integer(options$bootstrapSamples), 2000L)
    reps <- max(reps, 100L)
    formula <- stats::as.formula(paste(dep, "~", paste(rhs, collapse = " * ")))
    obs <- tryCatch(.je_effect_sizes({
        tab <- as.data.frame(stats::anova(stats::lm(formula, data = dat)))
        tab$Term <- rownames(tab)
        rownames(tab) <- NULL
        tab
    }), error = function(e) data.frame())
    if (nrow(obs) == 0)
        return("<p>Bootstrap effect-size CIs could not be calculated.</p>")
    boot_vals <- replicate(reps, {
        idx <- sample.int(nrow(dat), replace = TRUE)
        bdat <- dat[idx, , drop = FALSE]
        vals <- tryCatch({
            tab <- as.data.frame(stats::anova(stats::lm(formula, data = bdat)))
            tab$Term <- rownames(tab)
            rownames(tab) <- NULL
            eff <- .je_effect_sizes(tab)
            as.numeric(eff$`partial eta squared`)
        }, error = function(e) rep(NA_real_, nrow(obs)))
        if (length(vals) != nrow(obs)) rep(NA_real_, nrow(obs)) else vals
    })
    if (is.null(dim(boot_vals)))
        boot_vals <- matrix(boot_vals, nrow = nrow(obs))
    ci <- t(apply(boot_vals, 1, stats::quantile, probs = c(.025, .975), na.rm = TRUE))
    out <- data.frame(
        Term = obs$Term,
        `partial eta squared` = obs$`partial eta squared`,
        `bootstrap lower` = .je_fmt(ci[, 1]),
        `bootstrap upper` = .je_fmt(ci[, 2]),
        Samples = reps,
        check.names = FALSE
    )
    paste0("<h4>Bootstrap effect-size confidence intervals</h4>", .je_table_html(out))
}

.je_plot_status <- function(dat, dep, factors, options) {
    enabled <- c(
        if (isTRUE(options$qqPlotResiduals)) "Q-Q residual plot",
        if (isTRUE(options$residualPlots)) "residual plots",
        if (isTRUE(options$descriptivePlots)) "descriptive plots",
        if (isTRUE(options$barPlots)) "bar plots",
        if (isTRUE(options$raincloudPlots)) "raincloud plots"
    )
    if (length(enabled) == 0)
        return("<p>No plot options are enabled.</p>")
    html <- paste0("<p>Requested plot outputs: ", .je_escape(paste(enabled, collapse = ", ")), ".</p>")
    if (length(factors) > 0 && (isTRUE(options$descriptivePlots) || isTRUE(options$barPlots) || isTRUE(options$raincloudPlots))) {
        desc <- .je_descriptives(dat, dep, factors[1])
        html <- paste0(html, "<p>Plot data summary for ", .je_escape(factors[1]), ":</p>", .je_table_html(desc))
    }
    if (isTRUE(options$raincloudPlots)) {
        html <- paste0(html, "<p>Raincloud plot rendering is represented by the summary above in this build. Pixel-rendered raincloud graphics still need a jamovi Image renderer.</p>")
    }
    html
}

.je_publication_html <- function(options) {
    requested <- isTRUE(options$publicationMode) || isTRUE(options$publicationTables) ||
        isTRUE(options$publicationFigures) || isTRUE(options$exportWord) || isTRUE(options$exportPdf)
    if (!requested)
        return("<p>Publication mode is disabled.</p>")
    paste0(
        "<p>Publication mode is partially implemented: result tables are produced as copyable HTML. ",
        "Direct Word/PDF export from the analysis runner is not implemented because jamovi normally handles document export at the application/report level.</p>"
    )
}

.je_saved_columns_html <- function(options) {
    requested <- c(
        if (isTRUE(options$saveResiduals) && isTRUE(options$saveRawResiduals)) "raw residuals",
        if (isTRUE(options$saveResiduals) && isTRUE(options$saveStudentizedResiduals)) "studentized residuals",
        if (isTRUE(options$saveResiduals) && isTRUE(options$saveStandardizedResiduals)) "standardized residuals",
        if (isTRUE(options$savePredictions)) "predicted values"
    )
    if (length(requested) == 0)
        return("<p>No dataset output columns requested.</p>")
    paste0(
        "<p>Generated dataset columns requested: ",
        .je_escape(paste(requested, collapse = ", ")),
        ".</p><p>Use the Output controls in the options panel to choose/create the destination columns.</p>"
    )
}

.je_pairwise_html <- function(mat) {
    if (is.null(mat) || length(mat) == 0)
        return("<p>No pairwise p-values available.</p>")
    dat <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
    names(dat) <- c("Group 1", "Group 2", "p")
    dat <- dat[!is.na(dat$p), , drop = FALSE]
    if (nrow(dat) == 0)
        return("<p>No pairwise p-values available.</p>")
    dat$p <- .je_p(dat$p)
    .je_table_html(dat)
}

.je_posthoc_adjust <- function(options) {
    if (isTRUE(options$postHocBonferroni)) return("bonferroni")
    if (isTRUE(options$postHocHolm)) return("holm")
    if (isTRUE(options$postHocSidak)) return("bonferroni")
    "holm"
}

.je_adjust_method <- function(x) {
    x <- as.character(x %||% "none")
    switch(x,
        bonferroni = "bonferroni",
        holm = "holm",
        sidak = "bonferroni",
        tukey = "bonferroni",
        none = "none",
        "none"
    )
}

`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || is.na(x))
        y
    else
        x
}

.je_conditional_posthoc <- function(dat, dep, factors, options) {
    out <- character()
    for (moderator in factors) {
        simple <- setdiff(factors, moderator)[1]
        for (level in levels(dat[[moderator]])) {
            sub <- dat[dat[[moderator]] == level, , drop = FALSE]
            if (length(unique(sub[[simple]])) < 2)
                next
            pw <- tryCatch(stats::pairwise.t.test(sub[[dep]], sub[[simple]], p.adjust.method = .je_posthoc_adjust(options)), error = function(e) e)
            if (!inherits(pw, "error")) {
                out <- c(out, paste0("<h4>Conditional comparisons: ", .je_escape(simple), " at ", .je_escape(moderator), " = ", .je_escape(level), "</h4>", .je_pairwise_html(pw$p.value)))
            }
        }
    }
    paste(out, collapse = "")
}

.je_contr_diff <- function(n) {
    if (n < 2)
        return(matrix(numeric(), nrow = n, ncol = 0))
    mat <- matrix(0, nrow = n, ncol = n - 1)
    for (i in seq_len(n - 1)) {
        mat[i, i] <- -1
        mat[i + 1, i] <- 1
    }
    mat
}

.je_repro_html <- function(formula, options) {
    paste0(
        "<p>Analysis formula:</p><pre>", .je_escape(deparse(formula)), "</pre>",
        "<p>R syntax export is scaffolded. A full export should include model options, contrast settings, post hoc settings, and plot settings.</p>"
    )
}

.je_not_yet <- function(message) {
    paste0("<p><strong>Not implemented yet:</strong> ", .je_escape(message), "</p>")
}

.je_set_pending_outputs <- function(self, reason) {
    html <- paste0("<p>", .je_escape(reason), "</p>")
    self$results$modelSummary$setContent(html)
    self$results$descriptivesSection$setContent(html)
    self$results$effectSizesSection$setContent(html)
    self$results$assumptionsSection$setContent(html)
    self$results$contrasts$setContent(html)
    self$results$orderRestrictions$setContent(html)
    self$results$postHocSection$setContent(html)
    self$results$marginalMeansSection$setContent(html)
    self$results$simpleEffectsSection$setContent(html)
    self$results$nonparametric$setContent(html)
    self$results$plots$setContent(html)
    self$results$savedColumns$setContent(html)
    self$results$apa$setContent(html)
    self$results$teaching$setContent(html)
    self$results$publication$setContent(html)
    self$results$reproducibility$setContent(html)
}
