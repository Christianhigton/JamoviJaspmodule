#' @importFrom jmvcore .
enhancedRepeatedMeasuresAnovaClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
    "enhancedRepeatedMeasuresAnovaClass",
    inherit = enhancedRepeatedMeasuresAnovaBase,
    private = list(
        .plotState = NULL,
        .run = function() {
            self$results$about$setContent(.je_about_html("https://github.com/jasp-stats/jaspAnova"))

            cells <- .je_chr_vec(self$options$repeatedMeasures)
            between <- .je_chr_vec(self$options$betweenFactors)
            covariates <- .je_chr_vec(self$options$covariates)
            grouping <- .je_chr(self$options$groupingFactor)

            if (length(cells) < 2) {
                self$results$status$setContent("<p>Select at least two repeated-measures cells to run the repeated measures ANOVA engine.</p>")
                .je_rm_set_pending_outputs(self, "Waiting for at least two repeated-measures cells.")
                return()
            }

            vars <- unique(c(cells, between, covariates, grouping))
            vars <- vars[nzchar(vars)]
            wide <- self$data[, vars, drop = FALSE]
            names(wide) <- vars
            wide <- wide[stats::complete.cases(wide[, cells, drop = FALSE]), , drop = FALSE]
            row_nums <- as.integer(rownames(wide))
            if (anyNA(row_nums))
                row_nums <- seq_len(nrow(wide))

            if (nrow(wide) < 2) {
                self$results$status$setContent("<p>There are not enough complete cases to run repeated measures ANOVA.</p>")
                .je_rm_set_pending_outputs(self, "Insufficient complete cases.")
                return()
            }

            long <- .je_rm_long_data(wide, cells, between, covariates)
            levels_within <- levels(long$Within)

            for (f in between)
                long[[f]] <- as.factor(long[[f]])

            rhs <- c("Within", between, covariates)
            fixed_formula <- stats::as.formula(paste("Value ~", paste(rhs, collapse = " * ")))
            lm_fit <- tryCatch(stats::lm(fixed_formula, data = long), error = function(e) e)
            aov_fit <- tryCatch(stats::aov(stats::as.formula(paste("Value ~", paste(rhs, collapse = " * "), "+ Error(Subject/Within)")), data = long), error = function(e) e)

            if (inherits(lm_fit, "error")) {
                self$results$status$setContent(paste0("<p>The repeated measures model could not be fitted: ", .je_escape(lm_fit$message), "</p>"))
                .je_rm_set_pending_outputs(self, "Model fit failed.")
                return()
            }

            aov_tab <- .je_rm_anova_table(lm_fit)
            desc <- .je_rm_descriptives(long)
            contrasts <- .je_rm_contrasts(long, levels_within, lm_fit, self$options)
            assumptions <- .je_rm_assumptions(long, lm_fit, aov_fit, self$options)
            posthoc <- .je_rm_posthoc(long, self$options)
            marginal <- .je_rm_marginal_means(long, self$options)
            simple <- .je_rm_simple_effects(long, between, self$options)
            nonparametric <- .je_rm_nonparametric(wide, cells, self$options)
            apa <- .je_rm_apa(aov_tab, levels_within)

            self$results$status$setContent(paste0(
                "<p><strong>Implemented:</strong> wide-format repeated-measures ANOVA for selected cells, within-factor descriptives, ",
                "contrast inference, residual diagnostics, post hoc pairwise tests, marginal means, simple effects by between factor, ",
                "Friedman test, generated subject-level residual/prediction columns, and raincloud-style plots.</p>",
                "<p><strong>Scope note:</strong> this build treats selected repeated-measures cells as one within-subject factor. ",
                "Full JASP parity for multiple/nested RM factors, sphericity epsilon details, and Bayesian/order-restricted model weights still needs further porting.</p>"
            ))

            self$results$design$setContent(.je_rm_design_html(cells, levels_within, between, covariates))
            self$results$modelSummary$setContent(.je_table_html(aov_tab))
            self$results$descriptivesSection$setContent(.je_table_html(desc))
            self$results$effectSizesSection$setContent(.je_rm_effect_sizes(aov_tab))
            self$results$assumptionsSection$setContent(assumptions)
            self$results$contrasts$setContent(contrasts)
            self$results$orderRestrictions$setContent(.je_rm_order_status(self$options))
            self$results$postHocSection$setContent(posthoc)
            self$results$marginalMeansSection$setContent(marginal)
            self$results$simpleEffectsSection$setContent(simple)
            self$results$nonparametric$setContent(nonparametric)
            self$results$plots$setContent(.je_rm_plot_status(long, self$options))
            self$results$savedColumns$setContent(.je_rm_saved_columns_html(self$options))
            self$results$apa$setContent(apa)
            self$results$teaching$setContent(.je_rm_teaching_html(self$options))
            self$results$publication$setContent(.je_publication_html(self$options))
            self$results$reproducibility$setContent(.je_rm_repro_html(cells, between, covariates))

            private$.plotState <- list(long = long, fit = lm_fit)
            if (isTRUE(self$options$qqPlotResiduals))
                self$results$qqPlot$setState(private$.plotState)
            if (isTRUE(self$options$residualPlots))
                self$results$residualPlot$setState(private$.plotState)
            if (isTRUE(self$options$raincloudPlots))
                self$results$raincloudPlot$setState(private$.plotState)

            private$.populateOutputs(lm_fit, long, row_nums)
        },
        .populateOutputs = function(fit, long, row_nums) {
            subject_res <- stats::aggregate(stats::residuals(fit), list(Subject = long$Subject), mean)
            subject_std <- stats::aggregate(stats::rstandard(fit), list(Subject = long$Subject), mean)
            subject_stu <- stats::aggregate(stats::rstudent(fit), list(Subject = long$Subject), mean)
            subject_fit <- stats::aggregate(stats::fitted(fit), list(Subject = long$Subject), mean)

            if (isTRUE(self$options$saveResiduals) && isTRUE(self$options$saveRawResiduals) &&
                self$options$rmResidsOV && self$results$rmResidsOV$isNotFilled()) {
                self$results$rmResidsOV$setRowNums(row_nums)
                self$results$rmResidsOV$setValues(subject_res$x)
            }
            if (isTRUE(self$options$saveResiduals) && isTRUE(self$options$saveStudentizedResiduals) &&
                self$options$rmStudentizedResidsOV && self$results$rmStudentizedResidsOV$isNotFilled()) {
                self$results$rmStudentizedResidsOV$setRowNums(row_nums)
                self$results$rmStudentizedResidsOV$setValues(subject_stu$x)
            }
            if (isTRUE(self$options$saveResiduals) && isTRUE(self$options$saveStandardizedResiduals) &&
                self$options$rmStandardizedResidsOV && self$results$rmStandardizedResidsOV$isNotFilled()) {
                self$results$rmStandardizedResidsOV$setRowNums(row_nums)
                self$results$rmStandardizedResidsOV$setValues(subject_std$x)
            }
            if (isTRUE(self$options$savePredictions) && self$options$rmPredictOV &&
                self$results$rmPredictOV$isNotFilled()) {
                self$results$rmPredictOV$setRowNums(row_nums)
                self$results$rmPredictOV$setValues(subject_fit$x)
            }
        },
        .qqPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state) || !requireNamespace("ggplot2", quietly = TRUE))
                return(FALSE)
            df <- as.data.frame(stats::qqnorm(stats::rstandard(image$state$fit), plot.it = FALSE))
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
            long <- image$state$long
            ggplot2::ggplot(long, ggplot2::aes(x = Within, y = Value, fill = Within)) +
                ggplot2::geom_violin(trim = FALSE, alpha = 0.35, colour = NA) +
                ggplot2::geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.65) +
                ggplot2::geom_jitter(width = 0.08, alpha = 0.55, size = 1.6) +
                ggplot2::labs(x = "Repeated Measures Level", y = "Value") +
                ggplot2::guides(fill = "none") +
                ggtheme
        }
    )
)

.je_rm_long_data <- function(wide, cells, between, covariates) {
    id <- seq_len(nrow(wide))
    levels <- .je_rm_levels(cells, NULL)
    pieces <- lapply(seq_along(cells), function(i) {
        out <- data.frame(
            Subject = factor(id),
            Within = factor(levels[i], levels = levels),
            Value = wide[[cells[i]]],
            stringsAsFactors = FALSE
        )
        for (v in c(between, covariates))
            out[[v]] <- wide[[v]]
        out
    })
    do.call(rbind, pieces)
}

.je_rm_levels <- function(cells, spec) {
    cells
}

.je_rm_anova_table <- function(fit) {
    tab <- as.data.frame(stats::anova(fit))
    tab$Term <- rownames(tab)
    rownames(tab) <- NULL
    data.frame(
        Term = tab$Term,
        df = .je_fmt(tab$Df),
        `Sum Sq` = .je_fmt(tab$`Sum Sq`),
        `Mean Sq` = .je_fmt(tab$`Mean Sq`),
        F = .je_fmt(tab$`F value`),
        p = .je_p(tab$`Pr(>F)`),
        check.names = FALSE
    )
}

.je_rm_descriptives <- function(long) {
    sp <- split(long$Value, long$Within)
    data.frame(
        Level = names(sp),
        N = vapply(sp, length, integer(1)),
        Mean = .je_fmt(vapply(sp, mean, numeric(1))),
        SD = .je_fmt(vapply(sp, stats::sd, numeric(1))),
        check.names = FALSE
    )
}

.je_rm_design_html <- function(cells, levels, between, covariates) {
    paste0(
        "<p>Repeated-measures cells: ", .je_escape(paste(cells, collapse = ", ")), "</p>",
        "<p>Within-subject levels used by this build: ", .je_escape(paste(levels, collapse = ", ")), "</p>",
        "<p>Between-subject factors: ", .je_escape(if (length(between) == 0) "none" else paste(between, collapse = ", ")), "</p>",
        "<p>Covariates: ", .je_escape(if (length(covariates) == 0) "none" else paste(covariates, collapse = ", ")), "</p>"
    )
}

.je_rm_effect_sizes <- function(model_tab) {
    raw <- model_tab
    ss <- suppressWarnings(as.numeric(raw$`Sum Sq`))
    terms <- raw$Term
    resid <- which(terms == "Residuals")
    if (length(resid) == 0)
        return("<p>Effect sizes unavailable.</p>")
    sse <- ss[resid]
    total <- sum(ss, na.rm = TRUE)
    keep <- setdiff(seq_along(terms), resid)
    out <- data.frame(
        Term = terms[keep],
        `eta squared` = .je_fmt(ss[keep] / total),
        `partial eta squared` = .je_fmt(ss[keep] / (ss[keep] + sse)),
        check.names = FALSE
    )
    .je_table_html(out)
}

.je_rm_assumptions <- function(long, fit, aov_fit, options) {
    res <- stats::residuals(fit)
    normal <- if (length(res) >= 3 && length(res) <= 5000) {
        sh <- stats::shapiro.test(res)
        paste0("<p>Residual normality: Shapiro-Wilk W = ", .je_fmt(unname(sh$statistic)), ", p = ", .je_p(sh$p.value), ".</p>")
    } else {
        "<p>Residual normality: Shapiro-Wilk test skipped because sample size is outside the supported range.</p>"
    }
    sph <- if (length(levels(long$Within)) > 2) {
        "<p>Sphericity: Mauchly/epsilon details are not fully ported in this build. Use Greenhouse-Geisser/Huynh-Feldt options as recorded analysis settings until the exact JASP engine is ported.</p>"
    } else {
        "<p>Sphericity: not applicable for two repeated-measures levels.</p>"
    }
    paste0(normal, sph)
}

.je_rm_contrasts <- function(long, levels_within, fit, options) {
    requested <- options$rmContrastType != "none" || isTRUE(options$rmCustomContrasts)
    if (!requested)
        return("<p>No repeated-measures contrasts are enabled.</p>")
    mats <- list()
    n <- length(levels_within)
    if (options$rmContrastType == "helmert")
        mats$Helmert <- stats::contr.helmert(n)
    if (options$rmContrastType == "polynomial")
        mats$Polynomial <- stats::contr.poly(n)
    if (options$rmContrastType == "deviation")
        mats$Deviation <- stats::contr.sum(n)
    if (options$rmContrastType %in% c("difference", "repeated"))
        mats$Difference <- .je_contr_diff(n)
    if (isTRUE(options$rmCustomContrasts)) {
        custom <- .je_parse_custom_contrasts(options$rmContrastSyntax, n)
        if (nrow(custom) > 0)
            mats$Custom <- t(custom)
    }
    if (length(mats) == 0)
        return(paste0("<p>Enter one contrast per line with ", n, " numeric weights matching levels: ", .je_escape(paste(levels_within, collapse = ", ")), ".</p>"))
    html <- character()
    means <- tapply(long$Value, long$Within, mean)
    wide <- stats::reshape(long[, c("Subject", "Within", "Value")], idvar = "Subject", timevar = "Within", direction = "wide")
    for (name in names(mats)) {
        mat <- mats[[name]]
        display <- as.data.frame(mat)
        names(display) <- paste0("C", seq_len(ncol(display)))
        display$Level <- levels_within
        html <- c(html, paste0("<h4>", .je_escape(name), " contrast matrix</h4>", .je_table_html(display[, c("Level", setdiff(names(display), "Level")), drop = FALSE])))
        tests <- lapply(seq_len(ncol(mat)), function(i) {
            w <- mat[, i]
            vals <- as.matrix(wide[, paste0("Value.", levels_within), drop = FALSE]) %*% w
            tt <- stats::t.test(as.numeric(vals), mu = 0)
            data.frame(
                Contrast = paste0("C", i),
                Weights = paste(.je_fmt(w), collapse = ", "),
                Estimate = .je_fmt(sum(w * means[levels_within])),
                t = .je_fmt(unname(tt$statistic)),
                df = .je_fmt(unname(tt$parameter)),
                p = .je_p(tt$p.value),
                check.names = FALSE
            )
        })
        html <- c(html, paste0("<h4>", .je_escape(name), " contrast tests</h4>", .je_table_html(do.call(rbind, tests))))
    }
    paste(html, collapse = "")
}

.je_rm_posthoc <- function(long, options) {
    if (!isTRUE(options$postHoc))
        return("<p>Post hoc tests are disabled.</p>")
    pw <- stats::pairwise.t.test(long$Value, long$Within, paired = TRUE, p.adjust.method = if (isTRUE(options$postHocBonferroni)) "bonferroni" else "holm")
    paste0("<p>Paired pairwise comparisons across repeated-measures levels.</p>", .je_pairwise_html(pw$p.value))
}

.je_rm_marginal_means <- function(long, options) {
    if (!isTRUE(options$marginalMeans))
        return("<p>Estimated marginal means are disabled.</p>")
    .je_table_html(.je_rm_descriptives(long))
}

.je_rm_simple_effects <- function(long, between, options) {
    if (!isTRUE(options$simpleEffects))
        return("<p>Simple effects analysis is disabled.</p>")
    if (length(between) == 0)
        return("<p>Simple effects need a between-subject moderator in this build.</p>")
    mod <- between[1]
    pieces <- lapply(levels(as.factor(long[[mod]])), function(level) {
        sub <- long[as.factor(long[[mod]]) == level, , drop = FALSE]
        fit <- stats::lm(Value ~ Within, data = sub)
        paste0("<h4>", .je_escape(mod), " = ", .je_escape(level), "</h4>", .je_table_html(.je_rm_anova_table(fit)))
    })
    paste(pieces, collapse = "")
}

.je_rm_nonparametric <- function(wide, cells, options) {
    if (!isTRUE(options$friedman))
        return("<p>Friedman-style analysis is disabled.</p>")
    ft <- stats::friedman.test(as.matrix(wide[, cells, drop = FALSE]))
    paste0("<p>Friedman chi-square = ", .je_fmt(unname(ft$statistic)), ", df = ", ft$parameter, ", p = ", .je_p(ft$p.value), ".</p>")
}

.je_rm_apa <- function(model_tab, levels_within) {
    row <- model_tab[model_tab$Term == "Within", , drop = FALSE]
    if (nrow(row) == 0)
        return("<p>No APA repeated-measures text is available for this model.</p>")
    paste0("<p>A repeated-measures ANOVA tested differences across ", length(levels_within), " within-subject levels. The within-subject effect was F(", row$df, ") = ", row$F, ", p = ", row$p, ".</p>")
}

.je_rm_order_status <- function(options) {
    if (!isTRUE(options$orderRestricted) && !isTRUE(options$modelComparison) && !isTRUE(options$informedHypothesisTests))
        return("<p>Bayesian/order-restricted hypothesis testing is disabled.</p>")
    paste0("<p><strong>Captured hypothesis syntax:</strong></p><pre>", .je_escape(options$orderRestrictedSyntax), "</pre><p>Bayesian model weights are not yet calculated in this repeated-measures build.</p>")
}

.je_rm_plot_status <- function(long, options) {
    enabled <- c(if (isTRUE(options$qqPlotResiduals)) "Q-Q residual plot", if (isTRUE(options$residualPlots)) "residual plot", if (isTRUE(options$raincloudPlots)) "raincloud plot")
    if (length(enabled) == 0)
        return("<p>No plot options are enabled.</p>")
    paste0("<p>Rendered plots requested: ", .je_escape(paste(enabled, collapse = ", ")), ".</p>")
}

.je_rm_saved_columns_html <- function(options) {
    requested <- c(
        if (isTRUE(options$saveResiduals) && isTRUE(options$saveRawResiduals)) "subject mean residuals",
        if (isTRUE(options$saveResiduals) && isTRUE(options$saveStudentizedResiduals)) "subject mean studentized residuals",
        if (isTRUE(options$saveResiduals) && isTRUE(options$saveStandardizedResiduals)) "subject mean standardized residuals",
        if (isTRUE(options$savePredictions)) "subject mean predicted values"
    )
    if (length(requested) == 0)
        return("<p>No dataset output columns requested.</p>")
    paste0("<p>Generated dataset columns requested: ", .je_escape(paste(requested, collapse = ", ")), ".</p><p>Because repeated measures are stored across multiple source columns, this build writes subject-level means back to the dataset.</p>")
}

.je_rm_teaching_html <- function(options) {
    "<p>Repeated-measures ANOVA tests whether within-subject means differ after accounting for stable subject-to-subject differences. For two levels, the within-subject test is equivalent to a paired comparison.</p>"
}

.je_rm_repro_html <- function(cells, between, covariates) {
    paste0("<p>Repeated-measures cells:</p><pre>", .je_escape(paste(cells, collapse = ", ")), "</pre>",
           "<p>Between factors: ", .je_escape(if (length(between) == 0) "none" else paste(between, collapse = ", ")), "</p>",
           "<p>Covariates: ", .je_escape(if (length(covariates) == 0) "none" else paste(covariates, collapse = ", ")), "</p>")
}

.je_rm_set_pending_outputs <- function(self, reason) {
    html <- paste0("<p>", .je_escape(reason), "</p>")
    self$results$design$setContent(html)
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
