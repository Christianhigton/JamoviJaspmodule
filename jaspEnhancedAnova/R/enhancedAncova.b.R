#' @importFrom jmvcore .
enhancedAncovaClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
    "enhancedAncovaClass",
    inherit = enhancedAncovaBase,
    private = list(
        .plotState = NULL,
        .run = function() {
            self$results$about$setContent(.je_about_html("https://github.com/jasp-stats/jaspAnova"))

            dep     <- .je_chr(self$options$dependent)
            factors <- .je_chr_vec(self$options$fixedFactors)
            random  <- .je_chr_vec(self$options$randomFactors)
            covs    <- .je_chr_vec(self$options$covariates)
            weights <- .je_chr(self$options$wlsWeights)

            if (is.null(dep) || dep == "" || (length(factors) == 0 && length(covs) == 0 && length(random) == 0)) {
                self$results$status$setContent(
                    "<p>Select a dependent variable and at least one fixed factor, random factor, or covariate to run ANCOVA.</p>"
                )
                .je_set_pending_outputs(self, "Waiting for required variables.")
                return()
            }

            vars <- unique(c(dep, factors, random, covs, weights))
            vars <- vars[nzchar(vars)]
            dat  <- self$data[, vars, drop = FALSE]
            dat  <- dat[stats::complete.cases(dat), , drop = FALSE]
            row_nums <- as.integer(rownames(dat))
            if (anyNA(row_nums)) row_nums <- seq_len(nrow(dat))

            if (nrow(dat) < 3) {
                self$results$status$setContent("<p>Not enough complete cases to run ANCOVA (need at least 3).</p>")
                .je_set_pending_outputs(self, "Insufficient complete cases.")
                return()
            }

            all_factors <- unique(c(factors, random))
            for (f in all_factors) dat[[f]] <- as.factor(dat[[f]])

            rhs_terms <- .je_model_terms_from_options(self$options$modelTerms, all_factors, covs)
            formula <- stats::as.formula(paste(dep, "~", paste(rhs_terms, collapse = " + ")))

            fitted <- .je_fit_and_anova(formula, dat, all_factors, weights)
            if (inherits(fitted$fit, "error")) {
                self$results$status$setContent(paste0(
                    "<p>The ANCOVA model could not be fitted: ", .je_escape(fitted$fit$message), "</p>"
                ))
                .je_set_pending_outputs(self, "Model fit failed.")
                return()
            }

            fit      <- fitted$fit
            aov_tab  <- fitted$aov_tab
            ss_type  <- fitted$ss_type
            n        <- nrow(dat)

            effect_tab  <- .je_effect_sizes(aov_tab, n, self$options)
            desc_tab    <- .je_descriptives(dat, dep, all_factors, self$options)
            assumptions <- .je_ancova_assumptions(dat, dep, all_factors, covs, fit, self$options)
            posthoc     <- .je_posthoc(dat, dep, all_factors, fit, self$options)
            kruskal     <- if (length(.je_chr_vec(self$options$kruskalWallisFactors)) > 0) .je_kruskal(dat, dep, all_factors, self$options) else "<p>Assign factors to the Kruskal-Wallis Test box to enable nonparametric analysis.</p>"
            marginal    <- .je_marginal_means(dat, dep, all_factors, covs, fit, self$options)
            simple      <- .je_simple_effects(dat, dep, all_factors, fit, self$options)
            contrast    <- .je_contrasts(dat, dep, all_factors, fit, self$options)
            bootstrap   <- .je_bootstrap_effects(dat, dep, all_factors, covs, self$options)
            apa         <- .je_apa(aov_tab, effect_tab, dep, all_factors, ss_type)

            self$results$status$setContent(paste0(
                "<p><strong>Engine:</strong> Classical ANCOVA using Type ", ss_type,
                " Sum of Squares with covariates in the linear model.</p>",
                if (length(random) > 0)
                    "<p><em>Random factors are included as model factors in this jamovi-native port; mixed-effects random variance components are not estimated in this pass.</em></p>"
                else "",
                "<p><strong>Dependencies for full JASP parity:</strong> ",
                "<code>car</code> (Type III SS) - ", if (requireNamespace("car", quietly = TRUE)) "installed" else "<strong>missing</strong>",
                "; <code>emmeans</code> (post hoc/EMMs) - ", if (requireNamespace("emmeans", quietly = TRUE)) "installed" else "<strong>missing</strong>",
                ".</p>"
            ))

            self$results$modelSummary$setContent(
                .je_model_html(aov_tab, formula, n, weights, ss_type, NULL, self$options)
            )
            self$results$descriptivesSection$setContent(.je_table_html(desc_tab))
            self$results$effectSizesSection$setContent(.je_table_html(effect_tab))
            self$results$assumptionsSection$setContent(assumptions)
            self$results$contrasts$setContent(contrast)
            self$results$orderRestrictions$setContent(
                .je_order_restricted_full(dat, dep, all_factors, covs, fit, aov_tab, self$options)
            )
            self$results$postHocSection$setContent(posthoc)
            self$results$marginalMeansSection$setContent(marginal)
            self$results$simpleEffectsSection$setContent(simple)
            self$results$nonparametric$setContent(kruskal)
            self$results$plots$setContent(.je_plot_status(dat, dep, all_factors, self$options))
            self$results$savedColumns$setContent(.je_saved_columns_html(self$options))
            self$results$apa$setContent(paste0(apa, bootstrap))
            self$results$teaching$setContent(.je_teaching_html(effect_tab, assumptions))
            self$results$publication$setContent(.je_publication_html(self$options))
            self$results$reproducibility$setContent(.je_repro_html(formula, self$options))

            private$.plotState <- list(data = dat, dep = dep, factors = all_factors, fit = fit)
            if (isTRUE(self$options$qqPlot))
                self$results$qqPlot$setState(private$.plotState)
            if (!is.null(.je_chr(self$options$rainCloudHorizontalAxis)))
                self$results$raincloudPlot$setState(private$.plotState)

            private$.populateOutputs(fit, row_nums)
        },

        .populateOutputs = function(fit, row_nums) {
            if (isTRUE(self$options$residualsSavedToData) &&
                self$options$residualsSavedToDataType %in% c("raw", NULL) &&
                !is.null(self$options$residsOV) && self$results$residsOV$isNotFilled()) {
                self$results$residsOV$setRowNums(row_nums)
                self$results$residsOV$setValues(stats::residuals(fit))
            }
            if (isTRUE(self$options$residualsSavedToData) &&
                identical(self$options$residualsSavedToDataType, "student") &&
                !is.null(self$options$residsOV) && self$results$residsOV$isNotFilled()) {
                self$results$residsOV$setRowNums(row_nums)
                self$results$residsOV$setValues(stats::rstudent(fit))
            }
            if (isTRUE(self$options$residualsSavedToData) &&
                identical(self$options$residualsSavedToDataType, "standard") &&
                !is.null(self$options$residsOV) && self$results$residsOV$isNotFilled()) {
                self$results$residsOV$setRowNums(row_nums)
                self$results$residsOV$setValues(stats::rstandard(fit))
            }
            if (isTRUE(self$options$predictionsSavedToData) && self$options$predictOV &&
                self$results$predictOV$isNotFilled()) {
                self$results$predictOV$setRowNums(row_nums)
                self$results$predictOV$setValues(stats::fitted(fit))
            }
        },

        .qqPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state) || !requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)
            res <- stats::rstandard(image$state$fit)
            df  <- as.data.frame(stats::qqnorm(res, plot.it = FALSE))
            ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
                ggplot2::geom_abline(slope = 1, intercept = 0, colour = theme$color[1]) +
                ggplot2::geom_point(colour = theme$color[1]) +
                ggplot2::labs(x = "Theoretical Quantiles", y = "Standardized Residuals") +
                ggtheme
        },

        .residualPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state) || !requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)
            df <- data.frame(Fitted = stats::fitted(image$state$fit),
                             Residuals = stats::rstandard(image$state$fit))
            ggplot2::ggplot(df, ggplot2::aes(x = Fitted, y = Residuals)) +
                ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = theme$color[1]) +
                ggplot2::geom_point(colour = theme$color[1], alpha = 0.75) +
                ggplot2::labs(x = "Fitted Values", y = "Standardized Residuals") +
                ggtheme
        },

        .raincloudPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state) || !requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)
            dat     <- image$state$data
            dep     <- image$state$dep
            factors <- image$state$factors
            if (length(factors) == 0) return(FALSE)
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

.je_model_terms_from_options <- function(model_terms, factors, covariates) {
    default <- unique(c(factors, covariates))
    raw <- .je_chr(model_terms)
    if (is.null(raw) || !nzchar(trimws(raw)))
        return(default)

    terms <- unlist(strsplit(raw, "[,;\n\r]+"))
    terms <- trimws(terms)
    terms <- terms[nzchar(terms)]
    if (length(terms) == 0) default else unique(terms)
}

.je_ancova_assumptions <- function(dat, dep, factors, covariates, fit, options) {
    out <- .je_assumptions(dat, dep, factors, fit, options)

    if (length(factors) > 0 && length(covariates) > 0) {
        chunks <- character()
        for (cov in covariates) {
            rhs <- paste(factors, collapse = " + ")
            form <- stats::as.formula(paste(cov, "~", rhs))
            tab <- tryCatch(as.data.frame(stats::anova(stats::lm(form, data = dat))), error = function(e) NULL)
            if (is.null(tab)) next
            tab$Term <- rownames(tab)
            rownames(tab) <- NULL
            int <- tab[tab$Term != "Residuals", , drop = FALSE]
            if (nrow(int) == 0) next
            res <- data.frame(
                Covariate = cov,
                Term = int$Term,
                df = int$Df,
                F = .je_fmt(int$`F value`),
                p = .je_p(int$`Pr(>F)`),
                check.names = FALSE
            )
            chunks <- c(chunks, paste0("<h4>Factor-covariate independence: ", .je_escape(cov), "</h4>", .je_table_html(res)))
        }
        if (length(chunks) > 0)
            out <- paste0(out, paste(chunks, collapse = ""))
    }

    out
}
