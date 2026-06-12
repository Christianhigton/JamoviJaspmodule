#' @importFrom jmvcore .
enhancedAnovaClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
    "enhancedAnovaClass",
    inherit = enhancedAnovaBase,
    private = list(
        .plotState = NULL,
        .run = function() {
            self$results$about$setContent(.je_about_html("https://github.com/jasp-stats/jaspAnova"))

            dep     <- .je_chr(self$options$dependent)
            factors <- .je_chr_vec(self$options$fixedFactors)
            covs    <- .je_chr_vec(self$options$covariates)
            weights <- .je_chr(self$options$wlsWeights)

            if (is.null(dep) || dep == "" || (length(factors) == 0 && length(covs) == 0)) {
                self$results$status$setContent(
                    "<p>Select a dependent variable and at least one fixed factor or covariate to run the analysis.</p>"
                )
                .je_set_pending_outputs(self, "Waiting for required variables.")
                return()
            }

            vars <- unique(c(dep, factors, covs, weights))
            vars <- vars[nzchar(vars)]
            dat  <- self$data[, vars, drop = FALSE]
            dat  <- dat[stats::complete.cases(dat), , drop = FALSE]
            row_nums <- as.integer(rownames(dat))
            if (anyNA(row_nums)) row_nums <- seq_len(nrow(dat))

            if (nrow(dat) < 3) {
                self$results$status$setContent("<p>Not enough complete cases to run ANOVA (need at least 3).</p>")
                .je_set_pending_outputs(self, "Insufficient complete cases.")
                return()
            }

            for (f in factors) dat[[f]] <- as.factor(dat[[f]])

            rhs     <- c(factors, covs)
            formula <- stats::as.formula(paste(dep, "~", paste(rhs, collapse = " * ")))

            fitted <- .je_fit_and_anova(formula, dat, factors, weights)
            if (inherits(fitted$fit, "error")) {
                self$results$status$setContent(paste0(
                    "<p>The model could not be fitted: ", .je_escape(fitted$fit$message), "</p>"
                ))
                .je_set_pending_outputs(self, "Model fit failed.")
                return()
            }

            fit      <- fitted$fit
            aov_tab  <- fitted$aov_tab
            ss_type  <- fitted$ss_type
            n        <- nrow(dat)

            # Welch / Brown-Forsythe correction (one-way only)
            correction_result <- .je_welch_correction(dat, dep, factors, self$options)

            effect_tab  <- .je_effect_sizes(aov_tab, n, self$options)
            desc_tab    <- .je_descriptives(dat, dep, factors, self$options)
            assumptions <- .je_assumptions(dat, dep, factors, fit, self$options)
            posthoc     <- .je_posthoc(dat, dep, factors, fit, self$options)
            kruskal     <- if (length(.je_chr_vec(self$options$kruskalWallisFactors)) > 0) .je_kruskal(dat, dep, factors, self$options) else "<p>Assign factors to the Kruskal-Wallis Test box to enable nonparametric analysis.</p>"
            marginal    <- .je_marginal_means(dat, dep, factors, covs, fit, self$options)
            simple      <- .je_simple_effects(dat, dep, factors, fit, self$options)
            contrast    <- .je_contrasts(dat, dep, factors, fit, self$options)
            bootstrap   <- .je_bootstrap_effects(dat, dep, factors, covs, self$options)
            apa         <- .je_apa(aov_tab, effect_tab, dep, factors, ss_type)

            self$results$status$setContent(paste0(
                "<p><strong>Engine:</strong> Classical fixed-effects ANOVA using Type ", ss_type,
                " Sum of Squares",
                if (!is.null(correction_result)) paste0(" with ", correction_result$method, " correction") else "",
                ".</p>",
                "<p><strong>Dependencies for full JASP parity:</strong> ",
                "<code>car</code> (Type III SS) — ", if (requireNamespace("car", quietly = TRUE)) "installed" else "<strong>missing</strong>",
                "; <code>emmeans</code> (post hoc/EMMs) — ", if (requireNamespace("emmeans", quietly = TRUE)) "installed" else "<strong>missing</strong>",
                ".</p>"
            ))

            self$results$modelSummary$setContent(
                .je_model_html(aov_tab, formula, n, weights, ss_type, correction_result, self$options)
            )
            self$results$descriptivesSection$setContent(.je_table_html(desc_tab))
            self$results$effectSizesSection$setContent(.je_table_html(effect_tab))
            self$results$assumptionsSection$setContent(assumptions)
            self$results$contrasts$setContent(contrast)
            self$results$orderRestrictions$setContent(
                .je_order_restricted_full(dat, dep, factors, covs, fit, aov_tab, self$options)
            )
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

# ── Fitting ───────────────────────────────────────────────────────────────────

.je_fit_and_anova <- function(formula, dat, factors, weights) {
    # Fit with sum-to-zero contrasts so car::Anova(type="III") gives correct SS
    contrast_list <- if (length(factors) > 0)
        setNames(lapply(factors, function(f) stats::contr.sum(nlevels(dat[[f]]))), factors)
    else
        list()

    fit <- tryCatch({
        args <- list(formula = formula, data = dat)
        if (length(contrast_list) > 0) args$contrasts <- contrast_list
        if (!is.null(weights) && nzchar(weights) && weights %in% names(dat))
            args$weights <- dat[[weights]]
        do.call(stats::lm, args)
    }, error = function(e) e)

    if (inherits(fit, "error")) return(list(fit = fit, aov_tab = NULL, ss_type = NULL))

    # Attempt Type III via car
    aov_tab <- NULL
    ss_type <- "I"

    if (requireNamespace("car", quietly = TRUE)) {
        ct <- tryCatch(car::Anova(fit, type = "III"), error = function(e) NULL)
        if (!is.null(ct)) {
            tab <- as.data.frame(ct)
            tab$Term <- rownames(tab)
            rownames(tab) <- NULL
            tab <- tab[tab$Term != "(Intercept)", , drop = FALSE]
            tab$`Mean Sq` <- tab$`Sum Sq` / tab$Df
            # Normalise column names to match rest of code
            col_map <- c("Sum Sq" = "Sum Sq", "Df" = "Df", "F value" = "F value",
                         "Pr(>F)" = "Pr(>F)", "Mean Sq" = "Mean Sq", "Term" = "Term")
            aov_tab <- tab[, c("Term", "Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)"), drop = FALSE]
            ss_type <- "III"
        }
    }

    if (is.null(aov_tab)) {
        tab <- as.data.frame(stats::anova(fit))
        tab$Term <- rownames(tab)
        rownames(tab) <- NULL
        aov_tab <- tab[, c("Term", "Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)"), drop = FALSE]
        ss_type <- "I"
    }

    list(fit = fit, aov_tab = aov_tab, ss_type = ss_type)
}

# ── Welch / Brown-Forsythe omnibus correction ─────────────────────────────────

.je_welch_correction <- function(dat, dep, factors, options) {
    welch  <- isTRUE(options$homogeneityCorrectionWelch)
    brown  <- isTRUE(options$homogeneityCorrectionBrown)
    if ((!welch && !brown) || length(factors) != 1) return(NULL)

    f1   <- factors[1]
    form <- stats::as.formula(paste(dep, "~", f1))

    if (welch) {
        res <- tryCatch(stats::oneway.test(form, data = dat, var.equal = FALSE), error = function(e) NULL)
        if (is.null(res)) return(NULL)
        return(list(
            method = "Welch",
            F      = unname(res$statistic),
            df1    = unname(res$parameter[1]),
            df2    = unname(res$parameter[2]),
            p      = res$p.value
        ))
    }

    if (brown) {
        # BF test: ANOVA on absolute deviations from group medians
        group   <- dat[[f1]]
        medians <- tapply(dat[[dep]], group, stats::median, na.rm = TRUE)
        devs    <- abs(dat[[dep]] - medians[as.character(group)])
        bf_dat  <- data.frame(devs = devs, group = group)
        res <- tryCatch(
            stats::oneway.test(devs ~ group, data = bf_dat, var.equal = FALSE),
            error = function(e) NULL
        )
        if (is.null(res)) return(NULL)
        return(list(
            method = "Brown-Forsythe",
            F      = unname(res$statistic),
            df1    = unname(res$parameter[1]),
            df2    = unname(res$parameter[2]),
            p      = res$p.value
        ))
    }
    NULL
}

# ── Effect sizes ──────────────────────────────────────────────────────────────

.je_effect_sizes <- function(aov_tab, n, options) {
    resid_row <- which(aov_tab$Term == "Residuals")
    if (length(resid_row) == 0) return(data.frame())

    sse      <- aov_tab$`Sum Sq`[resid_row]
    df_err   <- aov_tab$Df[resid_row]
    mse      <- sse / df_err
    ss_total <- sum(aov_tab$`Sum Sq`[aov_tab$Term != "(Intercept)"], na.rm = TRUE)
    terms    <- aov_tab[setdiff(seq_len(nrow(aov_tab)), resid_row), , drop = FALSE]
    if (nrow(terms) == 0) return(data.frame())

    ss  <- terms$`Sum Sq`
    dft <- terms$Df
    fv  <- terms$`F value`

    eta         <- ss / ss_total
    partial_eta <- ss / (ss + sse)
    omega       <- pmax(0, (ss - dft * mse) / (ss_total + mse))
    partial_omega <- pmax(0, (ss - dft * mse) / (ss + sse + mse))

    out <- data.frame(Term = terms$Term, stringsAsFactors = FALSE)
    es_on <- isTRUE(options$effectSizeEstimates)
    if (es_on && isTRUE(options$effectSizeEtaSquared))          out$`eta squared`          <- .je_fmt(eta)
    if (es_on && isTRUE(options$effectSizePartialEtaSquared))   out$`partial eta squared`  <- .je_fmt(partial_eta)
    if (es_on && isTRUE(options$effectSizeOmegaSquared))        out$`omega squared`        <- .je_fmt(omega)
    if (es_on && isTRUE(options$effectSizePartialOmegaSquared)) out$`partial omega squared` <- .je_fmt(partial_omega)

    # Confidence intervals for partial eta squared via noncentral F
    if (isTRUE(options$effectSizeCi) && !is.null(fv)) {
        ci_level <- (options$effectSizeCiLevel %||% 95) / 100
        cis <- mapply(.je_partial_eta_ci, fv, dft, df_err, n,
                      MoreArgs = list(ci_level = ci_level), SIMPLIFY = FALSE)
        out$`peta2 CI lower` <- .je_fmt(vapply(cis, function(x) x["lower"], numeric(1)))
        out$`peta2 CI upper` <- .je_fmt(vapply(cis, function(x) x["upper"], numeric(1)))
    }

    out
}

.je_partial_eta_ci <- function(f_val, df1, df2, n, ci_level = 0.95) {
    if (is.na(f_val) || !is.finite(f_val) || f_val <= 0 || df1 <= 0 || df2 <= 0)
        return(c(lower = NA_real_, upper = NA_real_))
    alpha <- 1 - ci_level

    find_ncp <- function(target_p) {
        if (stats::pf(f_val, df1, df2, ncp = 0, lower.tail = FALSE) <= target_p)
            return(0)
        tryCatch(
            stats::uniroot(
                function(ncp) stats::pf(f_val, df1, df2, ncp = ncp, lower.tail = FALSE) - target_p,
                lower = 0, upper = f_val * df1 * 200 + 200, tol = 1e-6
            )$root,
            error = function(e) NA_real_
        )
    }

    ncp_lower <- find_ncp(1 - alpha / 2)
    ncp_upper <- find_ncp(alpha / 2)

    # eta²_p from ncp: ncp / (ncp + N) where N = total obs
    to_eta <- function(ncp) if (is.na(ncp) || ncp <= 0) 0 else ncp / (ncp + n)
    c(lower = max(0, to_eta(ncp_lower)),
      upper = min(1, to_eta(ncp_upper %||% NA_real_)))
}

# ── Vovk-Sellke maximum p-ratio ───────────────────────────────────────────────

.je_vovk_sellke <- function(p) {
    ifelse(is.na(p) | p >= 1 / exp(1), 1, -exp(1) * p * log(p))
}

# ── Levene / homogeneity ──────────────────────────────────────────────────────

.je_levene_test <- function(dat, dep, factors, center = "mean") {
    if (length(factors) == 0) return(NULL)
    group <- interaction(dat[, factors, drop = FALSE], drop = TRUE)
    y     <- dat[[dep]]
    z     <- if (center == "mean")
        abs(y - stats::ave(y, group, FUN = mean))
    else
        abs(y - stats::ave(y, group, FUN = stats::median))
    fit <- stats::lm(z ~ group)
    tab <- stats::anova(fit)
    list(F = tab$`F value`[1], df1 = tab$Df[1], df2 = tab$Df[2], p = tab$`Pr(>F)`[1],
         label = if (center == "mean") "Levene" else "Brown-Forsythe")
}

# ── Assumptions ───────────────────────────────────────────────────────────────

.je_assumptions <- function(dat, dep, factors, fit, options) {
    out <- character()

    # Residual normality
    res <- stats::residuals(fit)
    if (length(res) >= 3 && length(res) <= 5000) {
        sh  <- stats::shapiro.test(res)
        out <- c(out, paste0(
            "<p>Residual normality (Shapiro-Wilk): W = ", .je_fmt(unname(sh$statistic)),
            ", p = ", .je_p(sh$p.value), ".</p>"
        ))
    } else {
        out <- c(out, "<p>Shapiro-Wilk normality test skipped (n outside 3-5000).</p>")
    }

    # Homogeneity of variance
    if (length(factors) > 0) {
        lev <- .je_levene_test(dat, dep, factors, center = "mean")
        bf  <- .je_levene_test(dat, dep, factors, center = "median")
        if (!is.null(lev))
            out <- c(out, paste0(
                "<p>Levene's test (mean-based): F(", lev$df1, ", ", lev$df2, ") = ",
                .je_fmt(lev$F), ", p = ", .je_p(lev$p), ".</p>"
            ))
        if (!is.null(bf))
            out <- c(out, paste0(
                "<p>Brown-Forsythe test (median-based): F(", bf$df1, ", ", bf$df2, ") = ",
                .je_fmt(bf$F), ", p = ", .je_p(bf$p), ".</p>"
            ))
    } else {
        out <- c(out, "<p>Homogeneity test not applicable — no fixed factors.</p>")
    }

    paste(out, collapse = "")
}

# ── Descriptive statistics ────────────────────────────────────────────────────

.je_descriptives <- function(dat, dep, factors, options) {
    alpha <- 1 - ((options$effectSizeCiLevel %||% 95) / 100)
    if (length(factors) == 0) {
        x  <- dat[[dep]]
        se <- stats::sd(x) / sqrt(length(x))
        ci <- stats::qt(1 - alpha / 2, df = max(1, length(x) - 1)) * se
        return(data.frame(
            Group = "Overall", N = length(x),
            Mean = .je_fmt(mean(x)), SD = .je_fmt(stats::sd(x)),
            SE = .je_fmt(se),
            `CI lower` = .je_fmt(mean(x) - ci), `CI upper` = .je_fmt(mean(x) + ci),
            check.names = FALSE
        ))
    }

    group    <- interaction(dat[, factors, drop = FALSE], drop = TRUE, sep = " × ")
    split_y  <- split(dat[[dep]], group)
    do.call(rbind, lapply(names(split_y), function(g) {
        x  <- split_y[[g]]
        se <- stats::sd(x) / sqrt(length(x))
        ci <- stats::qt(1 - alpha / 2, df = max(1, length(x) - 1)) * se
        data.frame(
            Group = g, N = length(x),
            Mean = .je_fmt(mean(x)), SD = .je_fmt(stats::sd(x)),
            SE = .je_fmt(se),
            `CI lower` = .je_fmt(mean(x) - ci), `CI upper` = .je_fmt(mean(x) + ci),
            check.names = FALSE
        )
    }))
}

# ── Post hoc tests ────────────────────────────────────────────────────────────

.je_posthoc <- function(dat, dep, factors, fit, options) {
    if (length(.je_chr_vec(options$postHocTerms)) == 0)
        return("<p>Assign factors to Post Hoc Terms to run post hoc tests.</p>")
    if (length(factors) == 0) return("<p>Post hoc tests require at least one fixed factor.</p>")

    if (requireNamespace("emmeans", quietly = TRUE))
        return(.je_emmeans_posthoc(fit, factors, options))

    .je_fallback_posthoc(dat, dep, factors, options)
}

.je_emmeans_posthoc <- function(fit, factors, options) {
    chunks <- character()
    for (fac in factors) {
        emm <- tryCatch(emmeans::emmeans(fit, specs = fac), error = function(e) NULL)
        if (is.null(emm)) next

        methods <- character()
        if (isTRUE(options$postHocCorrectionTukey))    methods <- c(methods, "tukey")
        if (isTRUE(options$postHocCorrectionScheffe))  methods <- c(methods, "scheffe")
        if (isTRUE(options$postHocCorrectionBonferroni)) methods <- c(methods, "bonferroni")
        if (isTRUE(options$postHocCorrectionHolm))     methods <- c(methods, "holm")
        if (isTRUE(options$postHocCorrectionSidak))    methods <- c(methods, "sidak")
        if (length(methods) == 0) methods <- "tukey"

        for (adj in methods) {
            pw <- tryCatch(
                as.data.frame(emmeans::contrast(emm, method = "pairwise", adjust = adj)),
                error = function(e) NULL
            )
            if (is.null(pw)) next
            ci_level <- (options$effectSizeCiLevel %||% 95) / 100
            pw_ci <- tryCatch(
                as.data.frame(stats::confint(emmeans::contrast(emm, method = "pairwise", adjust = adj),
                                             level = ci_level)),
                error = function(e) NULL
            )
            out <- data.frame(
                Comparison = as.character(pw$contrast),
                Estimate   = .je_fmt(pw$estimate),
                SE         = .je_fmt(pw$SE),
                df         = .je_fmt(pw$df),
                t          = .je_fmt(pw$t.ratio),
                p          = .je_p(pw$p.value),
                check.names = FALSE
            )
            if (isTRUE(options$postHocCi) && !is.null(pw_ci)) {
                out$`lower CI` <- .je_fmt(pw_ci$lower.CL)
                out$`upper CI` <- .je_fmt(pw_ci$upper.CL)
            }
            if (isTRUE(options$vovkSellke))
                out$`VS-MPR` <- .je_fmt(.je_vovk_sellke(pw$p.value))

            chunks <- c(chunks, paste0(
                "<h4>", .je_escape(toupper(substr(adj, 1, 1))),
                .je_escape(substr(adj, 2, nchar(adj))), " — ", .je_escape(fac), "</h4>",
                .je_table_html(out)
            ))
        }

        # Games-Howell (unequal variances, requires manual implementation)
        if (isTRUE(options$postHocTypeGames)) {
            gh <- .je_games_howell(dat, dep, fac)
            if (!is.null(gh))
                chunks <- c(chunks, paste0("<h4>Games-Howell — ", .je_escape(fac), "</h4>", .je_table_html(gh)))
        }
    }

    if (length(chunks) == 0) return("<p>No post hoc output could be calculated.</p>")

    # Interaction post hoc (conditional comparisons)
    if (length(factors) >= 2 && isTRUE(options$postHocConditionalTable))
        chunks <- c(chunks, .je_emmeans_interaction_posthoc(fit, factors, options))

    paste(chunks, collapse = "")
}

.je_games_howell <- function(dat, dep, factor) {
    groups  <- levels(as.factor(dat[[factor]]))
    k       <- length(groups)
    if (k < 2) return(NULL)
    stats_list <- lapply(groups, function(g) {
        x <- dat[[dep]][dat[[factor]] == g]
        list(n = length(x), mean = mean(x), var = stats::var(x))
    })
    names(stats_list) <- groups

    rows <- list()
    for (i in seq_len(k - 1)) for (j in (i + 1):k) {
        gi <- stats_list[[groups[i]]]; gj <- stats_list[[groups[j]]]
        ni <- gi$n; nj <- gj$n; vi <- gi$var; vj <- gj$var
        if (ni < 2 || nj < 2) next
        se    <- sqrt(vi / ni + vj / nj)
        t_val <- (gi$mean - gj$mean) / se
        df_w  <- (vi / ni + vj / nj)^2 / ((vi / ni)^2 / (ni - 1) + (vj / nj)^2 / (nj - 1))
        p_val <- 1 - stats::ptukey(abs(t_val) * sqrt(2), nmeans = k, df = df_w)
        rows[[length(rows) + 1]] <- data.frame(
            Comparison = paste0(groups[i], " - ", groups[j]),
            Difference = .je_fmt(gi$mean - gj$mean),
            SE = .je_fmt(se),
            df = .je_fmt(df_w),
            p  = .je_p(p_val),
            check.names = FALSE
        )
    }
    if (length(rows) == 0) return(NULL)
    do.call(rbind, rows)
}

.je_emmeans_interaction_posthoc <- function(fit, factors, options) {
    if (length(factors) < 2 || !requireNamespace("emmeans", quietly = TRUE)) return("")
    adj <- .je_posthoc_adjust(options)
    chunks <- character()
    for (i in seq_len(length(factors))) {
        simple_fac <- factors[i]
        by_fac     <- factors[-i]
        emm <- tryCatch(
            emmeans::emmeans(fit, specs = simple_fac, by = by_fac),
            error = function(e) NULL
        )
        if (is.null(emm)) next
        pw <- tryCatch(
            as.data.frame(emmeans::contrast(emm, method = "pairwise", adjust = adj)),
            error = function(e) NULL
        )
        if (is.null(pw)) next
        out <- data.frame(
            Comparison = as.character(pw$contrast),
            By         = if ("by" %in% names(pw)) as.character(pw[[by_fac[1]]]) else "",
            Estimate   = .je_fmt(pw$estimate),
            SE         = .je_fmt(pw$SE),
            p          = .je_p(pw$p.value),
            check.names = FALSE
        )
        chunks <- c(chunks, paste0(
            "<h4>Conditional comparisons: ", .je_escape(simple_fac),
            " by ", .je_escape(paste(by_fac, collapse = " × ")), "</h4>",
            .je_table_html(out)
        ))
    }
    paste(chunks, collapse = "")
}

.je_fallback_posthoc <- function(dat, dep, factors, options) {
    chunks <- character()
    for (fac in factors) {
        formula <- stats::as.formula(paste(dep, "~", fac))
        fit     <- tryCatch(stats::aov(formula, data = dat), error = function(e) NULL)
        if (!is.null(fit) && isTRUE(options$postHocCorrectionTukey)) {
            tk <- tryCatch(stats::TukeyHSD(fit), error = function(e) NULL)
            if (!is.null(tk)) {
                tab <- as.data.frame(tk[[1]])
                tab$Comparison <- rownames(tab)
                rownames(tab) <- NULL
                out <- data.frame(
                    Factor = fac, Comparison = tab$Comparison,
                    Difference = .je_fmt(tab$diff),
                    Lower = .je_fmt(tab$lwr), Upper = .je_fmt(tab$upr),
                    p = .je_p(tab$`p adj`), check.names = FALSE
                )
                chunks <- c(chunks, paste0("<h4>Tukey HSD: ", .je_escape(fac), "</h4>", .je_table_html(out)))
            }
        }
        adj <- .je_posthoc_adjust(options)
        pw  <- tryCatch(
            stats::pairwise.t.test(dat[[dep]], dat[[fac]], p.adjust.method = adj),
            error = function(e) NULL
        )
        if (!is.null(pw))
            chunks <- c(chunks, paste0(
                "<h4>Pairwise t (", .je_escape(adj), "): ", .je_escape(fac), "</h4>",
                .je_pairwise_html(pw$p.value)
            ))
    }
    if (length(chunks) == 0) return("<p>No post hoc output could be calculated.</p>")
    paste0("<p><em>Install <code>emmeans</code> for full JASP-equivalent post hoc tests.</em></p>",
           paste(chunks, collapse = ""))
}

# ── Marginal means ────────────────────────────────────────────────────────────

.je_marginal_means <- function(dat, dep, factors, covariates, fit, options) {
    terms <- intersect(.je_chr_vec(options$marginalMeanTerms), factors)
    if (length(terms) == 0) return("<p>Assign factors to Marginal Means Terms to compute estimated marginal means.</p>")
    if (length(terms) == 0) terms <- factors
    if (length(terms) == 0) return("<p>Estimated marginal means require at least one fixed factor.</p>")

    if (requireNamespace("emmeans", quietly = TRUE)) {
        specs <- if (length(terms) == 1) terms else terms
        emm <- tryCatch(emmeans::emmeans(fit, specs = terms), error = function(e) NULL)
        if (!is.null(emm)) {
            ci_level <- (options$effectSizeCiLevel %||% 95) / 100
            em_df    <- as.data.frame(stats::confint(emm, level = ci_level))
            out <- data.frame(
                Term  = paste(terms, collapse = " × "),
                Level = do.call(paste, c(em_df[terms], sep = " × ")),
                EMM   = .je_fmt(em_df$emmean),
                SE    = .je_fmt(em_df$SE),
                df    = .je_fmt(em_df$df),
                Lower = .je_fmt(em_df$lower.CL),
                Upper = .je_fmt(em_df$upper.CL),
                check.names = FALSE
            )
            html <- paste0(
                "<p>Estimated marginal means (covariate-adjusted via <code>emmeans</code>).</p>",
                .je_table_html(out)
            )

            if (isTRUE(options$marginalMeansPairwise %||% FALSE) && length(terms) == 1) {
                adj <- .je_adjust_method(options$marginalMeanCiCorrection)
                pw  <- tryCatch(
                    as.data.frame(emmeans::contrast(emm, method = "pairwise", adjust = adj)),
                    error = function(e) NULL
                )
                if (!is.null(pw)) {
                    pw_out <- data.frame(
                        Comparison = as.character(pw$contrast),
                        Estimate   = .je_fmt(pw$estimate),
                        SE         = .je_fmt(pw$SE),
                        t          = .je_fmt(pw$t.ratio),
                        p          = .je_p(pw$p.value),
                        check.names = FALSE
                    )
                    if (isTRUE(options$vovkSellke))
                        pw_out$`VS-MPR` <- .je_fmt(.je_vovk_sellke(pw$p.value))
                    html <- paste0(html, "<h4>Pairwise comparisons (", .je_escape(adj), ")</h4>",
                                   .je_table_html(pw_out))
                }
            }

            if (isTRUE(options$marginalMeanComparedToZero) && length(terms) == 1) {
                cz <- tryCatch(
                    as.data.frame(emmeans::contrast(emm, method = "eff")),
                    error = function(e) NULL
                )
                if (!is.null(cz)) {
                    cz_out <- data.frame(
                        Level    = as.character(cz$contrast),
                        Estimate = .je_fmt(cz$estimate),
                        SE       = .je_fmt(cz$SE),
                        t        = .je_fmt(cz$t.ratio),
                        p        = .je_p(cz$p.value),
                        check.names = FALSE
                    )
                    html <- paste0(html, "<h4>Comparisons to zero</h4>", .je_table_html(cz_out))
                }
            }
            return(html)
        }
    }

    # Fallback: observed cell means
    group   <- interaction(dat[, terms, drop = FALSE], drop = TRUE, sep = " × ")
    split_y <- split(dat[[dep]], group)
    alpha   <- 1 - ((options$effectSizeCiLevel %||% 95) / 100)
    out <- do.call(rbind, lapply(names(split_y), function(g) {
        x  <- split_y[[g]]
        se <- stats::sd(x) / sqrt(length(x))
        ci <- stats::qt(1 - alpha / 2, df = max(1, length(x) - 1)) * se
        data.frame(
            Term = paste(terms, collapse = " × "), Level = g, N = length(x),
            Mean = .je_fmt(mean(x)), SE = .je_fmt(se),
            Lower = .je_fmt(mean(x) - ci), Upper = .je_fmt(mean(x) + ci),
            check.names = FALSE
        )
    }))
    paste0("<p><em>Showing observed cell means. Install <code>emmeans</code> for covariate-adjusted estimates.</em></p>",
           .je_table_html(out))
}

# ── Simple effects ────────────────────────────────────────────────────────────

.je_simple_effects <- function(dat, dep, factors, fit, options) {
    simple <- .je_chr(options$simpleMainEffectFactor)
    if (is.null(simple) || !nzchar(simple)) return("<p>Assign a factor to Simple Effect Factor to run simple effects analysis.</p>")
    mod1   <- .je_chr(options$simpleMainEffectModeratorFactorOne)
    mod2   <- .je_chr(options$simpleMainEffectModeratorFactorTwo)
    if (is.null(simple) || !(simple %in% factors))
        return("<p>Select a simple effect factor from the fixed factors.</p>")
    moderators <- intersect(c(mod1, mod2), factors)
    moderators <- setdiff(moderators[nzchar(moderators)], simple)
    if (length(moderators) == 0)
        return("<p>Select at least one moderator factor for simple effects.</p>")

    strata <- interaction(dat[, moderators, drop = FALSE], drop = TRUE, sep = " × ")
    pieces <- lapply(levels(strata), function(level) {
        sub <- dat[strata == level, , drop = FALSE]
        if (length(unique(sub[[simple]])) < 2 || nrow(sub) < 3)
            return(paste0("<h4>", .je_escape(level), "</h4><p>Insufficient data for this stratum.</p>"))
        f2 <- .je_fit_and_anova(
            stats::as.formula(paste(dep, "~", simple)),
            sub, simple, NULL
        )
        html <- paste0("<h4>", .je_escape(level), "</h4>",
                       .je_model_html(f2$aov_tab, stats::as.formula(paste(dep, "~", simple)),
                                      nrow(sub), NULL, f2$ss_type, NULL, options))
        if (isTRUE(options$simpleEffectsPostHoc %||% FALSE) && !inherits(f2$fit, "error"))
            html <- paste0(html, .je_posthoc(sub, dep, simple, f2$fit, options))
        html
    })
    paste0("<p>Simple effects are one-way ANOVAs within each stratum of the moderator(s).</p>",
           paste(pieces, collapse = ""))
}

# ── Contrasts ─────────────────────────────────────────────────────────────────

.je_contrasts <- function(dat, dep, factors, fit, options) {
    ct <- as.character(options$contrastType %||% "none")
    if (ct == "none") return("<p>Select a contrast type in the Contrasts section.</p>")
    if (length(factors) == 0) return("<p>Contrasts require at least one fixed factor.</p>")

    chunks <- character()
    for (fac in factors) {
        levs <- levels(dat[[fac]])
        mats <- list()
        if (ct == "helmert")    mats$Helmert    <- stats::contr.helmert(length(levs))
        if (ct == "polynomial") mats$Polynomial <- stats::contr.poly(length(levs))
        if (ct == "deviation")  mats$Deviation  <- stats::contr.sum(length(levs))
        if (ct %in% c("difference", "repeated"))
            mats$Difference <- .je_contr_diff(length(levs))
        if (ct == "simple") {
            # Simple contrasts: compare each level to the first
            mat <- matrix(0, nrow = length(levs), ncol = length(levs) - 1)
            for (j in seq_len(ncol(mat))) { mat[1, j] <- -1; mat[j + 1, j] <- 1 }
            mats$Simple <- mat
        }

        if (FALSE) {  # custom contrasts disabled in this UI
            custom <- .je_parse_custom_contrasts("", length(levs))
            if (nrow(custom) > 0) {
                mats$Custom <- t(custom)
                chunks <- c(chunks, .je_custom_contrast_tests(dat, dep, fac, custom, fit))
            } else {
                chunks <- c(chunks, paste0(
                    "<h4>Custom contrasts: ", .je_escape(fac), "</h4>",
                    "<p>Enter one contrast per line with ", length(levs),
                    " numeric weights matching levels: ",
                    .je_escape(paste(levs, collapse = ", ")), ".</p>"
                ))
            }
        }

        for (nm in names(mats)) {
            mat <- as.data.frame(mats[[nm]])
            names(mat) <- paste0("C", seq_len(ncol(mat)))
            mat$Level <- levs
            mat <- mat[, c("Level", setdiff(names(mat), "Level")), drop = FALSE]
            chunks <- c(chunks, paste0(
                "<h4>", .je_escape(nm), " contrasts: ", .je_escape(fac), "</h4>",
                .je_table_html(mat)
            ))
        }
    }
    paste0("<p>Contrast matrices generated for selected types. Custom one-factor inference uses the pooled residual MSE from the full model.</p>",
           paste(chunks, collapse = ""))
}

.je_parse_custom_contrasts <- function(text, n_levels) {
    if (is.null(text) || !nzchar(text))
        return(matrix(numeric(), nrow = 0, ncol = n_levels))
    lines <- unlist(strsplit(text, "\n", fixed = TRUE), use.names = FALSE)
    rows  <- lapply(lines, function(line) {
        line <- trimws(sub(".*:", "", line))
        if (!nzchar(line)) return(NULL)
        vals <- suppressWarnings(as.numeric(unlist(strsplit(gsub(",", " ", line), "\\s+"))))
        vals <- vals[!is.na(vals)]
        if (length(vals) != n_levels) return(NULL)
        vals
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) == 0)
        return(matrix(numeric(), nrow = 0, ncol = n_levels))
    do.call(rbind, rows)
}

.je_custom_contrast_tests <- function(dat, dep, factor, contrasts, fit) {
    levs     <- levels(dat[[factor]])
    means    <- tapply(dat[[dep]], dat[[factor]], mean)
    ns       <- tapply(dat[[dep]], dat[[factor]], length)
    aov_tab  <- as.data.frame(stats::anova(fit))
    resid_row <- which(rownames(aov_tab) == "Residuals")
    if (length(resid_row) == 0)
        return("<p>Custom contrasts could not be tested — residual error unavailable.</p>")
    mse  <- aov_tab$`Mean Sq`[resid_row]
    df_e <- aov_tab$Df[resid_row]
    rows <- lapply(seq_len(nrow(contrasts)), function(i) {
        w   <- contrasts[i, ]
        est <- sum(w * means[levs])
        se  <- sqrt(mse * sum((w^2) / ns[levs]))
        t   <- est / se
        p   <- 2 * stats::pt(abs(t), df = df_e, lower.tail = FALSE)
        data.frame(
            Contrast = paste0("C", i),
            Weights  = paste(w, collapse = ", "),
            Estimate = .je_fmt(est), SE = .je_fmt(se),
            t = .je_fmt(t), df = .je_fmt(df_e), p = .je_p(p),
            check.names = FALSE
        )
    })
    paste0("<h4>Custom contrast tests: ", .je_escape(factor), "</h4>",
           .je_table_html(do.call(rbind, rows)))
}

# ── Kruskal-Wallis ────────────────────────────────────────────────────────────

.je_kruskal <- function(dat, dep, factors, options) {
    if (length(factors) == 0)
        return("<p>Kruskal-Wallis requires at least one grouping factor.</p>")

    chunks <- character()
    for (fac in factors) {
        kt <- tryCatch(
            stats::kruskal.test(stats::as.formula(paste(dep, "~", fac)), data = dat),
            error = function(e) NULL
        )
        if (is.null(kt)) next
        H <- unname(kt$statistic)
        df <- kt$parameter
        p  <- kt$p.value
        n  <- nrow(dat)

        chunk <- paste0(
            "<h4>Kruskal-Wallis: ", .je_escape(fac), "</h4>",
            "<p>H(", df, ") = ", .je_fmt(H), ", p = ", .je_p(p), ", N = ", n, ".</p>"
        )

        if (isTRUE(options$kruskalEpsilon)) {
            eps <- H / (n - 1)
            chunk <- paste0(chunk, "<p>Epsilon squared = ", .je_fmt(eps), ".</p>")
        }
        if (isTRUE(options$kruskalEta)) {
            k   <- nlevels(as.factor(dat[[fac]]))
            eta <- (H - k + 1) / (n - k)
            chunk <- paste0(chunk, "<p>Eta squared (H-based) = ", .je_fmt(max(0, eta)), ".</p>")
        }
        if (isTRUE(options$postHocTypeDunn)) {
            chunk <- paste0(chunk, .je_dunn_posthoc(dat, dep, fac, options))
        }
        chunks <- c(chunks, chunk)
    }
    if (length(chunks) == 0) return("<p>Kruskal-Wallis could not be calculated.</p>")
    paste(chunks, collapse = "")
}

.je_dunn_posthoc <- function(dat, dep, fac, options) {
    groups <- split(dat[[dep]], dat[[fac]])
    if (length(groups) < 2) return("")
    n_all  <- nrow(dat)
    ranks  <- rank(dat[[dep]])
    ns     <- vapply(groups, length, integer(1))
    g_names <- names(groups)
    k       <- length(g_names)

    rank_means <- tapply(ranks, dat[[fac]], mean)
    rows <- list()
    for (i in seq_len(k - 1)) for (j in (i + 1):k) {
        ni <- ns[g_names[i]]; nj <- ns[g_names[j]]
        z  <- (rank_means[g_names[i]] - rank_means[g_names[j]]) /
              sqrt((n_all * (n_all + 1) / 12) * (1 / ni + 1 / nj))
        p  <- 2 * stats::pnorm(-abs(z))
        rows[[length(rows) + 1]] <- data.frame(
            Comparison = paste0(g_names[i], " - ", g_names[j]),
            z = .je_fmt(z), p = .je_p(p), check.names = FALSE
        )
    }
    if (length(rows) == 0) return("")
    paste0("<h4>Dunn post hoc (uncorrected)</h4>", .je_table_html(do.call(rbind, rows)))
}

# ── Bootstrap effect-size CIs ─────────────────────────────────────────────────

.je_bootstrap_effects <- function(dat, dep, factors, covariates, options) {
    if (!isTRUE(options$bootstrapCi %||% FALSE) && !isTRUE(options$postHocTypeStandardBootstrap %||% FALSE)) return("")
    rhs <- c(factors, covariates)
    if (length(rhs) == 0) return("<p>Bootstrap CIs require model terms.</p>")
    reps    <- max(100L, min(as.integer(options$bootstrapSamples %||% 1000), 2000L))
    formula <- stats::as.formula(paste(dep, "~", paste(rhs, collapse = " * ")))
    obs_eff <- tryCatch({
        f <- .je_fit_and_anova(formula, dat, factors, NULL)
        .je_effect_sizes(f$aov_tab, nrow(dat), options)
    }, error = function(e) data.frame())
    if (nrow(obs_eff) == 0) return("<p>Bootstrap could not start — model failed on observed data.</p>")

    boot_vals <- replicate(reps, {
        idx  <- sample.int(nrow(dat), replace = TRUE)
        bdat <- dat[idx, , drop = FALSE]
        for (f in factors) bdat[[f]] <- factor(bdat[[f]], levels = levels(dat[[f]]))
        tryCatch({
            bf  <- .je_fit_and_anova(formula, bdat, factors, NULL)
            eff <- .je_effect_sizes(bf$aov_tab, nrow(bdat), options)
            if ("partial eta squared" %in% names(eff))
                suppressWarnings(as.numeric(eff$`partial eta squared`))
            else
                rep(NA_real_, nrow(obs_eff))
        }, error = function(e) rep(NA_real_, nrow(obs_eff)))
    })
    if (is.null(dim(boot_vals)))
        boot_vals <- matrix(boot_vals, nrow = nrow(obs_eff))
    ci <- t(apply(boot_vals, 1, stats::quantile, probs = c(.025, .975), na.rm = TRUE))
    out <- data.frame(
        Term = obs_eff$Term,
        `partial eta squared` = if ("partial eta squared" %in% names(obs_eff)) obs_eff$`partial eta squared` else NA,
        `boot lower 2.5%` = .je_fmt(ci[, 1]),
        `boot upper 97.5%` = .je_fmt(ci[, 2]),
        Samples = reps,
        check.names = FALSE
    )
    paste0("<h4>Bootstrap effect-size confidence intervals (", reps, " samples)</h4>",
           .je_table_html(out))
}

# ── APA text ──────────────────────────────────────────────────────────────────

.je_apa <- function(aov_tab, effect_tab, dep, factors, ss_type) {
    test_rows <- aov_tab[aov_tab$Term != "Residuals", , drop = FALSE]
    if (nrow(test_rows) == 0) return("<p>No APA text available for this model.</p>")
    row <- test_rows[1, ]
    df_err  <- aov_tab$Df[aov_tab$Term == "Residuals"]
    eta_str <- if (nrow(effect_tab) > 0 && "partial eta squared" %in% names(effect_tab))
        paste0(", partial η² = ", effect_tab$`partial eta squared`[1]) else ""
    fac <- if (length(factors) > 0) paste(factors, collapse = ", ") else "the model term"
    paste0(
        "<p>A one-way ANOVA (Type ", .je_escape(ss_type), " SS) was conducted to examine the effect of ",
        .je_escape(fac), " on ", .je_escape(dep), ". ",
        "The omnibus test was <em>F</em>(", .je_fmt(row$Df), ", ",
        if (length(df_err) > 0) .je_fmt(df_err[1]) else "?",
        ") = ", .je_fmt(row$`F value`), ", <em>p</em> = ", .je_p(row$`Pr(>F)`), eta_str, ".</p>"
    )
}

# ── HTML helpers ──────────────────────────────────────────────────────────────

.je_model_html <- function(aov_tab, formula, n, weights, ss_type, correction, options) {
    display <- data.frame(
        Term       = aov_tab$Term,
        df         = .je_fmt(aov_tab$Df),
        `Sum Sq`   = .je_fmt(aov_tab$`Sum Sq`),
        `Mean Sq`  = .je_fmt(aov_tab$`Mean Sq`),
        F          = .je_fmt(aov_tab$`F value`),
        p          = .je_p(aov_tab$`Pr(>F)`),
        check.names = FALSE
    )

    # Vovk-Sellke column
    if (isTRUE(options$vovkSellke)) {
        vs <- .je_vovk_sellke(aov_tab$`Pr(>F)`)
        display$`VS-MPR` <- ifelse(is.na(vs), "", .je_fmt(vs))
    }

    wt_note <- if (!is.null(weights) && nzchar(weights))
        paste0("<p>WLS weights: <code>", .je_escape(weights), "</code></p>") else ""

    corr_html <- if (!is.null(correction)) paste0(
        "<table><thead><tr>",
        "<th>Correction</th><th>F</th><th>df1</th><th>df2</th><th>p</th>",
        "</tr></thead><tbody><tr><td>", .je_escape(correction$method), "</td>",
        "<td>", .je_fmt(correction$F), "</td><td>", .je_fmt(correction$df1), "</td>",
        "<td>", .je_fmt(correction$df2), "</td><td>", .je_p(correction$p), "</td>",
        "</tr></tbody></table>"
    ) else ""

    paste0(
        "<p>Model: <code>", .je_escape(deparse(formula)), "</code> &mdash; Type ", .je_escape(ss_type), " SS</p>",
        "<p>Complete cases: ", n, "</p>",
        wt_note,
        .je_table_html(display),
        if (nzchar(corr_html)) paste0("<h4>Homogeneity correction</h4>", corr_html) else ""
    )
}

# ── Order-restricted + Bayesian ANOVA dispatcher ─────────────────────────────

.je_order_restricted_full <- function(dat, dep, factors, covs, fit, aov_tab, options) {
    any_active <- isTRUE(options$orderRestricted %||% FALSE) ||
                  isTRUE(options$modelComparison %||% FALSE)  ||
                  isTRUE(options$informedHypothesisTests %||% FALSE) ||
                  nzchar(trimws(as.character(options$restrictedSyntax %||% "")))
    if (!any_active)
        return("<p>Bayesian / order-restricted hypothesis testing is disabled.</p>")

    chunks <- character()

    # ── Item 2: Order-restricted inference via bain ───────────────────────────
    if ((isTRUE(options$orderRestricted %||% FALSE) || isTRUE(options$informedHypothesisTests %||% FALSE)) &&
        nzchar(trimws(as.character(options$orderRestrictedSyntax %||% options$restrictedSyntax %||% "")))) {
        chunks <- c(chunks, .je_bain_analysis(fit, options))
    } else if (isTRUE(options$orderRestricted %||% FALSE) || isTRUE(options$informedHypothesisTests %||% FALSE)) {
        chunks <- c(chunks, paste0(
            "<p><strong>Order-restricted inference (bain):</strong> ",
            "enter hypothesis constraints in the syntax field, e.g. ",
            "<code>mu1 &gt; mu2 &gt; mu3</code> or <code>mu1 = mu2 &amp; mu3 &gt; mu1</code>.</p>"
        ))
    }

    # ── Item 1: Bayesian ANOVA via BayesFactor ────────────────────────────────
    if (isTRUE(options$modelComparison %||% FALSE))
        chunks <- c(chunks, .je_bayesian_anova(dat, dep, factors, covs, options))

    paste(chunks, collapse = "")
}

# ── Item 2: bain — order-restricted / informed hypothesis testing ─────────────

.je_bain_analysis <- function(fit, options) {
    h_raw <- trimws(as.character(options$orderRestrictedSyntax %||% options$restrictedSyntax %||% ""))
    if (!nzchar(h_raw))
        return("<p>No hypothesis syntax provided.</p>")

    if (!requireNamespace("bain", quietly = TRUE))
        return(paste0(
            "<p><strong>bain not installed.</strong> Install it with ",
            "<code>install.packages('bain')</code> to run order-restricted inference.</p>",
            "<p>Hypothesis syntax captured: <pre>", .je_escape(h_raw), "</pre></p>"
        ))

    result <- tryCatch({
        br <- bain::bain(fit, hypothesis = h_raw)
        br
    }, error = function(e) e)

    if (inherits(result, "error")) {
        return(paste0(
            "<p>bain could not process the hypothesis: <em>", .je_escape(result$message), "</em></p>",
            "<p>Syntax: <pre>", .je_escape(h_raw), "</pre></p>",
            "<p>Use <code>bain</code> coefficient names (visible in <code>coef(model)</code>). ",
            "Example: <code>groupB &gt; 0 &amp; groupC &gt; groupB</code></p>"
        ))
    }

    # Format the bain output
    .je_bain_html(result, h_raw)
}

.je_bain_html <- function(br, h_raw) {
    # bain result has: $fit (model fit table), $b (BF matrix), $BFmatrix
    fit_tab <- tryCatch(as.data.frame(br$fit), error = function(e) NULL)
    bf_tab  <- tryCatch({
        bfm <- br$BFmatrix
        if (!is.null(bfm)) {
            tab <- as.data.frame(round(bfm, 3))
            tab$Hypothesis <- rownames(tab)
            tab <- tab[, c("Hypothesis", setdiff(names(tab), "Hypothesis")), drop = FALSE]
            tab
        } else NULL
    }, error = function(e) NULL)

    chunks <- paste0("<h4>Order-restricted / Informed Hypothesis Test (bain)</h4>",
                     "<p>Hypothesis: <code>", .je_escape(h_raw), "</code></p>")

    if (!is.null(fit_tab)) {
        # bain::$fit contains: Fit_measure, Com_measure, BF_c, PMPb, PMPa
        display <- data.frame(check.names = FALSE)
        for (nm in names(fit_tab))
            display[[nm]] <- .je_fmt(as.numeric(fit_tab[[nm]]))
        display <- cbind(data.frame(Hypothesis = rownames(fit_tab), stringsAsFactors = FALSE),
                         display)
        chunks <- paste0(chunks, "<h5>Fit &amp; Complexity</h5>", .je_table_html(display))
    }

    if (!is.null(bf_tab))
        chunks <- paste0(chunks, "<h5>Bayes Factor Matrix</h5>", .je_table_html(bf_tab))

    # Posterior model probabilities
    pmp <- tryCatch({
        pm <- br$fit[, "PMPb", drop = TRUE]
        if (!is.null(pm)) {
            tab <- data.frame(
                Hypothesis = names(pm),
                PMP        = .je_fmt(as.numeric(pm)),
                check.names = FALSE
            )
            tab
        } else NULL
    }, error = function(e) NULL)

    if (!is.null(pmp))
        chunks <- paste0(chunks, "<h5>Posterior Model Probabilities</h5>",
                         .je_table_html(pmp))

    chunks
}

# ── Item 1: Bayesian ANOVA via BayesFactor ────────────────────────────────────

.je_bayesian_anova <- function(dat, dep, factors, covs, options) {
    if (!requireNamespace("BayesFactor", quietly = TRUE))
        return(paste0(
            "<p><strong>BayesFactor not installed.</strong> Install it with ",
            "<code>install.packages('BayesFactor')</code> to run Bayesian ANOVA.</p>"
        ))

    if (length(factors) == 0)
        return("<p>Bayesian ANOVA requires at least one fixed factor.</p>")

    rhs_terms <- c(factors, covs)
    formula   <- stats::as.formula(paste(dep, "~", paste(rhs_terms, collapse = " + ")))

    bf_result <- tryCatch(
        BayesFactor::anovaBF(
            formula       = formula,
            data          = dat,
            whichModels   = "withmain",
            progress      = FALSE
        ),
        error = function(e) e
    )

    if (inherits(bf_result, "error"))
        return(paste0(
            "<p>BayesFactor::anovaBF failed: <em>", .je_escape(bf_result$message), "</em></p>"
        ))

    .je_bayesian_anova_html(bf_result, dep)
}

.je_bayesian_anova_html <- function(bf_result, dep) {
    tab <- tryCatch({
        s  <- summary(bf_result)
        # BayesFactor summary returns a data.frame with rownames = model names
        df <- as.data.frame(s)
        df$Model <- rownames(df)
        rownames(df) <- NULL
        df
    }, error = function(e) NULL)

    if (is.null(tab)) {
        # Fallback: extract from the BFBayesFactor object directly
        tab <- tryCatch({
            bfs   <- as.vector(bf_result)
            nms   <- names(bf_result)
            errs  <- attr(bf_result, "error") %||% rep(NA_real_, length(bfs))
            data.frame(
                Model    = nms,
                `BF10`   = .je_fmt(bfs),
                `log(BF)`= .je_fmt(log(bfs)),
                `±Error%`= .je_fmt(errs * 100),
                check.names = FALSE
            )
        }, error = function(e) NULL)
    } else {
        # Normalise column names from summary output
        col_bf  <- intersect(c("bf", "BF", "Bayes factor", "bayes.factor"), names(tab))[1]
        col_err <- intersect(c("error", "Error", "error%"), names(tab))[1]
        out <- data.frame(Model = tab$Model, check.names = FALSE)
        if (!is.na(col_bf))  out$BF10     <- .je_fmt(as.numeric(tab[[col_bf]]))
        if (!is.na(col_bf))  out$`log(BF)`<- .je_fmt(log(pmax(1e-300, as.numeric(tab[[col_bf]]))))
        if (!is.na(col_err)) out$`±Error%`<- .je_fmt(as.numeric(tab[[col_err]]) * 100)
        tab <- out
    }

    # Best model BF against null
    best_note <- tryCatch({
        bfs <- as.vector(bf_result)
        best_idx <- which.max(bfs)
        paste0("<p>Best model: <strong>", .je_escape(names(bf_result)[best_idx]),
               "</strong> — BF<sub>10</sub> = ", .je_fmt(bfs[best_idx]),
               " (vs null model).</p>")
    }, error = function(e) "")

    # Inclusion BFs (averaging over models that include each predictor)
    inclusion_note <- tryCatch({
        inc <- BayesFactor::generalTestBF(bf_result)
        inc_df <- as.data.frame(inc)
        inc_df$Predictor <- rownames(inc_df)
        rownames(inc_df) <- NULL
        col_bf <- intersect(c("bf", "BF", "Bayes factor"), names(inc_df))[1]
        if (!is.na(col_bf)) {
            out <- data.frame(
                Predictor = inc_df$Predictor,
                `Incl BF` = .je_fmt(as.numeric(inc_df[[col_bf]])),
                check.names = FALSE
            )
            paste0("<h5>Inclusion Bayes Factors</h5>", .je_table_html(out))
        } else ""
    }, error = function(e) "")

    paste0(
        "<h4>Bayesian ANOVA (BayesFactor)</h4>",
        "<p>All BF values are relative to the null model (intercept only).</p>",
        best_note,
        if (!is.null(tab) && nrow(tab) > 0) .je_table_html(tab) else "",
        inclusion_note
    )
}

.je_teaching_html <- function(effect_tab, assumptions) {
    paste0(
        "<p>The ANOVA F-test asks whether the between-group variance is larger than would be expected from residual error alone. ",
        "Partial η² (partial eta squared) describes the fraction of explainable variance attributed to each term after removing other sources of variance.</p>",
        assumptions
    )
}

.je_plot_status <- function(dat, dep, factors, options) {
    enabled <- c(
        if (isTRUE(options$qqPlot))   "Q-Q residual plot",
        if (!is.null(.je_chr(options$descriptivePlotHorizontalAxis %||% NULL))) "descriptive plot",
        if (!is.null(.je_chr(options$barPlotHorizontalAxis %||% NULL)))         "bar plot",
        if (!is.null(.je_chr(options$rainCloudHorizontalAxis %||% NULL)))       "raincloud plot"
    )
    if (length(enabled) == 0) return("<p>No plot options are enabled.</p>")
    html <- paste0("<p>Requested: ", .je_escape(paste(enabled, collapse = ", ")), ".</p>")
    if (length(factors) > 0 && (!is.null(.je_chr(options$descriptivePlotHorizontalAxis %||% NULL)) ||
                                 !is.null(.je_chr(options$barPlotHorizontalAxis %||% NULL)) ||
                                 !is.null(.je_chr(options$rainCloudHorizontalAxis %||% NULL)))) {
        desc <- .je_descriptives(dat, dep, factors[1], options)
        html <- paste0(html, "<p>Summary for ", .je_escape(factors[1]), ":</p>", .je_table_html(desc))
    }
    html
}

.je_publication_html <- function(options) {
    requested <- isTRUE(options$publicationMode %||% FALSE) || isTRUE(options$publicationTables %||% FALSE) ||
        isTRUE(options$publicationFigures %||% FALSE) || isTRUE(options$exportWord %||% FALSE) || isTRUE(options$exportPdf %||% FALSE)
    if (!requested) return("<p>Publication mode is disabled.</p>")
    paste0("<p>Publication mode: result tables are produced as copyable HTML. Word/PDF export is handled at the jamovi application level.</p>")
}

.je_saved_columns_html <- function(options) {
    requested <- c(
        if (isTRUE(options$residualsSavedToData) && identical(as.character(options$residualsSavedToDataType %||% "raw"), "raw"))     "raw residuals",
        if (isTRUE(options$residualsSavedToData) && identical(as.character(options$residualsSavedToDataType), "student"))  "studentized residuals",
        if (isTRUE(options$residualsSavedToData) && identical(as.character(options$residualsSavedToDataType), "standard")) "standardized residuals",
        if (isTRUE(options$predictionsSavedToData))                                             "predicted values"
    )
    if (length(requested) == 0) return("<p>No dataset output columns requested.</p>")
    paste0("<p>Output columns: ", .je_escape(paste(requested, collapse = ", ")), ".</p>")
}

.je_repro_html <- function(formula, options) {
    paste0("<p>Model formula:</p><pre>", .je_escape(deparse(formula)), "</pre>",
           "<p>Full R syntax export (including all options) is planned for a future build.</p>")
}

# ── Shared utilities (also used by enhancedRepeatedMeasuresAnova.b.R) ─────────

.je_chr <- function(x) {
    if (is.null(x) || length(x) == 0) return(NULL)
    as.character(x[[1]])
}

.je_chr_vec <- function(x) {
    if (is.null(x) || length(x) == 0) return(character())
    x <- as.character(unlist(x, use.names = FALSE))
    x[nzchar(x)]
}

.je_escape <- function(x) {
    x <- as.character(x)
    x <- gsub("&",  "&amp;",  x, fixed = TRUE)
    x <- gsub("<",  "&lt;",   x, fixed = TRUE)
    x <- gsub(">",  "&gt;",   x, fixed = TRUE)
    x
}

.je_fmt <- function(x, digits = 3) {
    ifelse(is.na(x) | !is.finite(x), "", formatC(x, digits = digits, format = "f"))
}

.je_p <- function(p) {
    ifelse(is.na(p), "",
           ifelse(p < .001, "&lt; .001",
                  sub("^0", "", formatC(p, digits = 3, format = "f"))))
}

.je_table_html <- function(dat) {
    if (is.null(dat) || nrow(dat) == 0)
        return("<p>No results available for the current options.</p>")
    headers <- paste0("<th>", .je_escape(names(dat)), "</th>", collapse = "")
    rows <- apply(dat, 1, function(row) {
        paste0("<tr>", paste0("<td>", .je_escape(row), "</td>", collapse = ""), "</tr>")
    })
    paste0("<table><thead><tr>", headers, "</tr></thead><tbody>",
           paste(rows, collapse = ""), "</tbody></table>")
}

.je_pairwise_html <- function(mat) {
    if (is.null(mat) || length(mat) == 0) return("<p>No pairwise p-values available.</p>")
    dat <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
    names(dat) <- c("Group 1", "Group 2", "p")
    dat <- dat[!is.na(dat$p), , drop = FALSE]
    if (nrow(dat) == 0) return("<p>No pairwise p-values available.</p>")
    dat$p <- .je_p(as.numeric(dat$p))
    .je_table_html(dat)
}

.je_posthoc_adjust <- function(options) {
    if (isTRUE(options$postHocCorrectionTukey))      return("tukey")
    if (isTRUE(options$postHocCorrectionBonferroni)) return("bonferroni")
    if (isTRUE(options$postHocCorrectionHolm))       return("holm")
    if (isTRUE(options$postHocCorrectionSidak))      return("bonferroni")
    "holm"
}

.je_adjust_method <- function(x) {
    x <- as.character(x %||% "none")
    switch(x,
        bonferroni = "bonferroni", holm = "holm",
        sidak = "bonferroni", tukey = "tukey", none = "none", "none"
    )
}

.je_contr_diff <- function(n) {
    if (n < 2) return(matrix(numeric(), nrow = n, ncol = 0))
    mat <- matrix(0, nrow = n, ncol = n - 1)
    for (i in seq_len(n - 1)) { mat[i, i] <- -1; mat[i + 1, i] <- 1 }
    mat
}

`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
}

.je_set_pending_outputs <- function(self, reason) {
    html <- paste0("<p>", .je_escape(reason), "</p>")
    for (nm in c("modelSummary", "descriptivesSection", "effectSizesSection",
                 "assumptionsSection", "contrasts", "orderRestrictions",
                 "postHocSection", "marginalMeansSection", "simpleEffectsSection",
                 "nonparametric", "plots", "savedColumns", "apa", "teaching",
                 "publication", "reproducibility"))
        self$results[[nm]]$setContent(html)
}
