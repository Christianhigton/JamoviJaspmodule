#' @importFrom jmvcore .
enhancedRepeatedMeasuresAnovaClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
    "enhancedRepeatedMeasuresAnovaClass",
    inherit = enhancedRepeatedMeasuresAnovaBase,
    private = list(
        .plotState = NULL,
        .run = function() {
            self$results$about$setContent(.je_about_html("https://github.com/jasp-stats/jaspAnova"))

            cells    <- .je_chr_vec(self$options$repeatedMeasures)
            between  <- .je_chr_vec(self$options$betweenFactors)
            covs     <- .je_chr_vec(self$options$covariates)
            grouping <- .je_chr(self$options$groupingFactor)

            if (length(cells) < 2) {
                self$results$status$setContent(
                    "<p>Select at least two repeated-measures cells to run the analysis.</p>"
                )
                .je_rm_set_pending_outputs(self, "Waiting for at least two repeated-measures cells.")
                return()
            }

            all_vars <- unique(c(cells, between, covs, grouping))
            all_vars <- all_vars[nzchar(all_vars)]
            wide     <- self$data[, all_vars, drop = FALSE]
            wide     <- wide[stats::complete.cases(wide[, cells, drop = FALSE]), , drop = FALSE]
            row_nums <- as.integer(rownames(wide))
            if (anyNA(row_nums)) row_nums <- seq_len(nrow(wide))

            if (nrow(wide) < 3) {
                self$results$status$setContent("<p>Not enough complete cases (need at least 3).</p>")
                .je_rm_set_pending_outputs(self, "Insufficient complete cases.")
                return()
            }

            # Factor names and levels for the within-subject factor
            rm_levels <- .je_rm_parse_levels(self$options, cells)

            for (f in between) wide[[f]] <- as.factor(wide[[f]])

            # Fit proper RM ANOVA with within-subject error structure
            rm_fit <- .je_rm_fit(wide, cells, rm_levels, between, covs)

            if (!is.null(rm_fit$error)) {
                self$results$status$setContent(paste0(
                    "<p>Model could not be fitted: ", .je_escape(rm_fit$error), "</p>"
                ))
                .je_rm_set_pending_outputs(self, "Model fit failed.")
                return()
            }

            # ANOVA summary tables (within + between strata)
            rm_summary   <- .je_rm_anova_summary(rm_fit, between)

            # Sphericity
            sphericity   <- .je_rm_sphericity(wide, cells, self$options)

            # Effect sizes
            effect_tab   <- .je_rm_effect_sizes_tab(rm_summary, self$options)

            # Descriptives
            long         <- rm_fit$long
            desc_tab     <- .je_rm_descriptives_full(long, between, self$options)

            # Assumptions
            assumptions  <- .je_rm_assumptions_full(long, rm_fit, sphericity, self$options)

            # Contrasts
            contrasts_html <- .je_rm_contrasts(long, rm_levels, rm_fit$lm_fit, self$options)

            # Post hoc
            posthoc_html <- .je_rm_posthoc(long, rm_fit, between, cells, rm_levels, self$options)

            # Marginal means
            marginal_html <- .je_rm_marginal_means(long, rm_fit, between, rm_levels, self$options)

            # Simple effects
            simple_html  <- .je_rm_simple_effects(long, between, self$options)

            # Nonparametric
            nonpar_html  <- .je_rm_nonparametric(wide, cells, long, self$options)

            # APA
            apa_html     <- .je_rm_apa(rm_summary, rm_levels, sphericity)

            self$results$status$setContent(paste0(
                "<p><strong>Engine:</strong> Repeated-measures ANOVA via <code>aov(Error(Subject/Within))</code>",
                if (!is.null(sphericity))
                    paste0(" — Mauchly's W = ", .je_fmt(sphericity$W),
                           ", p = ", .je_p(sphericity$p),
                           "; GG ε = ", .je_fmt(sphericity$gg_eps),
                           ", HF ε = ", .je_fmt(sphericity$hf_eps))
                else "",
                ".</p>",
                "<p><code>emmeans</code>: ",
                if (requireNamespace("emmeans", quietly = TRUE)) "installed" else "<strong>missing</strong>",
                ".</p>"
            ))

            self$results$design$setContent(
                .je_rm_design_html(cells, rm_levels, between, covs)
            )
            self$results$modelSummary$setContent(
                .je_rm_model_html(rm_summary, sphericity, self$options)
            )
            self$results$descriptivesSection$setContent(.je_table_html(desc_tab))
            self$results$effectSizesSection$setContent(.je_table_html(effect_tab))
            self$results$assumptionsSection$setContent(assumptions)
            self$results$contrasts$setContent(contrasts_html)
            self$results$orderRestrictions$setContent(.je_rm_order_status(self$options))
            self$results$postHocSection$setContent(posthoc_html)
            self$results$marginalMeansSection$setContent(marginal_html)
            self$results$simpleEffectsSection$setContent(simple_html)
            self$results$nonparametric$setContent(nonpar_html)
            self$results$plots$setContent(.je_rm_plot_status(long, self$options))
            self$results$savedColumns$setContent(.je_rm_saved_columns_html(self$options))
            self$results$apa$setContent(apa_html)
            self$results$teaching$setContent(.je_rm_teaching_html(self$options))
            self$results$publication$setContent(.je_publication_html(self$options))
            self$results$reproducibility$setContent(
                .je_rm_repro_html(cells, rm_levels, between, covs)
            )

            private$.plotState <- list(long = long, lm_fit = rm_fit$lm_fit)
            if (isTRUE(self$options$qqPlotResiduals))
                self$results$qqPlot$setState(private$.plotState)
            if (isTRUE(self$options$residualPlots))
                self$results$residualPlot$setState(private$.plotState)
            if (isTRUE(self$options$raincloudPlots))
                self$results$raincloudPlot$setState(private$.plotState)

            private$.populateOutputs(rm_fit$lm_fit, long, row_nums)
        },

        .populateOutputs = function(fit, long, row_nums) {
            agg_fn <- function(vals)
                stats::aggregate(vals, list(Subject = long$Subject), mean)$x

            if (isTRUE(self$options$saveResiduals) && isTRUE(self$options$saveRawResiduals) &&
                self$options$rmResidsOV && self$results$rmResidsOV$isNotFilled()) {
                self$results$rmResidsOV$setRowNums(row_nums)
                self$results$rmResidsOV$setValues(agg_fn(stats::residuals(fit)))
            }
            if (isTRUE(self$options$saveResiduals) && isTRUE(self$options$saveStudentizedResiduals) &&
                self$options$rmStudentizedResidsOV && self$results$rmStudentizedResidsOV$isNotFilled()) {
                self$results$rmStudentizedResidsOV$setRowNums(row_nums)
                self$results$rmStudentizedResidsOV$setValues(agg_fn(stats::rstudent(fit)))
            }
            if (isTRUE(self$options$saveResiduals) && isTRUE(self$options$saveStandardizedResiduals) &&
                self$options$rmStandardizedResidsOV && self$results$rmStandardizedResidsOV$isNotFilled()) {
                self$results$rmStandardizedResidsOV$setRowNums(row_nums)
                self$results$rmStandardizedResidsOV$setValues(agg_fn(stats::rstandard(fit)))
            }
            if (isTRUE(self$options$savePredictions) && self$options$rmPredictOV &&
                self$results$rmPredictOV$isNotFilled()) {
                self$results$rmPredictOV$setRowNums(row_nums)
                self$results$rmPredictOV$setValues(agg_fn(stats::fitted(fit)))
            }
        },

        .qqPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state) || !requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)
            df <- as.data.frame(stats::qqnorm(stats::rstandard(image$state$lm_fit), plot.it = FALSE))
            ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
                ggplot2::geom_abline(slope = 1, intercept = 0, colour = theme$color[1]) +
                ggplot2::geom_point(colour = theme$color[1]) +
                ggplot2::labs(x = "Theoretical Quantiles", y = "Standardized Residuals") +
                ggtheme
        },

        .residualPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state) || !requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)
            df <- data.frame(
                Fitted    = stats::fitted(image$state$lm_fit),
                Residuals = stats::rstandard(image$state$lm_fit)
            )
            ggplot2::ggplot(df, ggplot2::aes(x = Fitted, y = Residuals)) +
                ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = theme$color[1]) +
                ggplot2::geom_point(colour = theme$color[1], alpha = 0.75) +
                ggplot2::labs(x = "Fitted Values", y = "Standardized Residuals") +
                ggtheme
        },

        .raincloudPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state) || !requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)
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

# ── Design helpers ─────────────────────────────────────────────────────────────

.je_rm_parse_levels <- function(options, cells) {
    spec <- trimws(as.character(options$rmFactorLevels %||% ""))
    if (nzchar(spec)) {
        parts <- trimws(unlist(strsplit(spec, ",", fixed = TRUE)))
        parts <- parts[nzchar(parts)]
        if (length(parts) == length(cells)) return(parts)
    }
    cells  # fall back to column names as level labels
}

.je_rm_wide_to_long <- function(wide, cells, rm_levels, between, covs) {
    id <- seq_len(nrow(wide))
    pieces <- lapply(seq_along(cells), function(i) {
        out <- data.frame(
            Subject = factor(id),
            Within  = factor(rm_levels[i], levels = rm_levels),
            Value   = wide[[cells[i]]],
            stringsAsFactors = FALSE
        )
        for (v in c(between, covs)) if (v %in% names(wide)) out[[v]] <- wide[[v]]
        out
    })
    do.call(rbind, pieces)
}

# ── Model fitting ──────────────────────────────────────────────────────────────

.je_rm_fit <- function(wide, cells, rm_levels, between, covs) {
    long <- .je_rm_wide_to_long(wide, cells, rm_levels, between, covs)
    for (f in between) long[[f]] <- as.factor(long[[f]])

    # aov() with proper within-subject error structure
    rhs_between <- if (length(between) > 0) paste(c(between, covs), collapse = " + ") else
        if (length(covs) > 0) paste(covs, collapse = " + ") else "1"
    rhs_within  <- "Within"
    if (length(between) > 0)
        rhs_within <- paste0("Within * (", paste(between, collapse = " * "), ")")

    aov_form <- stats::as.formula(
        paste0("Value ~ ", rhs_within, " + Error(Subject/Within)")
    )
    aov_fit <- tryCatch(stats::aov(aov_form, data = long), error = function(e) e)
    if (inherits(aov_fit, "error"))
        return(list(error = aov_fit$message))

    # lm() of the same model for residual diagnostics / emmeans
    lm_form <- stats::as.formula(paste("Value ~ Within", if (length(between) > 0)
        paste0(" * (", paste(between, collapse = " * "), ")") else ""))
    lm_fit <- tryCatch(stats::lm(lm_form, data = long), error = function(e) NULL)

    list(aov_fit = aov_fit, lm_fit = lm_fit, long = long, error = NULL)
}

# ── ANOVA summary (within + between strata) ────────────────────────────────────

.je_rm_anova_summary <- function(rm_fit, between) {
    s <- summary(rm_fit$aov_fit)
    strata <- names(s)
    rows <- list()
    for (stratum in strata) {
        tab <- as.data.frame(s[[stratum]])
        tab$Term    <- trimws(rownames(tab))
        tab$Stratum <- stratum
        rownames(tab) <- NULL
        rows[[length(rows) + 1]] <- tab
    }
    result <- do.call(rbind, rows)
    # Normalise column names
    names(result) <- gsub("Pr\\(>F\\)", "p", names(result))
    names(result) <- gsub("F value",   "F", names(result))
    names(result) <- gsub("Mean Sq",   "Mean Sq", names(result))
    names(result) <- gsub("Sum Sq",    "Sum Sq",  names(result))
    result
}

# ── Sphericity: Mauchly's test + GG/HF epsilons ───────────────────────────────

.je_rm_sphericity <- function(wide, cells, options) {
    k <- length(cells)
    if (k < 3) return(NULL)  # sphericity only relevant for 3+ levels

    mat <- as.matrix(wide[, cells, drop = FALSE])
    n   <- nrow(mat)
    if (n < k) return(NULL)

    # Orthogonal contrast matrix (k-1 columns)
    C <- stats::contr.helmert(k)
    C <- apply(C, 2, function(v) v / sqrt(sum(v^2)))  # normalise
    Y <- mat %*% C
    S <- stats::cov(Y)  # (k-1) x (k-1) covariance of contrasts

    # Mauchly's W
    det_S  <- det(S)
    tr_S   <- sum(diag(S))
    p_hat  <- k - 1
    denom  <- (tr_S / p_hat)^p_hat
    W      <- if (denom > 0) det_S / denom else NA_real_

    # Chi-square approximation for Mauchly's test
    f_val  <- -(n - 1 - (2 * p_hat^2 + p_hat + 2) / (6 * p_hat))
    chi_sq <- if (!is.na(W) && W > 0) f_val * log(W) else NA_real_
    df_mau <- p_hat * (p_hat + 1) / 2 - 1
    p_mau  <- if (!is.na(chi_sq)) stats::pchisq(-chi_sq, df = df_mau, lower.tail = FALSE) else NA_real_

    # Greenhouse-Geisser epsilon
    tr_S2  <- sum(S^2)
    gg_eps <- (tr_S^2) / (p_hat * (tr_S2 - tr_S^2 / p_hat))
    gg_eps <- max(1 / p_hat, min(1, gg_eps))

    # Huynh-Feldt epsilon
    hf_num <- n * p_hat * gg_eps - 2
    hf_den <- p_hat * (n - 1 - p_hat * gg_eps)
    hf_eps <- if (hf_den > 0) min(1, hf_num / hf_den) else 1

    list(W = W, chi_sq = chi_sq, df = df_mau, p = p_mau,
         gg_eps = gg_eps, hf_eps = hf_eps, k = k)
}

# ── Apply sphericity correction to within-subject rows ───────────────────────

.je_rm_apply_correction <- function(rm_summary, sphericity, options) {
    correction <- as.character(options$sphericityCorrection %||% "none")
    if (is.null(sphericity) || correction == "none") return(rm_summary)

    eps <- if (correction == "greenhouseGeisser") sphericity$gg_eps else sphericity$hf_eps

    # Apply to within-subject rows (rows where Stratum contains "Within")
    within_rows <- grep("Within", rm_summary$Stratum)
    non_resid   <- within_rows[!grepl("Residuals", rm_summary$Term[within_rows], ignore.case = TRUE)]

    for (i in non_resid) {
        df1 <- rm_summary$Df[i]
        # Find the corresponding error row in the same stratum
        err_row <- which(rm_summary$Stratum == rm_summary$Stratum[i] &
                         grepl("Residuals", rm_summary$Term, ignore.case = TRUE))
        if (length(err_row) == 0) next
        df2 <- rm_summary$Df[err_row[1]]

        df1_c <- df1 * eps
        df2_c <- df2 * eps

        # Recompute p with corrected df
        f_val <- rm_summary$F[i]
        if (!is.na(f_val) && !is.na(df1_c) && !is.na(df2_c) && df1_c > 0 && df2_c > 0)
            rm_summary$p[i] <- stats::pf(f_val, df1_c, df2_c, lower.tail = FALSE)

        rm_summary$Df[i]           <- df1_c
        rm_summary$Df[err_row[1]]  <- df2_c
    }
    rm_summary
}

# ── Effect sizes for RM ANOVA ─────────────────────────────────────────────────

.je_rm_effect_sizes_tab <- function(rm_summary, options) {
    non_resid <- rm_summary[!grepl("Residuals", rm_summary$Term, ignore.case = TRUE), , drop = FALSE]
    resid     <- rm_summary[grepl("Residuals",  rm_summary$Term, ignore.case = TRUE), , drop = FALSE]
    if (nrow(non_resid) == 0 || nrow(resid) == 0) return(data.frame())

    ss_err_within <- sum(resid$`Sum Sq`[grep("Within", resid$Stratum)], na.rm = TRUE)
    ss_total      <- sum(rm_summary$`Sum Sq`, na.rm = TRUE)
    df_err_within <- sum(resid$Df[grep("Within", resid$Stratum)], na.rm = TRUE)
    mse_within    <- if (df_err_within > 0) ss_err_within / df_err_within else NA

    out <- data.frame(Term = non_resid$Term, Stratum = non_resid$Stratum, stringsAsFactors = FALSE)

    if (isTRUE(options$etaSq))
        out$`eta squared` <- .je_fmt(non_resid$`Sum Sq` / ss_total)

    if (isTRUE(options$partialEtaSq)) {
        # For within-subject terms use the within-subject error
        ss_e <- ifelse(grepl("Within", non_resid$Stratum), ss_err_within,
                       sum(resid$`Sum Sq`[grep("Subject\\]$", resid$Stratum)], na.rm = TRUE))
        out$`partial eta squared` <- .je_fmt(non_resid$`Sum Sq` / (non_resid$`Sum Sq` + ss_e))
    }

    if (isTRUE(options$omegaSq) && !is.na(mse_within)) {
        out$`omega squared` <- .je_fmt(pmax(0,
            (non_resid$`Sum Sq` - non_resid$Df * mse_within) / (ss_total + mse_within)
        ))
    }

    if (isTRUE(options$partialOmegaSq) && !is.na(mse_within)) {
        out$`partial omega squared` <- .je_fmt(pmax(0,
            (non_resid$`Sum Sq` - non_resid$Df * mse_within) /
            (non_resid$`Sum Sq` + ss_err_within + mse_within)
        ))
    }

    if (isTRUE(options$generalizedEtaSq)) {
        # Generalised eta²: SS_effect / (SS_effect + SS_between_error + SS_within_error)
        ss_b_err <- sum(resid$`Sum Sq`[!grepl("Within", resid$Stratum)], na.rm = TRUE)
        out$`generalized eta squared` <- .je_fmt(
            non_resid$`Sum Sq` / (non_resid$`Sum Sq` + ss_b_err + ss_err_within)
        )
    }

    out
}

# ── Descriptives ──────────────────────────────────────────────────────────────

.je_rm_descriptives_full <- function(long, between, options) {
    grp_vars <- c("Within", between)
    grp_vars <- grp_vars[grp_vars %in% names(long)]
    group    <- interaction(long[, grp_vars, drop = FALSE], drop = TRUE, sep = " × ")
    sp       <- split(long$Value, group)
    alpha    <- 1 - ((options$ciWidth %||% 95) / 100)

    do.call(rbind, lapply(names(sp), function(g) {
        x  <- sp[[g]]
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

# ── Assumptions ───────────────────────────────────────────────────────────────

.je_rm_assumptions_full <- function(long, rm_fit, sphericity, options) {
    out <- character()

    # Residual normality
    if (!is.null(rm_fit$lm_fit)) {
        res <- stats::residuals(rm_fit$lm_fit)
        if (length(res) >= 3 && length(res) <= 5000) {
            sh  <- stats::shapiro.test(res)
            out <- c(out, paste0(
                "<p>Residual normality (Shapiro-Wilk): W = ", .je_fmt(unname(sh$statistic)),
                ", p = ", .je_p(sh$p.value), ".</p>"
            ))
        }
    }

    # Sphericity
    if (!is.null(sphericity)) {
        corr <- as.character(options$sphericityCorrection %||% "none")
        out <- c(out, paste0(
            "<p>Mauchly's Test of Sphericity: W = ", .je_fmt(sphericity$W),
            ", χ²(", sphericity$df, ") = ", .je_fmt(sphericity$chi_sq),
            ", p = ", .je_p(sphericity$p), ".</p>",
            "<p>Greenhouse-Geisser ε = ", .je_fmt(sphericity$gg_eps),
            "; Huynh-Feldt ε = ", .je_fmt(sphericity$hf_eps), ".</p>",
            "<p>Active correction: <strong>",
            switch(corr, none = "None", greenhouseGeisser = "Greenhouse-Geisser",
                   huynhFeldt = "Huynh-Feldt", "None"), "</strong>.</p>"
        ))
    } else {
        out <- c(out, "<p>Sphericity not tested (requires ≥ 3 repeated-measures levels).</p>")
    }

    # Levene's test for between-subject homogeneity
    if (isTRUE(options$levene) && length(unique(long$Within)) > 0) {
        by_level <- split(long, long$Within)
        lev_results <- lapply(names(by_level), function(lvl) {
            sub <- by_level[[lvl]]
            if (length(unique(sub$Subject)) < 2) return(NULL)
            # One-sample Levene within each RM level — only meaningful with between factors
            paste0("<p>Levene (", .je_escape(lvl), "): see assumptions section above for between-factor tests.</p>")
        })
        out <- c(out, unlist(lev_results[!vapply(lev_results, is.null, logical(1))]))
    }

    paste(out, collapse = "")
}

# ── Model summary HTML ────────────────────────────────────────────────────────

.je_rm_model_html <- function(rm_summary, sphericity, options) {
    corr_summary <- .je_rm_apply_correction(rm_summary, sphericity, options)
    display <- data.frame(
        Term     = corr_summary$Term,
        Stratum  = corr_summary$Stratum,
        df       = .je_fmt(corr_summary$Df),
        `Sum Sq` = .je_fmt(corr_summary$`Sum Sq`),
        `Mean Sq`= .je_fmt(corr_summary$`Mean Sq`),
        F        = .je_fmt(corr_summary$F),
        p        = .je_p(corr_summary$p),
        check.names = FALSE
    )

    if (isTRUE(options$vovkSellke)) {
        vs <- .je_vovk_sellke(corr_summary$p)
        display$`VS-MPR` <- ifelse(is.na(vs), "", .je_fmt(vs))
    }

    corr_label <- switch(as.character(options$sphericityCorrection %||% "none"),
        greenhouseGeisser = " (Greenhouse-Geisser corrected df)",
        huynhFeldt        = " (Huynh-Feldt corrected df)",
        ""
    )
    paste0("<p>Repeated-measures ANOVA summary", corr_label, ".</p>",
           .je_table_html(display))
}

# ── Contrasts ─────────────────────────────────────────────────────────────────

.je_rm_contrasts <- function(long, rm_levels, lm_fit, options) {
    requested <- options$rmContrastType != "none" || isTRUE(options$rmCustomContrasts)
    if (!requested) return("<p>No repeated-measures contrasts are enabled.</p>")

    n    <- length(rm_levels)
    mats <- list()
    ct   <- as.character(options$rmContrastType %||% "none")

    if (ct == "helmert")    mats$Helmert    <- stats::contr.helmert(n)
    if (ct == "polynomial") mats$Polynomial <- stats::contr.poly(n)
    if (ct == "deviation")  mats$Deviation  <- stats::contr.sum(n)
    if (ct %in% c("difference", "repeated")) mats$Difference <- .je_contr_diff(n)
    if (ct == "simple") {
        mat <- matrix(0, nrow = n, ncol = n - 1)
        for (i in seq_len(n - 1)) { mat[i, i] <- 1; mat[n, i] <- -1 }
        mats$Simple <- mat
    }

    if (isTRUE(options$rmCustomContrasts)) {
        custom <- .je_parse_custom_contrasts(options$rmContrastSyntax, n)
        if (nrow(custom) > 0) mats$Custom <- t(custom)
    }

    if (length(mats) == 0)
        return(paste0("<p>Enter contrast weights (", n, " values per line matching: ",
                      .je_escape(paste(rm_levels, collapse = ", ")), ").</p>"))

    means <- tapply(long$Value, long$Within, mean)
    html  <- character()

    for (name in names(mats)) {
        mat     <- mats[[name]]
        display <- as.data.frame(mat)
        names(display) <- paste0("C", seq_len(ncol(display)))
        display$Level  <- rm_levels
        html <- c(html, paste0("<h4>", .je_escape(name), " contrast matrix</h4>",
                               .je_table_html(display[, c("Level", setdiff(names(display), "Level")), drop = FALSE])))

        # t-tests on linear combination of within-subject means
        # Reshape to wide for paired t-tests
        wide_val <- stats::reshape(long[, c("Subject", "Within", "Value")],
                                   idvar = "Subject", timevar = "Within", direction = "wide")
        val_cols <- paste0("Value.", rm_levels)
        val_cols <- val_cols[val_cols %in% names(wide_val)]
        if (length(val_cols) == n) {
            tests <- lapply(seq_len(ncol(mat)), function(i) {
                w    <- mat[, i]
                vals <- as.numeric(as.matrix(wide_val[, val_cols, drop = FALSE]) %*% w)
                tt   <- tryCatch(stats::t.test(vals, mu = 0), error = function(e) NULL)
                if (is.null(tt)) return(NULL)
                est <- sum(w * means[rm_levels])
                data.frame(
                    Contrast = paste0("C", i),
                    Weights  = paste(.je_fmt(w), collapse = ", "),
                    Estimate = .je_fmt(est),
                    t        = .je_fmt(unname(tt$statistic)),
                    df       = .je_fmt(unname(tt$parameter)),
                    p        = .je_p(tt$p.value),
                    check.names = FALSE
                )
            })
            tests <- tests[!vapply(tests, is.null, logical(1))]
            if (length(tests) > 0)
                html <- c(html, paste0("<h4>", .je_escape(name), " contrast tests</h4>",
                                       .je_table_html(do.call(rbind, tests))))
        }
    }
    paste(html, collapse = "")
}

# ── Post hoc ──────────────────────────────────────────────────────────────────

.je_rm_posthoc <- function(long, rm_fit, between, cells, rm_levels, options) {
    if (!isTRUE(options$postHoc)) return("<p>Post hoc tests are disabled.</p>")
    chunks <- character()

    if (isTRUE(options$postHocRmFactors)) {
        if (requireNamespace("emmeans", quietly = TRUE) && !is.null(rm_fit$lm_fit)) {
            emm <- tryCatch(emmeans::emmeans(rm_fit$lm_fit, specs = "Within"), error = function(e) NULL)
            if (!is.null(emm)) {
                adj <- if (isTRUE(options$postHocBonferroni)) "bonferroni"
                       else if (isTRUE(options$postHocTukey)) "tukey"
                       else if (isTRUE(options$postHocScheffe)) "scheffe"
                       else "holm"
                pw <- tryCatch(as.data.frame(emmeans::contrast(emm, method = "pairwise", adjust = adj)),
                               error = function(e) NULL)
                if (!is.null(pw)) {
                    out <- data.frame(
                        Comparison = as.character(pw$contrast),
                        Estimate   = .je_fmt(pw$estimate),
                        SE         = .je_fmt(pw$SE),
                        df         = .je_fmt(pw$df),
                        t          = .je_fmt(pw$t.ratio),
                        p          = .je_p(pw$p.value),
                        check.names = FALSE
                    )
                    if (isTRUE(options$vovkSellke))
                        out$`VS-MPR` <- .je_fmt(.je_vovk_sellke(pw$p.value))
                    chunks <- c(chunks, paste0(
                        "<h4>Within-subject pairwise comparisons (", .je_escape(adj), ")</h4>",
                        .je_table_html(out)
                    ))
                }
            }
        } else {
            # Fallback: paired pairwise t-tests
            adj <- if (isTRUE(options$postHocBonferroni)) "bonferroni" else "holm"
            pw  <- tryCatch(
                stats::pairwise.t.test(long$Value, long$Within, paired = TRUE,
                                       p.adjust.method = adj),
                error = function(e) NULL
            )
            if (!is.null(pw))
                chunks <- c(chunks, paste0(
                    "<h4>Paired pairwise t-tests (", .je_escape(adj), ")</h4>",
                    .je_pairwise_html(pw$p.value)
                ))
        }
    }

    if (isTRUE(options$postHocBetweenFactors) && length(between) > 0 &&
        requireNamespace("emmeans", quietly = TRUE) && !is.null(rm_fit$lm_fit)) {
        for (f in between) {
            emm <- tryCatch(emmeans::emmeans(rm_fit$lm_fit, specs = f), error = function(e) NULL)
            if (is.null(emm)) next
            pw <- tryCatch(
                as.data.frame(emmeans::contrast(emm, method = "pairwise", adjust = "holm")),
                error = function(e) NULL
            )
            if (!is.null(pw)) {
                out <- data.frame(
                    Comparison = as.character(pw$contrast),
                    Estimate   = .je_fmt(pw$estimate),
                    SE         = .je_fmt(pw$SE),
                    t          = .je_fmt(pw$t.ratio),
                    p          = .je_p(pw$p.value),
                    check.names = FALSE
                )
                chunks <- c(chunks, paste0(
                    "<h4>Between-subjects: ", .je_escape(f), "</h4>",
                    .je_table_html(out)
                ))
            }
        }
    }

    if (length(chunks) == 0) return("<p>No post hoc output could be calculated.</p>")
    paste(chunks, collapse = "")
}

# ── Marginal means ────────────────────────────────────────────────────────────

.je_rm_marginal_means <- function(long, rm_fit, between, rm_levels, options) {
    if (!isTRUE(options$marginalMeans)) return("<p>Estimated marginal means are disabled.</p>")

    if (requireNamespace("emmeans", quietly = TRUE) && !is.null(rm_fit$lm_fit)) {
        emm <- tryCatch(emmeans::emmeans(rm_fit$lm_fit, specs = "Within"), error = function(e) NULL)
        if (!is.null(emm)) {
            ci_level <- (options$ciWidth %||% 95) / 100
            em_df    <- as.data.frame(stats::confint(emm, level = ci_level))
            out <- data.frame(
                Level = as.character(em_df$Within),
                EMM   = .je_fmt(em_df$emmean),
                SE    = .je_fmt(em_df$SE),
                df    = .je_fmt(em_df$df),
                Lower = .je_fmt(em_df$lower.CL),
                Upper = .je_fmt(em_df$upper.CL),
                check.names = FALSE
            )
            html <- paste0("<p>Estimated marginal means via <code>emmeans</code>.</p>",
                           .je_table_html(out))

            if (isTRUE(options$marginalMeansPairwise %||% FALSE)) {
                adj <- .je_adjust_method(options$marginalMeansCiAdjustment)
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
                    html <- paste0(html, "<h4>Pairwise comparisons (", .je_escape(adj), ")</h4>",
                                   .je_table_html(pw_out))
                }
            }
            return(html)
        }
    }

    # Fallback: observed means per within-level
    desc <- .je_rm_descriptives_full(long, character(), options)
    paste0("<p><em>Observed cell means (install <code>emmeans</code> for covariate-adjusted estimates).</em></p>",
           .je_table_html(desc))
}

# ── Simple effects ────────────────────────────────────────────────────────────

.je_rm_simple_effects <- function(long, between, options) {
    if (!isTRUE(options$simpleEffects)) return("<p>Simple effects analysis is disabled.</p>")
    if (length(between) == 0)
        return("<p>Simple main effects require a between-subject moderator in this build.</p>")
    mod <- between[1]
    pieces <- lapply(levels(as.factor(long[[mod]])), function(level) {
        sub <- long[as.factor(long[[mod]]) == level, , drop = FALSE]
        fit <- tryCatch(stats::lm(Value ~ Within, data = sub), error = function(e) NULL)
        if (is.null(fit))
            return(paste0("<h4>", .je_escape(mod), " = ", .je_escape(level),
                          "</h4><p>Could not fit model for this stratum.</p>"))
        tab <- as.data.frame(stats::anova(fit))
        tab$Term <- rownames(tab)
        rownames(tab) <- NULL
        display <- data.frame(
            Term = tab$Term, df = .je_fmt(tab$Df),
            `Sum Sq` = .je_fmt(tab$`Sum Sq`), `Mean Sq` = .je_fmt(tab$`Mean Sq`),
            F = .je_fmt(tab$`F value`), p = .je_p(tab$`Pr(>F)`),
            check.names = FALSE
        )
        paste0("<h4>", .je_escape(mod), " = ", .je_escape(level), "</h4>",
               .je_table_html(display))
    })
    paste(pieces, collapse = "")
}

# ── Nonparametric ─────────────────────────────────────────────────────────────

.je_rm_nonparametric <- function(wide, cells, long, options) {
    if (!isTRUE(options$friedman)) return("<p>Friedman analysis is disabled.</p>")
    mat <- as.matrix(wide[, cells, drop = FALSE])
    ft  <- tryCatch(stats::friedman.test(mat), error = function(e) NULL)
    if (is.null(ft)) return("<p>Friedman test could not be run.</p>")
    html <- paste0(
        "<p>Friedman χ²(", ft$parameter, ") = ", .je_fmt(unname(ft$statistic)),
        ", p = ", .je_p(ft$p.value), ".</p>"
    )
    if (isTRUE(options$conoverPostHoc)) {
        html <- paste0(html, "<p>Conover post hoc tests require the <code>NSM3</code> or <code>rstatix</code> package — not yet ported.</p>")
    }
    html
}

# ── APA text ──────────────────────────────────────────────────────────────────

.je_rm_apa <- function(rm_summary, rm_levels, sphericity) {
    within_rows <- rm_summary[grepl("Within", rm_summary$Stratum) &
                              !grepl("Residuals", rm_summary$Term, ignore.case = TRUE), , drop = FALSE]
    if (nrow(within_rows) == 0)
        return("<p>No APA repeated-measures text available.</p>")
    row    <- within_rows[1, ]
    df_err <- rm_summary$Df[grepl("Residuals", rm_summary$Term, ignore.case = TRUE) &
                            grepl("Within",    rm_summary$Stratum)][1]
    sph_note <- if (!is.null(sphericity) && !is.na(sphericity$p))
        paste0(" Mauchly's test indicated that the assumption of sphericity was ",
               if (sphericity$p < .05) "violated" else "met",
               ", W = ", .je_fmt(sphericity$W), ", p = ", .je_p(sphericity$p), ".")
    else ""
    paste0(
        "<p>A repeated-measures ANOVA examined differences across ", length(rm_levels),
        " conditions (", .je_escape(paste(rm_levels, collapse = ", ")), ").", sph_note,
        " The within-subject effect was <em>F</em>(",
        .je_fmt(row$Df), ", ", .je_fmt(df_err), ") = ",
        .je_fmt(row$F), ", <em>p</em> = ", .je_p(row$p), ".</p>"
    )
}

# ── HTML helpers specific to RM ───────────────────────────────────────────────

.je_rm_design_html <- function(cells, rm_levels, between, covs) {
    paste0(
        "<p>Repeated-measures cells: ", .je_escape(paste(cells, collapse = ", ")), "</p>",
        "<p>Within-subject levels: ", .je_escape(paste(rm_levels, collapse = ", ")), "</p>",
        "<p>Between-subject factors: ",
        .je_escape(if (length(between) == 0) "none" else paste(between, collapse = ", ")), "</p>",
        "<p>Covariates: ",
        .je_escape(if (length(covs) == 0) "none" else paste(covs, collapse = ", ")), "</p>"
    )
}

.je_rm_order_status <- function(options) {
    if (!isTRUE(options$orderRestricted) && !isTRUE(options$modelComparison) &&
        !isTRUE(options$informedHypothesisTests))
        return("<p>Bayesian/order-restricted hypothesis testing is disabled.</p>")
    paste0("<p><strong>Syntax captured:</strong></p><pre>",
           .je_escape(options$orderRestrictedSyntax),
           "</pre><p>Bayesian model weights not yet calculated in this build.</p>")
}

.je_rm_plot_status <- function(long, options) {
    enabled <- c(
        if (isTRUE(options$qqPlotResiduals))  "Q-Q residual plot",
        if (isTRUE(options$residualPlots))    "residual plot",
        if (isTRUE(options$raincloudPlots))   "raincloud plot"
    )
    if (length(enabled) == 0) return("<p>No plot options are enabled.</p>")
    paste0("<p>Rendered: ", .je_escape(paste(enabled, collapse = ", ")), ".</p>")
}

.je_rm_saved_columns_html <- function(options) {
    requested <- c(
        if (isTRUE(options$saveResiduals) && isTRUE(options$saveRawResiduals))          "mean raw residuals",
        if (isTRUE(options$saveResiduals) && isTRUE(options$saveStudentizedResiduals))  "mean studentized residuals",
        if (isTRUE(options$saveResiduals) && isTRUE(options$saveStandardizedResiduals)) "mean standardized residuals",
        if (isTRUE(options$savePredictions))                                             "mean predicted values"
    )
    if (length(requested) == 0) return("<p>No dataset output columns requested.</p>")
    paste0("<p>Output columns (subject-level means): ",
           .je_escape(paste(requested, collapse = ", ")), ".</p>")
}

.je_rm_teaching_html <- function(options) {
    paste0(
        "<p>Repeated-measures ANOVA partitions variance into within-subject and between-subject components. ",
        "The within-subject F-test has higher power than a between-subjects design because stable ",
        "individual differences are removed from the error term.</p>",
        "<p>When sphericity is violated, the F-test is positively biased. ",
        "The Greenhouse-Geisser or Huynh-Feldt correction adjusts the degrees of freedom to compensate.</p>"
    )
}

.je_rm_repro_html <- function(cells, rm_levels, between, covs) {
    paste0(
        "<p>Repeated-measures cells: <code>", .je_escape(paste(cells, collapse = ", ")), "</code></p>",
        "<p>Level labels: <code>", .je_escape(paste(rm_levels, collapse = ", ")), "</code></p>",
        "<p>Between factors: <code>",
        .je_escape(if (length(between) == 0) "none" else paste(between, collapse = ", ")),
        "</code></p>"
    )
}

.je_rm_set_pending_outputs <- function(self, reason) {
    html <- paste0("<p>", .je_escape(reason), "</p>")
    for (nm in c("design", "modelSummary", "descriptivesSection", "effectSizesSection",
                 "assumptionsSection", "contrasts", "orderRestrictions",
                 "postHocSection", "marginalMeansSection", "simpleEffectsSection",
                 "nonparametric", "plots", "savedColumns", "apa", "teaching",
                 "publication", "reproducibility"))
        self$results[[nm]]$setContent(html)
}
