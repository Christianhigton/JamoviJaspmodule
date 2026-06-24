#' @importFrom jmvcore .
enhancedRepeatedMeasuresAnovaClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
    "enhancedRepeatedMeasuresAnovaClass",
    inherit = enhancedRepeatedMeasuresAnovaBase,
    private = list(
        .plotState = NULL,
        .run = function() {
            self$results$about$setContent(.je_about_html("https://github.com/jasp-stats/jaspAnova"))

            cells    <- .je_chr_vec(self$options$repeatedMeasuresCells)
            between  <- .je_chr_vec(self$options$betweenSubjectFactors)
            covs     <- .je_chr_vec(self$options$covariates)
            grouping <- .je_chr(self$options$friedmanBetweenFactor)

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

            # Parse factor structure (single or multi-factor)
            factor_list <- .je_rm_parse_factor_spec(self$options, cells)
            factor_names <- vapply(factor_list, function(f) f$name, character(1))
            rm_levels    <- factor_list[[1]]$levels  # backward-compat for single-factor code
            within_model_terms <- .je_rm_parse_model_terms(
                self$options$withinModelTerms, factor_names, full_factorial = TRUE
            )
            between_model_terms <- .je_rm_parse_model_terms(
                self$options$betweenModelTerms, c(between, covs), full_factorial = FALSE
            )

            # Validate cell count
            n_cells_needed <- prod(vapply(factor_list, function(f) length(f$levels), integer(1)))
            if (n_cells_needed != length(cells)) {
                self$results$status$setContent(paste0(
                    "<p>Factor specification requires ", n_cells_needed,
                    " cells but ", length(cells), " are assigned. ",
                    "Check the repeated-measures factor specification.</p>"
                ))
                .je_rm_set_pending_outputs(self, "Cell count mismatch.")
                return()
            }

            for (f in between) wide[[f]] <- as.factor(wide[[f]])

            # Fit RM ANOVA with correct within-subject error structure
            rm_fit <- .je_rm_fit_multi(
                wide, cells, factor_list, between, covs,
                within_model_terms, between_model_terms, self$options
            )

            if (!is.null(rm_fit$error)) {
                self$results$status$setContent(paste0(
                    "<p>Model could not be fitted: ", .je_escape(rm_fit$error), "</p>"
                ))
                .je_rm_set_pending_outputs(self, "Model fit failed.")
                return()
            }

            # ANOVA summary tables (within + between strata)
            rm_summary   <- .je_rm_anova_summary(rm_fit, between)

            # Sphericity — per-effect for multi-factor designs
            sphericity_list <- .je_rm_sphericity_multi(wide, cells, factor_list, self$options)
            # For status bar, summarise the first within-subject effect
            sphericity <- if (length(sphericity_list) > 0) sphericity_list[[1]]$sphericity else NULL

            # Effect sizes
            effect_tab   <- .je_rm_effect_sizes_tab(rm_summary, self$options)

            # Descriptives
            long         <- rm_fit$long
            desc_tab     <- .je_rm_descriptives_full(long, between, self$options)

            # Assumptions
            assumptions  <- .je_rm_assumptions_full(long, rm_fit, sphericity_list, self$options)

            # Contrasts
            contrasts_html <- .je_rm_contrasts(long, rm_levels, rm_fit$lm_fit, self$options)

            # Post hoc
            posthoc_html <- .je_rm_posthoc(long, rm_fit, between, cells, factor_list, rm_levels, self$options)

            # Marginal means
            marginal_html <- .je_rm_marginal_means(long, rm_fit, between, factor_list, rm_levels, self$options)

            # Simple effects
            simple_html  <- .je_rm_simple_effects(long, between, self$options)

            # Nonparametric
            nonpar_html  <- .je_rm_nonparametric(wide, cells, long, self$options)

            # APA
            apa_html     <- .je_rm_apa(rm_summary, factor_list, sphericity_list)

            sph_summary <- if (!is.null(sphericity) && !is.na(sphericity$W))
                paste0(" — Mauchly's W = ", .je_fmt(sphericity$W),
                       ", p = ", .je_p(sphericity$p),
                       "; GG ε = ", .je_fmt(sphericity$gg_eps),
                       ", HF ε = ", .je_fmt(sphericity$hf_eps))
            else ""

            n_factors_str <- if (length(factor_list) == 1)
                paste0("one within-subject factor (", factor_names[1], ")")
            else
                paste0(length(factor_list), " within-subject factors (",
                       paste(factor_names, collapse = " × "), ")")

            self$results$status$setContent(paste0(
                "<p><strong>Engine:</strong> Repeated-measures ANOVA with ", n_factors_str,
                " via <code>aov(Error(Subject/...))</code>", sph_summary, ".</p>",
                "<p><code>emmeans</code>: ",
                if (requireNamespace("emmeans", quietly = TRUE)) "installed" else "<strong>missing</strong>",
                ".</p>"
            ))

            self$results$design$setContent(
                .je_rm_design_html(cells, factor_list, between, covs,
                                   within_model_terms, between_model_terms)
            )
            self$results$modelSummary$setContent(
                .je_rm_model_html(rm_summary, sphericity_list, self$options)
            )
            self$results$descriptivesSection$setContent(.je_table_html(desc_tab))
            self$results$effectSizesSection$setContent(.je_table_html(effect_tab))
            self$results$assumptionsSection$setContent(assumptions)
            self$results$contrasts$setContent(contrasts_html)
            self$results$orderRestrictions$setContent(
                .je_rm_order_restricted_full(long, rm_fit, factor_list, between, self$options)
            )
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
                .je_rm_repro_html(cells, factor_list, between, covs,
                                  within_model_terms, between_model_terms)
            )

            private$.plotState <- list(long = long, lm_fit = rm_fit$lm_fit)
            if (isTRUE(self$options$qqPlot))
                self$results$qqPlot$setState(private$.plotState)
            if (nzchar(as.character(self$options$rainCloudHorizontalAxis %||% "")))
                self$results$raincloudPlot$setState(private$.plotState)

            private$.populateOutputs(rm_fit$lm_fit, long, row_nums)
        },

        .populateOutputs = function(fit, long, row_nums) {
            agg_fn <- function(vals)
                stats::aggregate(vals, list(Subject = long$Subject), mean)$x

            if (isTRUE(self$options$residualsSavedToData) &&
                (is.null(self$options$residualsSavedToDataType) ||
                 identical(as.character(self$options$residualsSavedToDataType), "raw")) &&
                !is.null(self$options$residsOV) && self$results$residsOV$isNotFilled()) {
                self$results$residsOV$setRowNums(row_nums)
                self$results$residsOV$setValues(agg_fn(stats::residuals(fit)))
            }
            if (isTRUE(self$options$residualsSavedToData) &&
                identical(as.character(self$options$residualsSavedToDataType), "student") &&
                !is.null(self$options$residsOV) && self$results$residsOV$isNotFilled()) {
                self$results$residsOV$setRowNums(row_nums)
                self$results$residsOV$setValues(agg_fn(stats::rstudent(fit)))
            }
            if (isTRUE(self$options$residualsSavedToData) &&
                identical(as.character(self$options$residualsSavedToDataType), "standard") &&
                !is.null(self$options$residsOV) && self$results$residsOV$isNotFilled()) {
                self$results$residsOV$setRowNums(row_nums)
                self$results$residsOV$setValues(agg_fn(stats::rstandard(fit)))
            }
            if (isTRUE(self$options$predictionsSavedToData) && !is.null(self$options$predictOV) &&
                self$results$predictOV$isNotFilled()) {
                self$results$predictOV$setRowNums(row_nums)
                self$results$predictOV$setValues(agg_fn(stats::fitted(fit)))
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

# Parse the factor specification from options.
# Returns a list of lists: list(list(name="Time", levels=c("Pre","Post")), ...)
# Spec format: "FactorA:level1,level2;FactorB:level1,level2,level3"
# Falls back to rmFactorNames/rmFactorLevels for single-factor (or empty spec).
.je_rm_parse_factor_spec <- function(options, cells) {
    spec <- trimws(as.character(options$rmFactorSpec %||% ""))

    if (nzchar(spec) && grepl(":", spec, fixed = TRUE)) {
        factor_parts <- trimws(unlist(strsplit(spec, ";", fixed = TRUE)))
        factors <- lapply(factor_parts, function(fp) {
            pos <- regexpr(":", fp, fixed = TRUE)
            if (pos < 0) return(NULL)
            nm  <- trimws(substr(fp, 1, pos - 1))
            lvs <- trimws(unlist(strsplit(substr(fp, pos + 1, nchar(fp)), ",", fixed = TRUE)))
            lvs <- lvs[nzchar(lvs)]
            if (!nzchar(nm) || length(lvs) < 2) return(NULL)
            list(name = nm, levels = lvs)
        })
        factors <- factors[!vapply(factors, is.null, logical(1))]
        total <- prod(vapply(factors, function(f) length(f$levels), integer(1)))
        if (length(factors) > 0 && total == length(cells))
            return(factors)
    }

    # Single-factor fallback
    nm  <- trimws(as.character(options$rmFactorNames %||% "Within"))
    if (!nzchar(nm)) nm <- "Within"
    lvs_str <- trimws(as.character(options$rmFactorLevels %||% ""))
    lvs <- if (nzchar(lvs_str)) {
        l <- trimws(unlist(strsplit(lvs_str, ",", fixed = TRUE)))
        l[nzchar(l)]
    } else character()
    if (length(lvs) != length(cells)) lvs <- cells
    list(list(name = nm, levels = lvs))
}

.je_rm_split_term <- function(term, available) {
    parts <- trimws(unlist(strsplit(term, "\\*|:|×|x|X")))
    parts <- parts[nzchar(parts)]
    if (!length(parts))
        return(character())
    matched <- vapply(parts, function(part) {
        hit <- available[tolower(available) == tolower(part)]
        if (length(hit)) hit[[1]] else part
    }, character(1))
    matched[matched %in% available]
}

.je_rm_full_factorial_terms <- function(available) {
    available <- available[nzchar(available)]
    if (!length(available))
        return(character())
    terms <- character()
    for (k in seq_along(available)) {
        combos <- utils::combn(available, k, simplify = FALSE)
        terms <- c(terms, vapply(combos, function(x) paste(x, collapse = ":"), character(1)))
    }
    terms
}

.je_rm_parse_model_terms <- function(value, available, full_factorial = FALSE) {
    available <- unique(available[nzchar(available)])
    if (!length(available))
        return(character())

    raw <- trimws(as.character(value %||% ""))
    if (!nzchar(raw))
        return(if (isTRUE(full_factorial)) .je_rm_full_factorial_terms(available) else available)

    pieces <- trimws(unlist(strsplit(raw, "\\n|,|;")))
    pieces <- pieces[nzchar(pieces)]
    terms <- vapply(pieces, function(piece) {
        components <- .je_rm_split_term(piece, available)
        if (!length(components))
            return("")
        paste(unique(components), collapse = ":")
    }, character(1))
    terms <- unique(terms[nzchar(terms)])
    if (!length(terms))
        return(if (isTRUE(full_factorial)) .je_rm_full_factorial_terms(available) else available)
    terms
}

.je_rm_term_components <- function(terms) {
    unique(unlist(strsplit(terms, ":", fixed = TRUE), use.names = FALSE))
}

.je_rm_terms_formula <- function(terms) {
    terms <- terms[nzchar(terms)]
    if (!length(terms))
        return("1")
    paste(terms, collapse = " + ")
}

.je_rm_option <- function(options, name, default = NULL) {
    tryCatch(options[[name]], error = function(e) default)
}

.je_rm_term_factor <- function(long, term) {
    components <- trimws(unlist(strsplit(term, ":", fixed = TRUE)))
    components <- components[nzchar(components) & components %in% names(long)]
    if (!length(components))
        return(NULL)
    if (length(components) == 1)
        return(as.factor(long[[components]]))
    interaction(long[, components, drop = FALSE], sep = " × ", drop = TRUE)
}

.je_rm_term_label <- function(term) {
    gsub(":", " × ", term, fixed = TRUE)
}

.je_rm_emmeans_specs <- function(term) {
    components <- trimws(unlist(strsplit(term, ":", fixed = TRUE)))
    components <- components[nzchar(components)]
    if (!length(components))
        return(NULL)
    stats::as.formula(paste("~", paste(components, collapse = "*")))
}

# Build a long-format data frame for any number of within-subject factors.
# factor_list: output of .je_rm_parse_factor_spec()
# Cells are assumed to be in expand.grid order (last factor varies fastest).
.je_rm_wide_to_long_multi <- function(wide, cells, factor_list, between, covs) {
    id          <- seq_len(nrow(wide))
    factor_names <- vapply(factor_list, function(f) f$name, character(1))

    # All factor-level combinations in expand.grid order
    combo_grid <- expand.grid(
        lapply(rev(factor_list), function(f) f$levels),
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    combo_grid <- combo_grid[, rev(seq_len(ncol(combo_grid))), drop = FALSE]
    names(combo_grid) <- factor_names

    pieces <- lapply(seq_along(cells), function(i) {
        out <- data.frame(Subject = factor(id), Value = wide[[cells[i]]])
        for (fn in factor_names)
            out[[fn]] <- factor(combo_grid[i, fn],
                                levels = factor_list[[match(fn, factor_names)]]$levels)
        for (v in c(between, covs)) if (v %in% names(wide)) out[[v]] <- wide[[v]]
        out
    })
    long <- do.call(rbind, pieces)

    # Keep "Within" column for single-factor backward compat and fallback functions
    if (length(factor_names) == 1) {
        long$Within <- long[[factor_names]]
    } else {
        long$Within <- interaction(long[, factor_names, drop = FALSE],
                                   sep = " × ", drop = TRUE)
    }
    long
}

# ── Model fitting ──────────────────────────────────────────────────────────────

.je_rm_fit_multi <- function(wide, cells, factor_list, between, covs,
                             within_model_terms = NULL,
                             between_model_terms = NULL,
                             options = NULL) {
    long         <- .je_rm_wide_to_long_multi(wide, cells, factor_list, between, covs)
    factor_names <- vapply(factor_list, function(f) f$name, character(1))

    for (fn in factor_names) long[[fn]] <- as.factor(long[[fn]])
    for (f  in between)      long[[f]]  <- as.factor(long[[f]])

    within_model_terms <- within_model_terms %||% .je_rm_full_factorial_terms(factor_names)
    between_model_terms <- between_model_terms %||% c(between, covs)
    within_model_terms <- within_model_terms[nzchar(within_model_terms)]
    between_model_terms <- between_model_terms[nzchar(between_model_terms)]

    crossed_terms <- character()
    if (length(within_model_terms) && length(between_model_terms)) {
        crossed_terms <- as.vector(outer(
            within_model_terms, between_model_terms,
            FUN = function(a, b) paste(a, b, sep = ":")
        ))
    }
    fixed_terms <- unique(c(within_model_terms, between_model_terms, crossed_terms))
    fixed_rhs <- .je_rm_terms_formula(fixed_terms)

    # Error() term follows JASP's Subject/(selected within terms) structure.
    error_within <- .je_rm_terms_formula(within_model_terms)
    error_term <- paste0("Subject/(", error_within, ")")

    aov_form <- stats::as.formula(paste("Value ~", fixed_rhs, "+ Error(", error_term, ")"))
    aov_fit  <- tryCatch(stats::aov(aov_form, data = long), error = function(e) e)
    if (inherits(aov_fit, "error"))
        return(list(error = aov_fit$message))

    # Plain lm for residual diagnostics and emmeans (ignores repeated structure but sufficient)
    lm_form <- stats::as.formula(paste("Value ~", fixed_rhs))
    lm_fit  <- tryCatch(stats::lm(lm_form, data = long), error = function(e) NULL)

    list(aov_fit    = aov_fit,
         lm_fit     = lm_fit,
         long       = long,
         factor_list  = factor_list,
         factor_names = factor_names,
         within_model_terms = within_model_terms,
         between_model_terms = between_model_terms,
         fixed_terms = fixed_terms,
         formula = aov_form,
         error      = NULL)
}

# ── ANOVA summary (within + between strata) ────────────────────────────────────

.je_rm_anova_summary <- function(rm_fit, between) {
    s <- summary(rm_fit$aov_fit)
    strata <- names(s)
    rows <- list()
    for (stratum in strata) {
        tables <- s[[stratum]]
        if (inherits(tables, "summary.aov"))
            tables <- unclass(tables)
        for (tab in tables) {
            tab <- as.data.frame(tab)
            if (!nrow(tab))
                next
            tab$Term <- trimws(rownames(tab))
            tab$Stratum <- stratum
            rownames(tab) <- NULL
            rows[[length(rows) + 1]] <- tab
        }
    }
    if (!length(rows))
        return(data.frame())
    result <- do.call(rbind, rows)
    # Normalise column names
    names(result) <- gsub("Pr\\(>F\\)", "p", names(result))
    names(result) <- gsub("F value",   "F", names(result))
    names(result) <- gsub("Mean Sq",   "Mean Sq", names(result))
    names(result) <- gsub("Sum Sq",    "Sum Sq",  names(result))
    result
}

# ── Sphericity ────────────────────────────────────────────────────────────────
# Core Mauchly + GG/HF computation for a k-level wide matrix (k columns, n subjects)
.je_rm_sphericity_core <- function(mat) {
    k <- ncol(mat)
    n <- nrow(mat)
    if (k < 3 || n < k) return(NULL)

    C      <- stats::contr.helmert(k)
    C      <- apply(C, 2, function(v) v / sqrt(sum(v^2)))
    Y      <- mat %*% C
    S      <- stats::cov(Y)
    p_hat  <- k - 1

    det_S  <- det(S)
    tr_S   <- sum(diag(S))
    denom  <- (tr_S / p_hat)^p_hat
    W      <- if (denom > 0 && det_S > 0) det_S / denom else NA_real_

    f_val  <- -(n - 1 - (2 * p_hat^2 + p_hat + 2) / (6 * p_hat))
    chi_sq <- if (!is.na(W) && W > 0) f_val * log(W) else NA_real_
    df_mau <- p_hat * (p_hat + 1) / 2 - 1
    p_mau  <- if (!is.na(chi_sq) && !is.na(df_mau) && df_mau > 0)
        stats::pchisq(-chi_sq, df = df_mau, lower.tail = FALSE) else NA_real_

    tr_S2  <- sum(S^2)
    safe_d <- p_hat * (tr_S2 - tr_S^2 / p_hat)
    gg_eps <- if (safe_d > 0) max(1 / p_hat, min(1, tr_S^2 / safe_d)) else 1

    hf_num <- n * p_hat * gg_eps - 2
    hf_den <- p_hat * (n - 1 - p_hat * gg_eps)
    hf_eps <- if (hf_den > 0) min(1, hf_num / hf_den) else 1

    list(W = W, chi_sq = chi_sq, df = df_mau, p = p_mau,
         gg_eps = gg_eps, hf_eps = hf_eps, k = k)
}

# Per-effect sphericity for multi-factor designs.
# Returns a list of list(effect = <name>, sphericity = <result or NULL>).
.je_rm_sphericity_multi <- function(wide, cells, factor_list, options) {
    factor_names <- vapply(factor_list, function(f) f$name, character(1))
    n_factors    <- length(factor_names)

    # Combination grid matching cell order
    combo_grid <- expand.grid(
        lapply(rev(factor_list), function(f) f$levels),
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    combo_grid <- combo_grid[, rev(seq_len(ncol(combo_grid))), drop = FALSE]
    names(combo_grid) <- factor_names

    results <- list()

    # Test each main effect: marginalise over other factors by averaging
    for (fn in factor_names) {
        these_levels <- factor_list[[match(fn, factor_names)]]$levels
        if (length(these_levels) < 3) {
            results[[length(results) + 1]] <- list(effect = fn, sphericity = NULL)
            next
        }
        mat <- do.call(cbind, lapply(these_levels, function(lv) {
            idx <- which(combo_grid[[fn]] == lv)
            if (length(idx) == 1) wide[[cells[idx]]]
            else rowMeans(as.matrix(wide[, cells[idx], drop = FALSE]), na.rm = TRUE)
        }))
        colnames(mat) <- these_levels
        results[[length(results) + 1]] <- list(
            effect = fn, sphericity = .je_rm_sphericity_core(mat)
        )
    }

    # Interaction effect(s) — only for 2-factor case
    if (n_factors == 2) {
        ia_name <- paste(factor_names, collapse = " × ")
        mat     <- as.matrix(wide[, cells, drop = FALSE])
        results[[length(results) + 1]] <- list(
            effect = ia_name, sphericity = .je_rm_sphericity_core(mat)
        )
    }

    results
}

# ── Apply sphericity correction to within-subject rows ───────────────────────
# sphericity_list: output of .je_rm_sphericity_multi()

.je_rm_apply_correction <- function(rm_summary, sphericity_list, options) {
    use_gg <- isTRUE(options$sphericityCorrectionGreenhouseGeisser)
    use_hf <- isTRUE(options$sphericityCorrectionHuynhFeldt)
    if ((!use_gg && !use_hf) || length(sphericity_list) == 0) return(rm_summary)
    # GG takes precedence when both selected
    correction <- if (use_gg) "greenhouseGeisser" else "huynhFeldt"

    for (sph_entry in sphericity_list) {
        sph <- sph_entry$sphericity
        if (is.null(sph)) next
        eps <- if (correction == "greenhouseGeisser") sph$gg_eps else sph$hf_eps

        # Match strata that contain this effect name (Term column)
        effect_name <- sph_entry$effect
        hit_rows <- which(
            grepl(gsub(" × ", ":", effect_name, fixed = TRUE), rm_summary$Stratum, fixed = TRUE) |
            grepl(effect_name, rm_summary$Term, fixed = TRUE)
        )
        non_resid <- hit_rows[!grepl("Residuals", rm_summary$Term[hit_rows], ignore.case = TRUE)]

        for (i in non_resid) {
            err_row <- which(rm_summary$Stratum == rm_summary$Stratum[i] &
                             grepl("Residuals", rm_summary$Term, ignore.case = TRUE))
            if (length(err_row) == 0) next
            df1_c <- rm_summary$Df[i]           * eps
            df2_c <- rm_summary$Df[err_row[1]]  * eps
            f_val <- rm_summary$F[i]
            if (!is.na(f_val) && df1_c > 0 && df2_c > 0)
                rm_summary$p[i] <- stats::pf(f_val, df1_c, df2_c, lower.tail = FALSE)
            rm_summary$Df[i]          <- df1_c
            rm_summary$Df[err_row[1]] <- df2_c
        }
    }
    rm_summary
}

# ── Effect sizes for RM ANOVA ─────────────────────────────────────────────────

.je_rm_effect_sizes_tab <- function(rm_summary, options) {
    if (!isTRUE(options$effectSizeEstimates)) return(data.frame())
    non_resid <- rm_summary[!grepl("Residuals", rm_summary$Term, ignore.case = TRUE), , drop = FALSE]
    resid     <- rm_summary[grepl("Residuals",  rm_summary$Term, ignore.case = TRUE), , drop = FALSE]
    if (nrow(non_resid) == 0 || nrow(resid) == 0) return(data.frame())

    ss_err_within <- sum(resid$`Sum Sq`[grep("Within", resid$Stratum)], na.rm = TRUE)
    ss_total      <- sum(rm_summary$`Sum Sq`, na.rm = TRUE)
    df_err_within <- sum(resid$Df[grep("Within", resid$Stratum)], na.rm = TRUE)
    mse_within    <- if (df_err_within > 0) ss_err_within / df_err_within else NA

    out <- data.frame(Term = non_resid$Term, Stratum = non_resid$Stratum, stringsAsFactors = FALSE)

    if (isTRUE(options$effectSizeEtaSquared))
        out$`eta squared` <- .je_fmt(non_resid$`Sum Sq` / ss_total)

    if (isTRUE(options$effectSizePartialEtaSquared)) {
        # For within-subject terms use the within-subject error
        ss_e <- ifelse(grepl("Within", non_resid$Stratum), ss_err_within,
                       sum(resid$`Sum Sq`[grep("Subject\\]$", resid$Stratum)], na.rm = TRUE))
        out$`partial eta squared` <- .je_fmt(non_resid$`Sum Sq` / (non_resid$`Sum Sq` + ss_e))
    }

    if (isTRUE(options$effectSizeOmegaSquared) && !is.na(mse_within)) {
        out$`omega squared` <- .je_fmt(pmax(0,
            (non_resid$`Sum Sq` - non_resid$Df * mse_within) / (ss_total + mse_within)
        ))
    }

    if (isTRUE(options$effectSizePartialOmegaSquared) && !is.na(mse_within)) {
        out$`partial omega squared` <- .je_fmt(pmax(0,
            (non_resid$`Sum Sq` - non_resid$Df * mse_within) /
            (non_resid$`Sum Sq` + ss_err_within + mse_within)
        ))
    }

    if (isTRUE(options$effectSizeGeneralEtaSquared)) {
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
    alpha    <- 1 - ((options$effectSizeCiLevel %||% 95) / 100)

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

.je_rm_assumptions_full <- function(long, rm_fit, sphericity_list, options) {
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

    # Sphericity — one block per effect
    corr <- if (isTRUE(options$sphericityCorrectionGreenhouseGeisser)) "greenhouseGeisser"
            else if (isTRUE(options$sphericityCorrectionHuynhFeldt)) "huynhFeldt"
            else "none"
    corr_label <- switch(corr, none = "None",
                         greenhouseGeisser = "Greenhouse-Geisser",
                         huynhFeldt = "Huynh-Feldt", "None")

    if (length(sphericity_list) > 0) {
        tested <- FALSE
        for (e in sphericity_list) {
            sph <- e$sphericity
            if (is.null(sph)) {
                out <- c(out, paste0("<p>Sphericity (", .je_escape(e$effect),
                                     "): not applicable (≤ 2 levels).</p>"))
                next
            }
            tested <- TRUE
            out <- c(out, paste0(
                "<p>Mauchly's test (", .je_escape(e$effect), "): W = ",
                .je_fmt(sph$W), ", χ²(", sph$df, ") = ", .je_fmt(sph$chi_sq),
                ", p = ", .je_p(sph$p), ".</p>",
                "<p>GG ε = ", .je_fmt(sph$gg_eps),
                "; HF ε = ", .je_fmt(sph$hf_eps), ".</p>"
            ))
        }
        if (tested)
            out <- c(out, paste0("<p>Active correction: <strong>", corr_label, "</strong>.</p>"))
    } else {
        out <- c(out, "<p>Sphericity not tested.</p>")
    }

    paste(out, collapse = "")
}

# ── Model summary HTML ────────────────────────────────────────────────────────

.je_rm_model_html <- function(rm_summary, sphericity_list, options) {
    corr_summary <- .je_rm_apply_correction(rm_summary, sphericity_list, options)
    display <- data.frame(
        Term      = corr_summary$Term,
        Stratum   = corr_summary$Stratum,
        df        = .je_fmt(corr_summary$Df),
        `Sum Sq`  = .je_fmt(corr_summary$`Sum Sq`),
        `Mean Sq` = .je_fmt(corr_summary$`Mean Sq`),
        F         = .je_fmt(corr_summary$F),
        p         = .je_p(corr_summary$p),
        check.names = FALSE
    )

    if (isTRUE(options$vovkSellke)) {
        vs <- .je_vovk_sellke(corr_summary$p)
        display$`VS-MPR` <- ifelse(is.na(vs), "", .je_fmt(vs))
    }

    corr_label <- if (isTRUE(options$sphericityCorrectionGreenhouseGeisser))
        " (Greenhouse-Geisser corrected df)"
    else if (isTRUE(options$sphericityCorrectionHuynhFeldt))
        " (Huynh-Feldt corrected df)"
    else ""

    # Sphericity table
    sph_html <- if (length(sphericity_list) > 0) {
        sph_rows <- lapply(sphericity_list, function(e) {
            s <- e$sphericity
            if (is.null(s)) return(NULL)
            data.frame(
                Effect  = e$effect,
                W       = .je_fmt(s$W),
                `chi-sq`= .je_fmt(s$chi_sq),
                df      = .je_fmt(s$df),
                p       = .je_p(s$p),
                `GG-eps`= .je_fmt(s$gg_eps),
                `HF-eps`= .je_fmt(s$hf_eps),
                check.names = FALSE
            )
        })
        sph_rows <- sph_rows[!vapply(sph_rows, is.null, logical(1))]
        if (length(sph_rows) > 0)
            paste0("<h4>Mauchly's Test of Sphericity</h4>",
                   .je_table_html(do.call(rbind, sph_rows)))
        else ""
    } else ""

    paste0("<p>Repeated-measures ANOVA summary", corr_label, ".</p>",
           .je_table_html(display), sph_html)
}

# ── Contrasts ─────────────────────────────────────────────────────────────────

.je_rm_contrasts <- function(long, rm_levels, lm_fit, options) {
    requested <- options$contrastType != "none" ||
        isTRUE(.je_rm_option(options, "rmCustomContrasts", FALSE))
    if (!requested) return("<p>No repeated-measures contrasts are enabled.</p>")

    n    <- length(rm_levels)
    mats <- list()
    ct   <- as.character(options$contrastType %||% "none")

    if (ct == "helmert")    mats$Helmert    <- stats::contr.helmert(n)
    if (ct == "polynomial") mats$Polynomial <- stats::contr.poly(n)
    if (ct == "deviation")  mats$Deviation  <- stats::contr.sum(n)
    if (ct %in% c("difference", "repeated")) mats$Difference <- .je_contr_diff(n)
    if (ct == "simple") {
        mat <- matrix(0, nrow = n, ncol = n - 1)
        for (i in seq_len(n - 1)) { mat[i, i] <- 1; mat[n, i] <- -1 }
        mats$Simple <- mat
    }

    if (isTRUE(.je_rm_option(options, "rmCustomContrasts", FALSE))) {
        custom <- .je_parse_custom_contrasts(.je_rm_option(options, "rmContrastSyntax", ""), n)
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

.je_rm_posthoc <- function(long, rm_fit, between, cells, factor_list, rm_levels, options) {
    factor_names <- vapply(factor_list, function(f) f$name, character(1))
    rm_terms <- .je_rm_parse_model_terms(
        .je_rm_option(options, "postHocRmTerms", ""), factor_names, full_factorial = FALSE
    )
    between_terms <- .je_chr_vec(options$postHocTerms)
    between_terms <- between_terms[between_terms %in% between]
    if (length(rm_terms) == 0 && length(between_terms) == 0)
        return("<p>Post hoc tests are disabled.</p>")
    chunks  <- character()
    has_emm <- requireNamespace("emmeans", quietly = TRUE) && !is.null(rm_fit$lm_fit)

    adj <- if (isTRUE(options$postHocCorrectionBonferroni)) "bonferroni"
           else if (isTRUE(options$postHocCorrectionTukey))  "tukey"
           else if (isTRUE(options$postHocCorrectionScheffe)) "scheffe"
           else "holm"

    # ── Within-subject factor pairwise comparisons ────────────────────────────
    for (term in rm_terms) {
        term_label <- .je_rm_term_label(term)
        if (has_emm) {
            emm <- tryCatch(
                emmeans::emmeans(rm_fit$lm_fit, specs = .je_rm_emmeans_specs(term)),
                error = function(e) NULL
            )
            if (!is.null(emm)) {
                pw <- tryCatch(as.data.frame(emmeans::contrast(emm, method = "pairwise", adjust = adj)),
                               error = function(e) NULL)
                if (!is.null(pw)) {
                    out <- .je_rm_pw_table(pw, options)
                    chunks <- c(chunks, paste0(
                        "<h4>Within-subject pairwise comparisons - ", .je_escape(term_label),
                        " (", .je_escape(adj), ")</h4>",
                        .je_table_html(out)
                    ))
                    if (isTRUE(options$postHocEffectSize))
                        chunks <- c(chunks, .je_rm_posthoc_d(pw, long))
                }
            }
        } else {
            group <- .je_rm_term_factor(long, term)
            if (is.null(group))
                next
            pw <- tryCatch(
                stats::pairwise.t.test(long$Value, group, paired = TRUE, p.adjust.method = adj),
                error = function(e) NULL
            )
            if (!is.null(pw))
                chunks <- c(chunks, paste0(
                    "<h4>Paired pairwise t-tests - ", .je_escape(term_label),
                    " (", .je_escape(adj), ") - install emmeans for full output</h4>",
                    .je_pairwise_html(pw$p.value)
                ))
        }
    }

    # ── Between-subject factor pairwise comparisons ───────────────────────────
    if (length(between_terms) > 0 && has_emm) {
        for (f in between_terms) {
            emm <- tryCatch(emmeans::emmeans(rm_fit$lm_fit, specs = f), error = function(e) NULL)
            if (is.null(emm)) next
            pw <- tryCatch(
                as.data.frame(emmeans::contrast(emm, method = "pairwise", adjust = adj)),
                error = function(e) NULL
            )
            if (!is.null(pw)) {
                out <- .je_rm_pw_table(pw, options)
                chunks <- c(chunks, paste0(
                    "<h4>Between-subjects: ", .je_escape(f), " (", .je_escape(adj), ")</h4>",
                    .je_table_html(out)
                ))
            }
        }
    }

    # ── Interaction post hoc: Within at each level of Between ─────────────────
    if (length(between_terms) > 0 && has_emm) {
        for (f in between_terms) {
            # Within-subject comparisons at each level of the between factor
            emm_w_by_b <- tryCatch(
                emmeans::emmeans(rm_fit$lm_fit, specs = "Within", by = f),
                error = function(e) NULL
            )
            if (!is.null(emm_w_by_b)) {
                pw <- tryCatch(
                    as.data.frame(emmeans::contrast(emm_w_by_b, method = "pairwise", adjust = adj)),
                    error = function(e) NULL
                )
                if (!is.null(pw)) {
                    out <- .je_rm_pw_table_by(pw, f, options)
                    chunks <- c(chunks, paste0(
                        "<h4>Within-subject comparisons by ", .je_escape(f),
                        " (", .je_escape(adj), ")</h4>",
                        .je_table_html(out)
                    ))
                }
            }

            # Between-subject comparisons at each level of the within factor
            emm_b_by_w <- tryCatch(
                emmeans::emmeans(rm_fit$lm_fit, specs = f, by = "Within"),
                error = function(e) NULL
            )
            if (!is.null(emm_b_by_w)) {
                pw <- tryCatch(
                    as.data.frame(emmeans::contrast(emm_b_by_w, method = "pairwise", adjust = adj)),
                    error = function(e) NULL
                )
                if (!is.null(pw)) {
                    out <- .je_rm_pw_table_by(pw, "Within", options)
                    chunks <- c(chunks, paste0(
                        "<h4>Between-subject comparisons (", .je_escape(f),
                        ") by Within-level (", .je_escape(adj), ")</h4>",
                        .je_table_html(out)
                    ))
                }
            }
        }
    }

    if (length(chunks) == 0) return("<p>No post hoc output could be calculated.</p>")
    paste(chunks, collapse = "")
}

# Helper: format an emmeans pairwise contrast data frame into a display table
.je_rm_pw_table <- function(pw, options) {
    out <- data.frame(
        Comparison = as.character(pw$contrast),
        Estimate   = .je_fmt(pw$estimate),
        SE         = .je_fmt(pw$SE),
        df         = .je_fmt(pw$df),
        t          = .je_fmt(pw$t.ratio),
        p          = .je_p(pw$p.value),
        check.names = FALSE
    )
    if (isTRUE(options$postHocSignificanceFlag))
        out$sig <- ifelse(!is.na(pw$p.value) & pw$p.value < .05, "*", "")
    if (isTRUE(options$vovkSellke))
        out$`VS-MPR` <- .je_fmt(.je_vovk_sellke(pw$p.value))
    out
}

# Helper: format pairwise table when a "by" grouping column is present
.je_rm_pw_table_by <- function(pw, by_col, options) {
    by_vals <- if (by_col %in% names(pw)) as.character(pw[[by_col]]) else rep("", nrow(pw))
    out <- data.frame(
        By         = by_vals,
        Comparison = as.character(pw$contrast),
        Estimate   = .je_fmt(pw$estimate),
        SE         = .je_fmt(pw$SE),
        df         = .je_fmt(pw$df),
        t          = .je_fmt(pw$t.ratio),
        p          = .je_p(pw$p.value),
        check.names = FALSE
    )
    if (isTRUE(options$postHocSignificanceFlag))
        out$sig <- ifelse(!is.na(pw$p.value) & pw$p.value < .05, "*", "")
    if (isTRUE(options$vovkSellke))
        out$`VS-MPR` <- .je_fmt(.je_vovk_sellke(pw$p.value))
    out
}

# Cohen's d for within-subject contrasts (from paired differences)
.je_rm_posthoc_d <- function(pw, long) {
    wide_v <- tryCatch(
        stats::reshape(long[, c("Subject", "Within", "Value")],
                       idvar = "Subject", timevar = "Within", direction = "wide"),
        error = function(e) NULL
    )
    if (is.null(wide_v)) return("")
    rows <- lapply(seq_len(nrow(pw)), function(i) {
        parts <- trimws(strsplit(as.character(pw$contrast[i]), " - ")[[1]])
        if (length(parts) != 2) return(NULL)
        c1 <- paste0("Value.", parts[1])
        c2 <- paste0("Value.", parts[2])
        if (!(c1 %in% names(wide_v)) || !(c2 %in% names(wide_v))) return(NULL)
        diffs <- wide_v[[c1]] - wide_v[[c2]]
        d <- mean(diffs, na.rm = TRUE) / stats::sd(diffs, na.rm = TRUE)
        data.frame(Comparison = as.character(pw$contrast[i]),
                   `Cohen's d` = .je_fmt(d), check.names = FALSE)
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) == 0) return("")
    paste0("<h4>Within-subject effect sizes (Cohen's d<sub>z</sub>)</h4>",
           .je_table_html(do.call(rbind, rows)))
}

# ── Marginal means ────────────────────────────────────────────────────────────

.je_rm_marginal_means <- function(long, rm_fit, between, factor_list, rm_levels, options) {
    factor_names <- vapply(factor_list, function(f) f$name, character(1))
    rm_terms <- .je_rm_parse_model_terms(
        .je_rm_option(options, "marginalMeanRmTerms", ""), factor_names, full_factorial = FALSE
    )
    between_terms <- .je_chr_vec(options$marginalMeanTerms)
    between_terms <- between_terms[between_terms %in% between]
    terms <- unique(c(rm_terms, between_terms))
    if (length(terms) == 0) return("<p>Estimated marginal means are disabled.</p>")

    chunks <- character()

    if (requireNamespace("emmeans", quietly = TRUE) && !is.null(rm_fit$lm_fit)) {
        for (term in terms) {
        emm <- tryCatch(emmeans::emmeans(rm_fit$lm_fit, specs = .je_rm_emmeans_specs(term)), error = function(e) NULL)
        if (!is.null(emm)) {
            ci_level <- (options$effectSizeCiLevel %||% 95) / 100
            em_df    <- as.data.frame(stats::confint(emm, level = ci_level))
            level_col <- intersect(strsplit(term, ":", fixed = TRUE)[[1]], names(em_df))
            level_text <- if (length(level_col)) {
                apply(em_df[, level_col, drop = FALSE], 1, paste, collapse = " × ")
            } else seq_len(nrow(em_df))
            out <- data.frame(
                Level = as.character(level_text),
                EMM   = .je_fmt(em_df$emmean),
                SE    = .je_fmt(em_df$SE),
                df    = .je_fmt(em_df$df),
                Lower = .je_fmt(em_df$lower.CL),
                Upper = .je_fmt(em_df$upper.CL),
                check.names = FALSE
            )
            html <- paste0("<h4>", .je_escape(.je_rm_term_label(term)),
                           "</h4><p>Estimated marginal means via <code>emmeans</code>.</p>",
                           .je_table_html(out))

            if (isTRUE(options$marginalMeanComparedToZero) || length(strsplit(term, ":", fixed = TRUE)[[1]]) == 1) {
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
                    html <- paste0(html, "<h4>Pairwise comparisons (", .je_escape(adj), ")</h4>",
                                   .je_table_html(pw_out))
                }
            }
            chunks <- c(chunks, html)
        }
        }
        if (length(chunks))
            return(paste(chunks, collapse = ""))
    }

    # Fallback: observed means per within-level
    desc <- .je_rm_descriptives_full(long, character(), options)
    paste0("<p><em>Observed cell means (install <code>emmeans</code> for covariate-adjusted estimates).</em></p>",
           .je_table_html(desc))
}

# ── Simple effects ────────────────────────────────────────────────────────────

.je_rm_simple_effects <- function(long, between, options) {
    if (!nzchar(as.character(options$simpleMainEffectFactor %||% ""))) return("<p>Simple effects analysis is disabled.</p>")
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
    if (!nzchar(as.character(options$friedmanWithinFactor %||% ""))) return("<p>Friedman analysis is disabled.</p>")
    mat <- as.matrix(wide[, cells, drop = FALSE])
    ft  <- tryCatch(stats::friedman.test(mat), error = function(e) NULL)
    if (is.null(ft)) return("<p>Friedman test could not be run.</p>")
    html <- paste0(
        "<p>Friedman χ²(", ft$parameter, ") = ", .je_fmt(unname(ft$statistic)),
        ", p = ", .je_p(ft$p.value), ".</p>"
    )
    if (isTRUE(options$conoverTest)) {
        html <- paste0(html, "<p>Conover post hoc tests require the <code>NSM3</code> or <code>rstatix</code> package — not yet ported.</p>")
    }
    html
}

# ── APA text ──────────────────────────────────────────────────────────────────

.je_rm_apa <- function(rm_summary, factor_list, sphericity_list) {
    factor_names <- vapply(factor_list, function(f) f$name, character(1))
    design_str   <- paste(vapply(factor_list, function(f)
        paste0(f$name, " (", length(f$levels), " levels)"), character(1)), collapse = " × ")

    within_rows <- rm_summary[grepl(factor_names[1], rm_summary$Stratum) &
                              !grepl("Residuals", rm_summary$Term, ignore.case = TRUE), , drop = FALSE]
    if (nrow(within_rows) == 0)
        return("<p>No APA repeated-measures text available.</p>")

    first_within <- within_rows[1, ]
    df_err <- rm_summary$Df[grepl("Residuals", rm_summary$Term, ignore.case = TRUE) &
                            grepl(factor_names[1], rm_summary$Stratum)][1]

    sph_entry <- if (length(sphericity_list) > 0) sphericity_list[[1]] else NULL
    sph <- if (!is.null(sph_entry)) sph_entry$sphericity else NULL
    sph_note <- if (!is.null(sph) && !is.na(sph$p))
        paste0(" Mauchly's test of sphericity indicated the assumption was ",
               if (sph$p < .05) "violated" else "met",
               ", W = ", .je_fmt(sph$W), ", p = ", .je_p(sph$p), ".")
    else ""

    blocks <- paste0(
        "<p>A ", if (length(factor_list) == 1) "one-factor" else paste0(length(factor_list), "-factor"),
        " repeated-measures ANOVA was conducted: ", .je_escape(design_str), ".", sph_note,
        " The first within-subject effect (", .je_escape(first_within$Term),
        ") was <em>F</em>(", .je_fmt(first_within$Df), ", ", .je_fmt(df_err), ") = ",
        .je_fmt(first_within$F), ", <em>p</em> = ", .je_p(first_within$p), ".</p>"
    )

    # Add blocks for additional effects if multi-factor
    if (length(factor_list) > 1 && nrow(within_rows) > 1) {
        extra <- apply(within_rows[-1, , drop = FALSE], 1, function(r) {
            paste0("<p>", .je_escape(r["Term"]), ": <em>F</em> = ",
                   .je_fmt(as.numeric(r["F"])), ", <em>p</em> = ",
                   .je_p(as.numeric(r["p"])), ".</p>")
        })
        blocks <- paste0(blocks, paste(extra, collapse = ""))
    }
    blocks
}

# ── HTML helpers specific to RM ───────────────────────────────────────────────

.je_rm_design_html <- function(cells, factor_list, between, covs,
                               within_model_terms = character(),
                               between_model_terms = character()) {
    factor_names <- vapply(factor_list, function(f) f$name, character(1))
    factor_detail <- paste(vapply(factor_list, function(f)
        paste0(f$name, " (", paste(f$levels, collapse = ", "), ")"),
        character(1)), collapse = "; ")
    n_factors_str <- if (length(factor_list) == 1) "1 within-subject factor" else
        paste0(length(factor_list), " within-subject factors")
    paste0(
        "<p>Design: ", .je_escape(n_factors_str), " — ", .je_escape(factor_detail), "</p>",
        "<p>Cells assigned (", length(cells), "): ",
        .je_escape(paste(cells, collapse = ", ")), "</p>",
        "<p>Between-subject factors: ",
        .je_escape(if (length(between) == 0) "none" else paste(between, collapse = ", ")), "</p>",
        "<p>Covariates: ",
        .je_escape(if (length(covs) == 0) "none" else paste(covs, collapse = ", ")), "</p>",
        "<p>Within model terms: <code>",
        .je_escape(if (length(within_model_terms) == 0) "none" else paste(within_model_terms, collapse = ", ")),
        "</code></p>",
        "<p>Between model terms: <code>",
        .je_escape(if (length(between_model_terms) == 0) "none" else paste(between_model_terms, collapse = ", ")),
        "</code></p>"
    )
}

.je_rm_order_status <- function(options) {
    syntax <- .je_rm_option(options, "orderRestrictedSyntax",
                            .je_rm_option(options, "restrictedSyntax", ""))
    has_syntax <- isTRUE(nzchar(trimws(as.character(syntax %||% ""))))
    if (!isTRUE(.je_rm_option(options, "orderRestricted", FALSE)) &&
        !isTRUE(.je_rm_option(options, "modelComparison", FALSE)) &&
        !isTRUE(.je_rm_option(options, "informedHypothesisTests", FALSE)) &&
        !has_syntax)
        return("<p>Bayesian/order-restricted hypothesis testing is disabled.</p>")
    paste0("<p><strong>Syntax captured:</strong></p><pre>",
           .je_escape(syntax),
           "</pre>")
}

.je_rm_order_restricted_full <- function(long, rm_fit, factor_list, between, options) {
    syntax <- .je_rm_option(options, "orderRestrictedSyntax",
                            .je_rm_option(options, "restrictedSyntax", ""))
    has_syntax <- isTRUE(nzchar(trimws(as.character(syntax %||% ""))))
    any_active <- isTRUE(.je_rm_option(options, "orderRestricted", FALSE)) ||
                  isTRUE(.je_rm_option(options, "modelComparison", FALSE))  ||
                  isTRUE(.je_rm_option(options, "informedHypothesisTests", FALSE)) ||
                  has_syntax
    if (!any_active)
        return("<p>Bayesian / order-restricted hypothesis testing is disabled.</p>")

    chunks <- character()

    # Order-restricted inference via bain (uses the lm() fit)
    if ((isTRUE(.je_rm_option(options, "orderRestricted", FALSE)) ||
         isTRUE(.je_rm_option(options, "informedHypothesisTests", FALSE)) ||
         has_syntax) &&
        !is.null(rm_fit$lm_fit)) {
        h_raw <- trimws(as.character(syntax))
        if (nzchar(h_raw)) {
            chunks <- c(chunks, .je_bain_analysis(rm_fit$lm_fit, options))
        } else {
            chunks <- c(chunks, paste0(
                "<p><strong>Order-restricted inference (bain):</strong> enter hypothesis ",
                "constraints in the syntax field.</p>"
            ))
        }
    }

    # Bayesian repeated-measures ANOVA via BayesFactor (uses whichRandom = Subject)
    if (isTRUE(.je_rm_option(options, "modelComparison", FALSE)))
        chunks <- c(chunks, .je_bayesian_rm_anova(long, factor_list, between, options))

    paste(chunks, collapse = "")
}

.je_bayesian_rm_anova <- function(long, factor_list, between, options) {
    if (!requireNamespace("BayesFactor", quietly = TRUE))
        return(paste0(
            "<p><strong>BayesFactor not installed.</strong> ",
            "Install with <code>install.packages('BayesFactor')</code>.</p>"
        ))

    factor_names <- vapply(factor_list, function(f) f$name, character(1))
    rhs_terms    <- c(factor_names, between, "Subject")
    formula      <- stats::as.formula(
        paste("Value ~", paste(rhs_terms, collapse = " + "))
    )

    bf_result <- tryCatch(
        BayesFactor::anovaBF(
            formula       = formula,
            data          = long,
            whichRandom   = "Subject",
            whichModels   = "withmain",
            progress      = FALSE
        ),
        error = function(e) e
    )

    if (inherits(bf_result, "error"))
        return(paste0("<p>BayesFactor::anovaBF failed: <em>",
                      .je_escape(bf_result$message), "</em></p>"))

    .je_bayesian_anova_html(bf_result, "Value")
}

.je_rm_plot_status <- function(long, options) {
    enabled <- c(
        if (isTRUE(options$qqPlot)) "Q-Q residual plot",
        if (nzchar(as.character(options$rainCloudHorizontalAxis %||% ""))) "raincloud plot"
    )
    if (length(enabled) == 0) return("<p>No plot options are enabled.</p>")
    paste0("<p>Rendered: ", .je_escape(paste(enabled, collapse = ", ")), ".</p>")
}

.je_rm_saved_columns_html <- function(options) {
    requested <- c(
        if (isTRUE(options$residualsSavedToData) && identical(as.character(options$residualsSavedToDataType %||% "raw"), "raw"))      "mean raw residuals",
        if (isTRUE(options$residualsSavedToData) && identical(as.character(options$residualsSavedToDataType), "student"))   "mean studentized residuals",
        if (isTRUE(options$residualsSavedToData) && identical(as.character(options$residualsSavedToDataType), "standard"))  "mean standardized residuals",
        if (isTRUE(options$predictionsSavedToData))                                             "mean predicted values"
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

.je_rm_repro_html <- function(cells, factor_list, between, covs,
                              within_model_terms = character(),
                              between_model_terms = character()) {
    factor_str <- paste(vapply(factor_list, function(f)
        paste0(f$name, ": ", paste(f$levels, collapse = ", ")), character(1)), collapse = "; ")
    paste0(
        "<p>Cells assigned: <code>", .je_escape(paste(cells, collapse = ", ")), "</code></p>",
        "<p>Within-subject factor(s): <code>", .je_escape(factor_str), "</code></p>",
        "<p>Between factors: <code>",
        .je_escape(if (length(between) == 0) "none" else paste(between, collapse = ", ")),
        "</code></p>",
        "<p>Within model terms: <code>",
        .je_escape(if (length(within_model_terms) == 0) "none" else paste(within_model_terms, collapse = ", ")),
        "</code></p>",
        "<p>Between model terms: <code>",
        .je_escape(if (length(between_model_terms) == 0) "none" else paste(between_model_terms, collapse = ", ")),
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
