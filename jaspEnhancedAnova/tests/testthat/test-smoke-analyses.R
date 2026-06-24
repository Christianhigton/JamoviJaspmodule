test_that("classical ANOVA runs with optional boxes omitted", {
    set.seed(100)
    dat <- data.frame(
        y = rnorm(36),
        group = factor(rep(letters[1:3], each = 12)),
        x = rnorm(36)
    )

    expect_no_error(
        enhancedAnova(
            data = dat,
            dependent = "y",
            fixedFactors = "group",
            covariates = "x"
        )
    )
})

test_that("repeated measures ANOVA runs with optional boxes omitted", {
    set.seed(101)
    dat <- data.frame(
        id = seq_len(24),
        group = factor(rep(c("A", "B"), each = 12)),
        t1 = rnorm(24),
        t2 = rnorm(24),
        t3 = rnorm(24)
    )

    expect_no_error(
        enhancedRepeatedMeasuresAnova(
            data = dat,
            repeatedMeasuresCells = c("t1", "t2", "t3"),
            rmFactorNames = "time",
            rmFactorLevels = "pre,post,follow"
        )
    )

    expect_no_error(
        enhancedRepeatedMeasuresAnova(
            data = dat,
            repeatedMeasuresCells = c("t1", "t2", "t3"),
            betweenSubjectFactors = "group",
            rmFactorNames = "time",
            rmFactorLevels = "pre,post,follow",
            withinModelTerms = "time",
            betweenModelTerms = "group",
            postHocRmTerms = "time",
            marginalMeanRmTerms = "time"
        )
    )
})

test_that("ANCOVA runs with optional boxes omitted", {
    set.seed(102)
    dat <- data.frame(
        y = rnorm(40),
        group = factor(rep(c("A", "B"), each = 20)),
        x = rnorm(40),
        w = runif(40, 0.5, 2)
    )

    expect_no_error(
        enhancedAncova(
            data = dat,
            dependent = "y",
            fixedFactors = "group",
            covariates = "x"
        )
    )

    expect_no_error(
        enhancedAncova(
            data = dat,
            dependent = "y",
            fixedFactors = "group",
            covariates = "x",
            wlsWeights = "w",
            modelTerms = "group, x"
        )
    )
})

test_that("MANCOVA runs with optional boxes omitted", {
    set.seed(103)
    dat <- data.frame(
        y1 = rnorm(42),
        y2 = rnorm(42),
        group = factor(rep(c("A", "B", "C"), each = 14)),
        x = rnorm(42)
    )

    expect_no_error(
        enhancedManova(
            data = dat,
            dependent = c("y1", "y2"),
            fixedFactors = "group",
            covariates = "x"
        )
    )
})
