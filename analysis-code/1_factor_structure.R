# Exploratory and confirmatory factor-structure analyses.

.prepare_lavaan <- function() {
  configure_determinism()
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("The lavaan package is required for confirmatory factor analysis.")
  }
}

run_week_wlsmv_cfa <- function(analysis_set, instrument, weeks) {
  .prepare_lavaan()
  items <- if (instrument == "PHQ-8") PHQ8_ITEMS else PHQ9_ITEMS
  require_columns(analysis_set, c("studyWeek", items), "CFA input")
  model_list <- cfa_models(instrument)
  fit_keys <- c(
    "chisq.scaled", "df.scaled", "pvalue.scaled",
    "cfi.scaled", "tli.scaled", "rmsea.scaled",
    "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "srmr"
  )

  rows <- list()
  for (week in weeks) {
    df_week <- analysis_set[analysis_set$studyWeek == week, items, drop = FALSE]
    df_week <- df_week[stats::complete.cases(df_week), , drop = FALSE]

    for (model_name in names(model_list)) {
      fit <- lavaan::cfa(
        model_list[[model_name]],
        data = df_week,
        ordered = items,
        estimator = "WLSMV",
        parameterization = "theta",
        missing = "listwise"
      )
      fit_measures <- lavaan::fitMeasures(fit, fit_keys)
      rows[[paste(instrument, week, model_name, sep = "__")]] <- cbind(
        data.frame(
          instrument = instrument,
          Week = week,
          Model = model_name,
          N = nrow(df_week),
          stringsAsFactors = FALSE
        ),
        as.data.frame(as.list(fit_measures), check.names = FALSE)
      )
    }
  }

  fits <- do.call(rbind, rows)
  fits$Delta_CFI <- NA_real_
  fits$Delta_RMSEA <- NA_real_
  fits$Delta_SRMR <- NA_real_
  for (week in weeks) {
    primary <- fits[fits$Week == week & fits$Model == "2-Factor", ]
    week_rows <- fits$Week == week
    fits$Delta_CFI[week_rows] <- fits$cfi.scaled[week_rows] - primary$cfi.scaled[[1]]
    fits$Delta_RMSEA[week_rows] <- fits$rmsea.scaled[week_rows] - primary$rmsea.scaled[[1]]
    fits$Delta_SRMR[week_rows] <- fits$srmr[week_rows] - primary$srmr[[1]]
    fits[
      week_rows & fits$Model == "2-Factor",
      c("Delta_CFI", "Delta_RMSEA", "Delta_SRMR")
    ] <- NA_real_
  }
  rownames(fits) <- NULL
  fits
}

run_efa_phq9 <- function(
    analysis_set,
    weeks = c(0, 2, 4),
    n_factors_list = c(1, 2, 3)) {
  configure_determinism()
  if (!requireNamespace("psych", quietly = TRUE)) {
    stop("The psych package is required for exploratory factor analysis.")
  }
  require_columns(analysis_set, c("studyWeek", PHQ9_ITEMS), "EFA input")
  phq_cols <- PHQ9_ITEMS
  canon_sets <- lapply(
    PHQ_MODEL_SPECS[["PHQ-9"]][c("2-Factor", "3-Factor")],
    `[[`,
    "factors"
  )

  per_week <- lapply(weeks, function(w) {
    .efa_one_week(analysis_set, w, phq_cols, n_factors_list, canon_sets)
  })

  list(
    loadings = dplyr::bind_rows(lapply(per_week, `[[`, "loadings")),
    parallel = dplyr::bind_rows(lapply(per_week, `[[`, "parallel"))
  )
}

.efa_one_week <- function(analysis_set, week, phq_cols, n_factors_list, canon_sets) {
  df_week <- analysis_set[analysis_set$studyWeek == week, phq_cols, drop = FALSE]
  df_week <- df_week[stats::complete.cases(df_week), , drop = FALSE]
  n_used <- nrow(df_week)

  collect <- function(cnd) {
    invokeRestart("muffleWarning")
  }

  rho <- withCallingHandlers(
    psych::polychoric(df_week),
    warning = collect
  )$rho
  parallel_analysis <- withCallingHandlers(
    suppressMessages(psych::fa.parallel(
      df_week,
      cor = "poly",
      fm = "pa",
      fa = "fa",
      n.iter = 100,
      plot = FALSE
    )),
    warning = collect
  )

  loadings <- dplyr::bind_rows(lapply(n_factors_list, function(k) {
    fitted <- withCallingHandlers(
      suppressMessages(psych::fa(
        r = rho,
        n.obs = n_used,
        nfactors = k,
        rotate = "promax",
        fm = "pa"
      )),
      warning = collect
    )
    .canonicalize_efa_loadings(fitted, k, week, canon_sets)
  }))

  resampled <- if (!is.null(parallel_analysis$fa.simr)) {
    parallel_analysis$fa.simr
  } else {
    parallel_analysis$fa.sim
  }
  parallel <- data.frame(
    Week = as.numeric(week),
    Factor = seq_along(parallel_analysis$fa.values),
    Observed = unname(parallel_analysis$fa.values),
    Resampled = unname(resampled),
    stringsAsFactors = FALSE
  )

  list(loadings = loadings, parallel = parallel)
}

.canonicalize_efa_loadings <- function(fitted, n_factors, week, canon_sets) {
  loadings <- unclass(fitted$loadings)
  items <- rownames(loadings)
  model <- paste0(n_factors, "-Factor")

  for (j in seq_len(n_factors)) {
    if (loadings[which.max(abs(loadings[, j])), j] < 0) {
      loadings[, j] <- -loadings[, j]
    }
  }

  canonical_match <- TRUE
  n_matched <- NA_integer_
  if (n_factors > 1L) {
    canon <- canon_sets[[model]]
    permutations <- .perms(seq_len(n_factors))
    score <- vapply(permutations, function(permutation) {
      sum(vapply(seq_len(n_factors), function(canonical_factor) {
        mean(abs(loadings[
          intersect(canon[[canonical_factor]], items),
          permutation[canonical_factor]
        ]))
      }, numeric(1)))
    }, numeric(1))
    loadings <- loadings[, permutations[[which.max(score)]], drop = FALSE]
    dominant_factor <- max.col(abs(loadings), ties.method = "first")
    canonical_factor <- integer(length(items))
    names(canonical_factor) <- items
    for (j in seq_len(n_factors)) {
      canonical_factor[intersect(canon[[j]], items)] <- j
    }
    n_matched <- sum(dominant_factor == canonical_factor)
    canonical_match <- n_matched / length(items) >= 0.75
  }

  loading_matrix <- matrix(
    NA_real_,
    nrow = length(items),
    ncol = 3,
    dimnames = list(NULL, c("F1", "F2", "F3"))
  )
  loading_matrix[, seq_len(n_factors)] <- loadings
  data.frame(
    Week = as.numeric(week),
    Model = model,
    Item = items,
    Item_label = item_label(items),
    F1 = loading_matrix[, 1],
    F2 = loading_matrix[, 2],
    F3 = loading_matrix[, 3],
    canonical_match = canonical_match,
    canonical_n_matched = n_matched,
    stringsAsFactors = FALSE
  )
}

.perms <- function(x) {
  if (length(x) == 1L) {
    return(list(x))
  }
  do.call(c, lapply(seq_along(x), function(i) {
    lapply(.perms(x[-i]), function(permutation) c(x[i], permutation))
  }))
}

compute_phq9_2f_ordinal_omega <- function(analysis_set, weeks = c(1, 3)) {
  .prepare_lavaan()
  configure_determinism()
  if (!requireNamespace("semTools", quietly = TRUE)) {
    stop("The semTools package is required for factor-level ordinal omega.")
  }
  items <- PHQ9_ITEMS
  require_columns(analysis_set, c("studyWeek", items), "ordinal omega input")
  model <- cfa_models("PHQ-9")[["2-Factor"]]

  rows <- lapply(weeks, function(week) {
    df_week <- analysis_set[analysis_set$studyWeek == week, items, drop = FALSE]
    df_week <- df_week[stats::complete.cases(df_week), , drop = FALSE]
    fit <- lavaan::cfa(
      model,
      data = df_week,
      ordered = items,
      estimator = "WLSMV",
      parameterization = "theta",
      missing = "listwise"
    )
    reliability <- suppressWarnings(suppressMessages(semTools::reliability(fit)))
    omega <- reliability["omega3", ]
    data.frame(
      Week = week,
      N = nrow(df_week),
      Factor = names(omega),
      omega = unname(omega),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
