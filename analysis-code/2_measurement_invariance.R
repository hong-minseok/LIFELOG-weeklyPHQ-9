# Longitudinal and source measurement-invariance analyses.

nested_delta <- function(current, baseline) {
  current - baseline
}

make_wlsmv_model_specs_phq9 <- function() {
  PHQ_MODEL_SPECS[["PHQ-9"]][c("2-Factor", "3-Factor")]
}

make_wlsmv_model_specs_phq8 <- function() {
  PHQ_MODEL_SPECS[["PHQ-8"]]["2-Factor"]
}

build_wide_complete_case <- function(df_long, items, primary_weeks = ANCHOR_WEEKS) {
  df_wide <- df_long |>
    dplyr::filter(.data$studyWeek %in% primary_weeks) |>
    dplyr::select(dplyr::all_of(c("studyID", "studyWeek", items))) |>
    tidyr::pivot_wider(
      id_cols = dplyr::all_of("studyID"),
      names_from = dplyr::all_of("studyWeek"),
      values_from = dplyr::all_of(items),
      names_glue = "{.value}_w{studyWeek}"
    ) |>
    dplyr::arrange(.data$studyID)

  ordered_vars <- as.vector(outer(items, paste0("_w", primary_weeks), paste0))
  df_wide_cc <- df_wide |>
    dplyr::filter(stats::complete.cases(dplyr::across(dplyr::all_of(ordered_vars))))
  list(df_wide_cc = df_wide_cc, ordered_vars = ordered_vars)
}

is_weekly_cfa_estimable <- function(
    converged,
    post_check,
    max_abs_loading,
    cfi,
    rmsea,
    rmsea_low,
    rmsea_high) {
  estimable <- as.logical(converged) &
    as.logical(post_check) &
    !(is.finite(max_abs_loading) & max_abs_loading >= 0.999) &
    is.finite(cfi) &
    is.finite(rmsea) &
    rmsea > 0 &
    is.finite(rmsea_low) &
    is.finite(rmsea_high)
  estimable[is.na(estimable)] <- FALSE
  estimable
}

assess_weekly_cfa_estimability <- function(
    weekly_fit,
    weekly_loadings,
    model = "2-Factor") {
  loading_check <- weekly_loadings |>
    dplyr::filter(.data$Model == model) |>
    dplyr::group_by(.data$Week, .data$Model) |>
    dplyr::summarise(
      max_abs_loading = max(abs(.data$Loading), na.rm = TRUE),
      .groups = "drop"
    )

  weekly_fit |>
    dplyr::filter(.data$Model == model) |>
    dplyr::left_join(loading_check, by = c("Week", "Model")) |>
    dplyr::mutate(
      Estimable = is_weekly_cfa_estimable(
        .data$Converged,
        .data$Post_Check,
        .data$max_abs_loading,
        .data$CFI_report,
        .data$RMSEA_report,
        .data$RMSEA_low,
        .data$RMSEA_high
      )
    ) |>
    dplyr::arrange(.data$Week)
}

run_ordinal_wlsmv_longitudinal <- function(
    df_long,
    items,
    model_specs,
    primary_weeks = ANCHOR_WEEKS,
    shadow_weeks = 0:24) {
  configure_determinism()
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("The lavaan package is required for longitudinal invariance analysis.")
  }
  require_columns(
    df_long,
    c("studyID", "studyWeek", items),
    "longitudinal invariance input"
  )

  df_long <- df_long |>
    dplyr::filter(.data$studyWeek <= max(shadow_weeks)) |>
    dplyr::mutate(studyID = as.character(.data$studyID))

  fit_keys <- c(
    "chisq.scaled", "df.scaled", "pvalue.scaled",
    "cfi.scaled", "tli.scaled", "rmsea.scaled",
    "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "srmr",
    "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper"
  )

  extract_fit_row <- function(
      fit,
      analysis,
      model,
      level,
      n_used,
      warnings = character(),
      error = NA_character_) {
    fit_warnings <- character()
    with_fit_warnings <- function(expr, on_error) {
      tryCatch(
        withCallingHandlers(expr, warning = function(w) {
          fit_warnings <<- c(fit_warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }),
        error = on_error
      )
    }

    out <- data.frame(
      Analysis = analysis,
      Model = model,
      Level = level,
      N = n_used,
      Converged = FALSE,
      Post_Check = NA,
      Warning = paste_condition_messages(warnings),
      Error = error,
      stringsAsFactors = FALSE
    )
    for (key in fit_keys) {
      out[[key]] <- NA_real_
    }

    if (inherits(fit, "lavaan")) {
      fit_measures <- with_fit_warnings(
        lavaan::fitMeasures(fit, fit_keys),
        function(e) rep(NA_real_, length(fit_keys))
      )
      out$Converged <- isTRUE(with_fit_warnings(
        lavaan::lavInspect(fit, "converged"),
        function(e) FALSE
      ))
      out$Post_Check <- as.character(with_fit_warnings(
        lavaan::lavInspect(fit, "post.check"),
        function(e) NA
      ))
      for (key in fit_keys) {
        out[[key]] <- unname(fit_measures[key])
      }
      out$Warning <- paste_condition_messages(c(out$Warning, fit_warnings))
    }
    out
  }

  add_deltas <- function(df) {
    if (is.null(df) || nrow(df) == 0L) {
      return(df)
    }
    df |>
      dplyr::group_by(.data$Analysis, .data$Model) |>
      dplyr::arrange(match(.data$Level, c("Configural", "Threshold")), .by_group = TRUE) |>
      dplyr::mutate(
        CFI_report = coalesce_fit_measure(.data$cfi.scaled, .data$cfi),
        RMSEA_report = coalesce_fit_measure(.data$rmsea.scaled, .data$rmsea),
        RMSEA_low = coalesce_fit_measure(
          .data$rmsea.ci.lower.scaled,
          .data$rmsea.ci.lower
        ),
        RMSEA_high = coalesce_fit_measure(
          .data$rmsea.ci.upper.scaled,
          .data$rmsea.ci.upper
        ),
        Delta_CFI = c(NA_real_, nested_delta(.data$CFI_report[-1], .data$CFI_report[1])),
        Delta_RMSEA = c(
          NA_real_,
          nested_delta(.data$RMSEA_report[-1], .data$RMSEA_report[1])
        ),
        Delta_SRMR = c(NA_real_, diff(.data$srmr))
      ) |>
      dplyr::ungroup()
  }

  wide_cc <- build_wide_complete_case(df_long, items, primary_weeks)
  df_wide_cc <- wide_cc$df_wide_cc
  primary_ordered <- wide_cc$ordered_vars
  primary_n <- nrow(df_wide_cc)

  sparse_categories <- df_wide_cc |>
    dplyr::select(dplyr::all_of(primary_ordered)) |>
    tidyr::pivot_longer(
      dplyr::everything(),
      names_to = "Variable",
      values_to = "Score"
    ) |>
    dplyr::mutate(
      Item = sub("_w[0-9]+$", "", .data$Variable),
      Week = as.integer(sub("^.*_w", "", .data$Variable)),
      Score = factor(.data$Score, levels = 0:3)
    ) |>
    dplyr::count(.data$Item, .data$Week, .data$Score, .drop = FALSE) |>
    tidyr::pivot_wider(
      names_from = dplyr::all_of("Score"),
      values_from = dplyr::all_of("n"),
      names_prefix = "N_"
    ) |>
    dplyr::mutate(
      Min_Cell_N = pmin(.data$N_0, .data$N_1, .data$N_2, .data$N_3),
      Present_Categories = (.data$N_0 > 0) + (.data$N_1 > 0) +
        (.data$N_2 > 0) + (.data$N_3 > 0)
    ) |>
    dplyr::arrange(.data$Week, .data$Item)

  primary_sparse_vars <- sparse_categories |>
    dplyr::filter(.data$Present_Categories < 4) |>
    dplyr::mutate(Variable = paste0(.data$Item, "_w", .data$Week)) |>
    dplyr::pull("Variable")
  primary_sparse_error <- if (length(primary_sparse_vars) > 0L) {
    paste0(
      "Threshold invariance skipped because these primary item-weeks do not ",
      "contain all four response categories: ",
      paste(primary_sparse_vars, collapse = ", ")
    )
  } else {
    NA_character_
  }

  primary_rows <- list()
  for (model_name in names(model_specs)) {
    for (level in c("Configural", "Threshold")) {
      key <- paste(model_name, level, sep = "__")
      if (level == "Threshold" && !is.na(primary_sparse_error)) {
        primary_rows[[key]] <- extract_fit_row(
          NULL,
          "Primary 0-12 wide WLSMV",
          model_name,
          level,
          primary_n,
          error = primary_sparse_error
        )
      } else {
        fit_obj <- capture_with_warnings(lavaan::cfa(
          lav_wide_syntax(
            model_specs[[model_name]],
            items,
            primary_weeks,
            level,
            gsub("[^A-Za-z0-9]", "", model_name)
          ),
          data = df_wide_cc,
          ordered = primary_ordered,
          estimator = "WLSMV",
          parameterization = "theta",
          missing = "listwise"
        ))
        if (inherits(fit_obj$result, "error")) {
          primary_rows[[key]] <- extract_fit_row(
            NULL,
            "Primary 0-12 wide WLSMV",
            model_name,
            level,
            primary_n,
            fit_obj$warnings,
            conditionMessage(fit_obj$result)
          )
        } else {
          primary_rows[[key]] <- extract_fit_row(
            fit_obj$result,
            "Primary 0-12 wide WLSMV",
            model_name,
            level,
            primary_n,
            fit_obj$warnings
          )
        }
      }
      primary_rows[[key]]$Sparse_Category_Flag <- length(primary_sparse_vars) > 0L
    }
  }
  primary_fit <- add_deltas(dplyr::bind_rows(primary_rows))

  weekly_fit_rows <- list()
  weekly_loading_rows <- list()
  for (week in shadow_weeks) {
    df_week <- df_long |>
      dplyr::filter(.data$studyWeek == week) |>
      dplyr::select(dplyr::all_of(items)) |>
      dplyr::filter(stats::complete.cases(dplyr::across(dplyr::everything())))

    for (model_name in names(model_specs)) {
      fit_obj <- capture_with_warnings(lavaan::cfa(
        lav_long_syntax(model_specs[[model_name]]),
        data = df_week,
        ordered = items,
        estimator = "WLSMV",
        parameterization = "theta",
        missing = "listwise"
      ))
      key <- paste(week, model_name, sep = "__")
      if (inherits(fit_obj$result, "error")) {
        weekly_fit_rows[[key]] <- extract_fit_row(
          NULL,
          "Week-specific WLSMV CFA",
          model_name,
          paste0("Week_", week),
          nrow(df_week),
          fit_obj$warnings,
          conditionMessage(fit_obj$result)
        )
      } else {
        fit <- fit_obj$result
        fit_row <- extract_fit_row(
          fit,
          "Week-specific WLSMV CFA",
          model_name,
          paste0("Week_", week),
          nrow(df_week),
          fit_obj$warnings
        )
        parameter_obj <- capture_with_warnings(
          lavaan::parameterEstimates(fit, standardized = TRUE)
        )
        if (inherits(parameter_obj$result, "error")) {
          fit_row$Warning <- paste_condition_messages(c(
            fit_row$Warning,
            paste0("Parameter extraction failed: ", conditionMessage(parameter_obj$result))
          ))
        } else {
          fit_row$Warning <- paste_condition_messages(c(fit_row$Warning, parameter_obj$warnings))
          weekly_loading_rows[[key]] <- parameter_obj$result |>
            dplyr::filter(.data$op == "=~") |>
            dplyr::transmute(
              Week = week,
              Model = model_name,
              Factor = .data$lhs,
              Item = .data$rhs,
              Loading = .data$std.all
            )
        }
        weekly_fit_rows[[key]] <- fit_row
      }
    }
  }

  weekly_fit <- dplyr::bind_rows(weekly_fit_rows) |>
    dplyr::mutate(
      Week = as.integer(sub("Week_", "", .data$Level)),
      CFI_report = coalesce_fit_measure(.data$cfi.scaled, .data$cfi),
      RMSEA_report = coalesce_fit_measure(.data$rmsea.scaled, .data$rmsea),
      RMSEA_low = coalesce_fit_measure(
        .data$rmsea.ci.lower.scaled,
        .data$rmsea.ci.lower
      ),
      RMSEA_high = coalesce_fit_measure(
        .data$rmsea.ci.upper.scaled,
        .data$rmsea.ci.upper
      )
    )
  weekly_loadings <- dplyr::bind_rows(weekly_loading_rows)
  if (nrow(weekly_loadings) == 0L) {
    weekly_loadings <- data.frame(
      Week = integer(),
      Model = character(),
      Factor = character(),
      Item = character(),
      Loading = numeric(),
      stringsAsFactors = FALSE
    )
  }

  list(
    primary_fit = primary_fit,
    weekly_fit = weekly_fit,
    weekly_loadings = weekly_loadings,
    weekly_estimability = assess_weekly_cfa_estimability(weekly_fit, weekly_loadings)
  )
}

# Clinic-versus-community ordinal measurement-invariance analyses.

.GROUP_LEVELS <- c("Configural", "Threshold")

.group_equal_for <- function(level) {
  switch(
    level,
    Configural = character(0),
    # Threshold invariance for ordinal items also constrains factor loadings.
    Threshold = c("thresholds", "loadings"),
    stop("Unknown invariance level: ", level)
  )
}

.group_fit_indices <- function(fit) {
  keys <- c(
    "cfi.scaled", "tli.scaled", "rmsea.scaled",
    "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "srmr"
  )
  fit_measures <- tryCatch(
    lavaan::fitMeasures(fit, keys),
    error = function(e) rep(NA_real_, length(keys))
  )
  stats::setNames(
    as.numeric(fit_measures),
    c("CFI", "TLI", "RMSEA", "RMSEA_low", "RMSEA_high", "SRMR")
  )
}

.group_fit_level <- function(group_frame, level, items, use_cluster, model_spec) {
  args <- list(
    configural.model = lav_long_syntax(model_spec),
    data = group_frame,
    ordered = items,
    parameterization = "theta",
    ID.fac = "std.lv",
    ID.cat = "Wu.Estabrook.2016",
    group = "cohort_source",
    group.equal = .group_equal_for(level),
    estimator = "WLSMV",
    missing = "listwise",
    return.fit = TRUE
  )
  if (use_cluster) {
    args$cluster <- "studyID"
  }
  capture_with_warnings(do.call(semTools::measEq.syntax, args))
}

.group_fit_row <- function(fit_obj, level) {
  row <- data.frame(
    Level = level,
    Converged = FALSE,
    Post_Check = NA,
    Warning = paste_condition_messages(fit_obj$warnings),
    Error = NA_character_,
    CFI = NA_real_,
    TLI = NA_real_,
    RMSEA = NA_real_,
    RMSEA_low = NA_real_,
    RMSEA_high = NA_real_,
    SRMR = NA_real_,
    stringsAsFactors = FALSE
  )
  fit <- fit_obj$result
  if (inherits(fit, "error")) {
    row$Error <- conditionMessage(fit)
    return(row)
  }
  if (inherits(fit, "lavaan")) {
    row$Converged <- isTRUE(tryCatch(
      lavaan::lavInspect(fit, "converged"),
      error = function(e) FALSE
    ))
    row$Post_Check <- as.character(tryCatch(
      lavaan::lavInspect(fit, "post.check"),
      error = function(e) NA
    ))
    fit_indices <- .group_fit_indices(fit)
    row$CFI <- fit_indices[["CFI"]]
    row$TLI <- fit_indices[["TLI"]]
    row$RMSEA <- fit_indices[["RMSEA"]]
    row$RMSEA_low <- fit_indices[["RMSEA_low"]]
    row$RMSEA_high <- fit_indices[["RMSEA_high"]]
    row$SRMR <- fit_indices[["SRMR"]]
  }
  row
}

run_group_invariance_one <- function(group_frame, items, use_cluster, model_spec) {
  rows <- list()
  for (level in .GROUP_LEVELS) {
    fit_obj <- .group_fit_level(group_frame, level, items, use_cluster, model_spec)
    rows[[level]] <- .group_fit_row(fit_obj, level)
  }
  fit_table <- dplyr::bind_rows(rows) |>
    dplyr::arrange(match(.data$Level, .GROUP_LEVELS)) |>
    dplyr::mutate(
      Delta_CFI = c(NA_real_, nested_delta(.data$CFI[-1], .data$CFI[1])),
      Delta_RMSEA = c(NA_real_, nested_delta(.data$RMSEA[-1], .data$RMSEA[1]))
    )

  list(fit_table = fit_table)
}

.group_join_source <- function(frame, base_raw) {
  source <- base_raw |>
    dplyr::transmute(
      studyID = as.character(.data$studyID),
      cohort_source = .data$cohort_source
    )
  frame |>
    dplyr::mutate(studyID = as.character(.data$studyID)) |>
    dplyr::left_join(source, by = "studyID") |>
    dplyr::filter(!is.na(.data$cohort_source))
}

build_group_stacked_frame <- function(
    long_frame,
    base_raw,
    items = PHQ9_ITEMS,
    weeks = ANCHOR_WEEKS) {
  long_frame |>
    dplyr::filter(.data$studyWeek %in% weeks) |>
    dplyr::select(dplyr::all_of(c("studyID", "studyWeek", items))) |>
    .group_join_source(base_raw) |>
    dplyr::filter(stats::complete.cases(dplyr::across(dplyr::all_of(items)))) |>
    dplyr::arrange(.data$studyID, .data$studyWeek)
}

build_group_week_frame <- function(
    long_frame,
    base_raw,
    items = PHQ9_ITEMS,
    week = 1) {
  long_frame |>
    dplyr::filter(.data$studyWeek == week) |>
    dplyr::select(dplyr::all_of(c("studyID", "studyWeek", items))) |>
    .group_join_source(base_raw) |>
    dplyr::filter(stats::complete.cases(dplyr::across(dplyr::all_of(items)))) |>
    dplyr::arrange(.data$studyID)
}

.group_counts <- function(group_frame, analysis_label, cluster_robust) {
  id <- as.character(group_frame$studyID)
  source <- group_frame$cohort_source
  data.frame(
    Analysis = analysis_label,
    N_Participants = dplyr::n_distinct(id),
    N_Observations = nrow(group_frame),
    Clinic_Participants = dplyr::n_distinct(id[source == "Clinic"]),
    Community_Participants = dplyr::n_distinct(id[source == "Community"]),
    Clinic_Observations = sum(source == "Clinic"),
    Community_Observations = sum(source == "Community"),
    Cluster_Robust = cluster_robust,
    stringsAsFactors = FALSE
  )
}

.assert_group_inv_items <- function(model_spec, items, long_frame) {
  spec_items <- unique(c(
    unlist(model_spec$factors, use.names = FALSE),
    unlist(model_spec$residual_pairs, use.names = FALSE)
  ))
  missing_from_items <- setdiff(spec_items, items)
  if (length(missing_from_items) > 0L) {
    stop(
      "group invariance: model_spec references items not in `items`: ",
      paste(missing_from_items, collapse = ", ")
    )
  }
  extra_items <- setdiff(items, spec_items)
  if (length(extra_items) > 0L) {
    stop(
      "group invariance: `items` includes columns not referenced by model_spec: ",
      paste(extra_items, collapse = ", ")
    )
  }
  missing_from_data <- setdiff(items, names(long_frame))
  if (length(missing_from_data) > 0L) {
    stop(
      "group invariance: `items` not present in the data frame: ",
      paste(missing_from_data, collapse = ", ")
    )
  }
  invisible(TRUE)
}

run_group_invariance <- function(
    long_frame,
    base_raw,
    model_spec,
    items,
    weeks = ANCHOR_WEEKS,
    cross_check_weeks = c(1, 4)) {
  configure_determinism()
  require_columns(
    long_frame,
    c("studyID", "studyWeek", items),
    "source-invariance input"
  )
  .assert_group_inv_items(model_spec, items, long_frame)
  require_columns(base_raw, c("studyID", "cohort_source"), "source-invariance baseline input")
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("The lavaan package is required for source-invariance analysis.")
  }
  if (!requireNamespace("semTools", quietly = TRUE)) {
    stop("The semTools package is required for source-invariance analysis.")
  }

  stacked <- build_group_stacked_frame(long_frame, base_raw, items, weeks)
  cross_frames <- lapply(cross_check_weeks, function(week) {
    build_group_week_frame(long_frame, base_raw, items, week = week)
  })
  cross_labels <- sprintf("Week %d", cross_check_weeks)
  names(cross_frames) <- cross_labels

  # Categorical clustered WLSMV is unavailable in this environment.
  cluster_dropped <- TRUE

  stacked_result <- run_group_invariance_one(
    stacked,
    items,
    use_cluster = FALSE,
    model_spec
  )
  cross_results <- lapply(cross_frames, function(frame) {
    run_group_invariance_one(frame, items, use_cluster = FALSE, model_spec)
  })

  counts <- dplyr::bind_rows(
    .group_counts(stacked, STACKED_LABEL, cluster_robust = FALSE),
    dplyr::bind_rows(Map(
      function(frame, label) .group_counts(frame, label, cluster_robust = NA),
      cross_frames,
      cross_labels
    ))
  )

  tag_fit <- function(table, label, cluster_robust) {
    table |>
      dplyr::mutate(Analysis = label, Cluster_Robust = cluster_robust) |>
      dplyr::left_join(
        counts[
          counts$Analysis == label,
          c(
            "Analysis", "N_Participants", "N_Observations",
            "Clinic_Participants", "Community_Participants",
            "Clinic_Observations", "Community_Observations"
          )
        ],
        by = "Analysis"
      )
  }

  fit_table <- dplyr::bind_rows(
    tag_fit(stacked_result$fit_table, STACKED_LABEL, FALSE),
    dplyr::bind_rows(Map(
      function(result, label) tag_fit(result$fit_table, label, NA),
      cross_results,
      cross_labels
    ))
  )

  list(
    fit_table = fit_table,
    counts = counts,
    cluster_dropped = cluster_dropped
  )
}
