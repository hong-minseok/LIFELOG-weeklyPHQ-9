# Shared definitions for the public statistical-analysis code.

ANCHOR_WEEKS <- c(0, 4, 8, 12)
STACKED_LABEL <- paste("Stacked", paste(ANCHOR_WEEKS, collapse = "/"))
PHQ9_ITEMS <- paste0("phq9_", 1:9)
PHQ8_ITEMS <- PHQ9_ITEMS[1:8]
GAD2_ITEMS <- paste0("gad2_", 1:2)
CLUSTER_BOOTSTRAP_REPLICATES <- 2000L
CLUSTER_BOOTSTRAP_LABEL <- sprintf(
  "Participant-level cluster bootstrap (B=%d)",
  CLUSTER_BOOTSTRAP_REPLICATES
)

ANALYSIS_FILTERS <- list(
  min_phq_weeks = 4,
  min_paired_weeks = 2
)

MINIMUM_INPUT_COLUMNS <- list(
  weekly = c("studyID", "studyWeek", "phq9_sum", PHQ9_ITEMS),
  state = c("studyID", "studyWeek", "madrs_sum"),
  base = c("studyID", "cohort_source", "mde_base")
)

PHQ_MODEL_SPECS <- list(
  "PHQ-9" = list(
    "1-Factor" = list(
      factors = list(General = PHQ9_ITEMS),
      residual_pairs = list(c("phq9_6", "phq9_9"))
    ),
    "2-Factor" = list(
      factors = list(
        Affective = PHQ9_ITEMS[c(1, 2, 6, 7, 8, 9)],
        Somatic = PHQ9_ITEMS[c(3, 4, 5)]
      ),
      residual_pairs = list(PHQ9_ITEMS[c(6, 9)], PHQ9_ITEMS[c(7, 8)])
    ),
    "3-Factor" = list(
      factors = list(
        Affective = PHQ9_ITEMS[c(1, 2, 6, 9)],
        Somatic = PHQ9_ITEMS[c(3, 4, 5)],
        Cognitive = PHQ9_ITEMS[c(7, 8)]
      ),
      residual_pairs = list(PHQ9_ITEMS[c(6, 9)])
    )
  ),
  "PHQ-8" = list(
    "2-Factor" = list(
      factors = list(
        Affective = PHQ8_ITEMS[c(1, 2, 6, 7, 8)],
        Somatic = PHQ8_ITEMS[c(3, 4, 5)]
      ),
      residual_pairs = list(PHQ8_ITEMS[c(7, 8)])
    )
  )
)

ITEM_LABELS <- stats::setNames(paste("Item", 1:9), PHQ9_ITEMS)

require_columns <- function(data, columns, context = "input") {
  if (!is.data.frame(data)) {
    stop(context, " must be a data frame.")
  }
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(context, " is missing: ", paste(missing_columns, collapse = ", "))
  }
  invisible(data)
}

coalesce_fit_measure <- function(scaled, unscaled) {
  ifelse(!is.na(scaled), scaled, unscaled)
}

validate_analysis_input <- function(df, role = c("weekly", "state", "base")) {
  role <- match.arg(role)
  if (!is.data.frame(df)) {
    stop(role, " input must be a data frame.")
  }

  missing_columns <- setdiff(MINIMUM_INPUT_COLUMNS[[role]], names(df))
  if (length(missing_columns) > 0L) {
    stop(role, " input is missing: ", paste(missing_columns, collapse = ", "))
  }
  if (!is.character(df$studyID)) {
    stop(role, " input: studyID must be character.")
  }
  if (role != "base" && !is.numeric(df$studyWeek)) {
    stop(role, " input: studyWeek must be numeric.")
  }

  keys <- if (role == "base") "studyID" else c("studyID", "studyWeek")
  if (anyDuplicated(df[, keys, drop = FALSE])) {
    stop(role, " input contains duplicate rows for ", paste(keys, collapse = " + "), ".")
  }

  required <- MINIMUM_INPUT_COLUMNS[[role]]
  all_missing <- required[vapply(required, function(column) {
    all(is.na(df[[column]]))
  }, logical(1))]
  if (length(all_missing) > 0L) {
    stop(role, " input has entirely missing columns: ", paste(all_missing, collapse = ", "))
  }

  if (role == "weekly") {
    items <- PHQ9_ITEMS
    non_numeric <- items[!vapply(items, function(column) {
      is.numeric(df[[column]])
    }, logical(1))]
    if (length(non_numeric) > 0L) {
      stop("weekly input has nonnumeric PHQ items: ", paste(non_numeric, collapse = ", "))
    }
    non_integer <- items[vapply(items, function(column) {
      values <- df[[column]][!is.na(df[[column]])]
      any(values != floor(values))
    }, logical(1))]
    if (length(non_integer) > 0L) {
      stop("weekly input has noninteger PHQ items: ", paste(non_integer, collapse = ", "))
    }
    outside_item_range <- vapply(items, function(column) {
      any(df[[column]] < 0 | df[[column]] > 3, na.rm = TRUE)
    }, logical(1))
    if (any(outside_item_range)) {
      stop(
        "weekly input has PHQ items outside 0-3: ",
        paste(items[outside_item_range], collapse = ", ")
      )
    }
    if (any(df$phq9_sum < 0 | df$phq9_sum > 27, na.rm = TRUE)) {
      stop("weekly input has phq9_sum outside 0-27.")
    }
  }
  if (role == "base") {
    observed_sources <- unique(as.character(df$cohort_source[!is.na(df$cohort_source)]))
    unexpected_sources <- setdiff(observed_sources, c("Clinic", "Community"))
    if (length(unexpected_sources) > 0L) {
      stop(
        "base input has unexpected cohort_source values: ",
        paste(unexpected_sources, collapse = ", ")
      )
    }
  }
  invisible(df)
}

apply_eligibility_criteria <- function(
    weekly,
    state,
    min_phq = ANALYSIS_FILTERS$min_phq_weeks,
    min_paired = ANALYSIS_FILTERS$min_paired_weeks) {
  validate_analysis_input(weekly, "weekly")
  validate_analysis_input(state, "state")

  valid_phq <- weekly |>
    dplyr::filter(!is.na(.data$phq9_sum)) |>
    dplyr::group_by(.data$studyID) |>
    dplyr::summarise(
      n_phq = dplyr::n_distinct(.data$studyWeek),
      .groups = "drop"
    ) |>
    dplyr::filter(.data$n_phq >= min_phq)

  paired <- dplyr::inner_join(
    weekly |>
      dplyr::select(dplyr::all_of(c("studyID", "studyWeek", "phq9_sum"))) |>
      dplyr::filter(!is.na(.data$phq9_sum)),
    state |>
      dplyr::select(dplyr::all_of(c("studyID", "studyWeek", "madrs_sum"))) |>
      dplyr::filter(!is.na(.data$madrs_sum)),
    by = c("studyID", "studyWeek")
  ) |>
    dplyr::group_by(.data$studyID) |>
    dplyr::summarise(
      n_paired = dplyr::n_distinct(.data$studyWeek),
      .groups = "drop"
    ) |>
    dplyr::filter(.data$n_paired >= min_paired)

  valid_ids <- intersect(valid_phq$studyID, paired$studyID)
  weekly |>
    dplyr::filter(.data$studyID %in% valid_ids)
}

configure_determinism <- function(seed = 20250823L) {
  set.seed(seed)
  options(mc.cores = 1L)
  Sys.setenv(OMP_NUM_THREADS = "1")
  initialize_lavaan_cpu_fallback()
  invisible(TRUE)
}

initialize_lavaan_cpu_fallback <- function() {
  # lavaan 0.6-21 initializes ncpus through detectCores(); use one core if unavailable.
  detected <- suppressWarnings(parallel::detectCores())
  if (length(detected) == 1L && is.finite(detected) && detected >= 1L) {
    return(invisible(FALSE))
  }

  parallel_ns <- asNamespace("parallel")
  original <- get("detectCores", parallel_ns)
  unlockBinding("detectCores", parallel_ns)
  assign("detectCores", function(...) 1L, envir = parallel_ns)
  lockBinding("detectCores", parallel_ns)
  on.exit({
    unlockBinding("detectCores", parallel_ns)
    assign("detectCores", original, envir = parallel_ns)
    lockBinding("detectCores", parallel_ns)
  }, add = TRUE)
  get("lav_options_default", asNamespace("lavaan"))()
  invisible(TRUE)
}

item_label <- function(x) {
  key <- as.character(x)
  item_number <- suppressWarnings(as.integer(key))
  norm <- rep(NA_character_, length(key))
  norm[key %in% PHQ9_ITEMS] <- key[key %in% PHQ9_ITEMS]
  numeric_key <- !is.na(item_number) & item_number %in% seq_along(PHQ9_ITEMS)
  norm[numeric_key] <- PHQ9_ITEMS[item_number[numeric_key]]
  out <- unname(ITEM_LABELS[norm])
  if (any(is.na(out))) {
    stop("item_label(): unmapped item(s): ", paste(x[is.na(out)], collapse = ", "))
  }
  out
}

lav_wide_var <- function(item, week) {
  paste0(item, "_w", week)
}

lav_loading_line <- function(factor_name, item_vec, week, mode, model_tag) {
  lhs <- paste0(factor_name, "_w", week)
  rhs <- character(length(item_vec))
  for (i in seq_along(item_vec)) {
    item <- item_vec[[i]]
    variable <- lav_wide_var(item, week)
    if (i == 1L) {
      rhs[[i]] <- paste0("1*", variable)
    } else if (mode == "Configural") {
      rhs[[i]] <- variable
    } else {
      label <- paste0("l_", model_tag, "_", factor_name, "_", item)
      rhs[[i]] <- paste0(label, "*", variable)
    }
  }
  paste0(lhs, " =~ ", paste(rhs, collapse = " + "))
}

lav_within_residuals_wide <- function(residual_pairs, week) {
  vapply(residual_pairs, function(pair) {
    paste0(lav_wide_var(pair[[1]], week), " ~~ ", lav_wide_var(pair[[2]], week))
  }, character(1))
}

lav_time_residuals_wide <- function(items, weeks) {
  lines <- character()
  for (item in items) {
    pairs <- utils::combn(weeks, 2)
    for (i in seq_len(ncol(pairs))) {
      lines <- c(
        lines,
        paste0(lav_wide_var(item, pairs[1, i]), " ~~ ", lav_wide_var(item, pairs[2, i]))
      )
    }
  }
  lines
}

lav_thresholds_wide <- function(items, weeks) {
  lines <- character()
  for (item in items) {
    for (week in weeks) {
      variable <- lav_wide_var(item, week)
      lines <- c(lines, paste0(
        variable, " | ",
        "th_", item, "_1*t1 + ",
        "th_", item, "_2*t2 + ",
        "th_", item, "_3*t3"
      ))
    }
  }
  lines
}

lav_latent_means_wide <- function(factor_names, weeks) {
  lines <- character()
  for (factor_name in factor_names) {
    lines <- c(lines, paste0(factor_name, "_w", weeks[[1]], " ~ 0*1"))
    for (week in weeks[-1]) {
      lines <- c(lines, paste0(factor_name, "_w", week, " ~ NA*1"))
    }
  }
  lines
}

lav_wide_syntax <- function(model_spec, items, weeks, mode, model_tag) {
  # In Threshold mode, week-free labels impose equal loadings and thresholds over time.
  lines <- character()
  for (week in weeks) {
    for (factor_name in names(model_spec$factors)) {
      lines <- c(
        lines,
        lav_loading_line(
          factor_name,
          model_spec$factors[[factor_name]],
          week,
          mode,
          model_tag
        )
      )
    }
    lines <- c(lines, lav_within_residuals_wide(model_spec$residual_pairs, week))
  }
  lines <- c(lines, lav_time_residuals_wide(items, weeks))
  if (mode == "Threshold") {
    lines <- c(
      lines,
      lav_thresholds_wide(items, weeks),
      lav_latent_means_wide(names(model_spec$factors), weeks)
    )
  }
  paste(lines, collapse = "\n")
}

lav_long_syntax <- function(model_spec) {
  lines <- character()
  for (factor_name in names(model_spec$factors)) {
    lines <- c(
      lines,
      paste0(factor_name, " =~ ", paste(model_spec$factors[[factor_name]], collapse = " + "))
    )
  }
  lines <- c(lines, vapply(model_spec$residual_pairs, function(pair) {
    paste0(pair[[1]], " ~~ ", pair[[2]])
  }, character(1)))
  paste(lines, collapse = "\n")
}

lav_mlcfa_syntax <- function(model_spec) {
  paste0(
    "level: 1\n", lav_long_syntax(model_spec),
    "\nlevel: 2\n", lav_long_syntax(model_spec)
  )
}

cfa_models <- function(instrument) {
  lapply(PHQ_MODEL_SPECS[[instrument]], lav_long_syntax)
}

capture_with_warnings <- function(expr) {
  warnings <- character()
  result <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  list(result = result, warnings = unique(warnings))
}

paste_condition_messages <- function(msgs) {
  messages <- unlist(msgs, use.names = FALSE)
  messages <- messages[!is.na(messages) & nzchar(messages)]
  if (length(messages) == 0L) {
    return("")
  }
  messages <- gsub("[\r\n]+", " ", messages)
  messages <- gsub("[[:space:]]+", " ", messages)
  paste(unique(trimws(messages)), collapse = " | ")
}
