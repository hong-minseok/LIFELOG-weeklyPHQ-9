# External validity, trait-state correlations, MDE discrimination, and cutoffs.

# Harmonize baseline and follow-up flags into occasion-level concurrent MDE status.
attach_concurrent_mde <- function(state, base) {
  require_columns(
    state,
    c("studyID", "studyWeek", "madrs_sum", "mde_followup"),
    "MDE follow-up input"
  )
  require_columns(base, c("studyID", "mde_base"), "MDE baseline input")
  state_mde <- state |>
    dplyr::filter(.data$studyWeek <= 24) |>
    dplyr::transmute(
      studyID = as.character(.data$studyID),
      studyWeek = as.numeric(.data$studyWeek),
      madrs_sum = .data$madrs_sum,
      MDE = as.numeric(.data$mde_followup)
    )
  baseline_mde <- base |>
    dplyr::transmute(
      studyID = as.character(.data$studyID),
      studyWeek = 0,
      mde_base = as.numeric(.data$mde_base)
    )
  state_mde |>
    dplyr::full_join(baseline_mde, by = c("studyID", "studyWeek")) |>
    dplyr::mutate(MDE = ifelse(.data$studyWeek == 0, .data$mde_base, .data$MDE)) |>
    dplyr::select(dplyr::all_of(c("studyID", "studyWeek", "madrs_sum", "MDE")))
}

prepare_validity_frame <- function(weekly, state, base) {
  require_columns(weekly, c("studyID", "studyWeek", "phq9_sum"), "validity weekly input")
  weekly <- weekly |>
    dplyr::mutate(
      studyID = as.character(.data$studyID),
      studyWeek = as.numeric(.data$studyWeek)
    )
  mde_state <- attach_concurrent_mde(state, base)
  weekly |>
    dplyr::inner_join(mde_state, by = c("studyID", "studyWeek")) |>
    dplyr::select(dplyr::all_of(c(
      "studyID", "studyWeek", "phq9_sum", "madrs_sum", "MDE"
    )))
}

mde_auc_by_week <- function(validity_frame) {
  configure_determinism()
  require_columns(
    validity_frame,
    c("studyID", "studyWeek", "phq9_sum", "MDE"),
    "MDE AUC input"
  )
  data <- validity_frame |>
    dplyr::filter(!is.na(.data$MDE) & !is.na(.data$phq9_sum))
  weeks <- sort(unique(data$studyWeek))
  rows <- list()
  for (week in weeks) {
    subset <- data[data$studyWeek == week, ]
    if (nrow(subset) >= 10L && length(unique(subset$MDE)) > 1L) {
      roc <- pROC::roc(
        subset$MDE,
        subset$phq9_sum,
        quiet = TRUE,
        direction = "<"
      )
      ci <- .cluster_boot_auc_ci(subset$MDE, subset$phq9_sum, subset$studyID)
      rows[[as.character(week)]] <- data.frame(
        Week = week,
        AUC = as.numeric(roc$auc),
        CI_Lower = ci[1],
        CI_Upper = ci[2],
        N = nrow(subset),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

.scale_complete_sum <- function(data, variables) {
  complete <- stats::complete.cases(data[, variables, drop = FALSE])
  out <- rowSums(data[, variables, drop = FALSE])
  out[!complete] <- NA_real_
  out
}

.safe_cor_test <- function(x, y) {
  complete <- stats::complete.cases(x, y)
  n <- sum(complete)
  if (
    n < 4L ||
      stats::sd(x[complete]) == 0 ||
      stats::sd(y[complete]) == 0
  ) {
    return(list(
      estimate = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p = NA_real_
    ))
  }
  test <- tryCatch(
    stats::cor.test(x[complete], y[complete], method = "pearson"),
    error = function(e) NULL
  )
  if (is.null(test)) {
    return(list(
      estimate = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p = NA_real_
    ))
  }
  confidence_interval <- if (!is.null(test$conf.int)) {
    as.numeric(test$conf.int)
  } else {
    c(NA_real_, NA_real_)
  }
  list(
    estimate = unname(test$estimate),
    ci_low = confidence_interval[1],
    ci_high = confidence_interval[2],
    p = test$p.value
  )
}

.make_trait_state <- function(data, score_variables, min_valid = 3) {
  out <- data
  for (variable in score_variables) {
    by_id <- out |>
      dplyr::group_by(.data$studyID) |>
      dplyr::summarise(
        n_valid = sum(!is.na(.data[[variable]])),
        trait = ifelse(
          .data$n_valid >= min_valid,
          mean(.data[[variable]], na.rm = TRUE),
          NA_real_
        ),
        .groups = "drop"
      )
    count_name <- paste0(variable, "_n_valid")
    trait_name <- paste0(variable, "_trait")
    state_name <- paste0(variable, "_state")
    names(by_id)[names(by_id) == "n_valid"] <- count_name
    names(by_id)[names(by_id) == "trait"] <- trait_name
    out <- dplyr::left_join(out, by_id, by = "studyID")
    out[[state_name]] <- ifelse(
      !is.na(out[[trait_name]]),
      out[[variable]] - out[[trait_name]],
      NA_real_
    )
  }
  out
}

EXTERNAL_PHQ_VARS <- "PHQ9_Total"
EXTERNAL_VARS <- c(
  "MADRS_Total", "CGI_S", "HAMA_Total", "YMRS_Total",
  "GAD2_Total", "PSS_Total", "RULS_Total", "MHC_Total"
)

build_external_scores <- function(weekly, state, base) {
  require_columns(
    weekly,
    c(
      "studyID", "studyWeek", PHQ9_ITEMS, GAD2_ITEMS,
      paste0("pss", sprintf("%02d", 1:10)), paste0("ruls", 1:6),
      paste0("mhc", sprintf("%02d", 1:14))
    ),
    "external-validity weekly input"
  )
  require_columns(
    state,
    c(
      "studyID", "studyWeek", "cgis",
      paste0("madrs", sprintf("%02d", 1:10)),
      paste0("hama", sprintf("%02d", 1:14)),
      paste0("ymrs", sprintf("%02d", 1:10))
    ),
    "external-validity clinician input"
  )
  weekly <- weekly |>
    dplyr::mutate(
      studyID = as.character(.data$studyID),
      studyWeek = as.numeric(.data$studyWeek)
    ) |>
    dplyr::filter(.data$studyWeek <= 24)
  valid_ids <- unique(weekly$studyID)
  state_valid <- state |>
    dplyr::mutate(
      studyID = as.character(.data$studyID),
      studyWeek = as.numeric(.data$studyWeek)
    ) |>
    dplyr::filter(.data$studyID %in% valid_ids, .data$studyWeek <= 24)

  weekly_data <- as.data.frame(weekly)
  weekly_scores <- data.frame(
    studyID = weekly_data$studyID,
    studyWeek = weekly_data$studyWeek,
    PHQ9_Total = .scale_complete_sum(weekly_data, PHQ9_ITEMS),
    GAD2_Total = .scale_complete_sum(weekly_data, GAD2_ITEMS),
    PSS_Total = .scale_complete_sum(weekly_data, paste0("pss", sprintf("%02d", 1:10))),
    RULS_Total = .scale_complete_sum(weekly_data, paste0("ruls", 1:6)),
    MHC_Total = .scale_complete_sum(weekly_data, paste0("mhc", sprintf("%02d", 1:14))),
    stringsAsFactors = FALSE
  )

  state_data <- as.data.frame(state_valid)
  state_scores <- data.frame(
    studyID = state_data$studyID,
    studyWeek = state_data$studyWeek,
    MADRS_Total = .scale_complete_sum(state_data, paste0("madrs", sprintf("%02d", 1:10))),
    HAMA_Total = .scale_complete_sum(state_data, paste0("hama", sprintf("%02d", 1:14))),
    YMRS_Total = .scale_complete_sum(state_data, paste0("ymrs", sprintf("%02d", 1:10))),
    CGI_S = state_data$cgis,
    stringsAsFactors = FALSE
  )

  scores <- dplyr::left_join(
    weekly_scores,
    state_scores,
    by = c("studyID", "studyWeek")
  )
  .make_trait_state(scores, c(EXTERNAL_PHQ_VARS, EXTERNAL_VARS), min_valid = 3)
}

compute_correlation_grid <- function(scores) {
  configure_determinism()
  score_columns <- unlist(lapply(
    c(EXTERNAL_PHQ_VARS, EXTERNAL_VARS),
    function(variable) paste0(variable, c("_trait", "_state"))
  ))
  require_columns(scores, c("studyID", score_columns), "correlation-grid input")

  compute_correlation <- function(phq, external, level) {
    if (level == "between") {
      by_id <- data.frame(
        studyID = scores$studyID,
        x = scores[[paste0(phq, "_trait")]],
        y = scores[[paste0(external, "_trait")]]
      ) |>
        dplyr::distinct(.data$studyID, .keep_all = TRUE)
      x <- by_id$x
      y <- by_id$y
      ids <- by_id$studyID
    } else {
      x <- scores[[paste0(phq, "_state")]]
      y <- scores[[paste0(external, "_state")]]
      ids <- scores$studyID
    }

    complete <- stats::complete.cases(x, y)
    test <- .safe_cor_test(x, y)
    ci_method <- "Pearson cor.test"
    set.seed(20250823L)
    cluster_ci <- .cluster_boot_cor_ci(x, y, ids)
    if (all(is.finite(cluster_ci))) {
      test$ci_low <- cluster_ci[1]
      test$ci_high <- cluster_ci[2]
      ci_method <- CLUSTER_BOOTSTRAP_LABEL
    }
    data.frame(
      phq = phq,
      external = external,
      level = level,
      n_observations = sum(complete),
      n_subjects = length(unique(ids[complete])),
      estimate = test$estimate,
      ci_low = test$ci_low,
      ci_high = test$ci_high,
      ci_method = ci_method,
      p = test$p,
      stringsAsFactors = FALSE
    )
  }

  grid <- tidyr::expand_grid(
    phq = EXTERNAL_PHQ_VARS,
    external = EXTERNAL_VARS,
    level = c("between", "within")
  )
  do.call(rbind, Map(
    compute_correlation,
    grid$phq,
    grid$external,
    grid$level
  ))
}

weekly_validity_stability <- function(scores, weekly, state, base, mde_auc_by_week) {
  configure_determinism()
  require_columns(
    scores,
    c("studyID", "studyWeek", "PHQ9_Total", "MADRS_Total"),
    "weekly validity scores"
  )
  require_columns(weekly, c("studyID", "studyWeek", PHQ9_ITEMS), "weekly validity input")
  mde_state <- attach_concurrent_mde(state, base) |>
    dplyr::select(dplyr::all_of(c("studyID", "studyWeek", "MDE")))
  weekly <- weekly |>
    dplyr::mutate(
      studyID = as.character(.data$studyID),
      studyWeek = as.numeric(.data$studyWeek)
    )
  weekly_total <- weekly |>
    dplyr::transmute(
      studyID = .data$studyID,
      studyWeek = .data$studyWeek,
      PHQ9_Total = .scale_complete_sum(as.data.frame(weekly), PHQ9_ITEMS)
    )
  mde_week <- weekly_total |>
    dplyr::inner_join(mde_state, by = c("studyID", "studyWeek")) |>
    dplyr::filter(
      .data$studyWeek <= 24,
      stats::complete.cases(.data$PHQ9_Total, .data$MDE)
    )

  validity_weeks <- sort(unique(scores$studyWeek[
    !is.na(scores$PHQ9_Total) & !is.na(scores$MADRS_Total)
  ]))

  madrs_rows <- lapply(validity_weeks, function(week) {
    data <- scores |>
      dplyr::filter(.data$studyWeek == week) |>
      dplyr::select(dplyr::all_of(c(
        "studyID", "studyWeek", "PHQ9_Total", "MADRS_Total"
      ))) |>
      dplyr::filter(stats::complete.cases(.data$PHQ9_Total, .data$MADRS_Total))
    test <- .safe_cor_test(data$PHQ9_Total, data$MADRS_Total)
    cluster_ci <- .cluster_boot_cor_ci(data$PHQ9_Total, data$MADRS_Total, data$studyID)
    ci_method <- "Pearson cor.test"
    if (all(is.finite(cluster_ci))) {
      test$ci_low <- cluster_ci[1]
      test$ci_high <- cluster_ci[2]
      ci_method <- CLUSTER_BOOTSTRAP_LABEL
    }
    data.frame(
      panel = "MADRS",
      week = week,
      n = nrow(data),
      n_positive = NA_integer_,
      estimate = test$estimate,
      ci_low = test$ci_low,
      ci_high = test$ci_high,
      p = test$p,
      metric = "Pearson r",
      note = "Same-week PHQ-9 total vs MADRS total correlation",
      ci_method = ci_method,
      stringsAsFactors = FALSE
    )
  })

  required_auc_columns <- c("Week", "AUC", "CI_Lower", "CI_Upper", "N")
  if (!all(required_auc_columns %in% names(mde_auc_by_week))) {
    stop("Canonical MDE AUC table is missing required columns.")
  }
  mde_rows <- lapply(validity_weeks, function(week) {
    data <- mde_week[mde_week$studyWeek == week, ]
    canonical <- mde_auc_by_week[mde_auc_by_week$Week == week, , drop = FALSE]
    if (nrow(canonical) != 1L || nrow(data) != canonical$N[[1]]) {
      stop("Canonical MDE AUC row or N mismatch at Week ", week, ".")
    }
    data.frame(
      panel = "MDE",
      week = week,
      n = nrow(data),
      n_positive = sum(data$MDE == 1, na.rm = TRUE),
      estimate = canonical$AUC[[1]],
      ci_low = canonical$CI_Lower[[1]],
      ci_high = canonical$CI_Upper[[1]],
      p = NA_real_,
      metric = "AUC",
      note = "PHQ-9 total discrimination of major depressive episode status",
      ci_method = CLUSTER_BOOTSTRAP_LABEL,
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(madrs_rows, mde_rows)
}

compute_cluster_optimal_cutoff <- function(
    validity_frame,
    R = CLUSTER_BOOTSTRAP_REPLICATES,
    grid = 5:15) {
  configure_determinism()
  require_columns(
    validity_frame,
    c("studyID", "phq9_sum", "MDE"),
    "pooled cutoff input"
  )
  data <- validity_frame[
    !is.na(validity_frame$MDE) & !is.na(validity_frame$phq9_sum),
  ]

  roc <- pROC::roc(
    data$MDE,
    data$phq9_sum,
    quiet = TRUE,
    direction = "<",
    levels = c(0, 1)
  )
  auc <- as.numeric(roc$auc)
  optimal_cutoff <- .youden_cut(data, grid)
  optimal_performance <- .cutoff_perf(data, optimal_cutoff)
  standard_performance <- .cutoff_perf(data, 10)

  ids <- unique(data$studyID)
  bootstrap_auc <- numeric(R)
  bootstrap_cutoff <- numeric(R)
  bootstrap_optimal_sensitivity <- numeric(R)
  bootstrap_optimal_specificity <- numeric(R)
  bootstrap_standard_sensitivity <- numeric(R)
  bootstrap_standard_specificity <- numeric(R)
  for (iteration in seq_len(R)) {
    resampled <- .resample_clusters(data)
    resampled_roc <- tryCatch(
      pROC::roc(
        resampled$MDE,
        resampled$phq9_sum,
        quiet = TRUE,
        direction = "<",
        levels = c(0, 1)
      ),
      error = function(e) NULL
    )
    bootstrap_auc[iteration] <- if (is.null(resampled_roc)) {
      NA_real_
    } else {
      as.numeric(resampled_roc$auc)
    }
    bootstrap_cutoff[iteration] <- if (length(unique(resampled$MDE)) < 2L) {
      NA_integer_
    } else {
      .youden_cut(resampled, grid)
    }
    optimal <- .cutoff_perf(resampled, optimal_cutoff)
    standard <- .cutoff_perf(resampled, 10)
    bootstrap_optimal_sensitivity[iteration] <- optimal$sens
    bootstrap_optimal_specificity[iteration] <- optimal$spec
    bootstrap_standard_sensitivity[iteration] <- standard$sens
    bootstrap_standard_specificity[iteration] <- standard$spec
  }

  auc_ci <- .qci(bootstrap_auc)
  cutoff_ci <- .qci(bootstrap_cutoff)
  optimal_sensitivity_ci <- .qci(bootstrap_optimal_sensitivity)
  optimal_specificity_ci <- .qci(bootstrap_optimal_specificity)
  standard_sensitivity_ci <- .qci(bootstrap_standard_sensitivity)
  standard_specificity_ci <- .qci(bootstrap_standard_specificity)

  data.frame(
    rule = c("Youden-optimal", "Standard (>=10)"),
    cutoff = c(optimal_cutoff, 10L),
    cutoff_ci_low = c(round(cutoff_ci[1]), NA_integer_),
    cutoff_ci_high = c(round(cutoff_ci[2]), NA_integer_),
    sensitivity = c(optimal_performance$sens, standard_performance$sens),
    sens_ci_low = c(optimal_sensitivity_ci[1], standard_sensitivity_ci[1]),
    sens_ci_high = c(optimal_sensitivity_ci[2], standard_sensitivity_ci[2]),
    specificity = c(optimal_performance$spec, standard_performance$spec),
    spec_ci_low = c(optimal_specificity_ci[1], standard_specificity_ci[1]),
    spec_ci_high = c(optimal_specificity_ci[2], standard_specificity_ci[2]),
    youden_j = c(optimal_performance$youden, standard_performance$youden),
    auc = auc,
    auc_ci_low = auc_ci[1],
    auc_ci_high = auc_ci[2],
    n_observations = nrow(data),
    n_subjects = length(ids),
    n_mde_positive = sum(data$MDE == 1),
    stringsAsFactors = FALSE
  )
}

cutoff_stability_by_week <- function(
    validity_frame,
    mde_auc_by_week,
    B = CLUSTER_BOOTSTRAP_REPLICATES,
    grid = 5:15) {
  configure_determinism()
  required_auc_columns <- c("Week", "AUC", "CI_Lower", "CI_Upper", "N")
  if (!all(required_auc_columns %in% names(mde_auc_by_week))) {
    stop("Canonical MDE AUC table is missing required columns.")
  }
  require_columns(
    validity_frame,
    c("studyID", "studyWeek", "phq9_sum", "MDE"),
    "week-specific cutoff input"
  )
  data <- validity_frame[
    !is.na(validity_frame$MDE) & !is.na(validity_frame$phq9_sum),
  ]
  weeks <- sort(unique(data$studyWeek))
  rows <- list()

  for (week in weeks) {
    subset <- data[data$studyWeek == week, ]
    n <- nrow(subset)
    n_positive <- sum(subset$MDE == 1)
    if (n < 10L || length(unique(subset$MDE)) < 2L) {
      next
    }
    canonical <- mde_auc_by_week[mde_auc_by_week$Week == week, , drop = FALSE]
    if (nrow(canonical) != 1L || canonical$N[[1]] != n) {
      stop("Canonical MDE AUC row or N mismatch at Week ", week, ".")
    }

    optimal_cutoff <- .youden_cut(subset, grid)
    optimal_performance <- .cutoff_perf(subset, optimal_cutoff)
    standard_performance <- .cutoff_perf(subset, 10)
    bootstrap_cutoff <- integer(B)
    for (iteration in seq_len(B)) {
      resampled <- .resample_clusters(subset)
      if (
        length(unique(resampled$MDE)) < 2L ||
          stats::sd(resampled$phq9_sum) == 0
      ) {
        bootstrap_cutoff[iteration] <- NA_integer_
        next
      }
      bootstrap_cutoff[iteration] <- .youden_cut(resampled, grid)
    }
    cutoff_ci <- .qci(bootstrap_cutoff)

    rows[[as.character(week)]] <- data.frame(
      week = week,
      n = n,
      mde_positive = n_positive,
      auc = canonical$AUC[[1]],
      auc_ci_low = canonical$CI_Lower[[1]],
      auc_ci_high = canonical$CI_Upper[[1]],
      youden_cut = optimal_cutoff,
      cut_ci_low = round(cutoff_ci[1]),
      cut_ci_high = round(cutoff_ci[2]),
      youden_sensitivity = optimal_performance$sens,
      youden_specificity = optimal_performance$spec,
      sensitivity_at_10 = standard_performance$sens,
      specificity_at_10 = standard_performance$spec,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}
