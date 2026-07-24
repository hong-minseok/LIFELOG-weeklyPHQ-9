# Between- and within-person PHQ-9 structure.

mcfa_models_phq9 <- function() {
  list(
    "2-Factor" = lav_mlcfa_syntax(PHQ_MODEL_SPECS[["PHQ-9"]][["2-Factor"]])
  )
}

run_mcfa <- function(df_long, models_ml = mcfa_models_phq9()) {
  configure_determinism()
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("The lavaan package is required for multilevel CFA.")
  }
  require_columns(df_long, c("studyID", PHQ9_ITEMS), "multilevel CFA input")
  df_long <- df_long |>
    dplyr::mutate(studyID = as.character(.data$studyID))

  icc_list <- list()
  isomorphism_list <- list()
  loading_list <- list()
  sample_size_list <- list()

  for (model_name in names(models_ml)) {
    fit <- lavaan::cfa(
      models_ml[[model_name]],
      data = df_long,
      cluster = "studyID",
      estimator = "MLR"
    )
    sample_size_list[[model_name]] <- data.frame(
      Model = model_name,
      N_Participants = lavaan::lavInspect(fit, "nclusters"),
      N_Observations = lavaan::lavInspect(fit, "ntotal")
    )

    parameters <- lavaan::parameterEstimates(fit, standardized = TRUE) |>
      dplyr::filter(.data$op == "=~")
    level_1 <- parameters |>
      dplyr::filter(.data$level == 1) |>
      dplyr::transmute(lhs = .data$lhs, rhs = .data$rhs, L1 = .data$std.all)
    level_2 <- parameters |>
      dplyr::filter(.data$level == 2) |>
      dplyr::transmute(lhs = .data$lhs, rhs = .data$rhs, L2 = .data$std.all)
    loading_comparison <- dplyr::inner_join(level_1, level_2, by = c("lhs", "rhs"))
    loading_comparison$Model <- model_name
    loading_list[[model_name]] <- loading_comparison
    isomorphism_list[[model_name]] <- data.frame(
      Model = model_name,
      Correlation = stats::cor(loading_comparison$L1, loading_comparison$L2)
    )

    implied <- lavaan::lavInspect(fit, "implied")
    covariance_within <- diag(implied[[1]]$cov)
    covariance_between <- diag(implied[[2]]$cov)
    items <- intersect(PHQ9_ITEMS, names(covariance_within))
    icc_list[[model_name]] <- data.frame(
      Model = model_name,
      Item = items,
      Var_Within = covariance_within[items],
      Var_Between = covariance_between[items]
    ) |>
      dplyr::mutate(
        Total = .data$Var_Within + .data$Var_Between,
        ICC = .data$Var_Between / .data$Total,
        State_Pct = 1 - .data$ICC
      )
  }

  list(
    icc = dplyr::bind_rows(icc_list),
    iso_stats = dplyr::bind_rows(isomorphism_list),
    iso_data = dplyr::bind_rows(loading_list),
    sample_sizes = dplyr::bind_rows(sample_size_list)
  )
}

build_trait_state_summary <- function(mcfa) {
  isomorphism <- mcfa$iso_stats[mcfa$iso_stats$Model == "2-Factor", ]
  isomorphism$instrument <- "PHQ-9"
  isomorphism$status <- "primary"

  item_icc <- mcfa$icc[mcfa$icc$Model == "2-Factor", ]
  item_icc$instrument <- "PHQ-9"
  item_icc$status <- "primary"
  mean_icc <- stats::aggregate(
    cbind(ICC, State_Pct) ~ instrument + Model + status,
    item_icc,
    function(x) mean(x, na.rm = TRUE)
  )

  list(
    iso_primary = isomorphism,
    icc_primary = mean_icc
  )
}
