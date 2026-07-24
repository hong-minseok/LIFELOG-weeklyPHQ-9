# Participant-level cluster bootstrap and cutoff helpers.

.cluster_boot_ci <- function(data, statistic, R) {
  rows_by_id <- split(seq_len(nrow(data)), data$studyID)
  n_ids <- length(rows_by_id)
  bootstrap_values <- vapply(seq_len(R), function(iteration) {
    sampled_ids <- sample.int(n_ids, n_ids, replace = TRUE)
    sampled_rows <- unlist(rows_by_id[sampled_ids], use.names = FALSE)
    statistic(data[sampled_rows, , drop = FALSE])
  }, numeric(1))
  bootstrap_values <- bootstrap_values[is.finite(bootstrap_values)]
  if (length(bootstrap_values) < 2L) {
    return(c(NA_real_, NA_real_))
  }
  stats::quantile(
    bootstrap_values,
    c(0.025, 0.975),
    names = FALSE
  )
}

.cluster_boot_cor_ci <- function(
    x,
    y,
    ids,
    R = CLUSTER_BOOTSTRAP_REPLICATES) {
  complete <- stats::complete.cases(x, y, ids)
  data <- data.frame(
    studyID = as.character(ids[complete]),
    x = x[complete],
    y = y[complete],
    stringsAsFactors = FALSE
  )
  if (
    nrow(data) < 4L ||
      length(unique(data$studyID)) < 4L ||
      stats::sd(data$x) == 0 ||
      stats::sd(data$y) == 0
  ) {
    return(c(NA_real_, NA_real_))
  }
  .cluster_boot_ci(data, function(resampled) {
    if (
      nrow(resampled) < 4L ||
        stats::sd(resampled$x) == 0 ||
        stats::sd(resampled$y) == 0
    ) {
      return(NA_real_)
    }
    stats::cor(resampled$x, resampled$y, method = "pearson")
  }, R)
}

.cluster_boot_auc_ci <- function(
    outcome,
    score,
    ids,
    R = CLUSTER_BOOTSTRAP_REPLICATES) {
  complete <- stats::complete.cases(outcome, score, ids)
  data <- data.frame(
    studyID = as.character(ids[complete]),
    outcome = outcome[complete],
    score = score[complete],
    stringsAsFactors = FALSE
  )
  if (
    nrow(data) < 10L ||
      length(unique(data$studyID)) < 4L ||
      length(unique(data$outcome)) < 2L ||
      stats::sd(data$score) == 0
  ) {
    return(c(NA_real_, NA_real_))
  }
  .cluster_boot_ci(data, function(resampled) {
    if (
      nrow(resampled) < 10L ||
        length(unique(resampled$outcome)) < 2L ||
        stats::sd(resampled$score) == 0
    ) {
      return(NA_real_)
    }
    roc <- pROC::roc(
      resampled$outcome,
      resampled$score,
      quiet = TRUE,
      direction = "<"
    )
    as.numeric(roc$auc)
  }, R)
}

.cutoff_perf <- function(data, cutoff) {
  positive <- data$phq9_sum >= cutoff
  true_positive <- sum(positive & data$MDE == 1)
  false_negative <- sum(!positive & data$MDE == 1)
  false_positive <- sum(positive & data$MDE == 0)
  true_negative <- sum(!positive & data$MDE == 0)
  sensitivity <- true_positive / (true_positive + false_negative)
  specificity <- true_negative / (true_negative + false_positive)
  list(
    sens = sensitivity,
    spec = specificity,
    youden = sensitivity + specificity - 1
  )
}

.youden_cut <- function(data, grid) {
  youden <- vapply(grid, function(cutoff) {
    .cutoff_perf(data, cutoff)$youden
  }, numeric(1))
  grid[which.max(youden)]
}

.qci <- function(x) {
  stats::quantile(
    x,
    c(0.025, 0.975),
    na.rm = TRUE,
    names = FALSE
  )
}

.resample_clusters <- function(data) {
  ids <- unique(as.character(data$studyID))
  rows_by_id <- split(data, as.character(data$studyID))
  do.call(rbind, rows_by_id[sample(ids, length(ids), replace = TRUE)])
}
