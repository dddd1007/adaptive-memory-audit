brm_validate_public_data <- function(project_dir, write_report = TRUE) {
  data_path <- file.path(
    project_dir, "data", "analysis_transitions_deidentified.csv"
  )
  validation_dir <- file.path(project_dir, "validation")
  if (!dir.exists(validation_dir)) {
    dir.create(validation_dir, recursive = TRUE)
  }

  data <- brm_read_csv(data_path)
  required_fields <- c(
    "participant_code", "card_code", "study_day", "daily_review_index",
    "first_ease", "last_ease", "last_scheduled_interval_days",
    "next_gap_days", "next_first_ease", "scheduler_encoding"
  )
  if (!identical(names(data), required_fields)) {
    stop(
      "Unexpected public-data schema. Expected: ",
      paste(required_fields, collapse = ", ")
    )
  }
  data$next_failure <- as.integer(data$next_first_ease == 1)

  checks <- list()
  add_check <- function(name, passed, observed, expected) {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = name,
      passed = isTRUE(passed),
      observed = as.character(observed),
      expected = as.character(expected),
      stringsAsFactors = FALSE
    )
  }

  composite_key <- paste(
    data$participant_code, data$card_code, data$study_day, sep = "\r"
  )
  forbidden_column_pattern <- paste(
    c(
      "name", "email", "phone", "address", "device", "filename", "path",
      "timestamp", "clock", "midnight", "absolute_date", "event_id",
      "original_id", "linkage"
    ),
    collapse = "|"
  )

  add_check(
    "public transition rows and fields",
    nrow(data) == 2045L && ncol(data) == 11L,
    sprintf("%d rows; %d released fields", nrow(data), length(required_fields)),
    "2,045 rows; 10 released fields"
  )
  add_check(
    "pseudonymous participant and card code formats",
    all(grepl("^P[0-9]{2}$", data$participant_code)) &&
      all(grepl("^C[0-9]{3}$", data$card_code)) &&
      length(unique(data$participant_code)) == 12L &&
      length(unique(data$card_code)) == 100L,
    sprintf(
      "%d participant codes; %d card codes",
      length(unique(data$participant_code)),
      length(unique(data$card_code))
    ),
    "12 participant codes; 100 card codes"
  )
  add_check(
    "no direct-identifier or absolute-time fields",
    !any(grepl(forbidden_column_pattern, required_fields, ignore.case = TRUE)),
    paste(required_fields, collapse = ", "),
    "no forbidden field names"
  )
  add_check(
    "participant-card-day key is unique",
    !anyDuplicated(composite_key),
    sum(duplicated(composite_key)),
    0L
  )
  add_check(
    "required fields are complete",
    !anyNA(data[, required_fields, drop = FALSE]),
    sum(is.na(data[, required_fields, drop = FALSE])),
    0L
  )
  add_check(
    "ratings and relative study days are valid",
    all(data$first_ease %in% 1:4) &&
      all(data$last_ease %in% 1:4) &&
      all(data$next_first_ease %in% 1:4) &&
      all(data$study_day %in% 1:13),
    sprintf(
      "rating range %d--%d; study-day range %d--%d",
      min(c(data$first_ease, data$last_ease, data$next_first_ease)),
      max(c(data$first_ease, data$last_ease, data$next_first_ease)),
      min(data$study_day), max(data$study_day)
    ),
    "ratings 1--4; relative study days 1--13"
  )
  add_check(
    "delay fields are finite and positive",
    all(is.finite(data$last_scheduled_interval_days)) &&
      all(is.finite(data$next_gap_days)) &&
      all(data$last_scheduled_interval_days > 0) &&
      all(data$next_gap_days > 0),
    sprintf(
      "planned %.6f--%.3f days; actual %.6f--%.3f days",
      min(data$last_scheduled_interval_days),
      max(data$last_scheduled_interval_days),
      min(data$next_gap_days), max(data$next_gap_days)
    ),
    "finite positive values"
  )
  add_check(
    "Again outcome count",
    sum(data$next_failure) == 199L,
    sum(data$next_failure),
    199L
  )

  validation <- do.call(rbind, checks)
  if (write_report) {
    write.csv(
      validation,
      file.path(validation_dir, "public_data_validation.csv"),
      row.names = FALSE
    )
  }
  if (!all(validation$passed)) {
    stop(
      "Public data validation failed: ",
      paste(validation$check[!validation$passed], collapse = "; ")
    )
  }

  list(
    data = data,
    checks = validation,
    required_fields = required_fields,
    path = data_path
  )
}
