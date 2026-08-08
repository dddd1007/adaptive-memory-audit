project_dir <- Sys.getenv("BRM_PROJECT_DIR", unset = normalizePath(getwd(), mustWork = TRUE))
source_dir <- file.path(project_dir, "src", "R")
sources <- list.files(source_dir, pattern = "\\.R$", full.names = TRUE)
results <- lapply(sources, function(path) {
  error <- NULL
  tryCatch(parse(file = path, keep.source = TRUE), error = function(condition) error <<- conditionMessage(condition))
  data.frame(
    file = basename(path),
    parsed = is.null(error),
    error = if (is.null(error)) "" else error,
    stringsAsFactors = FALSE
  )
})
results <- do.call(rbind, results)
validation_dir <- file.path(project_dir, "validation")
if (!dir.exists(validation_dir)) dir.create(validation_dir, recursive = TRUE)
write.csv(results, file.path(validation_dir, "r_source_parse.csv"), row.names = FALSE)
if (!all(results$parsed)) stop("R source parsing failed: ", paste(results$file[!results$parsed], collapse = ", "))
message("R parsed all ", nrow(results), " source files.")
