#!/usr/bin/env Rscript

arguments_full <- commandArgs(trailingOnly = FALSE)
script_argument <- grep("^--file=", arguments_full, value = TRUE)
script_path <- if (length(script_argument)) sub("^--file=", "", script_argument[[1]]) else "run_empirical.R"
project_dir <- normalizePath(dirname(script_path), mustWork = TRUE)
Sys.setenv(BRM_PROJECT_DIR = project_dir)

dir.create(file.path(project_dir, "outputs"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(project_dir, "validation"), showWarnings = FALSE, recursive = TRUE)

source(file.path(project_dir, "src", "R", "analyse_empirical_case.R"), echo = FALSE)
source(file.path(project_dir, "src", "R", "audit_empirical_sensitivity.R"), echo = FALSE)

message("Descriptive empirical reconstruction completed in ", project_dir)
