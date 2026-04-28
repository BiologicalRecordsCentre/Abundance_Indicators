# 01_inventory.R
# Purpose: Map all existing scripts and functions on DataLabs for the England Species Abundance Indicator codebase.
# Output: "code_inventory.csv" and "package_usage.csv" written to the working directory.
# Run: In RStudio/terminal, setwd("<your-project-root>"); source("01_inventory.R")

suppressPackageStartupMessages({
  pkgs <- c("dplyr","purrr","stringr","tibble","readr","tools")
  to_install <- setdiff(pkgs, rownames(installed.packages()))
  if(length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
  lapply(pkgs, library, character.only = TRUE)
})

# ---------- Configuration ----------
# Root directory to scan (change if needed)
root_dir <- getwd()

# File patterns to include
include_ext <- c(".R", ".r", ".Rmd", ".qmd", ".Rmarkdown")
# -----------------------------------

# Helper: list files recursively with allowed extensions
list_code_files <- function(root = ".", exts = include_ext){
  all <- list.files(root, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
  all[tolower(tools::file_ext(all)) %in% tolower(sub("^\\.", "", exts))]
}

# Helper: read file safely
read_file_safe <- function(path){
  tryCatch(paste(readLines(path, warn = FALSE), collapse = "\n"),
           error = function(e) "")
}

# Extract function names via regex (lightweight; robust enough for inventory)
# Matches "name <- function(" at line starts or after whitespace
find_functions <- function(text){
  # Remove roxygen comment lines for cleaner matching
  txt <- gsub("^#'\\s?.*$", "", text, perl = TRUE)
  m <- str_match_all(txt, "(?m)^\\s*([a-zA-Z0-9._]+)\\s*<-\\s*function\\s*\\(")[[1]]
  unique(m[,2])
}

# Extract roxygen @title or first roxygen line above a function
extract_roxygen_title <- function(text, fname){
  # Find the function definition line
  pattern <- paste0("(?m)^\\s*", fname, "\\s*<-\\s*function\\s*\\(")
  loc <- str_locate(text, pattern)
  if(any(is.na(loc))) return(NA_character_)
  # Get text above the function and scan for roxygen
  up_to <- substr(text, 1, loc[1]-1)
  # Grab last roxygen block
  roxy <- str_match(up_to, "(?s)(#'.*?)(?:\\n[^#']|\\Z)")
  if(is.na(roxy[1])) return(NA_character_)
  # First non-empty roxygen line without #' prefix
  lines <- str_split(roxy[2], "\n")[[1]]
  lines <- gsub("^#'\\s?", "", lines)
  lines <- lines[nzchar(trimws(lines))]
  if(length(lines)==0) return(NA_character_)
  lines[[1]]
}

# Detect package usage: library(), require(), and pkg::fun tokens
detect_packages <- function(text){
  libs <- str_match_all(text, "(?m)\\b(library|require)\\s*\\(\\s*([A-Za-z0-9._]+)\\s*\\)")[[1]]
  p1 <- if(nrow(libs)) unique(libs[,3]) else character(0)
  colons <- str_match_all(text, "\\b([A-Za-z0-9._]+)\\s*::\\s*[A-Za-z0-9._]+")[[1]]
  p2 <- if(nrow(colons)) unique(colons[,2]) else character(0)
  unique(sort(c(p1,p2)))
}

# Detect IO patterns: reads and writes (basic)
detect_io <- function(text){
  reads <- c("read_csv","read_csv2","read_tsv","read_delim","fread","readRDS","readr::read_",
             "vroom","readxl::read_","sf::st_read","arrow::read_","qs::qread","jsonlite::fromJSON")
  writes <- c("write_csv","write_csv2","write_tsv","write_delim","fwrite","saveRDS","readr::write_",
              "writexl::write_","sf::st_write","arrow::write_","qs::qsave","jsonlite::toJSON",
              "openxlsx::write")
  r <- unique(unlist(str_extract_all(text, paste0("\\b(", paste(reads, collapse="|"), ")\\w*"))))
  w <- unique(unlist(str_extract_all(text, paste0("\\b(", paste(writes, collapse="|"), ")\\w*"))))
  list(reads = sort(r), writes = sort(w))
}

# Detect sourcing of other scripts
detect_sources <- function(text){
  s <- str_match_all(text, "(?m)\\bsource\\s*\\(\\s*['\"]([^'\"]+)['\"]")[[1]]
  if(nrow(s)) unique(s[,2]) else character(0)
}

# Main scan
files <- list_code_files(root_dir, include_ext)
if(length(files)==0){
  message("No code files found under: ", root_dir)
}

records <- purrr::map_dfr(files, function(f){
  txt <- read_file_safe(f)
  fns <- find_functions(txt)
  pkgs <- detect_packages(txt)
  io   <- detect_io(txt)
  srcs <- detect_sources(txt)

  if(length(fns)==0){
    tibble::tibble(
      file = f,
      `function` = NA_character_,
      title = NA_character_,
      packages = paste(pkgs, collapse=";"),
      reads = paste(io$reads, collapse=";"),
      writes = paste(io$writes, collapse=";"),
      sources = paste(srcs, collapse=";")
    )
  } else {
    tibble::tibble(
      file = f,
      `function` = fns,
      title = purrr::map_chr(fns, ~extract_roxygen_title(txt, .x)),
      packages = paste(pkgs, collapse=";"),
      reads = paste(io$reads, collapse=";"),
      writes = paste(io$writes, collapse=";"),
      sources = paste(srcs, collapse=";")
    )
  }
})

# Summarise package usage across repo
pkg_usage <- records |>
  tidyr::separate_rows(packages, sep=";") |>
  dplyr::filter(nzchar(packages)) |>
  dplyr::count(packages, sort = TRUE)

# Write outputs
readr::write_csv(records, "code_inventory.csv")
readr::write_csv(pkg_usage, "package_usage.csv")

message("Wrote: code_inventory.csv (function-level map)")
message("Wrote: package_usage.csv (packages referenced)")

# (Optional) Print a quick preview
print(head(records, 20))
print(head(pkg_usage, 20))
