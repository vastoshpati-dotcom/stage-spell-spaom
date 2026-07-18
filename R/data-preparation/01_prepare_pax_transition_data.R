# =============================================================================
# PA-X TRANSITION DATASET GENERATOR
# =============================================================================
#
# PURPOSE
# -------
# This program creates:
#
#   01_primary_transition_dataset.csv
#
# directly from the original PA-X agreement-level workbook, normally named:
#
#   pax_data_2257_agreements_13-04-26(1).xlsx
#
# The resulting CSV is the input required by:
#
#   PA_X_SPAOM_Placeholder_Value_Generator.R
#
# PRIMARY SAMPLE RULE
# -------------------
# Agreements are ordered within peace process by:
#   1. agreement date;
#   2. agreement identifier;
#   3. version number.
#
# Consecutive agreements are converted into transitions. Self-transitions are
# retained. The primary sample retains processes contributing at least two
# observed transitions, matching the revised PA-X application.
#
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Packages
# -----------------------------------------------------------------------------

REQUIRED_PACKAGES <- c(
  "readxl",
  "dplyr"
)

missing_packages <- REQUIRED_PACKAGES[
  !vapply(
    REQUIRED_PACKAGES,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  install.packages(missing_packages)
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
})


# -----------------------------------------------------------------------------
# 2. User settings
# -----------------------------------------------------------------------------

# Leave as NA to search Downloads automatically. When the workbook is not found,
# a file-selection window will open. Select the ORIGINAL PA-X .xlsx workbook.
WORKBOOK_FILE <- NA_character_

MINIMUM_TRANSITIONS_PER_PROCESS <- 2L


# -----------------------------------------------------------------------------
# 3. Locate Downloads and the original PA-X workbook
# -----------------------------------------------------------------------------

find_downloads_directory <- function() {
  candidates <- unique(
    c(
      file.path(
        Sys.getenv("USERPROFILE"),
        "Downloads"
      ),
      file.path(
        Sys.getenv("HOME"),
        "Downloads"
      ),
      file.path(
        path.expand("~"),
        "Downloads"
      ),
      getwd()
    )
  )

  candidates <- candidates[
    nzchar(candidates)
  ]

  for (candidate in candidates) {
    if (dir.exists(candidate)) {
      return(
        normalizePath(
          candidate,
          winslash = "/",
          mustWork = FALSE
        )
      )
    }
  }

  normalizePath(
    getwd(),
    winslash = "/",
    mustWork = FALSE
  )
}

DOWNLOADS_DIR <- find_downloads_directory()

locate_workbook <- function() {
  if (
    is.character(WORKBOOK_FILE) &&
    length(WORKBOOK_FILE) == 1L &&
    !is.na(WORKBOOK_FILE) &&
    file.exists(WORKBOOK_FILE)
  ) {
    return(
      normalizePath(
        WORKBOOK_FILE,
        winslash = "/",
        mustWork = TRUE
      )
    )
  }

  exact_candidates <- c(
    file.path(
      DOWNLOADS_DIR,
      "pax_data_2257_agreements_13-04-26(1).xlsx"
    ),
    file.path(
      DOWNLOADS_DIR,
      "pax_data_2257_agreements_13-04-26.xlsx"
    ),
    file.path(
      getwd(),
      "pax_data_2257_agreements_13-04-26(1).xlsx"
    ),
    file.path(
      getwd(),
      "pax_data_2257_agreements_13-04-26.xlsx"
    )
  )

  existing <- exact_candidates[
    file.exists(exact_candidates)
  ]

  if (length(existing) > 0L) {
    return(
      normalizePath(
        existing[1L],
        winslash = "/",
        mustWork = TRUE
      )
    )
  }

  search_directories <- unique(
    c(
      DOWNLOADS_DIR,
      getwd()
    )
  )

  for (directory in search_directories) {
    matches <- list.files(
      path = directory,
      pattern = "^pax_data_.*agreements.*\\.xlsx$",
      full.names = TRUE,
      ignore.case = TRUE
    )

    if (length(matches) > 0L) {
      return(
        normalizePath(
          matches[1L],
          winslash = "/",
          mustWork = TRUE
        )
      )
    }
  }

  if (interactive()) {
    message(
      "Select the ORIGINAL PA-X agreement-level Excel workbook."
    )

    selected <- file.choose()

    if (
      !grepl(
        "\\.xlsx?$",
        selected,
        ignore.case = TRUE
      )
    ) {
      stop(
        "The selected file is not an Excel workbook."
      )
    }

    return(
      normalizePath(
        selected,
        winslash = "/",
        mustWork = TRUE
      )
    )
  }

  stop(
    "The original PA-X Excel workbook could not be found."
  )
}

DATA_FILE <- locate_workbook()

OUTPUT_FILE <- file.path(
  DOWNLOADS_DIR,
  "01_primary_transition_dataset.csv"
)

SUMMARY_FILE <- file.path(
  DOWNLOADS_DIR,
  "01_primary_transition_dataset_summary.txt"
)

cat("Original PA-X workbook:\n")
cat(DATA_FILE, "\n\n")

cat("Transition dataset will be saved as:\n")
cat(OUTPUT_FILE, "\n\n")


# -----------------------------------------------------------------------------
# 4. General helpers
# -----------------------------------------------------------------------------

safe_numeric <- function(x) {
  suppressWarnings(
    as.numeric(x)
  )
}

convert_pax_date <- function(x) {
  if (inherits(x, "Date")) {
    return(x)
  }

  if (inherits(x, "POSIXt")) {
    return(
      as.Date(x)
    )
  }

  if (is.numeric(x)) {
    return(
      as.Date(
        x,
        origin = "1899-12-30"
      )
    )
  }

  character_x <- trimws(
    as.character(x)
  )

  parsed <- suppressWarnings(
    as.Date(character_x)
  )

  unresolved <- is.na(parsed) &
    nzchar(character_x)

  if (any(unresolved)) {
    formats <- c(
      "%d-%m-%Y",
      "%d/%m/%Y",
      "%m/%d/%Y",
      "%d-%b-%Y",
      "%d %B %Y"
    )

    for (format_value in formats) {
      attempt <- suppressWarnings(
        as.Date(
          character_x[unresolved],
          format = format_value
        )
      )

      replaceable <- unresolved
      replaceable[unresolved] <- !is.na(attempt)

      parsed[replaceable] <- attempt[
        !is.na(attempt)
      ]

      unresolved <- is.na(parsed) &
        nzchar(character_x)

      if (!any(unresolved)) {
        break
      }
    }
  }

  parsed
}

standardize_stage <- function(x) {
  original <- trimws(
    as.character(x)
  )

  normalized <- tolower(
    gsub(
      "[^a-z0-9]+",
      "",
      original
    )
  )

  output <- rep(
    NA_character_,
    length(normalized)
  )

  output[
    normalized %in% c(
      "pre",
      "prenegotiation",
      "prenegotiationprocess"
    )
  ] <- "Pre-negotiation"

  output[
    normalized %in% c(
      "cea",
      "ceasefire"
    )
  ] <- "Ceasefire"

  output[
    normalized %in% c(
      "subpar",
      "substantivepartial",
      "partial",
      "partialsubstantive"
    )
  ] <- "Substantive Partial"

  output[
    normalized %in% c(
      "subcomp",
      "substantivecomprehensive",
      "comprehensive",
      "comprehensivesubstantive"
    )
  ] <- "Substantive Comprehensive"

  output[
    normalized %in% c(
      "imp",
      "implementation"
    )
  ] <- "Implementation"

  output[
    normalized %in% c(
      "ren",
      "renegotiation",
      "oth",
      "other",
      "renegotiationother"
    )
  ] <- "Renegotiation/Other"

  output
}


# -----------------------------------------------------------------------------
# 5. Read and validate the workbook
# -----------------------------------------------------------------------------

sheet_names <- readxl::excel_sheets(
  DATA_FILE
)

data_sheet <- if (
  "in" %in% sheet_names
) {
  "in"
} else {
  sheet_names[1L]
}

raw_data <- readxl::read_excel(
  DATA_FILE,
  sheet = data_sheet
)

essential_columns <- c(
  "PP",
  "Dat",
  "Stage"
)

missing_essential <- setdiff(
  essential_columns,
  names(raw_data)
)

if (length(missing_essential) > 0L) {
  stop(
    "The workbook is missing essential columns: ",
    paste(
      missing_essential,
      collapse = ", "
    )
  )
}

# Optional ordering and descriptive columns are created when absent.
if (!("AgtId" %in% names(raw_data))) {
  raw_data$AgtId <- seq_len(
    nrow(raw_data)
  )
}

if (!("Ver" %in% names(raw_data))) {
  raw_data$Ver <- 1
}

if (!("PPName" %in% names(raw_data))) {
  raw_data$PPName <- NA_character_
}

if (!("Agt" %in% names(raw_data))) {
  raw_data$Agt <- NA_character_
}

raw_data$.original_row <- seq_len(
  nrow(raw_data)
)


# -----------------------------------------------------------------------------
# 6. Harmonize agreements and construct consecutive transitions
# -----------------------------------------------------------------------------

agreement_data <- raw_data %>%
  mutate(
    PP = trimws(
      as.character(PP)
    ),
    PPName = as.character(PPName),
    Dat = convert_pax_date(Dat),
    AgtId_numeric = safe_numeric(AgtId),
    Ver_numeric = safe_numeric(Ver),
    Stage6 = standardize_stage(Stage)
  )

unknown_stage_values <- agreement_data %>%
  filter(
    !is.na(Stage),
    is.na(Stage6)
  ) %>%
  distinct(Stage) %>%
  pull(Stage)

if (length(unknown_stage_values) > 0L) {
  warning(
    "The following unrecognized stage values were excluded: ",
    paste(
      unknown_stage_values,
      collapse = ", "
    )
  )
}

agreement_data <- agreement_data %>%
  filter(
    !is.na(PP),
    nzchar(PP),
    !is.na(Dat),
    !is.na(Stage6)
  ) %>%
  mutate(
    AgtId_numeric = ifelse(
      is.finite(AgtId_numeric),
      AgtId_numeric,
      .original_row
    ),
    Ver_numeric = ifelse(
      is.finite(Ver_numeric),
      Ver_numeric,
      1
    )
  ) %>%
  arrange(
    PP,
    Dat,
    AgtId_numeric,
    Ver_numeric,
    .original_row
  ) %>%
  group_by(PP) %>%
  mutate(
    agreement_order = row_number(),
    process_agreements = n(),
    origin_stage = lag(Stage6),
    destination_stage = Stage6,
    origin_date = lag(Dat),
    destination_date = Dat,
    origin_agreement_id = lag(AgtId_numeric),
    destination_agreement_id = AgtId_numeric,
    origin_agreement_title = lag(
      as.character(Agt)
    ),
    destination_agreement_title = as.character(Agt),
    transition_index = agreement_order - 1L
  ) %>%
  ungroup()

all_transitions <- agreement_data %>%
  filter(
    transition_index >= 1L,
    !is.na(origin_stage),
    !is.na(destination_stage)
  ) %>%
  mutate(
    D = as.integer(
      origin_stage != destination_stage
    ),
    transition_id = paste(
      PP,
      transition_index,
      sep = "::"
    ),
    same_date_transition = as.integer(
      origin_date ==
        destination_date
    )
  )

process_transition_counts <- all_transitions %>%
  count(
    PP,
    name = "process_transitions"
  )

primary_processes <- process_transition_counts %>%
  filter(
    process_transitions >=
      MINIMUM_TRANSITIONS_PER_PROCESS
  ) %>%
  pull(PP)

primary_transition_data <- all_transitions %>%
  filter(
    PP %in% primary_processes
  ) %>%
  left_join(
    process_transition_counts,
    by = "PP"
  ) %>%
  transmute(
    PP = as.character(PP),
    PPName = as.character(PPName),
    origin_stage = as.character(origin_stage),
    destination_stage = as.character(
      destination_stage
    ),
    D = as.integer(D),
    transition_index = as.integer(
      transition_index
    ),
    transition_id = as.character(
      transition_id
    ),
    origin_date = as.character(
      origin_date
    ),
    destination_date = as.character(
      destination_date
    ),
    same_date_transition = as.integer(
      same_date_transition
    ),
    origin_agreement_id =
      origin_agreement_id,
    destination_agreement_id =
      destination_agreement_id,
    origin_agreement_title =
      origin_agreement_title,
    destination_agreement_title =
      destination_agreement_title,
    process_transitions =
      process_transitions
  ) %>%
  arrange(
    PP,
    transition_index
  )


# -----------------------------------------------------------------------------
# 7. Validate and save
# -----------------------------------------------------------------------------

required_output_columns <- c(
  "PP",
  "origin_stage",
  "destination_stage",
  "D",
  "transition_index",
  "transition_id"
)

missing_output_columns <- setdiff(
  required_output_columns,
  names(primary_transition_data)
)

if (length(missing_output_columns) > 0L) {
  stop(
    "Internal error: required output columns are missing."
  )
}

if (nrow(primary_transition_data) == 0L) {
  stop(
    "No primary-sample transitions were constructed."
  )
}

write.csv(
  primary_transition_data,
  OUTPUT_FILE,
  row.names = FALSE,
  na = ""
)

summary_lines <- c(
  "PA-X PRIMARY TRANSITION DATASET",
  "==============================",
  "",
  paste0(
    "Original workbook: ",
    DATA_FILE
  ),
  paste0(
    "Workbook sheet: ",
    data_sheet
  ),
  paste0(
    "Agreement records read: ",
    nrow(raw_data)
  ),
  paste0(
    "Usable agreement records: ",
    nrow(agreement_data)
  ),
  paste0(
    "Peace processes with at least one transition: ",
    dplyr::n_distinct(
      all_transitions$PP
    )
  ),
  paste0(
    "Primary minimum transitions per process: ",
    MINIMUM_TRANSITIONS_PER_PROCESS
  ),
  paste0(
    "Primary peace processes: ",
    dplyr::n_distinct(
      primary_transition_data$PP
    )
  ),
  paste0(
    "Primary transitions: ",
    nrow(primary_transition_data)
  ),
  paste0(
    "Self-transitions: ",
    sum(
      primary_transition_data$D == 0L
    )
  ),
  paste0(
    "Departures: ",
    sum(
      primary_transition_data$D == 1L
    )
  ),
  paste0(
    "Departure proportion: ",
    sprintf(
      "%.6f",
      mean(
        primary_transition_data$D
      )
    )
  ),
  paste0(
    "Same-date transitions: ",
    sum(
      primary_transition_data$
        same_date_transition,
      na.rm = TRUE
    )
  ),
  "",
  paste0(
    "CSV saved at: ",
    OUTPUT_FILE
  )
)

writeLines(
  summary_lines,
  SUMMARY_FILE
)

cat("\n============================================================\n")
cat("TRANSITION DATASET CREATED SUCCESSFULLY\n")
cat("============================================================\n")
cat("Primary peace processes: ")
cat(
  dplyr::n_distinct(
    primary_transition_data$PP
  ),
  "\n"
)
cat("Primary transitions: ")
cat(
  nrow(
    primary_transition_data
  ),
  "\n"
)
cat("Departure proportion: ")
cat(
  sprintf(
    "%.6f",
    mean(
      primary_transition_data$D
    )
  ),
  "\n\n"
)
cat("Created file:\n")
cat(OUTPUT_FILE, "\n\n")
cat(
  "Now run PA_X_SPAOM_Placeholder_Value_Generator.R. ",
  "It will locate this CSV automatically.\n"
)
