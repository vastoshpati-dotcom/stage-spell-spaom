# =============================================================================
# STAGE-SPELL AGREEMENT OPPORTUNITY MODEL (SS-SPAOM)
# Complete empirical validation program for the PA-X transition data
# =============================================================================
#
# PRIMARY DEPARTURE MODEL
# -----------------------
# Let s_ir be the number of consecutive documented transitions for which
# peace process i has occupied its current origin stage before transition r.
#
#   N_iar ~ Poisson(mu_iar)
#
#   log(mu_iar)
#     = log(eta_a) - beta_s * log(s_ir),
#
#   q_ia(r)
#     = 1 - exp[-mu_iar],
#
# where eta_a > 0 is an origin-stage-specific baseline opportunity intensity
# and beta_s >= 0 is the common stage-spell persistence/fatigue parameter.
#
# The model reduces to the origin-specific homogeneous model when beta_s = 0.
#
# IMPORTANT CONSTRUCTION RULE
# ---------------------------
# The spell age is computed only from the sequence of origin stages already
# occupied before the current destination is observed. For
#
#   A -> A -> A -> B
#
# the spell ages are 1, 2, and 3. The clock resets to one after entry into B.
#
# PRIMARY COMPARATORS
# -------------------
# 1. Origin-specific homogeneous model
# 2. Earlier global-index kernel SPAOM
# 3. Global alpha-beta model
# 4. Stage-Spell SPAOM
# 5. Dual-clock model using global agreement position and stage-spell age
# 6. Flexible natural spline in log stage-spell age
# 7. Quadratic global-index model
#
# SENSITIVITY MODEL
# -----------------
# Origin-stage-specific spell coefficients beta_a.
#
# INFERENCE AND VALIDATION
# ------------------------
# - constrained maximum likelihood;
# - parametric-bootstrap LR test of H0: beta_s = 0;
# - process-bootstrap uncertainty for beta_s and eta_a;
# - repeated grouped process-level cross-validation;
# - process-bootstrap paired predictive differences;
# - departure and complete-transition metrics;
# - same-date, stage aggregation, process-length, truncation, weighting,
#   and leave-one-process-out robustness analyses;
# - publication-ready graphs, CSV tables, Excel workbook, and LaTeX block.
#
# =============================================================================


# -----------------------------------------------------------------------------
# 0. USER SETTINGS
# -----------------------------------------------------------------------------

BASE_DIR <- "C:/Users/Hp/Downloads"

# Leave as NA to locate a valid transition file recursively.
TRANSITION_FILE <- NA_character_

ANALYSIS_MODE <- "FINAL"  # Frozen publication run

OUTPUT_DIR <- file.path(
  BASE_DIR,
  "PA_X_SS_SPAOM_FINAL_Analysis"
)

CHECKPOINT_DIR <- file.path(
  OUTPUT_DIR,
  "Checkpoints"
)

GRAPH_DIR <- file.path(
  OUTPUT_DIR,
  "Graphs"
)

if (ANALYSIS_MODE == "QUICK") {
  NULL_BOOTSTRAP_B <- 49L
  PARAMETER_BOOTSTRAP_B <- 99L
  CV_REPEATS <- 3L
  CV_FOLDS <- 5L
  PAIRED_BOOTSTRAP_B <- 500L
  ROBUST_CV_REPEATS <- 1L
  ROBUST_CV_FOLDS <- 5L
  MULTISTARTS <- 12L
  OLD_TAU_GRID <- seq(0.50, 30, by = 1.00)
} else {
  NULL_BOOTSTRAP_B <- 499L
  PARAMETER_BOOTSTRAP_B <- 499L
  CV_REPEATS <- 10L
  CV_FOLDS <- 5L
  PAIRED_BOOTSTRAP_B <- 5000L
  ROBUST_CV_REPEATS <- 3L
  ROBUST_CV_FOLDS <- 5L
  MULTISTARTS <- 30L
  OLD_TAU_GRID <- seq(0.25, 30, by = 0.25)
}

BASE_SEED <- 20260714L
PROBABILITY_FLOOR <- 1e-12
DESTINATION_SMOOTHING <- 0.50
SPLINE_DF <- 4L

KAPPA_LOWER <- -12
KAPPA_UPPER <- 5
DYNAMIC_UPPER <- 8

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  CHECKPOINT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  GRAPH_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

options(
  stringsAsFactors = FALSE,
  warn = 1
)



# -----------------------------------------------------------------------------
# FROZEN PRIMARY SPECIFICATION
# -----------------------------------------------------------------------------
# This FINAL script does not search for a new model.  It estimates the
# pre-specified Stage-Spell SPAOM
#
#   log(mu_iar) = log(eta_a) - beta_s log(s_ir),  beta_s >= 0,
#
# using 499 null-bootstrap replications, 499 process-bootstrap parameter
# replications, and 10 repetitions of 5-fold grouped process-level
# cross-validation.  A new output/checkpoint directory is used so QUICK-run
# checkpoints cannot contaminate the publication run.

# -----------------------------------------------------------------------------
# 1. PACKAGES
# -----------------------------------------------------------------------------

REQUIRED_PACKAGES <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "writexl"
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
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})


# -----------------------------------------------------------------------------
# 2. GENERAL HELPERS
# -----------------------------------------------------------------------------

clip_probability <- function(probability) {
  pmin(
    pmax(probability, PROBABILITY_FLOOR),
    1 - PROBABILITY_FLOOR
  )
}

stable_log_q <- function(mu) {
  output <- numeric(length(mu))
  small <- mu <= log(2)

  output[small] <- log(
    -expm1(-mu[small])
  )

  output[!small] <- log1p(
    -exp(-mu[!small])
  )

  output
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    NA_real_
  } else {
    mean(x)
  }
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) <= 1L) {
    NA_real_
  } else {
    sd(x)
  }
}

format_number <- function(
  value,
  digits = 4L
) {
  if (
    length(value) == 0L ||
    !is.finite(value)
  ) {
    "NA"
  } else {
    formatC(
      value,
      format = "f",
      digits = digits
    )
  }
}

save_csv <- function(
  object,
  filename
) {
  write.csv(
    object,
    file.path(
      OUTPUT_DIR,
      filename
    ),
    row.names = FALSE,
    na = ""
  )
}

save_plot <- function(
  plot_object,
  filename,
  width = 10,
  height = 7
) {
  ggsave(
    filename = file.path(
      GRAPH_DIR,
      filename
    ),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300
  )
}

append_error <- function(
  stage,
  identifier,
  message
) {
  error_file <- file.path(
    OUTPUT_DIR,
    "analysis_errors.csv"
  )

  row <- data.frame(
    time = as.character(Sys.time()),
    stage = as.character(stage),
    identifier = as.character(identifier),
    message = as.character(message),
    stringsAsFactors = FALSE
  )

  write.table(
    row,
    file = error_file,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(error_file),
    append = file.exists(error_file),
    quote = TRUE
  )
}

as_date_safely <- function(x) {
  if (inherits(x, "Date")) {
    return(x)
  }

  if (inherits(x, "POSIXt")) {
    return(as.Date(x))
  }

  if (is.numeric(x)) {
    return(
      as.Date(
        x,
        origin = "1899-12-30"
      )
    )
  }

  suppressWarnings(
    as.Date(
      as.character(x)
    )
  )
}


# -----------------------------------------------------------------------------
# 3. LOCATE AND VALIDATE THE TRANSITION DATASET
# -----------------------------------------------------------------------------

locate_transition_file <- function() {

  required_columns_local <- c(
    "PP",
    "origin_stage",
    "destination_stage",
    "D",
    "transition_index",
    "transition_id"
  )

  validate_candidate <- function(candidate) {

    if (
      !is.character(candidate) ||
      length(candidate) != 1L ||
      is.na(candidate) ||
      !file.exists(candidate)
    ) {
      return(
        list(
          valid = FALSE,
          reason = "file does not exist",
          data = NULL
        )
      )
    }

    file_size <- file.info(candidate)$size

    if (
      !is.finite(file_size) ||
      file_size <= 0
    ) {
      return(
        list(
          valid = FALSE,
          reason = "file is empty",
          data = NULL
        )
      )
    }

    read_result <- tryCatch(
      list(
        data = read.csv(
          candidate,
          stringsAsFactors = FALSE,
          check.names = FALSE
        ),
        error = NULL
      ),
      error = function(e) {
        list(
          data = NULL,
          error = conditionMessage(e)
        )
      }
    )

    if (!is.null(read_result$error)) {
      return(
        list(
          valid = FALSE,
          reason = paste0(
            "read error: ",
            read_result$error
          ),
          data = NULL
        )
      )
    }

    test_data <- read_result$data

    if (nrow(test_data) == 0L) {
      return(
        list(
          valid = FALSE,
          reason = "CSV contains zero data rows",
          data = test_data
        )
      )
    }

    missing_columns_local <- setdiff(
      required_columns_local,
      names(test_data)
    )

    if (length(missing_columns_local) > 0L) {
      return(
        list(
          valid = FALSE,
          reason = paste0(
            "missing columns: ",
            paste(
              missing_columns_local,
              collapse = ", "
            )
          ),
          data = test_data
        )
      )
    }

    list(
      valid = TRUE,
      reason = "valid",
      data = test_data
    )
  }

  if (
    is.character(TRANSITION_FILE) &&
    length(TRANSITION_FILE) == 1L &&
    !is.na(TRANSITION_FILE) &&
    nzchar(TRANSITION_FILE)
  ) {
    direct_check <- validate_candidate(
      TRANSITION_FILE
    )

    if (isTRUE(direct_check$valid)) {
      selected_file <- normalizePath(
        TRANSITION_FILE,
        winslash = "/",
        mustWork = TRUE
      )

      cat(
        "Using explicitly supplied transition file:\n",
        selected_file,
        "\n\n",
        sep = ""
      )

      return(selected_file)
    }

    warning(
      "The supplied TRANSITION_FILE was rejected: ",
      direct_check$reason
    )
  }

  preferred_candidates <- c(
    file.path(
      BASE_DIR,
      "01_primary_transition_dataset.csv"
    ),
    file.path(
      BASE_DIR,
      "PA_X_SPAOM_Complete_Final_Analysis_V2",
      "01_primary_transition_dataset.csv"
    ),
    file.path(
      BASE_DIR,
      "PA_X_SPAOM_Complete_Final_Analysis",
      "01_primary_transition_dataset.csv"
    )
  )

  recursive_candidates <- list.files(
    BASE_DIR,
    pattern = "^01_primary_transition_dataset.*\\.csv$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )

  candidates <- unique(
    c(
      preferred_candidates,
      recursive_candidates
    )
  )

  candidates <- candidates[
    file.exists(candidates)
  ]

  if (length(candidates) == 0L) {
    stop(
      paste0(
        "No file beginning with '01_primary_transition_dataset' ",
        "was found under ",
        BASE_DIR,
        "."
      )
    )
  }

  candidate_information <- file.info(
    candidates
  )

  candidates <- candidates[
    order(
      candidate_information$mtime,
      decreasing = TRUE,
      na.last = TRUE
    )
  ]

  diagnostics <- list()
  diagnostic_index <- 1L

  for (candidate in candidates) {

    check <- validate_candidate(candidate)

    diagnostics[[diagnostic_index]] <- data.frame(
      candidate = candidate,
      size_bytes = file.info(candidate)$size,
      status = check$reason,
      stringsAsFactors = FALSE
    )

    diagnostic_index <- diagnostic_index + 1L

    if (isTRUE(check$valid)) {

      selected_file <- normalizePath(
        candidate,
        winslash = "/",
        mustWork = TRUE
      )

      cat(
        "Valid transition dataset selected:\n",
        selected_file,
        "\nRows: ",
        nrow(check$data),
        "\nColumns: ",
        ncol(check$data),
        "\n\n",
        sep = ""
      )

      return(selected_file)
    }
  }

  diagnostic_table <- dplyr::bind_rows(
    diagnostics
  )

  diagnostic_file <- file.path(
    OUTPUT_DIR,
    "transition_file_diagnostics.csv"
  )

  write.csv(
    diagnostic_table,
    diagnostic_file,
    row.names = FALSE
  )

  stop(
    paste0(
      "Transition files were found but none was valid. See ",
      diagnostic_file
    )
  )
}

DATA_FILE <- locate_transition_file()

RAW_TRANSITIONS <- read.csv(
  DATA_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# -----------------------------------------------------------------------------
# 4. DATA PREPARATION AND STAGE-SPELL AGE
# -----------------------------------------------------------------------------

STAGE_LEVELS_SIX <- c(
  "Pre-negotiation",
  "Ceasefire",
  "Substantive Partial",
  "Substantive Comprehensive",
  "Implementation",
  "Renegotiation/Other"
)

add_stage_spell_age <- function(
  data,
  stage_levels
) {
  output <- data %>%
    mutate(
      PP = as.character(PP),
      origin_stage = factor(
        as.character(origin_stage),
        levels = stage_levels
      ),
      destination_stage = factor(
        as.character(destination_stage),
        levels = stage_levels
      ),
      D = as.integer(D),
      transition_index = as.numeric(
        transition_index
      ),
      transition_id = as.character(
        transition_id
      )
    ) %>%
    filter(
      !is.na(PP),
      nzchar(PP),
      !is.na(origin_stage),
      !is.na(destination_stage),
      D %in% c(0L, 1L),
      is.finite(transition_index),
      transition_index >= 1
    ) %>%
    arrange(
      PP,
      transition_index
    ) %>%
    group_by(PP) %>%
    mutate(
      origin_character = as.character(
        origin_stage
      ),
      spell_start = row_number() == 1L |
        origin_character != dplyr::lag(
          origin_character,
          default = "__START__"
        ),
      spell_id = cumsum(
        spell_start
      )
    ) %>%
    group_by(
      PP,
      spell_id
    ) %>%
    mutate(
      spell_age = row_number()
    ) %>%
    ungroup() %>%
    group_by(PP) %>%
    mutate(
      process_length = n()
    ) %>%
    ungroup() %>%
    mutate(
      log_spell_age = log(
        spell_age
      ),
      log_transition_index = log(
        transition_index
      ),
      process_weight_equal = 1 /
        process_length
    ) %>%
    select(
      -origin_character,
      -spell_start
    )

  output
}

TRANSITIONS <- RAW_TRANSITIONS %>%
  mutate(
    PP = as.character(PP),
    origin_stage = as.character(
      origin_stage
    ),
    destination_stage = as.character(
      destination_stage
    )
  )

observed_stages <- unique(
  c(
    TRANSITIONS$origin_stage,
    TRANSITIONS$destination_stage
  )
)

additional_stages <- setdiff(
  observed_stages,
  STAGE_LEVELS_SIX
)

STAGE_LEVELS <- c(
  STAGE_LEVELS_SIX[
    STAGE_LEVELS_SIX %in%
      observed_stages
  ],
  sort(additional_stages)
)

TRANSITIONS <- add_stage_spell_age(
  TRANSITIONS,
  STAGE_LEVELS
)

if (nrow(TRANSITIONS) == 0L) {
  stop(
    "No valid transitions remain after preparation."
  )
}

cat(
  "Transitions: ",
  nrow(TRANSITIONS),
  "\nPeace processes: ",
  dplyr::n_distinct(
    TRANSITIONS$PP
  ),
  "\nMaximum stage-spell age: ",
  max(TRANSITIONS$spell_age),
  "\n\n",
  sep = ""
)

DATA_SUMMARY <- data.frame(
  input_file = DATA_FILE,
  transitions = nrow(TRANSITIONS),
  peace_processes = dplyr::n_distinct(
    TRANSITIONS$PP
  ),
  maximum_transition_index = max(
    TRANSITIONS$transition_index
  ),
  maximum_spell_age = max(
    TRANSITIONS$spell_age
  ),
  same_date_transitions = if (
    "same_date_transition" %in%
      names(TRANSITIONS)
  ) {
    sum(
      as.numeric(
        TRANSITIONS$same_date_transition
      ) == 1,
      na.rm = TRUE
    )
  } else {
    NA_integer_
  },
  analysis_mode = ANALYSIS_MODE,
  stringsAsFactors = FALSE
)

save_csv(
  DATA_SUMMARY,
  "00_data_summary.csv"
)


# -----------------------------------------------------------------------------
# 5. CONSTRAINED OPPORTUNITY MODELS
# -----------------------------------------------------------------------------

model_parameter_count <- function(
  model,
  stage_count
) {
  switch(
    model,
    homogeneous = stage_count,
    stage_spell = stage_count + 1L,
    global_alpha_beta = stage_count + 2L,
    dual_clock = stage_count + 2L,
    stage_specific_spell = 2L * stage_count,
    stop(
      "Unknown constrained model: ",
      model
    )
  )
}

decode_parameters <- function(
  parameters,
  model,
  stage_levels
) {
  stage_count <- length(stage_levels)

  kappa <- setNames(
    parameters[
      seq_len(stage_count)
    ],
    stage_levels
  )

  output <- list(
    kappa = kappa,
    eta = exp(kappa),
    alpha = 0,
    beta_spell = 0,
    beta_global = 0,
    beta_by_stage = setNames(
      rep(0, stage_count),
      stage_levels
    )
  )

  if (model == "stage_spell") {
    output$beta_spell <- parameters[
      stage_count + 1L
    ]
  }

  if (model == "global_alpha_beta") {
    output$alpha <- parameters[
      stage_count + 1L
    ]

    output$beta_global <- parameters[
      stage_count + 2L
    ]
  }

  if (model == "dual_clock") {
    output$alpha <- parameters[
      stage_count + 1L
    ]

    output$beta_spell <- parameters[
      stage_count + 2L
    ]
  }

  if (model == "stage_specific_spell") {
    output$beta_by_stage <- setNames(
      parameters[
        (stage_count + 1L):
          (2L * stage_count)
      ],
      stage_levels
    )
  }

  output
}

constrained_negative_loglikelihood <- function(
  parameters,
  data,
  model,
  stage_levels,
  weights
) {
  decoded <- decode_parameters(
    parameters,
    model,
    stage_levels
  )

  origin <- as.character(
    data$origin_stage
  )

  linear_predictor <- decoded$kappa[
    origin
  ]

  if (model == "stage_spell") {
    linear_predictor <- linear_predictor -
      decoded$beta_spell *
      data$log_spell_age
  }

  if (model == "global_alpha_beta") {
    linear_predictor <- linear_predictor +
      decoded$alpha *
      data$log_transition_index -
      decoded$beta_global *
      (
        data$transition_index - 1
      )
  }

  if (model == "dual_clock") {
    linear_predictor <- linear_predictor +
      decoded$alpha *
      data$log_transition_index -
      decoded$beta_spell *
      data$log_spell_age
  }

  if (model == "stage_specific_spell") {
    linear_predictor <- linear_predictor -
      decoded$beta_by_stage[
        origin
      ] *
      data$log_spell_age
  }

  mu <- exp(
    pmin(
      linear_predictor,
      700
    )
  )

  loglikelihood <- sum(
    weights * (
      data$D *
        stable_log_q(mu) +
      (1 - data$D) *
        (-mu)
    )
  )

  if (!is.finite(loglikelihood)) {
    1e100
  } else {
    -loglikelihood
  }
}

empirical_kappa_start <- function(
  data,
  stage_levels
) {
  overall_probability <- clip_probability(
    mean(data$D)
  )

  probabilities <- vapply(
    stage_levels,
    function(stage) {
      selected <- as.character(
        data$origin_stage
      ) == stage

      if (sum(selected) == 0L) {
        overall_probability
      } else {
        clip_probability(
          mean(
            data$D[selected]
          )
        )
      }
    },
    numeric(1)
  )

  log(
    -log(
      1 - probabilities
    )
  )
}

make_start_values <- function(
  data,
  stage_levels,
  model,
  embedded_fits,
  number_of_starts
) {
  kappa_start <- empirical_kappa_start(
    data,
    stage_levels
  )

  starts <- list()

  if (model == "homogeneous") {
    starts[[1L]] <- kappa_start
  }

  if (model == "stage_spell") {
    beta_starts <- c(
      0,
      0.10,
      0.25,
      0.50,
      0.80,
      1.20,
      2.00
    )

    starts <- lapply(
      beta_starts,
      function(beta_value) {
        c(
          kappa_start,
          beta_value
        )
      }
    )

    if (
      !is.null(
        embedded_fits$homogeneous
      )
    ) {
      starts <- c(
        list(
          c(
            embedded_fits$homogeneous$kappa,
            0
          )
        ),
        starts
      )
    }
  }

  if (model == "global_alpha_beta") {
    alpha_starts <- c(
      0,
      0.05,
      0.20,
      0.50,
      1.00
    )

    beta_starts <- c(
      0,
      0.005,
      0.02,
      0.08,
      0.25
    )

    start_index <- 1L

    for (alpha_value in alpha_starts) {
      for (beta_value in beta_starts) {
        starts[[start_index]] <- c(
          kappa_start,
          alpha_value,
          beta_value
        )

        start_index <- start_index + 1L
      }
    }

    if (
      !is.null(
        embedded_fits$homogeneous
      )
    ) {
      starts <- c(
        list(
          c(
            embedded_fits$homogeneous$kappa,
            0,
            0
          )
        ),
        starts
      )
    }
  }

  if (model == "dual_clock") {
    alpha_starts <- c(
      0,
      0.02,
      0.10,
      0.30,
      0.75
    )

    beta_starts <- c(
      0,
      0.10,
      0.30,
      0.60,
      1.00
    )

    start_index <- 1L

    for (alpha_value in alpha_starts) {
      for (beta_value in beta_starts) {
        starts[[start_index]] <- c(
          kappa_start,
          alpha_value,
          beta_value
        )

        start_index <- start_index + 1L
      }
    }

    if (
      !is.null(
        embedded_fits$stage_spell
      )
    ) {
      starts <- c(
        list(
          c(
            embedded_fits$stage_spell$kappa,
            0,
            embedded_fits$stage_spell$beta_spell
          )
        ),
        starts
      )
    }

    if (
      !is.null(
        embedded_fits$homogeneous
      )
    ) {
      starts <- c(
        list(
          c(
            embedded_fits$homogeneous$kappa,
            0,
            0
          )
        ),
        starts
      )
    }
  }

  if (model == "stage_specific_spell") {
    common_beta <- if (
      !is.null(
        embedded_fits$stage_spell
      )
    ) {
      embedded_fits$stage_spell$beta_spell
    } else {
      0.50
    }

    starts <- list(
      c(
        kappa_start,
        rep(
          common_beta,
          length(stage_levels)
        )
      ),
      c(
        kappa_start,
        rep(
          0,
          length(stage_levels)
        )
      ),
      c(
        kappa_start,
        seq(
          0.10,
          0.90,
          length.out =
            length(stage_levels)
        )
      )
    )

    if (
      !is.null(
        embedded_fits$stage_spell
      )
    ) {
      starts <- c(
        list(
          c(
            embedded_fits$stage_spell$kappa,
            rep(
              embedded_fits$stage_spell$beta_spell,
              length(stage_levels)
            )
          )
        ),
        starts
      )
    }
  }

  if (length(starts) > number_of_starts) {
    selected <- unique(
      round(
        seq(
          1,
          length(starts),
          length.out =
            number_of_starts
        )
      )
    )

    starts <- starts[selected]
  }

  starts
}

fit_constrained_model <- function(
  data,
  stage_levels,
  model = c(
    "homogeneous",
    "stage_spell",
    "global_alpha_beta",
    "dual_clock",
    "stage_specific_spell"
  ),
  weights = rep(1, nrow(data)),
  embedded_fits = list(),
  number_of_starts = MULTISTARTS
) {
  model <- match.arg(model)

  stage_count <- length(stage_levels)
  parameter_count <- model_parameter_count(
    model,
    stage_count
  )

  lower <- rep(
    KAPPA_LOWER,
    stage_count
  )

  upper <- rep(
    KAPPA_UPPER,
    stage_count
  )

  if (model == "stage_spell") {
    lower <- c(
      lower,
      0
    )

    upper <- c(
      upper,
      DYNAMIC_UPPER
    )
  }

  if (
    model %in% c(
      "global_alpha_beta",
      "dual_clock"
    )
  ) {
    lower <- c(
      lower,
      0,
      0
    )

    upper <- c(
      upper,
      DYNAMIC_UPPER,
      DYNAMIC_UPPER
    )
  }

  if (model == "stage_specific_spell") {
    lower <- c(
      lower,
      rep(
        0,
        stage_count
      )
    )

    upper <- c(
      upper,
      rep(
        DYNAMIC_UPPER,
        stage_count
      )
    )
  }

  starts <- make_start_values(
    data,
    stage_levels,
    model,
    embedded_fits,
    number_of_starts
  )

  candidate_rows <- list()
  candidate_index <- 1L

  for (start in starts) {

    bounded_start <- pmin(
      pmax(
        start,
        lower
      ),
      upper
    )

    start_value <- constrained_negative_loglikelihood(
      bounded_start,
      data,
      model,
      stage_levels,
      weights
    )

    candidate_rows[[candidate_index]] <- list(
      parameters = bounded_start,
      value = start_value,
      convergence = 0L,
      source = "embedded_or_start"
    )

    candidate_index <- candidate_index + 1L

    fitted <- tryCatch(
      optim(
        par = bounded_start,
        fn = constrained_negative_loglikelihood,
        data = data,
        model = model,
        stage_levels = stage_levels,
        weights = weights,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper,
        control = list(
          maxit = 6000,
          factr = 1e7,
          pgtol = 1e-8
        )
      ),
      error = function(e) {
        list(
          par = bounded_start,
          value = Inf,
          convergence = 999L,
          message = conditionMessage(e)
        )
      }
    )

    candidate_rows[[candidate_index]] <- list(
      parameters = fitted$par,
      value = fitted$value,
      convergence = fitted$convergence,
      source = "optimized"
    )

    candidate_index <- candidate_index + 1L
  }

  values <- vapply(
    candidate_rows,
    function(candidate) {
      if (
        is.null(candidate$value) ||
        !is.finite(candidate$value)
      ) {
        Inf
      } else {
        candidate$value
      }
    },
    numeric(1)
  )

  if (all(!is.finite(values))) {
    stop(
      "All optimizations failed for ",
      model,
      "."
    )
  }

  best <- candidate_rows[[which.min(values)]]

  decoded <- decode_parameters(
    best$parameters,
    model,
    stage_levels
  )

  list(
    converged = is.finite(best$value),
    model = model,
    parameters = best$parameters,
    kappa = decoded$kappa,
    eta = decoded$eta,
    alpha = decoded$alpha,
    beta_spell = decoded$beta_spell,
    beta_global = decoded$beta_global,
    beta_by_stage = decoded$beta_by_stage,
    loglik = -best$value,
    parameter_count = parameter_count,
    convergence_code = best$convergence,
    solution_source = best$source,
    stage_levels = stage_levels
  )
}

predict_constrained_model <- function(
  fit,
  new_data
) {
  origin <- as.character(
    new_data$origin_stage
  )

  linear_predictor <- fit$kappa[
    origin
  ]

  if (fit$model == "stage_spell") {
    linear_predictor <- linear_predictor -
      fit$beta_spell *
      new_data$log_spell_age
  }

  if (fit$model == "global_alpha_beta") {
    linear_predictor <- linear_predictor +
      fit$alpha *
      new_data$log_transition_index -
      fit$beta_global *
      (
        new_data$transition_index - 1
      )
  }

  if (fit$model == "dual_clock") {
    linear_predictor <- linear_predictor +
      fit$alpha *
      new_data$log_transition_index -
      fit$beta_spell *
      new_data$log_spell_age
  }

  if (fit$model == "stage_specific_spell") {
    linear_predictor <- linear_predictor -
      fit$beta_by_stage[
        origin
      ] *
      new_data$log_spell_age
  }

  mu <- exp(
    pmin(
      linear_predictor,
      700
    )
  )

  clip_probability(
    -expm1(-mu)
  )
}


# -----------------------------------------------------------------------------
# 6. EARLIER GLOBAL-INDEX KERNEL SPAOM
# -----------------------------------------------------------------------------

old_kernel <- function(
  r,
  tau
) {
  x <- r / tau

  x * exp(
    1 - x
  )
}

old_kernel_nll <- function(
  parameters,
  tau,
  data,
  stage_levels,
  weights
) {
  stage_count <- length(stage_levels)

  kappa <- setNames(
    parameters[
      seq_len(stage_count)
    ],
    stage_levels
  )

  theta <- parameters[
    stage_count + 1L
  ]

  linear_predictor <- kappa[
    as.character(
      data$origin_stage
    )
  ] +
    theta *
    old_kernel(
      data$transition_index,
      tau
    )

  mu <- exp(
    pmin(
      linear_predictor,
      700
    )
  )

  loglikelihood <- sum(
    weights * (
      data$D *
        stable_log_q(mu) +
      (1 - data$D) *
        (-mu)
    )
  )

  if (!is.finite(loglikelihood)) {
    1e100
  } else {
    -loglikelihood
  }
}

fit_old_kernel_at_tau <- function(
  tau,
  data,
  stage_levels,
  homogeneous_fit,
  weights,
  number_of_starts = 7L
) {
  theta_starts <- c(
    0,
    0.05,
    0.20,
    0.50,
    1.00,
    2.00,
    4.00
  )

  starts <- lapply(
    theta_starts,
    function(theta_value) {
      c(
        homogeneous_fit$kappa,
        theta_value
      )
    }
  )

  if (length(starts) > number_of_starts) {
    starts <- starts[
      unique(
        round(
          seq(
            1,
            length(starts),
            length.out =
              number_of_starts
          )
        )
      )
    ]
  }

  lower <- c(
    rep(
      KAPPA_LOWER,
      length(stage_levels)
    ),
    0
  )

  upper <- c(
    rep(
      KAPPA_UPPER,
      length(stage_levels)
    ),
    DYNAMIC_UPPER
  )

  candidates <- list()
  candidate_index <- 1L

  for (start in starts) {

    candidates[[candidate_index]] <- list(
      par = start,
      value = old_kernel_nll(
        start,
        tau,
        data,
        stage_levels,
        weights
      )
    )

    candidate_index <- candidate_index + 1L

    fit <- tryCatch(
      optim(
        par = start,
        fn = old_kernel_nll,
        tau = tau,
        data = data,
        stage_levels = stage_levels,
        weights = weights,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper,
        control = list(
          maxit = 5000,
          factr = 1e7,
          pgtol = 1e-8
        )
      ),
      error = function(e) {
        list(
          par = start,
          value = Inf
        )
      }
    )

    candidates[[candidate_index]] <- fit
    candidate_index <- candidate_index + 1L
  }

  values <- vapply(
    candidates,
    function(candidate) {
      if (
        is.null(candidate$value) ||
        !is.finite(candidate$value)
      ) {
        Inf
      } else {
        candidate$value
      }
    },
    numeric(1)
  )

  if (all(!is.finite(values))) {
    return(
      list(
        converged = FALSE,
        loglik = -Inf
      )
    )
  }

  best <- candidates[[which.min(values)]]

  stage_count <- length(stage_levels)

  list(
    converged = TRUE,
    loglik = -best$value,
    kappa = setNames(
      best$par[
        seq_len(stage_count)
      ],
      stage_levels
    ),
    eta = exp(
      setNames(
        best$par[
          seq_len(stage_count)
        ],
        stage_levels
      )
    ),
    theta = best$par[
      stage_count + 1L
    ]
  )
}

fit_old_kernel_model <- function(
  data,
  stage_levels,
  homogeneous_fit,
  tau_grid = OLD_TAU_GRID,
  weights = rep(1, nrow(data))
) {
  profile_rows <- lapply(
    tau_grid,
    function(tau_value) {
      fit <- fit_old_kernel_at_tau(
        tau = tau_value,
        data = data,
        stage_levels = stage_levels,
        homogeneous_fit = homogeneous_fit,
        weights = weights
      )

      data.frame(
        tau = tau_value,
        loglik = fit$loglik,
        converged = isTRUE(
          fit$converged
        ),
        stringsAsFactors = FALSE
      )
    }
  )

  profile <- bind_rows(
    profile_rows
  )

  valid <- which(
    profile$converged &
      is.finite(profile$loglik)
  )

  if (length(valid) == 0L) {
    stop(
      "Earlier kernel SPAOM failed at every tau value."
    )
  }

  best_index <- valid[
    which.max(
      profile$loglik[valid]
    )
  ]

  tau_hat <- profile$tau[
    best_index
  ]

  final_fit <- fit_old_kernel_at_tau(
    tau = tau_hat,
    data = data,
    stage_levels = stage_levels,
    homogeneous_fit = homogeneous_fit,
    weights = weights
  )

  list(
    converged = final_fit$converged,
    model = "old_kernel",
    kappa = final_fit$kappa,
    eta = final_fit$eta,
    theta = final_fit$theta,
    tau = tau_hat,
    loglik = final_fit$loglik,
    parameter_count =
      length(stage_levels) + 2L,
    profile = profile,
    stage_levels = stage_levels
  )
}

predict_old_kernel_model <- function(
  fit,
  new_data
) {
  linear_predictor <- fit$kappa[
    as.character(
      new_data$origin_stage
    )
  ] +
    fit$theta *
    old_kernel(
      new_data$transition_index,
      fit$tau
    )

  mu <- exp(
    pmin(
      linear_predictor,
      700
    )
  )

  clip_probability(
    -expm1(-mu)
  )
}


# -----------------------------------------------------------------------------
# 7. FLEXIBLE GLM COMPARATORS
# -----------------------------------------------------------------------------

prepare_glm_data <- function(
  data,
  r_center = NULL,
  r_scale = NULL
) {
  if (is.null(r_center)) {
    r_center <- mean(
      data$transition_index
    )
  }

  if (is.null(r_scale)) {
    r_scale <- sd(
      data$transition_index
    )

    if (
      !is.finite(r_scale) ||
      r_scale <= 0
    ) {
      r_scale <- 1
    }
  }

  output <- data %>%
    mutate(
      origin_stage = factor(
        origin_stage,
        levels = levels(
          data$origin_stage
        )
      ),
      r_scaled = (
        transition_index -
          r_center
      ) /
        r_scale
    )

  list(
    data = output,
    r_center = r_center,
    r_scale = r_scale
  )
}

fit_glm_comparator <- function(
  data,
  comparator = c(
    "spline_spell",
    "quadratic_r"
  ),
  weights = rep(1, nrow(data))
) {
  comparator <- match.arg(
    comparator
  )

  prepared <- prepare_glm_data(
    data
  )

  formula_object <- if (
    comparator == "spline_spell"
  ) {
    as.formula(
      paste0(
        "D ~ 0 + origin_stage + ",
        "splines::ns(log_spell_age, df = ",
        SPLINE_DF,
        ")"
      )
    )
  } else {
    D ~ 0 + origin_stage +
      r_scaled +
      I(r_scaled^2)
  }

  fit <- tryCatch(
    suppressWarnings(
      glm(
        formula_object,
        data = prepared$data,
        family = binomial(
          link = "cloglog"
        ),
        weights = weights,
        control = glm.control(
          maxit = 100
        )
      )
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    stop(
      comparator,
      " comparator failed."
    )
  }

  coefficients <- coef(fit)

  coefficients[
    !is.finite(coefficients)
  ] <- 0

  fit$coefficients <- coefficients

  fitted_probability <- clip_probability(
    fitted(fit)
  )

  loglikelihood <- sum(
    weights * (
      data$D *
        log(
          fitted_probability
        ) +
      (1 - data$D) *
        log(
          1 -
            fitted_probability
        )
    )
  )

  list(
    fit = fit,
    model = comparator,
    loglik = loglikelihood,
    parameter_count = sum(
      is.finite(coefficients)
    ),
    r_center = prepared$r_center,
    r_scale = prepared$r_scale,
    stage_levels = levels(
      data$origin_stage
    )
  )
}

predict_glm_comparator <- function(
  fit,
  new_data
) {
  prepared <- prepare_glm_data(
    new_data,
    r_center = fit$r_center,
    r_scale = fit$r_scale
  )

  prediction <- suppressWarnings(
    predict(
      fit$fit,
      newdata = prepared$data,
      type = "response"
    )
  )

  clip_probability(
    prediction
  )
}


# -----------------------------------------------------------------------------
# 8. DESTINATION COMPONENT AND PERFORMANCE METRICS
# -----------------------------------------------------------------------------

fit_destination_matrix <- function(
  data,
  stage_levels,
  smoothing =
    DESTINATION_SMOOTHING
) {
  probability_matrix <- matrix(
    0,
    nrow = length(stage_levels),
    ncol = length(stage_levels),
    dimnames = list(
      stage_levels,
      stage_levels
    )
  )

  for (origin in stage_levels) {
    admissible <- setdiff(
      stage_levels,
      origin
    )

    counts <- setNames(
      rep(
        smoothing,
        length(admissible)
      ),
      admissible
    )

    selected <- data$D == 1L &
      as.character(
        data$origin_stage
      ) == origin

    if (sum(selected) > 0L) {
      observed_counts <- table(
        as.character(
          data$destination_stage[
            selected
          ]
        )
      )

      common <- intersect(
        names(observed_counts),
        admissible
      )

      counts[common] <- counts[common] +
        as.numeric(
          observed_counts[common]
        )
    }

    probability_matrix[
      origin,
      admissible
    ] <- counts /
      sum(counts)
  }

  probability_matrix
}

combine_transition_probabilities <- function(
  data,
  departure_probability,
  destination_matrix,
  stage_levels
) {
  output <- matrix(
    0,
    nrow = nrow(data),
    ncol = length(stage_levels),
    dimnames = list(
      NULL,
      stage_levels
    )
  )

  for (
    row_index in
    seq_len(
      nrow(data)
    )
  ) {
    origin <- as.character(
      data$origin_stage[
        row_index
      ]
    )

    output[
      row_index,
      origin
    ] <- 1 -
      departure_probability[
        row_index
      ]

    admissible <- setdiff(
      stage_levels,
      origin
    )

    output[
      row_index,
      admissible
    ] <- departure_probability[
      row_index
    ] *
      destination_matrix[
        origin,
        admissible
      ]
  }

  output /
    rowSums(output)
}

calibration_intercept_slope <- function(
  actual,
  predicted
) {
  predicted <- clip_probability(
    predicted
  )

  fit <- tryCatch(
    suppressWarnings(
      glm(
        actual ~ qlogis(predicted),
        family = binomial()
      )
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(
      c(
        intercept = NA_real_,
        slope = NA_real_
      )
    )
  }

  coefficients <- unname(
    coef(fit)
  )

  if (length(coefficients) < 2L) {
    return(
      c(
        intercept = NA_real_,
        slope = NA_real_
      )
    )
  }

  c(
    intercept = coefficients[1L],
    slope = coefficients[2L]
  )
}

evaluate_predictions <- function(
  data,
  departure_probability,
  transition_probability_matrix,
  stage_levels
) {
  departure_probability <- clip_probability(
    departure_probability
  )

  calibration <- calibration_intercept_slope(
    data$D,
    departure_probability
  )

  row_index <- seq_len(
    nrow(data)
  )

  actual_stage_index <- match(
    as.character(
      data$destination_stage
    ),
    stage_levels
  )

  observed_probability <- transition_probability_matrix[
    cbind(
      row_index,
      actual_stage_index
    )
  ]

  one_hot <- matrix(
    0,
    nrow = nrow(data),
    ncol = length(stage_levels)
  )

  one_hot[
    cbind(
      row_index,
      actual_stage_index
    )
  ] <- 1

  contributions <- data.frame(
    PP = as.character(data$PP),
    transition_id = as.character(
      data$transition_id
    ),
    D = data$D,
    q = departure_probability,
    departure_logloss = -(
      data$D *
        log(
          departure_probability
        ) +
      (1 - data$D) *
        log(
          1 -
            departure_probability
        )
    ),
    departure_brier = (
      data$D -
        departure_probability
    )^2,
    full_logloss = -log(
      pmax(
        observed_probability,
        PROBABILITY_FLOOR
      )
    ),
    full_brier = rowSums(
      (
        one_hot -
          transition_probability_matrix
      )^2
    ),
    stringsAsFactors = FALSE
  )

  metrics <- c(
    departure_logloss = mean(
      contributions$departure_logloss
    ),
    departure_brier = mean(
      contributions$departure_brier
    ),
    calibration_intercept = unname(
      calibration["intercept"]
    ),
    calibration_slope = unname(
      calibration["slope"]
    ),
    full_logloss = mean(
      contributions$full_logloss
    ),
    full_brier = mean(
      contributions$full_brier
    )
  )

  list(
    metrics = metrics,
    contributions = contributions
  )
}


# -----------------------------------------------------------------------------
# 9. FULL-DATA FITS
# -----------------------------------------------------------------------------

cat(
  "Fitting full-data SS-SPAOM and comparators...\n"
)

FIT_HOMOGENEOUS <- fit_constrained_model(
  TRANSITIONS,
  STAGE_LEVELS,
  model = "homogeneous"
)

FIT_STAGE_SPELL <- fit_constrained_model(
  TRANSITIONS,
  STAGE_LEVELS,
  model = "stage_spell",
  embedded_fits = list(
    homogeneous =
      FIT_HOMOGENEOUS
  )
)

FIT_GLOBAL_ALPHA_BETA <- fit_constrained_model(
  TRANSITIONS,
  STAGE_LEVELS,
  model = "global_alpha_beta",
  embedded_fits = list(
    homogeneous =
      FIT_HOMOGENEOUS
  )
)

FIT_DUAL_CLOCK <- fit_constrained_model(
  TRANSITIONS,
  STAGE_LEVELS,
  model = "dual_clock",
  embedded_fits = list(
    homogeneous =
      FIT_HOMOGENEOUS,
    stage_spell =
      FIT_STAGE_SPELL
  )
)

FIT_STAGE_SPECIFIC <- fit_constrained_model(
  TRANSITIONS,
  STAGE_LEVELS,
  model =
    "stage_specific_spell",
  embedded_fits = list(
    homogeneous =
      FIT_HOMOGENEOUS,
    stage_spell =
      FIT_STAGE_SPELL
  )
)

FIT_OLD_KERNEL <- fit_old_kernel_model(
  TRANSITIONS,
  STAGE_LEVELS,
  homogeneous_fit =
    FIT_HOMOGENEOUS
)

FIT_SPLINE_SPELL <- fit_glm_comparator(
  TRANSITIONS,
  comparator =
    "spline_spell"
)

FIT_QUADRATIC_R <- fit_glm_comparator(
  TRANSITIONS,
  comparator =
    "quadratic_r"
)

FULL_MODEL_COMPARISON <- data.frame(
  model = c(
    "Stage-Spell SPAOM",
    "Origin-specific homogeneous",
    "Dual-clock model",
    "Global alpha-beta model",
    "Earlier kernel SPAOM",
    "Spline spell-age model",
    "Quadratic global-r model",
    "Stage-specific spell sensitivity"
  ),
  loglik = c(
    FIT_STAGE_SPELL$loglik,
    FIT_HOMOGENEOUS$loglik,
    FIT_DUAL_CLOCK$loglik,
    FIT_GLOBAL_ALPHA_BETA$loglik,
    FIT_OLD_KERNEL$loglik,
    FIT_SPLINE_SPELL$loglik,
    FIT_QUADRATIC_R$loglik,
    FIT_STAGE_SPECIFIC$loglik
  ),
  parameters = c(
    FIT_STAGE_SPELL$parameter_count,
    FIT_HOMOGENEOUS$parameter_count,
    FIT_DUAL_CLOCK$parameter_count,
    FIT_GLOBAL_ALPHA_BETA$parameter_count,
    FIT_OLD_KERNEL$parameter_count,
    FIT_SPLINE_SPELL$parameter_count,
    FIT_QUADRATIC_R$parameter_count,
    FIT_STAGE_SPECIFIC$parameter_count
  ),
  stringsAsFactors = FALSE
) %>%
  mutate(
    AIC = -2 * loglik +
      2 * parameters,
    BIC_cluster = -2 * loglik +
      parameters *
      log(
        dplyr::n_distinct(
          TRANSITIONS$PP
        )
      ),
    delta_AIC = AIC -
      min(AIC),
    delta_BIC_cluster =
      BIC_cluster -
      min(BIC_cluster)
  ) %>%
  arrange(AIC)

PRIMARY_ESTIMATES <- data.frame(
  beta_spell_hat =
    FIT_STAGE_SPELL$beta_spell,
  LR_vs_homogeneous = max(
    0,
    2 * (
      FIT_STAGE_SPELL$loglik -
        FIT_HOMOGENEOUS$loglik
    )
  ),
  relative_intensity_s2 =
    2^(
      -FIT_STAGE_SPELL$beta_spell
    ),
  relative_intensity_s5 =
    5^(
      -FIT_STAGE_SPELL$beta_spell
    ),
  relative_intensity_s10 =
    10^(
      -FIT_STAGE_SPELL$beta_spell
    ),
  loglik =
    FIT_STAGE_SPELL$loglik,
  transitions =
    nrow(TRANSITIONS),
  peace_processes =
    dplyr::n_distinct(
      TRANSITIONS$PP
    ),
  analysis_mode =
    ANALYSIS_MODE,
  stringsAsFactors = FALSE
)

ETA_ESTIMATES <- data.frame(
  origin_stage =
    STAGE_LEVELS,
  eta_hat = as.numeric(
    FIT_STAGE_SPELL$eta[
      STAGE_LEVELS
    ]
  ),
  departure_at_spell_age_1 =
    1 -
    exp(
      -as.numeric(
        FIT_STAGE_SPELL$eta[
          STAGE_LEVELS
        ]
      )
    ),
  stringsAsFactors = FALSE
)

STAGE_SPECIFIC_BETA <- data.frame(
  origin_stage =
    STAGE_LEVELS,
  beta_stage_hat =
    as.numeric(
      FIT_STAGE_SPECIFIC$
        beta_by_stage[
          STAGE_LEVELS
        ]
    ),
  stringsAsFactors = FALSE
)

save_csv(
  PRIMARY_ESTIMATES,
  "01_primary_stage_spell_estimates.csv"
)

save_csv(
  ETA_ESTIMATES,
  "02_origin_specific_eta_estimates.csv"
)

save_csv(
  STAGE_SPECIFIC_BETA,
  "03_stage_specific_beta_sensitivity.csv"
)

save_csv(
  FULL_MODEL_COMPARISON,
  "04_full_data_model_comparison.csv"
)

save_csv(
  FIT_OLD_KERNEL$profile,
  "05_earlier_kernel_profile.csv"
)


# -----------------------------------------------------------------------------
# 10. PARAMETRIC-BOOTSTRAP TEST OF H0: BETA_S = 0
# -----------------------------------------------------------------------------

NULL_BOOTSTRAP_FILE <- file.path(
  CHECKPOINT_DIR,
  paste0(
    "SS_Null_Bootstrap_",
    ANALYSIS_MODE,
    ".csv"
  )
)

run_null_bootstrap <- function() {

  if (file.exists(NULL_BOOTSTRAP_FILE)) {
    results <- read.csv(
      NULL_BOOTSTRAP_FILE,
      stringsAsFactors = FALSE
    )
  } else {
    results <- data.frame(
      replication_id = integer(0),
      LR = numeric(0),
      beta_spell_hat = numeric(0),
      converged = logical(0),
      stringsAsFactors = FALSE
    )
  }

  completed <- results$replication_id

  null_probability <- predict_constrained_model(
    FIT_HOMOGENEOUS,
    TRANSITIONS
  )

  set.seed(
    BASE_SEED + 10000L
  )

  bootstrap_seeds <- sample.int(
    .Machine$integer.max,
    NULL_BOOTSTRAP_B
  )

  for (
    replication_id in
    seq_len(
      NULL_BOOTSTRAP_B
    )
  ) {
    if (
      replication_id %in%
        completed
    ) {
      next
    }

    cat(
      "Null bootstrap ",
      replication_id,
      " of ",
      NULL_BOOTSTRAP_B,
      "\n",
      sep = ""
    )

    set.seed(
      bootstrap_seeds[
        replication_id
      ]
    )

    bootstrap_data <- TRANSITIONS

    bootstrap_data$D <- rbinom(
      nrow(bootstrap_data),
      size = 1L,
      prob = null_probability
    )

    null_fit <- tryCatch(
      fit_constrained_model(
        bootstrap_data,
        STAGE_LEVELS,
        model = "homogeneous",
        number_of_starts =
          max(
            4L,
            MULTISTARTS %/% 2L
          )
      ),
      error = function(e) NULL
    )

    spell_fit <- tryCatch(
      fit_constrained_model(
        bootstrap_data,
        STAGE_LEVELS,
        model = "stage_spell",
        embedded_fits = list(
          homogeneous =
            null_fit
        ),
        number_of_starts =
          MULTISTARTS
      ),
      error = function(e) {
        append_error(
          "null_bootstrap",
          replication_id,
          conditionMessage(e)
        )
        NULL
      }
    )

    converged <-
      !is.null(null_fit) &&
      !is.null(spell_fit) &&
      is.finite(null_fit$loglik) &&
      is.finite(spell_fit$loglik)

    LR <- if (converged) {
      max(
        0,
        2 * (
          spell_fit$loglik -
            null_fit$loglik
        )
      )
    } else {
      NA_real_
    }

    results <- bind_rows(
      results,
      data.frame(
        replication_id =
          replication_id,
        LR = LR,
        beta_spell_hat = if (
          converged
        ) {
          spell_fit$beta_spell
        } else {
          NA_real_
        },
        converged =
          converged,
        stringsAsFactors = FALSE
      )
    )

    write.csv(
      results,
      NULL_BOOTSTRAP_FILE,
      row.names = FALSE
    )
  }

  observed_LR <-
    PRIMARY_ESTIMATES$
      LR_vs_homogeneous[1L]

  successful_LR <- results$LR[
    is.finite(results$LR)
  ]

  p_value <- (
    1 +
      sum(
        successful_LR >=
          observed_LR
      )
  ) /
    (
      1 +
        length(
          successful_LR
        )
    )

  list(
    replications = results,
    summary = data.frame(
      observed_LR =
        observed_LR,
      bootstrap_p_value =
        p_value,
      successful_replications =
        length(successful_LR),
      requested_replications =
        NULL_BOOTSTRAP_B,
      stringsAsFactors = FALSE
    )
  )
}

NULL_BOOTSTRAP_RESULT <-
  run_null_bootstrap()

save_csv(
  NULL_BOOTSTRAP_RESULT$
    replications,
  "06_dynamic_test_bootstrap_replications.csv"
)

save_csv(
  NULL_BOOTSTRAP_RESULT$
    summary,
  "07_dynamic_test_bootstrap_summary.csv"
)


# -----------------------------------------------------------------------------
# 11. PROCESS-BOOTSTRAP PARAMETER UNCERTAINTY
# -----------------------------------------------------------------------------

PARAMETER_BOOTSTRAP_FILE <- file.path(
  CHECKPOINT_DIR,
  paste0(
    "SS_Parameter_Bootstrap_",
    ANALYSIS_MODE,
    ".csv"
  )
)

run_parameter_bootstrap <- function() {

  if (
    file.exists(
      PARAMETER_BOOTSTRAP_FILE
    )
  ) {
    results <- read.csv(
      PARAMETER_BOOTSTRAP_FILE,
      stringsAsFactors = FALSE
    )
  } else {
    results <- data.frame(
      replication_id = integer(0),
      beta_spell_hat = numeric(0),
      relative_intensity_s2 = numeric(0),
      relative_intensity_s5 = numeric(0),
      relative_intensity_s10 = numeric(0),
      converged = logical(0),
      stringsAsFactors = FALSE
    )
  }

  completed <- results$replication_id

  process_ids <- unique(
    as.character(
      TRANSITIONS$PP
    )
  )

  set.seed(
    BASE_SEED + 20000L
  )

  bootstrap_seeds <- sample.int(
    .Machine$integer.max,
    PARAMETER_BOOTSTRAP_B
  )

  for (
    replication_id in
    seq_len(
      PARAMETER_BOOTSTRAP_B
    )
  ) {
    if (
      replication_id %in%
        completed
    ) {
      next
    }

    cat(
      "Parameter bootstrap ",
      replication_id,
      " of ",
      PARAMETER_BOOTSTRAP_B,
      "\n",
      sep = ""
    )

    set.seed(
      bootstrap_seeds[
        replication_id
      ]
    )

    sampled_processes <- sample(
      process_ids,
      size = length(process_ids),
      replace = TRUE
    )

    bootstrap_parts <- vector(
      "list",
      length(sampled_processes)
    )

    for (
      sample_index in
      seq_along(
        sampled_processes
      )
    ) {
      process_data <- TRANSITIONS %>%
        filter(
          PP ==
            sampled_processes[
              sample_index
            ]
        )

      process_data$PP <- paste0(
        "bootstrap_",
        replication_id,
        "_",
        sample_index
      )

      bootstrap_parts[[sample_index]] <-
        process_data
    }

    bootstrap_data <- bind_rows(
      bootstrap_parts
    )

    null_fit <- tryCatch(
      fit_constrained_model(
        bootstrap_data,
        STAGE_LEVELS,
        model = "homogeneous",
        number_of_starts = 4L
      ),
      error = function(e) NULL
    )

    bootstrap_fit <- tryCatch(
      fit_constrained_model(
        bootstrap_data,
        STAGE_LEVELS,
        model = "stage_spell",
        embedded_fits = list(
          homogeneous =
            null_fit
        ),
        number_of_starts =
          MULTISTARTS
      ),
      error = function(e) {
        append_error(
          "parameter_bootstrap",
          replication_id,
          conditionMessage(e)
        )
        NULL
      }
    )

    converged <-
      !is.null(bootstrap_fit) &&
      is.finite(
        bootstrap_fit$loglik
      )

    beta_hat <- if (
      converged
    ) {
      bootstrap_fit$beta_spell
    } else {
      NA_real_
    }

    results <- bind_rows(
      results,
      data.frame(
        replication_id =
          replication_id,
        beta_spell_hat =
          beta_hat,
        relative_intensity_s2 = if (
          converged
        ) {
          2^(-beta_hat)
        } else {
          NA_real_
        },
        relative_intensity_s5 = if (
          converged
        ) {
          5^(-beta_hat)
        } else {
          NA_real_
        },
        relative_intensity_s10 = if (
          converged
        ) {
          10^(-beta_hat)
        } else {
          NA_real_
        },
        converged =
          converged,
        stringsAsFactors = FALSE
      )
    )

    write.csv(
      results,
      PARAMETER_BOOTSTRAP_FILE,
      row.names = FALSE
    )
  }

  results
}

PARAMETER_BOOTSTRAP_RESULTS <-
  run_parameter_bootstrap()

summarize_bootstrap_parameter <- function(
  values,
  estimate,
  parameter_name
) {
  finite_values <- values[
    is.finite(values)
  ]

  data.frame(
    parameter =
      parameter_name,
    estimate =
      estimate,
    bootstrap_mean =
      safe_mean(finite_values),
    bootstrap_standard_error =
      safe_sd(finite_values),
    ci_lower = if (
      length(finite_values) > 0L
    ) {
      unname(
        quantile(
          finite_values,
          0.025
        )
      )
    } else {
      NA_real_
    },
    ci_upper = if (
      length(finite_values) > 0L
    ) {
      unname(
        quantile(
          finite_values,
          0.975
        )
      )
    } else {
      NA_real_
    },
    successful_replications =
      length(finite_values),
    stringsAsFactors = FALSE
  )
}

PARAMETER_BOOTSTRAP_SUMMARY <- bind_rows(
  summarize_bootstrap_parameter(
    PARAMETER_BOOTSTRAP_RESULTS$
      beta_spell_hat,
    FIT_STAGE_SPELL$beta_spell,
    "beta_spell"
  ),
  summarize_bootstrap_parameter(
    PARAMETER_BOOTSTRAP_RESULTS$
      relative_intensity_s2,
    2^(
      -FIT_STAGE_SPELL$beta_spell
    ),
    "relative_intensity_s2"
  ),
  summarize_bootstrap_parameter(
    PARAMETER_BOOTSTRAP_RESULTS$
      relative_intensity_s5,
    5^(
      -FIT_STAGE_SPELL$beta_spell
    ),
    "relative_intensity_s5"
  ),
  summarize_bootstrap_parameter(
    PARAMETER_BOOTSTRAP_RESULTS$
      relative_intensity_s10,
    10^(
      -FIT_STAGE_SPELL$beta_spell
    ),
    "relative_intensity_s10"
  )
)

save_csv(
  PARAMETER_BOOTSTRAP_RESULTS,
  "08_parameter_bootstrap_replications.csv"
)

save_csv(
  PARAMETER_BOOTSTRAP_SUMMARY,
  "09_parameter_bootstrap_summary.csv"
)


# -----------------------------------------------------------------------------
# 12. GROUPED PROCESS-LEVEL CROSS-VALIDATION
# -----------------------------------------------------------------------------

make_grouped_folds <- function(
  data,
  fold_count,
  seed
) {
  set.seed(seed)

  process_sizes <- data %>%
    count(
      PP,
      name = "transitions"
    ) %>%
    mutate(
      random_tie =
        runif(n())
    ) %>%
    arrange(
      desc(transitions),
      random_tie
    )

  fold_load <- rep(
    0,
    fold_count
  )

  folds <- vector(
    "list",
    fold_count
  )

  for (
    process_index in
    seq_len(
      nrow(process_sizes)
    )
  ) {
    candidates <- which(
      fold_load ==
        min(fold_load)
    )

    selected_fold <- sample(
      candidates,
      1L
    )

    folds[[selected_fold]] <- c(
      folds[[selected_fold]],
      process_sizes$PP[
        process_index
      ]
    )

    fold_load[selected_fold] <-
      fold_load[selected_fold] +
      process_sizes$transitions[
        process_index
      ]
  }

  folds
}

fit_cv_models <- function(
  training_data,
  stage_levels
) {
  homogeneous <- fit_constrained_model(
    training_data,
    stage_levels,
    model = "homogeneous"
  )

  stage_spell <- fit_constrained_model(
    training_data,
    stage_levels,
    model = "stage_spell",
    embedded_fits = list(
      homogeneous =
        homogeneous
    )
  )

  global_alpha_beta <- fit_constrained_model(
    training_data,
    stage_levels,
    model =
      "global_alpha_beta",
    embedded_fits = list(
      homogeneous =
        homogeneous
    )
  )

  dual_clock <- fit_constrained_model(
    training_data,
    stage_levels,
    model = "dual_clock",
    embedded_fits = list(
      homogeneous =
        homogeneous,
      stage_spell =
        stage_spell
    )
  )

  list(
    homogeneous =
      homogeneous,
    stage_spell =
      stage_spell,
    global_alpha_beta =
      global_alpha_beta,
    dual_clock =
      dual_clock,
    old_kernel =
      fit_old_kernel_model(
        training_data,
        stage_levels,
        homogeneous_fit =
          homogeneous
      ),
    spline_spell =
      fit_glm_comparator(
        training_data,
        comparator =
          "spline_spell"
      ),
    quadratic_r =
      fit_glm_comparator(
        training_data,
        comparator =
          "quadratic_r"
      ),
    destination =
      fit_destination_matrix(
        training_data,
        stage_levels
      )
  )
}

predict_cv_models <- function(
  fits,
  test_data
) {
  list(
    `Stage-Spell SPAOM` =
      predict_constrained_model(
        fits$stage_spell,
        test_data
      ),
    `Origin-specific homogeneous` =
      predict_constrained_model(
        fits$homogeneous,
        test_data
      ),
    `Dual-clock model` =
      predict_constrained_model(
        fits$dual_clock,
        test_data
      ),
    `Global alpha-beta model` =
      predict_constrained_model(
        fits$global_alpha_beta,
        test_data
      ),
    `Earlier kernel SPAOM` =
      predict_old_kernel_model(
        fits$old_kernel,
        test_data
      ),
    `Spline spell-age model` =
      predict_glm_comparator(
        fits$spline_spell,
        test_data
      ),
    `Quadratic global-r model` =
      predict_glm_comparator(
        fits$quadratic_r,
        test_data
      )
  )
}

CV_DIRECTORY <- file.path(
  CHECKPOINT_DIR,
  paste0(
    "SS_CV_",
    ANALYSIS_MODE
  )
)

dir.create(
  CV_DIRECTORY,
  recursive = TRUE,
  showWarnings = FALSE
)

run_cross_validation <- function() {

  metric_files <- character(0)
  contribution_files <- character(0)

  for (
    repeat_id in
    seq_len(
      CV_REPEATS
    )
  ) {
    folds <- make_grouped_folds(
      TRANSITIONS,
      CV_FOLDS,
      BASE_SEED +
        30000L +
        repeat_id
    )

    for (
      fold_id in
      seq_len(
        CV_FOLDS
      )
    ) {
      prefix <- paste0(
        "repeat_",
        sprintf(
          "%02d",
          repeat_id
        ),
        "_fold_",
        sprintf(
          "%02d",
          fold_id
        )
      )

      metric_file <- file.path(
        CV_DIRECTORY,
        paste0(
          prefix,
          "_metrics.csv"
        )
      )

      contribution_file <- file.path(
        CV_DIRECTORY,
        paste0(
          prefix,
          "_contributions.csv"
        )
      )

      metric_files <- c(
        metric_files,
        metric_file
      )

      contribution_files <- c(
        contribution_files,
        contribution_file
      )

      if (
        file.exists(metric_file) &&
        file.exists(
          contribution_file
        )
      ) {
        next
      }

      cat(
        "Cross-validation repeat ",
        repeat_id,
        "/",
        CV_REPEATS,
        ", fold ",
        fold_id,
        "/",
        CV_FOLDS,
        "\n",
        sep = ""
      )

      test_processes <- folds[[fold_id]]

      training_data <- TRANSITIONS %>%
        filter(
          !(PP %in%
              test_processes)
        )

      test_data <- TRANSITIONS %>%
        filter(
          PP %in%
            test_processes
        )

      fits <- tryCatch(
        fit_cv_models(
          training_data,
          STAGE_LEVELS
        ),
        error = function(e) {
          append_error(
            "cross_validation_fit",
            prefix,
            conditionMessage(e)
          )
          NULL
        }
      )

      if (is.null(fits)) {
        next
      }

      predictions <- tryCatch(
        predict_cv_models(
          fits,
          test_data
        ),
        error = function(e) {
          append_error(
            "cross_validation_prediction",
            prefix,
            conditionMessage(e)
          )
          NULL
        }
      )

      if (is.null(predictions)) {
        next
      }

      metric_rows <- list()
      contribution_rows <- list()
      result_index <- 1L

      for (
        model_name in
        names(predictions)
      ) {
        q <- predictions[[model_name]]

        transition_probability_matrix <-
          combine_transition_probabilities(
            test_data,
            q,
            fits$destination,
            STAGE_LEVELS
          )

        evaluation <- evaluate_predictions(
          test_data,
          q,
          transition_probability_matrix,
          STAGE_LEVELS
        )

        metric_rows[[result_index]] <- data.frame(
          repeat_id =
            repeat_id,
          fold_id =
            fold_id,
          model =
            model_name,
          departure_logloss =
            unname(
              evaluation$metrics[
                "departure_logloss"
              ]
            ),
          departure_brier =
            unname(
              evaluation$metrics[
                "departure_brier"
              ]
            ),
          calibration_intercept =
            unname(
              evaluation$metrics[
                "calibration_intercept"
              ]
            ),
          calibration_slope =
            unname(
              evaluation$metrics[
                "calibration_slope"
              ]
            ),
          full_logloss =
            unname(
              evaluation$metrics[
                "full_logloss"
              ]
            ),
          full_brier =
            unname(
              evaluation$metrics[
                "full_brier"
              ]
            ),
          stringsAsFactors = FALSE
        )

        contribution <- evaluation$
          contributions

        contribution$repeat_id <-
          repeat_id

        contribution$fold_id <-
          fold_id

        contribution$model <-
          model_name

        contribution_rows[[result_index]] <-
          contribution

        result_index <- result_index + 1L
      }

      write.csv(
        bind_rows(metric_rows),
        metric_file,
        row.names = FALSE
      )

      write.csv(
        bind_rows(
          contribution_rows
        ),
        contribution_file,
        row.names = FALSE
      )
    }
  }

  existing_metric_files <- metric_files[
    file.exists(metric_files)
  ]

  existing_contribution_files <-
    contribution_files[
      file.exists(
        contribution_files
      )
    ]

  if (
    length(existing_metric_files) == 0L ||
    length(existing_contribution_files) == 0L
  ) {
    stop(
      "No cross-validation files were produced."
    )
  }

  list(
    fold_metrics = bind_rows(
      lapply(
        existing_metric_files,
        read.csv,
        stringsAsFactors = FALSE
      )
    ),
    contributions = bind_rows(
      lapply(
        existing_contribution_files,
        read.csv,
        stringsAsFactors = FALSE
      )
    )
  )
}

CV_RESULT <- run_cross_validation()

CV_REPEAT_METRICS <- CV_RESULT$
  contributions %>%
  group_by(
    repeat_id,
    model
  ) %>%
  summarise(
    departure_logloss =
      mean(departure_logloss),
    departure_brier =
      mean(departure_brier),
    full_logloss =
      mean(full_logloss),
    full_brier =
      mean(full_brier),
    calibration_intercept =
      unname(
        calibration_intercept_slope(
          D,
          q
        )[
          "intercept"
        ]
      ),
    calibration_slope =
      unname(
        calibration_intercept_slope(
          D,
          q
        )[
          "slope"
        ]
      ),
    .groups = "drop"
  )

CV_MODEL_SUMMARY <- CV_REPEAT_METRICS %>%
  group_by(model) %>%
  summarise(
    repetitions = n(),
    departure_logloss_mean =
      mean(departure_logloss),
    departure_logloss_sd =
      safe_sd(departure_logloss),
    departure_brier_mean =
      mean(departure_brier),
    departure_brier_sd =
      safe_sd(departure_brier),
    calibration_intercept_mean =
      mean(
        calibration_intercept,
        na.rm = TRUE
      ),
    calibration_slope_mean =
      mean(
        calibration_slope,
        na.rm = TRUE
      ),
    full_logloss_mean =
      mean(full_logloss),
    full_brier_mean =
      mean(full_brier),
    .groups = "drop"
  ) %>%
  arrange(
    departure_logloss_mean
  )

save_csv(
  CV_RESULT$fold_metrics,
  "10_cv_fold_metrics.csv"
)

save_csv(
  CV_REPEAT_METRICS,
  "11_cv_repeat_metrics.csv"
)

save_csv(
  CV_MODEL_SUMMARY,
  "12_cv_model_summary.csv"
)


# -----------------------------------------------------------------------------
# 13. PAIRED PROCESS-BOOTSTRAP PREDICTIVE ADVANTAGES
# -----------------------------------------------------------------------------

PRIMARY_MODEL_NAME <- "Stage-Spell SPAOM"

comparison_models <- setdiff(
  unique(
    CV_RESULT$contributions$model
  ),
  PRIMARY_MODEL_NAME
)

metrics_to_compare <- c(
  "departure_logloss",
  "departure_brier",
  "full_logloss",
  "full_brier"
)

calculate_process_bootstrap_interval <- function(
  contribution_data,
  comparison_model,
  metric_name
) {
  primary <- contribution_data %>%
    filter(
      model ==
        PRIMARY_MODEL_NAME
    ) %>%
    select(
      repeat_id,
      transition_id,
      PP,
      primary_value =
        all_of(metric_name)
    )

  comparison <- contribution_data %>%
    filter(
      model ==
        comparison_model
    ) %>%
    select(
      repeat_id,
      transition_id,
      comparison_value =
        all_of(metric_name)
    )

  paired <- inner_join(
    primary,
    comparison,
    by = c(
      "repeat_id",
      "transition_id"
    )
  ) %>%
    mutate(
      advantage =
        comparison_value -
        primary_value
    ) %>%
    group_by(PP) %>%
    summarise(
      advantage_sum =
        sum(advantage),
      observation_count =
        n(),
      .groups = "drop"
    )

  observed_advantage <- sum(
    paired$advantage_sum
  ) /
    sum(
      paired$observation_count
    )

  set.seed(
    BASE_SEED +
      40000L +
      sum(
        utf8ToInt(
          paste0(
            comparison_model,
            metric_name
          )
        )
      )
  )

  bootstrap_values <- replicate(
    PAIRED_BOOTSTRAP_B,
    {
      sampled_indices <- sample(
        seq_len(
          nrow(paired)
        ),
        size = nrow(paired),
        replace = TRUE
      )

      sum(
        paired$advantage_sum[
          sampled_indices
        ]
      ) /
        sum(
          paired$observation_count[
            sampled_indices
          ]
        )
    }
  )

  interval <- quantile(
    bootstrap_values,
    probs = c(
      0.025,
      0.975
    ),
    na.rm = TRUE
  )

  c(
    advantage =
      observed_advantage,
    lower =
      unname(interval[1L]),
    upper =
      unname(interval[2L])
  )
}

paired_rows <- list()
paired_index <- 1L

for (
  comparison_model in
  comparison_models
) {
  primary_repeat <- CV_REPEAT_METRICS %>%
    filter(
      model ==
        PRIMARY_MODEL_NAME
    )

  comparison_repeat <- CV_REPEAT_METRICS %>%
    filter(
      model ==
        comparison_model
    )

  paired_repeat <- inner_join(
    primary_repeat,
    comparison_repeat,
    by = "repeat_id",
    suffix = c(
      "_primary",
      "_comparison"
    )
  )

  for (
    metric_name in
    metrics_to_compare
  ) {
    differences <- paired_repeat[[
      paste0(
        metric_name,
        "_comparison"
      )
    ]] -
      paired_repeat[[
        paste0(
          metric_name,
          "_primary"
        )
      ]]

    differences <- differences[
      is.finite(differences)
    ]

    process_interval <-
      calculate_process_bootstrap_interval(
        CV_RESULT$contributions,
        comparison_model,
        metric_name
      )

    paired_rows[[paired_index]] <- data.frame(
      comparison_model =
        comparison_model,
      metric =
        metric_name,
      mean_repeat_advantage =
        safe_mean(differences),
      repeat_standard_error = if (
        length(differences) > 1L
      ) {
        safe_sd(differences) /
          sqrt(length(differences))
      } else {
        NA_real_
      },
      repeat_win_proportion =
        mean(differences > 0),
      process_bootstrap_advantage =
        unname(
          process_interval[
            "advantage"
          ]
        ),
      process_bootstrap_lower =
        unname(
          process_interval[
            "lower"
          ]
        ),
      process_bootstrap_upper =
        unname(
          process_interval[
            "upper"
          ]
        ),
      stringsAsFactors = FALSE
    )

    paired_index <- paired_index + 1L
  }
}

PAIRED_ADVANTAGES <- bind_rows(
  paired_rows
)

save_csv(
  PAIRED_ADVANTAGES,
  "13_paired_predictive_advantages.csv"
)


# -----------------------------------------------------------------------------
# 14. CALIBRATION TABLE
# -----------------------------------------------------------------------------

make_calibration_bins <- function(
  actual,
  predicted,
  bin_count = 10L
) {
  breaks <- unique(
    quantile(
      predicted,
      probs = seq(
        0,
        1,
        length.out =
          bin_count + 1L
      ),
      na.rm = TRUE
    )
  )

  if (length(breaks) < 3L) {
    breaks <- seq(
      0,
      1,
      length.out =
        bin_count + 1L
    )
  }

  bin <- cut(
    predicted,
    breaks = breaks,
    include.lowest = TRUE,
    labels = FALSE
  )

  data.frame(
    actual = actual,
    predicted = predicted,
    bin = bin
  ) %>%
    filter(
      !is.na(bin)
    ) %>%
    group_by(bin) %>%
    summarise(
      mean_predicted =
        mean(predicted),
      observed_departure =
        mean(actual),
      transitions = n(),
      .groups = "drop"
    )
}

CALIBRATION_TABLE <- CV_RESULT$
  contributions %>%
  group_by(
    model,
    transition_id
  ) %>%
  summarise(
    D = first(D),
    q = mean(q),
    .groups = "drop"
  ) %>%
  group_by(model) %>%
  group_modify(
    ~ make_calibration_bins(
      .x$D,
      .x$q
    )
  ) %>%
  ungroup()

save_csv(
  CALIBRATION_TABLE,
  "14_departure_calibration_table.csv"
)


# -----------------------------------------------------------------------------
# 15. AGREEMENT-SEQUENCE RECONSTRUCTION FOR ROBUSTNESS
# -----------------------------------------------------------------------------

can_reconstruct_agreements <- all(
  c(
    "origin_agreement_id",
    "destination_agreement_id",
    "origin_date",
    "destination_date"
  ) %in%
    names(TRANSITIONS)
)

reconstruct_agreements <- function(
  transition_data
) {
  first_rows <- transition_data %>%
    arrange(
      PP,
      transition_index
    ) %>%
    group_by(PP) %>%
    slice_head(n = 1L) %>%
    ungroup()

  initial_agreements <- first_rows %>%
    transmute(
      PP = as.character(PP),
      agreement_id =
        as.character(
          origin_agreement_id
        ),
      agreement_date =
        as_date_safely(
          origin_date
        ),
      stage =
        as.character(
          origin_stage
        ),
      original_order = 0
    )

  destination_agreements <- transition_data %>%
    transmute(
      PP = as.character(PP),
      agreement_id =
        as.character(
          destination_agreement_id
        ),
      agreement_date =
        as_date_safely(
          destination_date
        ),
      stage =
        as.character(
          destination_stage
        ),
      original_order =
        as.numeric(
          transition_index
        )
    )

  bind_rows(
    initial_agreements,
    destination_agreements
  ) %>%
    filter(
      !is.na(PP),
      nzchar(PP),
      !is.na(agreement_date),
      !is.na(stage),
      nzchar(stage)
    ) %>%
    arrange(
      PP,
      original_order
    ) %>%
    distinct(
      PP,
      agreement_id,
      .keep_all = TRUE
    )
}

build_transitions_from_agreements <- function(
  agreements,
  order_mode = c(
    "current",
    "reversed_within_date"
  ),
  one_per_date = FALSE,
  exclude_ambiguous_dates = FALSE,
  four_stage = FALSE
) {
  order_mode <- match.arg(
    order_mode
  )

  data <- agreements

  if (four_stage) {
    data <- data %>%
      mutate(
        stage = case_when(
          stage ==
            "Pre-negotiation" ~
            "Early",
          stage ==
            "Ceasefire" ~
            "Ceasefire",
          stage %in% c(
            "Substantive Partial",
            "Substantive Comprehensive"
          ) ~
            "Substantive",
          stage %in% c(
            "Implementation",
            "Renegotiation/Other"
          ) ~
            "Implementation/Renegotiation",
          TRUE ~
            stage
        )
      )
  }

  if (exclude_ambiguous_dates) {
    data <- data %>%
      group_by(
        PP,
        agreement_date
      ) %>%
      filter(
        n() == 1L
      ) %>%
      ungroup()
  }

  if (one_per_date) {
    data <- data %>%
      arrange(
        PP,
        agreement_date,
        original_order
      ) %>%
      group_by(
        PP,
        agreement_date
      ) %>%
      slice_tail(n = 1L) %>%
      ungroup()
  }

  if (order_mode ==
      "current") {
    data <- data %>%
      arrange(
        PP,
        agreement_date,
        original_order
      )
  } else {
    data <- data %>%
      arrange(
        PP,
        agreement_date,
        desc(original_order)
      )
  }

  transitions <- data %>%
    group_by(PP) %>%
    mutate(
      origin_stage =
        dplyr::lag(stage),
      destination_stage =
        stage,
      origin_date =
        dplyr::lag(
          agreement_date
        ),
      destination_date =
        agreement_date,
      transition_index =
        row_number() - 1L
    ) %>%
    filter(
      transition_index >= 1L
    ) %>%
    ungroup() %>%
    mutate(
      D = as.integer(
        origin_stage !=
          destination_stage
      ),
      same_date_transition =
        as.integer(
          origin_date ==
            destination_date
        ),
      transition_id = paste0(
        PP,
        "::",
        transition_index
      )
    )

  stage_levels <- unique(
    c(
      as.character(
        transitions$origin_stage
      ),
      as.character(
        transitions$destination_stage
      )
    )
  )

  preferred <- if (four_stage) {
    c(
      "Early",
      "Ceasefire",
      "Substantive",
      "Implementation/Renegotiation"
    )
  } else {
    STAGE_LEVELS
  }

  stage_levels <- c(
    preferred[
      preferred %in%
        stage_levels
    ],
    sort(
      setdiff(
        stage_levels,
        preferred
      )
    )
  )

  list(
    data = add_stage_spell_age(
      transitions,
      stage_levels
    ),
    stage_levels =
      stage_levels
  )
}


# -----------------------------------------------------------------------------
# 16. ROBUSTNESS ANALYSES
# -----------------------------------------------------------------------------

PRIMARY_AGREEMENTS <- if (
  can_reconstruct_agreements
) {
  reconstruct_agreements(
    TRANSITIONS
  )
} else {
  NULL
}

make_robustness_scenarios <- function() {
  scenarios <- list(
    `Current ordering` = list(
      data = TRANSITIONS,
      stage_levels =
        STAGE_LEVELS,
      weights =
        rep(
          1,
          nrow(TRANSITIONS)
        )
    )
  )

  process_lengths <- TRANSITIONS %>%
    count(
      PP,
      name = "R"
    )

  for (minimum_length in c(
    2L,
    5L,
    10L
  )) {
    eligible <- process_lengths %>%
      filter(
        R >= minimum_length
      ) %>%
      pull(PP)

    scenario_data <- TRANSITIONS %>%
      filter(
        PP %in%
          eligible
      )

    scenarios[[
      paste0(
        "Processes with R at least ",
        minimum_length
      )
    ]] <- list(
      data = scenario_data,
      stage_levels =
        STAGE_LEVELS,
      weights =
        rep(
          1,
          nrow(scenario_data)
        )
    )
  }

  truncated <- TRANSITIONS %>%
    filter(
      transition_index <= 30
    )

  scenarios[["Truncated at r=30"]] <- list(
    data = truncated,
    stage_levels =
      STAGE_LEVELS,
    weights =
      rep(
        1,
        nrow(truncated)
      )
  )

  scenarios[["Equal process weighting"]] <- list(
    data = TRANSITIONS,
    stage_levels =
      STAGE_LEVELS,
    weights =
      TRANSITIONS$
        process_weight_equal
  )

  if (
    "same_date_transition" %in%
      names(TRANSITIONS)
  ) {
    no_same_date <- TRANSITIONS %>%
      filter(
        as.numeric(
          same_date_transition
        ) != 1
      )

    no_same_date <- add_stage_spell_age(
      no_same_date,
      STAGE_LEVELS
    )

    scenarios[["Exclude same-date transitions"]] <- list(
      data = no_same_date,
      stage_levels =
        STAGE_LEVELS,
      weights =
        rep(
          1,
          nrow(no_same_date)
        )
    )
  }

  if (!is.null(PRIMARY_AGREEMENTS)) {

    reversed <- build_transitions_from_agreements(
      PRIMARY_AGREEMENTS,
      order_mode =
        "reversed_within_date"
    )

    scenarios[["Reversed within-date ordering"]] <- list(
      data = reversed$data,
      stage_levels =
        reversed$stage_levels,
      weights =
        rep(
          1,
          nrow(reversed$data)
        )
    )

    one_per_date <- build_transitions_from_agreements(
      PRIMARY_AGREEMENTS,
      one_per_date = TRUE
    )

    scenarios[["One agreement per process-date"]] <- list(
      data = one_per_date$data,
      stage_levels =
        one_per_date$stage_levels,
      weights =
        rep(
          1,
          nrow(one_per_date$data)
        )
    )

    unambiguous <- build_transitions_from_agreements(
      PRIMARY_AGREEMENTS,
      exclude_ambiguous_dates =
        TRUE
    )

    scenarios[["Exclude ambiguous process-dates"]] <- list(
      data = unambiguous$data,
      stage_levels =
        unambiguous$stage_levels,
      weights =
        rep(
          1,
          nrow(unambiguous$data)
        )
    )

    four_stage <- build_transitions_from_agreements(
      PRIMARY_AGREEMENTS,
      four_stage = TRUE
    )

    scenarios[["Four-stage aggregation"]] <- list(
      data = four_stage$data,
      stage_levels =
        four_stage$stage_levels,
      weights =
        rep(
          1,
          nrow(four_stage$data)
        )
    )
  }

  scenarios
}

ROBUSTNESS_SCENARIOS <-
  make_robustness_scenarios()

robustness_rows <- list()
robustness_index <- 1L

for (
  scenario_name in
  names(
    ROBUSTNESS_SCENARIOS
  )
) {
  cat(
    "Robustness fit: ",
    scenario_name,
    "\n",
    sep = ""
  )

  scenario <- ROBUSTNESS_SCENARIOS[[scenario_name]]

  scenario_data <- scenario$data
  scenario_stages <-
    scenario$stage_levels
  scenario_weights <-
    scenario$weights

  if (
    nrow(scenario_data) < 20L ||
    dplyr::n_distinct(
      scenario_data$PP
    ) < 5L
  ) {
    next
  }

  homogeneous <- tryCatch(
    fit_constrained_model(
      scenario_data,
      scenario_stages,
      model = "homogeneous",
      weights =
        scenario_weights,
      number_of_starts = 5L
    ),
    error = function(e) NULL
  )

  spell <- tryCatch(
    fit_constrained_model(
      scenario_data,
      scenario_stages,
      model = "stage_spell",
      weights =
        scenario_weights,
      embedded_fits = list(
        homogeneous =
          homogeneous
      ),
      number_of_starts =
        max(
          6L,
          MULTISTARTS %/% 2L
        )
    ),
    error = function(e) {
      append_error(
        "robustness_fit",
        scenario_name,
        conditionMessage(e)
      )
      NULL
    }
  )

  if (
    is.null(homogeneous) ||
    is.null(spell)
  ) {
    next
  }

  robustness_rows[[robustness_index]] <-
    data.frame(
      scenario =
        scenario_name,
      transitions =
        nrow(scenario_data),
      peace_processes =
        dplyr::n_distinct(
          scenario_data$PP
        ),
      beta_spell_hat =
        spell$beta_spell,
      relative_intensity_s5 =
        5^(
          -spell$beta_spell
        ),
      LR_vs_homogeneous =
        max(
          0,
          2 * (
            spell$loglik -
              homogeneous$loglik
          )
        ),
      delta_AIC_homogeneous_minus_spell =
        (
          -2 *
            homogeneous$loglik +
            2 *
            homogeneous$parameter_count
        ) -
        (
          -2 *
            spell$loglik +
            2 *
            spell$parameter_count
        ),
      stringsAsFactors = FALSE
    )

  robustness_index <- robustness_index + 1L
}

ROBUSTNESS_RESULTS <- bind_rows(
  robustness_rows
)

save_csv(
  ROBUSTNESS_RESULTS,
  "15_robustness_parameter_results.csv"
)


# -----------------------------------------------------------------------------
# 17. COMPACT ROBUSTNESS CROSS-VALIDATION
# -----------------------------------------------------------------------------

run_two_model_grouped_cv <- function(
  scenario_data,
  stage_levels,
  repeats,
  folds,
  scenario_name
) {
  rows <- list()
  row_index <- 1L

  for (
    repeat_id in seq_len(repeats)
  ) {
    fold_list <- make_grouped_folds(
      scenario_data,
      folds,
      BASE_SEED +
        50000L +
        repeat_id +
        sum(
          utf8ToInt(
            scenario_name
          )
        )
    )

    for (
      fold_id in seq_len(folds)
    ) {
      test_processes <- fold_list[[fold_id]]

      training_data <- scenario_data %>%
        filter(
          !(PP %in%
              test_processes)
        )

      test_data <- scenario_data %>%
        filter(
          PP %in%
            test_processes
        )

      homogeneous <- tryCatch(
        fit_constrained_model(
          training_data,
          stage_levels,
          model = "homogeneous",
          number_of_starts = 4L
        ),
        error = function(e) NULL
      )

      spell <- tryCatch(
        fit_constrained_model(
          training_data,
          stage_levels,
          model = "stage_spell",
          embedded_fits = list(
            homogeneous =
              homogeneous
          ),
          number_of_starts = 6L
        ),
        error = function(e) NULL
      )

      if (
        is.null(homogeneous) ||
        is.null(spell)
      ) {
        next
      }

      q_homogeneous <-
        predict_constrained_model(
          homogeneous,
          test_data
        )

      q_spell <-
        predict_constrained_model(
          spell,
          test_data
        )

      homogeneous_logloss <- mean(
        -(
          test_data$D *
            log(q_homogeneous) +
          (1 - test_data$D) *
            log(
              1 -
                q_homogeneous
            )
        )
      )

      spell_logloss <- mean(
        -(
          test_data$D *
            log(q_spell) +
          (1 - test_data$D) *
            log(
              1 -
                q_spell
            )
        )
      )

      homogeneous_brier <- mean(
        (
          test_data$D -
            q_homogeneous
        )^2
      )

      spell_brier <- mean(
        (
          test_data$D -
            q_spell
        )^2
      )

      rows[[row_index]] <- data.frame(
        scenario =
          scenario_name,
        repeat_id =
          repeat_id,
        fold_id =
          fold_id,
        logloss_advantage =
          homogeneous_logloss -
          spell_logloss,
        brier_advantage =
          homogeneous_brier -
          spell_brier,
        stringsAsFactors = FALSE
      )

      row_index <- row_index + 1L
    }
  }

  bind_rows(rows)
}

robust_cv_rows <- list()
robust_cv_index <- 1L

for (
  scenario_name in
  names(
    ROBUSTNESS_SCENARIOS
  )
) {
  scenario <- ROBUSTNESS_SCENARIOS[[scenario_name]]

  scenario_data <- scenario$data

  if (
    nrow(scenario_data) < 50L ||
    dplyr::n_distinct(
      scenario_data$PP
    ) < 10L
  ) {
    next
  }

  cat(
    "Robustness CV: ",
    scenario_name,
    "\n",
    sep = ""
  )

  robust_cv_rows[[robust_cv_index]] <-
    run_two_model_grouped_cv(
      scenario_data =
        scenario_data,
      stage_levels =
        scenario$stage_levels,
      repeats =
        ROBUST_CV_REPEATS,
      folds =
        ROBUST_CV_FOLDS,
      scenario_name =
        scenario_name
    )

  robust_cv_index <- robust_cv_index + 1L
}

ROBUSTNESS_CV_FOLDS <- bind_rows(
  robust_cv_rows
)

ROBUSTNESS_CV_SUMMARY <- ROBUSTNESS_CV_FOLDS %>%
  group_by(scenario) %>%
  summarise(
    folds = n(),
    mean_logloss_advantage =
      mean(logloss_advantage),
    logloss_win_proportion =
      mean(
        logloss_advantage > 0
      ),
    mean_brier_advantage =
      mean(brier_advantage),
    brier_win_proportion =
      mean(
        brier_advantage > 0
      ),
    .groups = "drop"
  )

save_csv(
  ROBUSTNESS_CV_FOLDS,
  "16_robustness_cv_fold_results.csv"
)

save_csv(
  ROBUSTNESS_CV_SUMMARY,
  "17_robustness_cv_summary.csv"
)


# -----------------------------------------------------------------------------
# 18. LEAVE-ONE-PROCESS-OUT INFLUENCE
# -----------------------------------------------------------------------------

process_ids <- unique(
  as.character(
    TRANSITIONS$PP
  )
)

influence_rows <- list()

for (
  process_index in
  seq_along(process_ids)
) {
  process_id <- process_ids[
    process_index
  ]

  cat(
    "Influence fit ",
    process_index,
    " of ",
    length(process_ids),
    "\n",
    sep = ""
  )

  reduced_data <- TRANSITIONS %>%
    filter(
      PP !=
        process_id
    )

  homogeneous <- tryCatch(
    fit_constrained_model(
      reduced_data,
      STAGE_LEVELS,
      model = "homogeneous",
      number_of_starts = 3L
    ),
    error = function(e) NULL
  )

  spell <- tryCatch(
    fit_constrained_model(
      reduced_data,
      STAGE_LEVELS,
      model = "stage_spell",
      embedded_fits = list(
        homogeneous =
          homogeneous
      ),
      number_of_starts = 5L
    ),
    error = function(e) NULL
  )

  if (
    is.null(spell) ||
    !is.finite(
      spell$beta_spell
    )
  ) {
    influence_rows[[process_index]] <- data.frame(
      PP = process_id,
      beta_leave_one_out = NA_real_,
      beta_change = NA_real_,
      process_transitions = sum(
        TRANSITIONS$PP ==
          process_id
      ),
      stringsAsFactors = FALSE
    )
  } else {
    influence_rows[[process_index]] <- data.frame(
      PP = process_id,
      beta_leave_one_out =
        spell$beta_spell,
      beta_change =
        spell$beta_spell -
        FIT_STAGE_SPELL$beta_spell,
      process_transitions = sum(
        TRANSITIONS$PP ==
          process_id
      ),
      stringsAsFactors = FALSE
    )
  }
}

INFLUENCE_RESULTS <- bind_rows(
  influence_rows
) %>%
  arrange(
    desc(
      abs(beta_change)
    )
  )

save_csv(
  INFLUENCE_RESULTS,
  "18_leave_one_process_out_influence.csv"
)


# -----------------------------------------------------------------------------
# 19. DESCRIPTIVE SPELL-AGE TABLE
# -----------------------------------------------------------------------------

SPELL_AGE_SUMMARY <- TRANSITIONS %>%
  mutate(
    spell_age_group = ifelse(
      spell_age <= 10,
      as.character(spell_age),
      "11+"
    ),
    spell_age_group = factor(
      spell_age_group,
      levels = c(
        as.character(
          1:10
        ),
        "11+"
      )
    )
  ) %>%
  group_by(
    spell_age_group
  ) %>%
  summarise(
    transitions = n(),
    departures = sum(D),
    departure_probability =
      mean(D),
    standard_error = sqrt(
      departure_probability *
        (
          1 -
            departure_probability
        ) /
        transitions
    ),
    ci_lower = pmax(
      0,
      departure_probability -
        qnorm(0.975) *
        standard_error
    ),
    ci_upper = pmin(
      1,
      departure_probability +
        qnorm(0.975) *
        standard_error
    ),
    .groups = "drop"
  )

save_csv(
  SPELL_AGE_SUMMARY,
  "19_empirical_spell_age_departure_rates.csv"
)


# -----------------------------------------------------------------------------
# 20. GRAPHS
# -----------------------------------------------------------------------------

plot_spell_empirical <- ggplot(
  SPELL_AGE_SUMMARY,
  aes(
    x = spell_age_group,
    y = departure_probability
  )
) +
  geom_point(
    size = 3
  ) +
  geom_errorbar(
    aes(
      ymin = ci_lower,
      ymax = ci_upper
    ),
    width = 0.15
  ) +
  labs(
    title =
      "Empirical departure probability by stage-spell age",
    subtitle =
      "Spell age counts consecutive documented transitions in the current stage",
    x =
      "Stage-spell age",
    y =
      "Departure probability"
  ) +
  theme_bw()

save_plot(
  plot_spell_empirical,
  "01_empirical_spell_age_departure_rates.png"
)

curve_data <- expand.grid(
  spell_age = seq(
    1,
    min(
      20,
      max(
        TRANSITIONS$spell_age
      )
    )
  ),
  origin_stage =
    STAGE_LEVELS,
  stringsAsFactors = FALSE
)

curve_data$origin_stage <- factor(
  curve_data$origin_stage,
  levels = STAGE_LEVELS
)

curve_data$transition_index <- 1
curve_data$log_spell_age <- log(
  curve_data$spell_age
)
curve_data$log_transition_index <- 0

curve_data$fitted_probability <-
  predict_constrained_model(
    FIT_STAGE_SPELL,
    curve_data
  )

empirical_stage_spell <- TRANSITIONS %>%
  filter(
    spell_age <= 20
  ) %>%
  group_by(
    origin_stage,
    spell_age
  ) %>%
  summarise(
    empirical_probability =
      mean(D),
    transitions = n(),
    .groups = "drop"
  )

plot_fitted_curves <- ggplot() +
  geom_line(
    data = curve_data,
    aes(
      x = spell_age,
      y =
        fitted_probability
    ),
    linewidth = 0.9
  ) +
  geom_point(
    data =
      empirical_stage_spell,
    aes(
      x = spell_age,
      y =
        empirical_probability,
      size = transitions
    ),
    alpha = 0.55
  ) +
  facet_wrap(
    ~ origin_stage,
    scales = "free_y"
  ) +
  labs(
    title =
      "Stage-specific fitted departure profiles under SS-SPAOM",
    subtitle = paste0(
      "Common spell-persistence estimate beta = ",
      format_number(
        FIT_STAGE_SPELL$beta_spell
      )
    ),
    x =
      "Stage-spell age",
    y =
      "Departure probability",
    size =
      "Transitions"
  ) +
  theme_bw()

save_plot(
  plot_fitted_curves,
  "02_stage_specific_spell_departure_curves.png",
  width = 12,
  height = 9
)

relative_decay <- data.frame(
  spell_age = 1:20
) %>%
  mutate(
    relative_opportunity_intensity =
      spell_age^(
        -FIT_STAGE_SPELL$beta_spell
      )
  )

plot_relative_decay <- ggplot(
  relative_decay,
  aes(
    x = spell_age,
    y =
      relative_opportunity_intensity
  )
) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_point() +
  labs(
    title =
      "Estimated decline in opportunity intensity across a stage spell",
    subtitle =
      "Values are relative to spell age one",
    x =
      "Stage-spell age",
    y =
      expression(
        mu[a](s) / mu[a](1)
      )
  ) +
  theme_bw()

save_plot(
  plot_relative_decay,
  "03_relative_opportunity_intensity.png"
)

plot_delta_aic <- FULL_MODEL_COMPARISON %>%
  mutate(
    model = reorder(
      model,
      delta_AIC
    )
  ) %>%
  ggplot(
    aes(
      x = model,
      y = delta_AIC
    )
  ) +
  geom_col() +
  coord_flip() +
  labs(
    title =
      "Departure-model comparison by delta AIC",
    x = "Model",
    y = expression(
      Delta * AIC
    )
  ) +
  theme_bw()

save_plot(
  plot_delta_aic,
  "04_delta_AIC_model_comparison.png"
)

plot_cv_logloss <- ggplot(
  CV_REPEAT_METRICS,
  aes(
    x = reorder(
      model,
      departure_logloss,
      FUN = median
    ),
    y =
      departure_logloss
  )
) +
  geom_boxplot() +
  coord_flip() +
  labs(
    title =
      "Repeated process-level cross-validated departure log loss",
    x = "Model",
    y =
      "Departure log loss"
  ) +
  theme_bw()

save_plot(
  plot_cv_logloss,
  "05_cv_departure_logloss.png"
)

plot_cv_brier <- ggplot(
  CV_REPEAT_METRICS,
  aes(
    x = reorder(
      model,
      departure_brier,
      FUN = median
    ),
    y =
      departure_brier
  )
) +
  geom_boxplot() +
  coord_flip() +
  labs(
    title =
      "Repeated process-level cross-validated departure Brier score",
    x = "Model",
    y =
      "Departure Brier score"
  ) +
  theme_bw()

save_plot(
  plot_cv_brier,
  "06_cv_departure_brier.png"
)

principal_calibration_models <- c(
  "Stage-Spell SPAOM",
  "Origin-specific homogeneous",
  "Dual-clock model",
  "Spline spell-age model"
)

plot_calibration <- CALIBRATION_TABLE %>%
  filter(
    model %in%
      principal_calibration_models
  ) %>%
  ggplot(
    aes(
      x = mean_predicted,
      y =
        observed_departure
    )
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = 2
  ) +
  geom_line() +
  geom_point(
    aes(
      size = transitions
    )
  ) +
  facet_wrap(
    ~ model
  ) +
  coord_equal(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  labs(
    title =
      "Cross-validated departure calibration",
    x =
      "Mean predicted departure probability",
    y =
      "Observed departure proportion",
    size =
      "Transitions"
  ) +
  theme_bw()

save_plot(
  plot_calibration,
  "07_departure_calibration.png",
  width = 11,
  height = 9
)

advantage_plot_data <- PAIRED_ADVANTAGES %>%
  filter(
    metric %in% c(
      "departure_logloss",
      "departure_brier"
    )
  ) %>%
  mutate(
    label = paste(
      comparison_model,
      metric,
      sep = " | "
    ),
    label = reorder(
      label,
      process_bootstrap_advantage
    )
  )

plot_advantages <- ggplot(
  advantage_plot_data,
  aes(
    x = label,
    y =
      process_bootstrap_advantage
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  geom_point(
    size = 3
  ) +
  geom_errorbar(
    aes(
      ymin =
        process_bootstrap_lower,
      ymax =
        process_bootstrap_upper
    ),
    width = 0.15
  ) +
  coord_flip() +
  labs(
    title =
      "Paired predictive advantages of SS-SPAOM",
    subtitle =
      "Positive values favour SS-SPAOM; bars are process-bootstrap 95% intervals",
    x =
      "Comparator and metric",
    y =
      "Predictive advantage"
  ) +
  theme_bw()

save_plot(
  plot_advantages,
  "08_paired_predictive_advantages.png",
  width = 11,
  height = 8
)

plot_beta_bootstrap <- PARAMETER_BOOTSTRAP_RESULTS %>%
  filter(
    converged,
    is.finite(
      beta_spell_hat
    )
  ) %>%
  ggplot(
    aes(
      x = beta_spell_hat
    )
  ) +
  geom_histogram(
    bins = 30
  ) +
  geom_vline(
    xintercept =
      FIT_STAGE_SPELL$beta_spell,
    linetype = 2
  ) +
  labs(
    title =
      "Process-bootstrap distribution of the stage-spell effect",
    subtitle =
      "Dashed line marks the full-data estimate",
    x = expression(
      hat(beta)[s]
    ),
    y = "Frequency"
  ) +
  theme_bw()

save_plot(
  plot_beta_bootstrap,
  "09_beta_process_bootstrap.png"
)

plot_null_bootstrap <- NULL_BOOTSTRAP_RESULT$
  replications %>%
  filter(
    is.finite(LR)
  ) %>%
  ggplot(
    aes(
      x = LR
    )
  ) +
  geom_histogram(
    bins = 30
  ) +
  geom_vline(
    xintercept =
      NULL_BOOTSTRAP_RESULT$
        summary$
        observed_LR[1L],
    linetype = 2
  ) +
  labs(
    title =
      "Parametric-bootstrap null distribution for the stage-spell LR statistic",
    subtitle =
      "Dashed line marks the observed statistic",
    x =
      "Likelihood-ratio statistic",
    y = "Frequency"
  ) +
  theme_bw()

save_plot(
  plot_null_bootstrap,
  "10_dynamic_test_bootstrap_distribution.png"
)

plot_robustness_beta <- ROBUSTNESS_RESULTS %>%
  mutate(
    scenario = reorder(
      scenario,
      beta_spell_hat
    )
  ) %>%
  ggplot(
    aes(
      x = scenario,
      y = beta_spell_hat
    )
  ) +
  geom_hline(
    yintercept =
      FIT_STAGE_SPELL$beta_spell,
    linetype = 2
  ) +
  geom_point(
    size = 3
  ) +
  coord_flip() +
  labs(
    title =
      "Robustness of the estimated stage-spell effect",
    subtitle =
      "Dashed line marks the primary estimate",
    x =
      "Robustness analysis",
    y = expression(
      hat(beta)[s]
    )
  ) +
  theme_bw()

save_plot(
  plot_robustness_beta,
  "11_robustness_beta.png"
)

top_influence <- INFLUENCE_RESULTS %>%
  filter(
    is.finite(beta_change)
  ) %>%
  slice_head(n = 20L) %>%
  mutate(
    PP = reorder(
      PP,
      beta_change
    )
  )

plot_influence <- ggplot(
  top_influence,
  aes(
    x = PP,
    y = beta_change,
    size =
      process_transitions
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  geom_point() +
  coord_flip() +
  labs(
    title =
      "Most influential peace processes for the stage-spell estimate",
    x =
      "Peace-process identifier",
    y =
      "Change in beta after omission",
    size =
      "Transitions"
  ) +
  theme_bw()

save_plot(
  plot_influence,
  "12_leave_one_process_out_influence.png"
)

all_plots <- list(
  plot_spell_empirical,
  plot_fitted_curves,
  plot_relative_decay,
  plot_delta_aic,
  plot_cv_logloss,
  plot_cv_brier,
  plot_calibration,
  plot_advantages,
  plot_beta_bootstrap,
  plot_null_bootstrap,
  plot_robustness_beta,
  plot_influence
)

pdf(
  file.path(
    OUTPUT_DIR,
    "SS_SPAOM_Mandatory_Graphs.pdf"
  ),
  width = 11,
  height = 8.5
)

for (plot_object in all_plots) {
  print(plot_object)
}

dev.off()


# -----------------------------------------------------------------------------
# 21. AUTOMATIC SCIENTIFIC VERDICT
# -----------------------------------------------------------------------------

homogeneous_logloss <- PAIRED_ADVANTAGES %>%
  filter(
    comparison_model ==
      "Origin-specific homogeneous",
    metric ==
      "departure_logloss"
  )

homogeneous_brier <- PAIRED_ADVANTAGES %>%
  filter(
    comparison_model ==
      "Origin-specific homogeneous",
    metric ==
      "departure_brier"
  )

spline_logloss <- PAIRED_ADVANTAGES %>%
  filter(
    comparison_model ==
      "Spline spell-age model",
    metric ==
      "departure_logloss"
  )

dynamic_p_value <- NULL_BOOTSTRAP_RESULT$
  summary$
  bootstrap_p_value[1L]

strong_dynamic_evidence <-
  is.finite(dynamic_p_value) &&
  dynamic_p_value < 0.05

predictive_logloss_evidence <-
  nrow(homogeneous_logloss) > 0L &&
  homogeneous_logloss$
    process_bootstrap_lower[1L] > 0

predictive_brier_evidence <-
  nrow(homogeneous_brier) > 0L &&
  homogeneous_brier$
    process_bootstrap_lower[1L] > 0

competitive_with_spline <-
  nrow(spline_logloss) > 0L &&
  spline_logloss$
    process_bootstrap_lower[1L] > -0.002

verdict <- if (
  strong_dynamic_evidence &&
  predictive_logloss_evidence &&
  predictive_brier_evidence
) {
  "STRONG: SS-SPAOM shows inferential and predictive evidence over the origin-specific homogeneous model."
} else if (
  predictive_logloss_evidence ||
  predictive_brier_evidence
) {
  "PROMISING: SS-SPAOM shows predictive improvement, but all inferential criteria are not yet decisive."
} else if (
  strong_dynamic_evidence
) {
  "MIXED: the stage-spell effect is statistically supported, but predictive superiority is not established."
} else {
  "NOT YET DECISIVE: the stage-spell effect or its predictive advantage requires stronger evidence."
}

VERDICT_TABLE <- data.frame(
  criterion = c(
    "Bootstrap p-value below 0.05",
    "Log-loss paired CI above zero",
    "Brier paired CI above zero",
    "Competitive with spell-age spline"
  ),
  satisfied = c(
    strong_dynamic_evidence,
    predictive_logloss_evidence,
    predictive_brier_evidence,
    competitive_with_spline
  ),
  stringsAsFactors = FALSE
)

save_csv(
  VERDICT_TABLE,
  "20_model_verdict_criteria.csv"
)

writeLines(
  c(
    verdict,
    "",
    paste0(
      "Estimated beta_spell: ",
      FIT_STAGE_SPELL$beta_spell
    ),
    paste0(
      "Parametric-bootstrap p-value: ",
      dynamic_p_value
    ),
    paste0(
      "Relative opportunity intensity at spell age 5: ",
      5^(
        -FIT_STAGE_SPELL$beta_spell
      )
    ),
    paste0(
      "Output directory: ",
      OUTPUT_DIR
    )
  ),
  file.path(
    OUTPUT_DIR,
    "MODEL_VERDICT.txt"
  )
)


# -----------------------------------------------------------------------------
# 22. LATEX RESULTS BLOCK
# -----------------------------------------------------------------------------

beta_summary <- PARAMETER_BOOTSTRAP_SUMMARY %>%
  filter(
    parameter ==
      "beta_spell"
  )

latex_lines <- c(
  "\\subsection{Stage-Spell Opportunity Dynamics}",
  "",
  "The primary opportunity intensity was specified as",
  "\\[",
  "\\log \\mu_{iar}=\\log\\eta_a-\\beta_s\\log s_{ir},",
  "\\]",
  "where \\(s_{ir}\\) is the number of consecutive documented transitions for which process \\(i\\) has occupied its current origin stage.",
  "",
  "The estimated common stage-spell parameter was",
  "\\[",
  paste0(
    "\\widehat{\\beta}_s=",
    format_number(
      FIT_STAGE_SPELL$beta_spell
    ),
    "."
  ),
  "\\]",
  "The process-bootstrap confidence interval was",
  "\\[",
  paste0(
    "\\left[",
    format_number(
      beta_summary$ci_lower[1L]
    ),
    ",\\,",
    format_number(
      beta_summary$ci_upper[1L]
    ),
    "\\right]."
  ),
  "\\]",
  paste0(
    "At spell age five, the estimated opportunity intensity relative to spell age one was \\(",
    format_number(
      5^(
        -FIT_STAGE_SPELL$beta_spell
      )
    ),
    "\\)."
  ),
  "",
  paste0(
    "The parametric-bootstrap likelihood-ratio statistic for \\(H_0:\\beta_s=0\\) was \\(",
    format_number(
      NULL_BOOTSTRAP_RESULT$
        summary$
        observed_LR[1L]
    ),
    "\\), with bootstrap \\(p=",
    format_number(
      dynamic_p_value
    ),
    "\\)."
  ),
  "",
  paste0(
    "\\textit{Automated diagnostic verdict:} ",
    verdict
  )
)

writeLines(
  latex_lines,
  file.path(
    OUTPUT_DIR,
    "Filled_SS_SPAOM_Results.tex"
  )
)


# -----------------------------------------------------------------------------
# 23. EXCEL WORKBOOK, R OBJECT, STATUS, AND SESSION INFORMATION
# -----------------------------------------------------------------------------

writexl::write_xlsx(
  list(
    Data_Summary =
      DATA_SUMMARY,
    Primary_Estimates =
      PRIMARY_ESTIMATES,
    Origin_Eta =
      ETA_ESTIMATES,
    Stage_Beta_Sensitivity =
      STAGE_SPECIFIC_BETA,
    Full_Model_Comparison =
      FULL_MODEL_COMPARISON,
    Dynamic_Test =
      NULL_BOOTSTRAP_RESULT$
        summary,
    Parameter_Bootstrap =
      PARAMETER_BOOTSTRAP_SUMMARY,
    CV_Model_Summary =
      CV_MODEL_SUMMARY,
    Paired_Advantages =
      PAIRED_ADVANTAGES,
    Calibration =
      CALIBRATION_TABLE,
    Robustness =
      ROBUSTNESS_RESULTS,
    Robustness_CV =
      ROBUSTNESS_CV_SUMMARY,
    Influence =
      INFLUENCE_RESULTS,
    Spell_Age_Description =
      SPELL_AGE_SUMMARY,
    Verdict =
      VERDICT_TABLE
  ),
  path = file.path(
    OUTPUT_DIR,
    "SS_SPAOM_Complete_Results.xlsx"
  )
)

saveRDS(
  list(
    data_file =
      DATA_FILE,
    analysis_mode =
      ANALYSIS_MODE,
    transitions =
      TRANSITIONS,
    homogeneous_fit =
      FIT_HOMOGENEOUS,
    stage_spell_fit =
      FIT_STAGE_SPELL,
    global_alpha_beta_fit =
      FIT_GLOBAL_ALPHA_BETA,
    dual_clock_fit =
      FIT_DUAL_CLOCK,
    old_kernel_fit =
      FIT_OLD_KERNEL,
    spline_spell_fit =
      FIT_SPLINE_SPELL,
    quadratic_r_fit =
      FIT_QUADRATIC_R,
    stage_specific_fit =
      FIT_STAGE_SPECIFIC,
    dynamic_test =
      NULL_BOOTSTRAP_RESULT,
    parameter_bootstrap =
      PARAMETER_BOOTSTRAP_RESULTS,
    cross_validation =
      CV_RESULT,
    paired_advantages =
      PAIRED_ADVANTAGES,
    robustness =
      ROBUSTNESS_RESULTS,
    robustness_cv =
      ROBUSTNESS_CV_SUMMARY,
    influence =
      INFLUENCE_RESULTS
  ),
  file.path(
    OUTPUT_DIR,
    "SS_SPAOM_Complete_Results.rds"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUTPUT_DIR,
    "R_session_info.txt"
  )
)

status_lines <- c(
  paste0(
    "Analysis mode: ",
    ANALYSIS_MODE
  ),
  paste0(
    "Input file: ",
    DATA_FILE
  ),
  paste0(
    "Transitions: ",
    nrow(TRANSITIONS)
  ),
  paste0(
    "Peace processes: ",
    dplyr::n_distinct(
      TRANSITIONS$PP
    )
  ),
  paste0(
    "Beta spell estimate: ",
    FIT_STAGE_SPELL$beta_spell
  ),
  paste0(
    "Bootstrap p-value: ",
    dynamic_p_value
  ),
  paste0(
    "Completed primary CV folds: ",
    nrow(
      unique(
        CV_RESULT$
          fold_metrics[
            ,
            c(
              "repeat_id",
              "fold_id"
            )
          ]
      )
    ),
    " / ",
    CV_REPEATS *
      CV_FOLDS
  ),
  paste0(
    "Successful null bootstraps: ",
    NULL_BOOTSTRAP_RESULT$
      summary$
      successful_replications[1L],
    " / ",
    NULL_BOOTSTRAP_RESULT$
      summary$
      requested_replications[1L]
  ),
  paste0(
    "Successful parameter bootstraps: ",
    sum(
      PARAMETER_BOOTSTRAP_RESULTS$
        converged
    ),
    " / ",
    PARAMETER_BOOTSTRAP_B
  ),
  paste0(
    "Verdict: ",
    verdict
  ),
  paste0(
    "Output directory: ",
    OUTPUT_DIR
  )
)

writeLines(
  status_lines,
  file.path(
    OUTPUT_DIR,
    "RUN_STATUS.txt"
  )
)

cat("\n============================================================\n")
cat("STAGE-SPELL SPAOM ANALYSIS COMPLETED\n")
cat("============================================================\n")
cat("Results are saved in:\n")
cat(OUTPUT_DIR, "\n\n")
cat("Open these files first:\n")
cat(
  file.path(
    OUTPUT_DIR,
    "MODEL_VERDICT.txt"
  ),
  "\n"
)
cat(
  file.path(
    OUTPUT_DIR,
    "SS_SPAOM_Complete_Results.xlsx"
  ),
  "\n"
)
cat(
  file.path(
    OUTPUT_DIR,
    "SS_SPAOM_Mandatory_Graphs.pdf"
  ),
  "\n"
)
cat(
  file.path(
    OUTPUT_DIR,
    "Filled_SS_SPAOM_Results.tex"
  ),
  "\n"
)
cat(
  file.path(
    OUTPUT_DIR,
    "RUN_STATUS.txt"
  ),
  "\n"
)
