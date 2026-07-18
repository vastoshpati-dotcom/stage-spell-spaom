# =============================================================================
# REVISED MONTE CARLO STUDY FOR THE STAGE-SPELL SPAOM
# =============================================================================
#
# Data-generating departure model:
#
#   log(mu_iar) = log(eta_a) - beta_s log(s_ir)
#   q_ia(r)      = 1 - exp[-mu_iar]
#
# The study separates origin-stage persistence (eta_a) from within-stage
# duration dependence (beta_s).  It compares:
#
#   1. Origin-specific homogeneous model
#   2. Stage-Spell SPAOM (primary)
#   3. Dual-clock model: spell age plus global agreement position
#   4. Global alpha-beta model
#   5. Natural spline in log stage-spell age
#
# Primary outcomes:
#   - recovery, bias, RMSE, coverage, and boundary rate for beta_s;
#   - recovery of eta_a;
#   - test-sample departure log loss, Brier score, and calibration;
#   - complete-transition log loss and Brier score;
#   - paired predictive advantages of SS-SPAOM.
#
# The destination matrix is estimated identically for every departure model.
#
# =============================================================================


# -----------------------------------------------------------------------------
# 0. USER SETTINGS
# -----------------------------------------------------------------------------

BASE_DIR <- "C:/Users/Hp/Downloads"

TRANSITION_FILE <- NA_character_

SIMULATION_MODE <- "FINAL"  # Publication simulation

OUTPUT_DIR <- file.path(
  BASE_DIR,
  paste0(
    "SS_SPAOM_Revised_Simulation_",
    SIMULATION_MODE
  )
)

CHECKPOINT_DIR <- file.path(
  OUTPUT_DIR,
  "Scenario_Checkpoints"
)

GRAPH_DIR <- file.path(
  OUTPUT_DIR,
  "Graphs"
)

if (SIMULATION_MODE == "PILOT") {
  B_REPLICATIONS <- 10L
  N_VALUES <- c(100L, 184L)
  MEAN_LENGTH_VALUES <- c(8L, 20L)
  BETA_VALUES <- c(0, 0.52)
} else {
  B_REPLICATIONS <- 300L
  N_VALUES <- c(100L, 184L, 300L, 500L)
  MEAN_LENGTH_VALUES <- c(8L, 12L, 20L)
  BETA_VALUES <- c(0, 0.25, 0.52, 0.80)
}

BASELINE_DESIGNS <- c(
  "common",
  "unequal"
)

USE_PARALLEL <- TRUE

N_WORKERS <- max(
  1L,
  min(
    4L,
    parallel::detectCores() - 1L
  )
)

BASE_SEED <- 20260714L

DESTINATION_SMOOTHING <- 0.50
PROBABILITY_FLOOR <- 1e-12
SPLINE_DF <- 4L

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
# 1. PACKAGES
# -----------------------------------------------------------------------------

REQUIRED_PACKAGES <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "writexl",
  "future",
  "future.apply"
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
# 2. HELPERS
# -----------------------------------------------------------------------------

clip_probability <- function(probability) {
  pmin(
    pmax(
      probability,
      PROBABILITY_FLOOR
    ),
    1 - PROBABILITY_FLOOR
  )
}

inverse_cloglog <- function(eta) {
  clip_probability(
    -expm1(
      -exp(
        pmin(
          eta,
          700
        )
      )
    )
  )
}

safe_mean <- function(x) {
  x <- x[
    is.finite(x)
  ]

  if (length(x) == 0L) {
    NA_real_
  } else {
    mean(x)
  }
}

safe_sd <- function(x) {
  x <- x[
    is.finite(x)
  ]

  if (length(x) <= 1L) {
    NA_real_
  } else {
    sd(x)
  }
}

append_error <- function(
  stage,
  identifier,
  message
) {
  error_file <- file.path(
    OUTPUT_DIR,
    "simulation_errors.csv"
  )

  row <- data.frame(
    time = as.character(
      Sys.time()
    ),
    stage = as.character(stage),
    identifier = as.character(
      identifier
    ),
    message = as.character(message),
    stringsAsFactors = FALSE
  )

  write.table(
    row,
    file = error_file,
    sep = ",",
    row.names = FALSE,
    col.names =
      !file.exists(
        error_file
      ),
    append =
      file.exists(
        error_file
      ),
    quote = TRUE
  )
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


# -----------------------------------------------------------------------------
# 3. LOCATE AND PREPARE THE EMPIRICAL CALIBRATION DATA
# -----------------------------------------------------------------------------

locate_transition_file <- function() {
  required_columns <- c(
    "PP",
    "origin_stage",
    "destination_stage",
    "D",
    "transition_index",
    "transition_id"
  )

  validate_file <- function(path) {
    if (
      !file.exists(path) ||
      file.info(path)$size <= 0
    ) {
      return(FALSE)
    }

    test <- tryCatch(
      read.csv(
        path,
        nrows = 10L,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      error = function(e) NULL
    )

    !is.null(test) &&
      all(
        required_columns %in%
          names(test)
      )
  }

  if (
    is.character(TRANSITION_FILE) &&
    length(TRANSITION_FILE) == 1L &&
    !is.na(TRANSITION_FILE) &&
    nzchar(TRANSITION_FILE) &&
    validate_file(
      TRANSITION_FILE
    )
  ) {
    return(
      normalizePath(
        TRANSITION_FILE,
        winslash = "/",
        mustWork = TRUE
      )
    )
  }

  candidates <- list.files(
    BASE_DIR,
    pattern =
      "^01_primary_transition_dataset.*\\.csv$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )

  candidates <- candidates[
    file.exists(candidates)
  ]

  if (length(candidates) == 0L) {
    stop(
      "No primary transition dataset was found under ",
      BASE_DIR,
      "."
    )
  }

  candidate_info <- file.info(
    candidates
  )

  candidates <- candidates[
    order(
      candidate_info$mtime,
      decreasing = TRUE,
      na.last = TRUE
    )
  ]

  for (candidate in candidates) {
    if (validate_file(candidate)) {
      return(
        normalizePath(
          candidate,
          winslash = "/",
          mustWork = TRUE
        )
      )
    }
  }

  stop(
    "Transition candidates were found but none was valid."
  )
}

add_spell_age <- function(
  data,
  stage_levels
) {
  data %>%
    mutate(
      PP = as.character(PP),
      origin_stage = factor(
        as.character(
          origin_stage
        ),
        levels = stage_levels
      ),
      destination_stage = factor(
        as.character(
          destination_stage
        ),
        levels = stage_levels
      ),
      D = as.integer(D),
      transition_index =
        as.numeric(
          transition_index
        )
    ) %>%
    filter(
      !is.na(PP),
      !is.na(origin_stage),
      !is.na(destination_stage),
      D %in% c(0L, 1L),
      is.finite(
        transition_index
      )
    ) %>%
    arrange(
      PP,
      transition_index
    ) %>%
    group_by(PP) %>%
    mutate(
      origin_character =
        as.character(
          origin_stage
        ),
      spell_start =
        row_number() == 1L |
        origin_character !=
          dplyr::lag(
            origin_character,
            default =
              "__START__"
          ),
      spell_id =
        cumsum(
          spell_start
        )
    ) %>%
    group_by(
      PP,
      spell_id
    ) %>%
    mutate(
      spell_age =
        row_number()
    ) %>%
    ungroup() %>%
    mutate(
      log_spell_age =
        log(spell_age),
      log_transition_index =
        log(
          transition_index
        )
    ) %>%
    select(
      -origin_character,
      -spell_start
    )
}

DATA_FILE <- locate_transition_file()

EMPIRICAL_RAW <- read.csv(
  DATA_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

PREFERRED_STAGES <- c(
  "Pre-negotiation",
  "Ceasefire",
  "Substantive Partial",
  "Substantive Comprehensive",
  "Implementation",
  "Renegotiation/Other"
)

observed_stages <- unique(
  c(
    as.character(
      EMPIRICAL_RAW$origin_stage
    ),
    as.character(
      EMPIRICAL_RAW$destination_stage
    )
  )
)

STAGE_LEVELS <- c(
  PREFERRED_STAGES[
    PREFERRED_STAGES %in%
      observed_stages
  ],
  sort(
    setdiff(
      observed_stages,
      PREFERRED_STAGES
    )
  )
)

EMPIRICAL_DATA <- add_spell_age(
  EMPIRICAL_RAW,
  STAGE_LEVELS
)


# -----------------------------------------------------------------------------
# 4. CALIBRATION OF ETA, INITIAL STAGES, AND DESTINATIONS
# -----------------------------------------------------------------------------

estimate_eta_at_spell_one <- function(
  data,
  stage_levels
) {
  overall_probability <- clip_probability(
    mean(
      data$D[
        data$spell_age == 1
      ]
    )
  )

  stage_probability <- vapply(
    stage_levels,
    function(stage) {
      selected <-
        as.character(
          data$origin_stage
        ) == stage &
        data$spell_age == 1

      successes <- sum(
        data$D[selected]
      )

      total <- sum(selected)

      if (total == 0L) {
        overall_probability
      } else {
        clip_probability(
          (
            successes + 0.5
          ) /
            (
              total + 1
            )
        )
      }
    },
    numeric(1)
  )

  setNames(
    -log(
      1 - stage_probability
    ),
    stage_levels
  )
}

estimate_initial_stage_distribution <- function(
  data,
  stage_levels
) {
  initial <- data %>%
    arrange(
      PP,
      transition_index
    ) %>%
    group_by(PP) %>%
    slice_head(n = 1L) %>%
    ungroup()

  counts <- table(
    factor(
      initial$origin_stage,
      levels = stage_levels
    )
  )

  probabilities <- (
    as.numeric(counts) +
      0.5
  ) /
    sum(
      as.numeric(counts) +
        0.5
    )

  setNames(
    probabilities,
    stage_levels
  )
}

estimate_destination_matrix <- function(
  data,
  stage_levels,
  smoothing =
    DESTINATION_SMOOTHING
) {
  output <- matrix(
    0,
    nrow =
      length(stage_levels),
    ncol =
      length(stage_levels),
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

    selected <-
      data$D == 1L &
      as.character(
        data$origin_stage
      ) == origin

    observed <- table(
      as.character(
        data$destination_stage[
          selected
        ]
      )
    )

    common <- intersect(
      names(observed),
      admissible
    )

    counts[common] <-
      counts[common] +
      as.numeric(
        observed[common]
      )

    output[
      origin,
      admissible
    ] <- counts /
      sum(counts)
  }

  output
}

ETA_UNEQUAL <- estimate_eta_at_spell_one(
  EMPIRICAL_DATA,
  STAGE_LEVELS
)

ETA_COMMON_VALUE <- exp(
  mean(
    log(
      ETA_UNEQUAL
    )
  )
)

ETA_COMMON <- setNames(
  rep(
    ETA_COMMON_VALUE,
    length(STAGE_LEVELS)
  ),
  STAGE_LEVELS
)

INITIAL_STAGE_PROBABILITY <-
  estimate_initial_stage_distribution(
    EMPIRICAL_DATA,
    STAGE_LEVELS
  )

TRUE_DESTINATION_MATRIX <-
  estimate_destination_matrix(
    EMPIRICAL_DATA,
    STAGE_LEVELS
  )

CALIBRATION_TABLE <- data.frame(
  origin_stage =
    STAGE_LEVELS,
  eta_unequal =
    as.numeric(
      ETA_UNEQUAL[
        STAGE_LEVELS
      ]
    ),
  eta_common =
    as.numeric(
      ETA_COMMON[
        STAGE_LEVELS
      ]
    ),
  initial_stage_probability =
    as.numeric(
      INITIAL_STAGE_PROBABILITY[
        STAGE_LEVELS
      ]
    ),
  stringsAsFactors = FALSE
)

save_csv(
  CALIBRATION_TABLE,
  "00_simulation_calibration.csv"
)


# -----------------------------------------------------------------------------
# 5. SIMULATE MULTISTATE PEACE-PROCESS HISTORIES
# -----------------------------------------------------------------------------

simulate_process_length <- function(
  mean_length
) {
  max(
    2L,
    rpois(
      1L,
      lambda =
        max(
          0.1,
          mean_length - 2
        )
    ) +
      2L
  )
}

simulate_transition_dataset <- function(
  n_processes,
  mean_length,
  beta_true,
  eta_true,
  stage_levels,
  initial_stage_probability,
  destination_matrix,
  dataset_label
) {
  process_rows <- vector(
    "list",
    n_processes
  )

  for (
    process_index in
    seq_len(n_processes)
  ) {
    process_length <-
      simulate_process_length(
        mean_length
      )

    current_stage <- sample(
      stage_levels,
      size = 1L,
      prob =
        initial_stage_probability[
          stage_levels
        ]
    )

    spell_age <- 1L

    rows <- vector(
      "list",
      process_length
    )

    for (
      transition_index in
      seq_len(
        process_length
      )
    ) {
      mu <- eta_true[
        current_stage
      ] *
        spell_age^(
          -beta_true
        )

      departure_probability <-
        clip_probability(
          1 -
            exp(-mu)
        )

      departure <- rbinom(
        1L,
        size = 1L,
        prob =
          departure_probability
      )

      if (departure == 1L) {
        admissible <- setdiff(
          stage_levels,
          current_stage
        )

        destination_stage <- sample(
          admissible,
          size = 1L,
          prob =
            destination_matrix[
              current_stage,
              admissible
            ]
        )
      } else {
        destination_stage <-
          current_stage
      }

      rows[[transition_index]] <-
        data.frame(
          PP = paste0(
            dataset_label,
            "_",
            process_index
          ),
          transition_id = paste0(
            dataset_label,
            "_",
            process_index,
            "::",
            transition_index
          ),
          transition_index =
            transition_index,
          origin_stage =
            current_stage,
          destination_stage =
            destination_stage,
          D = departure,
          spell_age =
            spell_age,
          log_spell_age =
            log(spell_age),
          log_transition_index =
            log(
              transition_index
            ),
          true_departure_probability =
            departure_probability,
          stringsAsFactors = FALSE
        )

      if (
        destination_stage ==
          current_stage
      ) {
        spell_age <- spell_age + 1L
      } else {
        current_stage <-
          destination_stage

        spell_age <- 1L
      }
    }

    process_rows[[process_index]] <-
      bind_rows(rows)
  }

  bind_rows(
    process_rows
  ) %>%
    mutate(
      origin_stage = factor(
        origin_stage,
        levels = stage_levels
      ),
      destination_stage = factor(
        destination_stage,
        levels = stage_levels
      )
    )
}


# -----------------------------------------------------------------------------
# 6. FIT THE COMPETING DEPARTURE MODELS
# -----------------------------------------------------------------------------

fit_glm_safely <- function(
  formula_object,
  data
) {
  tryCatch(
    suppressWarnings(
      glm(
        formula_object,
        data = data,
        family =
          binomial(
            link =
              "cloglog"
          ),
        control =
          glm.control(
            maxit = 100
          )
      )
    ),
    error = function(e) NULL
  )
}

fit_homogeneous_model <- function(
  data
) {
  fit <- fit_glm_safely(
    D ~ 0 + origin_stage,
    data
  )

  list(
    model =
      "Origin-specific homogeneous",
    fit = fit,
    converged =
      !is.null(fit) &&
      isTRUE(
        fit$converged
      )
  )
}

fit_stage_spell_model <- function(
  data,
  homogeneous_fit
) {
  full_fit <- fit_glm_safely(
    D ~ 0 +
      origin_stage +
      log_spell_age,
    data
  )

  if (is.null(full_fit)) {
    return(
      list(
        model =
          "Stage-Spell SPAOM",
        fit =
          homogeneous_fit$fit,
        raw_fit =
          NULL,
        beta_hat =
          0,
        beta_se =
          NA_real_,
        beta_ci_lower =
          0,
        beta_ci_upper =
          NA_real_,
        boundary =
          TRUE,
        converged =
          FALSE
      )
    )
  }

  coefficient_name <-
    "log_spell_age"

  coefficient <- coef(
    full_fit
  )[
    coefficient_name
  ]

  standard_error <- tryCatch(
    sqrt(
      diag(
        vcov(full_fit)
      )
    )[
      coefficient_name
    ],
    error = function(e) {
      NA_real_
    }
  )

  beta_unconstrained <-
    -unname(coefficient)

  if (
    !is.finite(
      beta_unconstrained
    ) ||
    beta_unconstrained <= 0
  ) {
    beta_hat <- 0

    selected_fit <-
      homogeneous_fit$fit

    boundary <- TRUE
  } else {
    beta_hat <-
      beta_unconstrained

    selected_fit <-
      full_fit

    boundary <- FALSE
  }

  beta_ci_lower <- max(
    0,
    beta_hat -
      1.96 *
      standard_error
  )

  beta_ci_upper <- max(
    0,
    beta_hat +
      1.96 *
      standard_error
  )

  list(
    model =
      "Stage-Spell SPAOM",
    fit =
      selected_fit,
    raw_fit =
      full_fit,
    beta_hat =
      beta_hat,
    beta_se =
      standard_error,
    beta_ci_lower =
      beta_ci_lower,
    beta_ci_upper =
      beta_ci_upper,
    boundary =
      boundary,
    converged =
      isTRUE(
        full_fit$converged
      )
  )
}

fit_dual_clock_model <- function(
  data
) {
  fit <- fit_glm_safely(
    D ~ 0 +
      origin_stage +
      log_spell_age +
      log_transition_index,
    data
  )

  list(
    model =
      "Dual-clock model",
    fit = fit,
    converged =
      !is.null(fit) &&
      isTRUE(
        fit$converged
      )
  )
}

fit_global_alpha_beta_model <- function(
  data
) {
  fit <- fit_glm_safely(
    D ~ 0 +
      origin_stage +
      log_transition_index +
      I(
        transition_index - 1
      ),
    data
  )

  list(
    model =
      "Global alpha-beta model",
    fit = fit,
    converged =
      !is.null(fit) &&
      isTRUE(
        fit$converged
      )
  )
}

fit_spline_spell_model <- function(
  data
) {
  fit <- fit_glm_safely(
    as.formula(
      paste0(
        "D ~ 0 + origin_stage + ",
        "splines::ns(log_spell_age, df = ",
        SPLINE_DF,
        ")"
      )
    ),
    data
  )

  list(
    model =
      "Spline spell-age model",
    fit = fit,
    converged =
      !is.null(fit) &&
      isTRUE(
        fit$converged
      )
  )
}

fit_all_departure_models <- function(
  training_data
) {
  homogeneous <-
    fit_homogeneous_model(
      training_data
    )

  list(
    homogeneous =
      homogeneous,
    stage_spell =
      fit_stage_spell_model(
        training_data,
        homogeneous
      ),
    dual_clock =
      fit_dual_clock_model(
        training_data
      ),
    global_alpha_beta =
      fit_global_alpha_beta_model(
        training_data
      ),
    spline_spell =
      fit_spline_spell_model(
        training_data
      )
  )
}

predict_model <- function(
  fit_object,
  new_data
) {
  if (
    is.null(fit_object$fit)
  ) {
    return(
      rep(
        NA_real_,
        nrow(new_data)
      )
    )
  }

  prediction <- tryCatch(
    suppressWarnings(
      predict(
        fit_object$fit,
        newdata = new_data,
        type = "response"
      )
    ),
    error = function(e) {
      rep(
        NA_real_,
        nrow(new_data)
      )
    }
  )

  clip_probability(
    prediction
  )
}


# -----------------------------------------------------------------------------
# 7. ETA RECOVERY AND CONFIDENCE INTERVALS
# -----------------------------------------------------------------------------

extract_eta_recovery <- function(
  stage_spell_fit,
  stage_levels,
  eta_true,
  replication_id,
  scenario
) {
  selected_fit <-
    stage_spell_fit$fit

  if (is.null(selected_fit)) {
    return(
      data.frame()
    )
  }

  new_data <- data.frame(
    origin_stage = factor(
      stage_levels,
      levels = stage_levels
    ),
    spell_age = 1,
    log_spell_age = 0,
    transition_index = 1,
    log_transition_index = 0
  )

  prediction <- tryCatch(
    predict(
      selected_fit,
      newdata = new_data,
      type = "link",
      se.fit = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(prediction)) {
    return(
      data.frame()
    )
  }

  eta_hat <- exp(
    as.numeric(
      prediction$fit
    )
  )

  eta_lower <- exp(
    as.numeric(
      prediction$fit -
        1.96 *
        prediction$se.fit
    )
  )

  eta_upper <- exp(
    as.numeric(
      prediction$fit +
        1.96 *
        prediction$se.fit
    )
  )

  data.frame(
    scenario_id =
      scenario$scenario_id,
    replication_id =
      replication_id,
    baseline_design =
      scenario$baseline_design,
    n_processes =
      scenario$n_processes,
    mean_length =
      scenario$mean_length,
    beta_true =
      scenario$beta_true,
    origin_stage =
      stage_levels,
    eta_true =
      as.numeric(
        eta_true[
          stage_levels
        ]
      ),
    eta_hat =
      eta_hat,
    eta_ci_lower =
      eta_lower,
    eta_ci_upper =
      eta_upper,
    eta_covered =
      as.numeric(
        eta_true[
          stage_levels
        ]
      ) >=
        eta_lower &
      as.numeric(
        eta_true[
          stage_levels
        ]
      ) <=
        eta_upper,
    stringsAsFactors = FALSE
  )
}


# -----------------------------------------------------------------------------
# 8. DESTINATION MODEL AND PREDICTIVE METRICS
# -----------------------------------------------------------------------------

fit_destination_matrix <- function(
  data,
  stage_levels,
  smoothing =
    DESTINATION_SMOOTHING
) {
  output <- matrix(
    0,
    nrow =
      length(stage_levels),
    ncol =
      length(stage_levels),
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

    selected <-
      data$D == 1L &
      as.character(
        data$origin_stage
      ) == origin

    observed <- table(
      as.character(
        data$destination_stage[
          selected
        ]
      )
    )

    common <- intersect(
      names(observed),
      admissible
    )

    counts[common] <-
      counts[common] +
      as.numeric(
        observed[common]
      )

    output[
      origin,
      admissible
    ] <- counts /
      sum(counts)
  }

  output
}

make_transition_probabilities <- function(
  data,
  departure_probability,
  destination_matrix,
  stage_levels
) {
  probability_matrix <- matrix(
    0,
    nrow = nrow(data),
    ncol =
      length(stage_levels),
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

    probability_matrix[
      row_index,
      origin
    ] <-
      1 -
      departure_probability[
        row_index
      ]

    admissible <- setdiff(
      stage_levels,
      origin
    )

    probability_matrix[
      row_index,
      admissible
    ] <-
      departure_probability[
        row_index
      ] *
      destination_matrix[
        origin,
        admissible
      ]
  }

  probability_matrix /
    rowSums(
      probability_matrix
    )
}

calibration_intercept_slope <- function(
  actual,
  predicted
) {
  valid <- is.finite(actual) &
    is.finite(predicted)

  actual <- actual[valid]
  predicted <- predicted[valid]

  if (
    length(actual) < 10L ||
    length(unique(actual)) < 2L ||
    sd(predicted) < 1e-8
  ) {
    return(
      c(
        intercept = NA_real_,
        slope = NA_real_
      )
    )
  }

  fit <- tryCatch(
    suppressWarnings(
      glm(
        actual ~ qlogis(
          clip_probability(predicted)
        ),
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

  coefficients <- coef(fit)

  if (
    length(coefficients) < 2L ||
    any(!is.finite(coefficients))
  ) {
    return(
      c(
        intercept = NA_real_,
        slope = NA_real_
      )
    )
  }

  c(
    intercept =
      unname(coefficients[1L]),
    slope =
      unname(coefficients[2L])
  )
}

evaluate_model <- function(
  test_data,
  departure_probability,
  destination_matrix,
  stage_levels
) {
  if (
    any(
      !is.finite(
        departure_probability
      )
    )
  ) {
    return(
      c(
        departure_logloss =
          NA_real_,
        departure_brier =
          NA_real_,
        calibration_intercept =
          NA_real_,
        calibration_slope =
          NA_real_,
        full_logloss =
          NA_real_,
        full_brier =
          NA_real_
      )
    )
  }

  departure_probability <-
    clip_probability(
      departure_probability
    )

  transition_probability_matrix <-
    make_transition_probabilities(
      test_data,
      departure_probability,
      destination_matrix,
      stage_levels
    )

  actual_index <- match(
    as.character(
      test_data$destination_stage
    ),
    stage_levels
  )

  row_index <- seq_len(
    nrow(test_data)
  )

  actual_probability <-
    transition_probability_matrix[
      cbind(
        row_index,
        actual_index
      )
    ]

  one_hot <- matrix(
    0,
    nrow =
      nrow(test_data),
    ncol =
      length(stage_levels)
  )

  one_hot[
    cbind(
      row_index,
      actual_index
    )
  ] <- 1

  calibration <-
    calibration_intercept_slope(
      test_data$D,
      departure_probability
    )

  c(
    departure_logloss =
      mean(
        -(
          test_data$D *
            log(
              departure_probability
            ) +
          (1 -
             test_data$D) *
            log(
              1 -
                departure_probability
            )
        )
      ),
    departure_brier =
      mean(
        (
          test_data$D -
            departure_probability
        )^2
      ),
    calibration_intercept =
      unname(
        calibration[
          "intercept"
        ]
      ),
    calibration_slope =
      unname(
        calibration[
          "slope"
        ]
      ),
    full_logloss =
      mean(
        -log(
          pmax(
            actual_probability,
            PROBABILITY_FLOOR
          )
        )
      ),
    full_brier =
      mean(
        rowSums(
          (
            one_hot -
              transition_probability_matrix
          )^2
        )
      )
  )
}


# -----------------------------------------------------------------------------
# 9. ONE MONTE CARLO REPLICATION
# -----------------------------------------------------------------------------

run_one_replication <- function(
  replication_id,
  scenario
) {
  tryCatch(
    {
      set.seed(
        scenario$scenario_seed +
          replication_id
      )

      eta_true <- if (
        scenario$baseline_design ==
          "common"
      ) {
        ETA_COMMON
      } else {
        ETA_UNEQUAL
      }

      training_data <-
        simulate_transition_dataset(
          n_processes =
            scenario$n_processes,
          mean_length =
            scenario$mean_length,
          beta_true =
            scenario$beta_true,
          eta_true =
            eta_true,
          stage_levels =
            STAGE_LEVELS,
          initial_stage_probability =
            INITIAL_STAGE_PROBABILITY,
          destination_matrix =
            TRUE_DESTINATION_MATRIX,
          dataset_label = paste0(
            "train_",
            scenario$scenario_id,
            "_",
            replication_id
          )
        )

      test_data <-
        simulate_transition_dataset(
          n_processes =
            scenario$n_processes,
          mean_length =
            scenario$mean_length,
          beta_true =
            scenario$beta_true,
          eta_true =
            eta_true,
          stage_levels =
            STAGE_LEVELS,
          initial_stage_probability =
            INITIAL_STAGE_PROBABILITY,
          destination_matrix =
            TRUE_DESTINATION_MATRIX,
          dataset_label = paste0(
            "test_",
            scenario$scenario_id,
            "_",
            replication_id
          )
        )

      fitted_models <-
        fit_all_departure_models(
          training_data
        )

      estimated_destination_matrix <-
        fit_destination_matrix(
          training_data,
          STAGE_LEVELS
        )

      prediction_list <- list(
        `Origin-specific homogeneous` =
          predict_model(
            fitted_models$homogeneous,
            test_data
          ),
        `Stage-Spell SPAOM` =
          predict_model(
            fitted_models$stage_spell,
            test_data
          ),
        `Dual-clock model` =
          predict_model(
            fitted_models$dual_clock,
            test_data
          ),
        `Global alpha-beta model` =
          predict_model(
            fitted_models$global_alpha_beta,
            test_data
          ),
        `Spline spell-age model` =
          predict_model(
            fitted_models$spline_spell,
            test_data
          )
      )

      model_rows <- lapply(
        names(
          prediction_list
        ),
        function(model_name) {
          metrics <- evaluate_model(
            test_data,
            prediction_list[[model_name]],
            estimated_destination_matrix,
            STAGE_LEVELS
          )

          data.frame(
            scenario_id =
              scenario$scenario_id,
            replication_id =
              replication_id,
            baseline_design =
              scenario$baseline_design,
            n_processes =
              scenario$n_processes,
            mean_length =
              scenario$mean_length,
            beta_true =
              scenario$beta_true,
            model =
              model_name,
            departure_logloss =
              unname(
                metrics[
                  "departure_logloss"
                ]
              ),
            departure_brier =
              unname(
                metrics[
                  "departure_brier"
                ]
              ),
            calibration_intercept =
              unname(
                metrics[
                  "calibration_intercept"
                ]
              ),
            calibration_slope =
              unname(
                metrics[
                  "calibration_slope"
                ]
              ),
            full_logloss =
              unname(
                metrics[
                  "full_logloss"
                ]
              ),
            full_brier =
              unname(
                metrics[
                  "full_brier"
                ]
              ),
            stringsAsFactors = FALSE
          )
        }
      )

      stage_spell_fit <-
        fitted_models$stage_spell

      parameter_row <- data.frame(
        scenario_id =
          scenario$scenario_id,
        replication_id =
          replication_id,
        baseline_design =
          scenario$baseline_design,
        n_processes =
          scenario$n_processes,
        mean_length =
          scenario$mean_length,
        beta_true =
          scenario$beta_true,
        beta_hat =
          stage_spell_fit$beta_hat,
        beta_se =
          stage_spell_fit$beta_se,
        beta_ci_lower =
          stage_spell_fit$
            beta_ci_lower,
        beta_ci_upper =
          stage_spell_fit$
            beta_ci_upper,
        beta_covered =
          scenario$beta_true >=
            stage_spell_fit$
              beta_ci_lower &&
          scenario$beta_true <=
            stage_spell_fit$
              beta_ci_upper,
        boundary_estimate =
          stage_spell_fit$boundary,
        converged =
          stage_spell_fit$converged,
        training_transitions =
          nrow(training_data),
        test_transitions =
          nrow(test_data),
        stringsAsFactors = FALSE
      )

      eta_rows <-
        extract_eta_recovery(
          stage_spell_fit,
          STAGE_LEVELS,
          eta_true,
          replication_id,
          scenario
        )

      list(
        model_metrics =
          bind_rows(
            model_rows
          ),
        parameter_results =
          parameter_row,
        eta_results =
          eta_rows,
        error = NULL
      )
    },
    error = function(e) {
      list(
        model_metrics =
          data.frame(),
        parameter_results =
          data.frame(),
        eta_results =
          data.frame(),
        error = data.frame(
          scenario_id =
            scenario$scenario_id,
          replication_id =
            replication_id,
          message =
            conditionMessage(e),
          stringsAsFactors = FALSE
        )
      )
    }
  )
}


# -----------------------------------------------------------------------------
# 10. SCENARIO GRID AND RESUMABLE EXECUTION
# -----------------------------------------------------------------------------

SCENARIO_GRID <- expand.grid(
  baseline_design =
    BASELINE_DESIGNS,
  n_processes =
    N_VALUES,
  mean_length =
    MEAN_LENGTH_VALUES,
  beta_true =
    BETA_VALUES,
  stringsAsFactors = FALSE
) %>%
  arrange(
    baseline_design,
    n_processes,
    mean_length,
    beta_true
  ) %>%
  mutate(
    scenario_id =
      row_number(),
    scenario_seed =
      BASE_SEED +
      scenario_id *
      100000L
  ) %>%
  select(
    scenario_id,
    everything()
  )

save_csv(
  SCENARIO_GRID,
  "01_scenario_grid.csv"
)

if (
  USE_PARALLEL &&
  N_WORKERS > 1L
) {
  future::plan(
    future::multisession,
    workers = N_WORKERS
  )
} else {
  future::plan(
    future::sequential
  )
}

for (
  scenario_index in
  seq_len(
    nrow(
      SCENARIO_GRID
    )
  )
) {
  scenario <- as.list(
    SCENARIO_GRID[
      scenario_index,
      ,
      drop = FALSE
    ]
  )

  scenario_id <-
    scenario$scenario_id

  checkpoint_file <- file.path(
    CHECKPOINT_DIR,
    paste0(
      "scenario_",
      sprintf(
        "%03d",
        scenario_id
      ),
      ".rds"
    )
  )

  if (file.exists(checkpoint_file)) {
    cat(
      "Scenario ",
      scenario_id,
      " already completed.\n",
      sep = ""
    )

    next
  }

  cat(
    "\nScenario ",
    scenario_id,
    " of ",
    nrow(SCENARIO_GRID),
    ": design=",
    scenario$baseline_design,
    ", n=",
    scenario$n_processes,
    ", mean R=",
    scenario$mean_length,
    ", beta=",
    scenario$beta_true,
    "\n",
    sep = ""
  )

  replication_results <-
    future.apply::future_lapply(
      seq_len(
        B_REPLICATIONS
      ),
      function(replication_id) {
        run_one_replication(
          replication_id,
          scenario
        )
      },
      future.seed = TRUE
    )

  scenario_result <- list(
    model_metrics = bind_rows(
      lapply(
        replication_results,
        `[[`,
        "model_metrics"
      )
    ),
    parameter_results = bind_rows(
      lapply(
        replication_results,
        `[[`,
        "parameter_results"
      )
    ),
    eta_results = bind_rows(
      lapply(
        replication_results,
        `[[`,
        "eta_results"
      )
    ),
    errors = bind_rows(
      lapply(
        replication_results,
        `[[`,
        "error"
      )
    )
  )

  saveRDS(
    scenario_result,
    checkpoint_file
  )

  if (
    nrow(
      scenario_result$errors
    ) > 0L
  ) {
    write.table(
      scenario_result$errors,
      file = file.path(
        OUTPUT_DIR,
        "simulation_errors.csv"
      ),
      sep = ",",
      row.names = FALSE,
      col.names =
        !file.exists(
          file.path(
            OUTPUT_DIR,
            "simulation_errors.csv"
          )
        ),
      append =
        file.exists(
          file.path(
            OUTPUT_DIR,
            "simulation_errors.csv"
          )
        ),
      quote = TRUE
    )
  }
}

future::plan(
  future::sequential
)

checkpoint_files <- list.files(
  CHECKPOINT_DIR,
  pattern =
    "^scenario_[0-9]+\\.rds$",
  full.names = TRUE
)

if (length(checkpoint_files) == 0L) {
  stop(
    "No scenario checkpoint files were produced."
  )
}

all_scenarios <- lapply(
  checkpoint_files,
  readRDS
)

MODEL_RESULTS <- bind_rows(
  lapply(
    all_scenarios,
    `[[`,
    "model_metrics"
  )
)

PARAMETER_RESULTS <- bind_rows(
  lapply(
    all_scenarios,
    `[[`,
    "parameter_results"
  )
)

ETA_RESULTS <- bind_rows(
  lapply(
    all_scenarios,
    `[[`,
    "eta_results"
  )
)

ERROR_RESULTS <- bind_rows(
  lapply(
    all_scenarios,
    `[[`,
    "errors"
  )
)


# -----------------------------------------------------------------------------
# 11. MONTE CARLO SUMMARIES
# -----------------------------------------------------------------------------

PARAMETER_SUMMARY <- PARAMETER_RESULTS %>%
  group_by(
    baseline_design,
    n_processes,
    mean_length,
    beta_true
  ) %>%
  summarise(
    successful_replications =
      sum(
        converged,
        na.rm = TRUE
      ),
    mean_beta_hat =
      mean(
        beta_hat,
        na.rm = TRUE
      ),
    bias_beta =
      mean(
        beta_hat -
          beta_true,
        na.rm = TRUE
      ),
    rmse_beta =
      sqrt(
        mean(
          (
            beta_hat -
              beta_true
          )^2,
          na.rm = TRUE
        )
      ),
    coverage_beta =
      mean(
        beta_covered,
        na.rm = TRUE
      ),
    boundary_rate =
      mean(
        boundary_estimate,
        na.rm = TRUE
      ),
    mean_training_transitions =
      mean(
        training_transitions,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

ETA_SUMMARY <- ETA_RESULTS %>%
  group_by(
    baseline_design,
    n_processes,
    mean_length,
    beta_true,
    origin_stage
  ) %>%
  summarise(
    mean_eta_hat =
      mean(
        eta_hat,
        na.rm = TRUE
      ),
    eta_true =
      first(eta_true),
    bias_eta =
      mean(
        eta_hat -
          eta_true,
        na.rm = TRUE
      ),
    rmse_eta =
      sqrt(
        mean(
          (
            eta_hat -
              eta_true
          )^2,
          na.rm = TRUE
        )
      ),
    coverage_eta =
      mean(
        eta_covered,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

MODEL_SUMMARY <- MODEL_RESULTS %>%
  group_by(
    baseline_design,
    n_processes,
    mean_length,
    beta_true,
    model
  ) %>%
  summarise(
    successful_replications =
      sum(
        is.finite(
          departure_logloss
        )
      ),
    departure_logloss_mean =
      mean(
        departure_logloss,
        na.rm = TRUE
      ),
    departure_logloss_sd =
      safe_sd(
        departure_logloss
      ),
    departure_brier_mean =
      mean(
        departure_brier,
        na.rm = TRUE
      ),
    departure_brier_sd =
      safe_sd(
        departure_brier
      ),
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
      mean(
        full_logloss,
        na.rm = TRUE
      ),
    full_brier_mean =
      mean(
        full_brier,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

PRIMARY_REPLICATION <- MODEL_RESULTS %>%
  filter(
    model ==
      "Stage-Spell SPAOM"
  ) %>%
  select(
    scenario_id,
    replication_id,
    primary_departure_logloss =
      departure_logloss,
    primary_departure_brier =
      departure_brier,
    primary_full_logloss =
      full_logloss,
    primary_full_brier =
      full_brier
  )

PAIRED_REPLICATIONS <- MODEL_RESULTS %>%
  filter(
    model !=
      "Stage-Spell SPAOM"
  ) %>%
  inner_join(
    PRIMARY_REPLICATION,
    by = c(
      "scenario_id",
      "replication_id"
    )
  ) %>%
  mutate(
    delta_departure_logloss =
      departure_logloss -
      primary_departure_logloss,
    delta_departure_brier =
      departure_brier -
      primary_departure_brier,
    delta_full_logloss =
      full_logloss -
      primary_full_logloss,
    delta_full_brier =
      full_brier -
      primary_full_brier
  )

summarise_difference <- function(
  values
) {
  values <- values[
    is.finite(values)
  ]

  standard_error <- if (
    length(values) > 1L
  ) {
    sd(values) /
      sqrt(
        length(values)
      )
  } else {
    NA_real_
  }

  c(
    mean =
      safe_mean(values),
    se =
      standard_error,
    lower =
      safe_mean(values) -
      1.96 *
      standard_error,
    upper =
      safe_mean(values) +
      1.96 *
      standard_error,
    win_rate =
      mean(
        values > 0
      )
  )
}

PAIRED_SUMMARY_ROWS <- list()
paired_index <- 1L

group_variables <- c(
  "baseline_design",
  "n_processes",
  "mean_length",
  "beta_true",
  "model"
)

split_data <- split(
  PAIRED_REPLICATIONS,
  interaction(
    PAIRED_REPLICATIONS[
      group_variables
    ],
    drop = TRUE
  )
)

for (group_data in split_data) {
  if (nrow(group_data) == 0L) {
    next
  }

  for (
    metric_name in c(
      "delta_departure_logloss",
      "delta_departure_brier",
      "delta_full_logloss",
      "delta_full_brier"
    )
  ) {
    summary_values <-
      summarise_difference(
        group_data[[metric_name]]
      )

    PAIRED_SUMMARY_ROWS[[paired_index]] <- data.frame(
      baseline_design =
        group_data$
          baseline_design[1L],
      n_processes =
        group_data$
          n_processes[1L],
      mean_length =
        group_data$
          mean_length[1L],
      beta_true =
        group_data$
          beta_true[1L],
      comparator =
        group_data$
          model[1L],
      metric =
        sub(
          "^delta_",
          "",
          metric_name
        ),
      mean_advantage =
        unname(
          summary_values[
            "mean"
          ]
        ),
      standard_error =
        unname(
          summary_values[
            "se"
          ]
        ),
      ci_lower =
        unname(
          summary_values[
            "lower"
          ]
        ),
      ci_upper =
        unname(
          summary_values[
            "upper"
          ]
        ),
      win_rate =
        unname(
          summary_values[
            "win_rate"
          ]
        ),
      stringsAsFactors = FALSE
    )

    paired_index <-
      paired_index + 1L
  }
}

PAIRED_SUMMARY <- bind_rows(
  PAIRED_SUMMARY_ROWS
)

BEST_MODEL_REGIONS <- MODEL_SUMMARY %>%
  group_by(
    baseline_design,
    n_processes,
    mean_length,
    beta_true
  ) %>%
  slice_min(
    departure_logloss_mean,
    n = 1L,
    with_ties = TRUE
  ) %>%
  ungroup()

save_csv(
  MODEL_RESULTS,
  "02_replication_model_results.csv"
)

save_csv(
  PARAMETER_RESULTS,
  "03_replication_parameter_results.csv"
)

save_csv(
  ETA_RESULTS,
  "04_replication_eta_results.csv"
)

save_csv(
  PARAMETER_SUMMARY,
  "05_parameter_recovery_summary.csv"
)

save_csv(
  ETA_SUMMARY,
  "06_eta_recovery_summary.csv"
)

save_csv(
  MODEL_SUMMARY,
  "07_model_performance_summary.csv"
)

save_csv(
  PAIRED_REPLICATIONS,
  "08_paired_replication_advantages.csv"
)

save_csv(
  PAIRED_SUMMARY,
  "09_paired_advantage_summary.csv"
)

save_csv(
  BEST_MODEL_REGIONS,
  "10_best_model_regions.csv"
)

if (nrow(ERROR_RESULTS) > 0L) {
  save_csv(
    ERROR_RESULTS,
    "simulation_errors_combined.csv"
  )
}


# -----------------------------------------------------------------------------
# 12. MANDATORY GRAPHS
# -----------------------------------------------------------------------------

plot_beta_recovery <- PARAMETER_SUMMARY %>%
  ggplot(
    aes(
      x = beta_true,
      y = mean_beta_hat,
      group =
        factor(
          mean_length
        ),
      linetype =
        factor(
          mean_length
        )
    )
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = 2
  ) +
  geom_line() +
  geom_point() +
  facet_grid(
    baseline_design ~
      n_processes
  ) +
  labs(
    title =
      "Recovery of the stage-spell parameter",
    x =
      "True beta",
    y =
      "Mean estimated beta",
    linetype =
      "Mean process length"
  ) +
  theme_bw()

save_plot(
  plot_beta_recovery,
  "01_beta_recovery.png",
  width = 13,
  height = 8
)

plot_beta_rmse <- PARAMETER_SUMMARY %>%
  ggplot(
    aes(
      x =
        factor(
          n_processes
        ),
      y = rmse_beta,
      group =
        factor(
          mean_length
        ),
      linetype =
        factor(
          mean_length
        )
    )
  ) +
  geom_line() +
  geom_point() +
  facet_grid(
    baseline_design ~
      beta_true,
    labeller =
      label_both
  ) +
  labs(
    title =
      "RMSE of the stage-spell parameter",
    x =
      "Number of peace processes",
    y =
      "RMSE",
    linetype =
      "Mean process length"
  ) +
  theme_bw()

save_plot(
  plot_beta_rmse,
  "02_beta_RMSE.png",
  width = 13,
  height = 8
)

plot_beta_coverage <- PARAMETER_SUMMARY %>%
  ggplot(
    aes(
      x =
        factor(
          n_processes
        ),
      y =
        coverage_beta,
      group =
        factor(
          mean_length
        ),
      linetype =
        factor(
          mean_length
        )
    )
  ) +
  geom_hline(
    yintercept = 0.95,
    linetype = 2
  ) +
  geom_line() +
  geom_point() +
  facet_grid(
    baseline_design ~
      beta_true,
    labeller =
      label_both
  ) +
  labs(
    title =
      "Approximate confidence-interval coverage for beta",
    x =
      "Number of peace processes",
    y =
      "Coverage",
    linetype =
      "Mean process length"
  ) +
  theme_bw()

save_plot(
  plot_beta_coverage,
  "03_beta_coverage.png",
  width = 13,
  height = 8
)

heatmap_data <- PAIRED_SUMMARY %>%
  filter(
    comparator ==
      "Origin-specific homogeneous",
    metric ==
      "departure_logloss"
  )

plot_advantage_heatmap <- ggplot(
  heatmap_data,
  aes(
    x =
      factor(
        n_processes
      ),
    y =
      factor(
        mean_length
      ),
    fill =
      mean_advantage
  )
) +
  geom_tile() +
  geom_text(
    aes(
      label =
        sprintf(
          "%.3f",
          mean_advantage
        )
    )
  ) +
  facet_grid(
    baseline_design ~
      beta_true,
    labeller =
      label_both
  ) +
  labs(
    title =
      "Departure-log-loss advantage of SS-SPAOM over the homogeneous model",
    subtitle =
      "Positive values favour SS-SPAOM",
    x =
      "Number of peace processes",
    y =
      "Mean process length",
    fill =
      "Advantage"
  ) +
  theme_bw()

save_plot(
  plot_advantage_heatmap,
  "04_logloss_advantage_heatmap.png",
  width = 14,
  height = 8
)

plot_brier_heatmap <- PAIRED_SUMMARY %>%
  filter(
    comparator ==
      "Origin-specific homogeneous",
    metric ==
      "departure_brier"
  ) %>%
  ggplot(
    aes(
      x =
        factor(
          n_processes
        ),
      y =
        factor(
          mean_length
        ),
      fill =
        mean_advantage
    )
  ) +
  geom_tile() +
  geom_text(
    aes(
      label =
        sprintf(
          "%.3f",
          mean_advantage
        )
    )
  ) +
  facet_grid(
    baseline_design ~
      beta_true,
    labeller =
      label_both
  ) +
  labs(
    title =
      "Departure-Brier advantage of SS-SPAOM over the homogeneous model",
    subtitle =
      "Positive values favour SS-SPAOM",
    x =
      "Number of peace processes",
    y =
      "Mean process length",
    fill =
      "Advantage"
  ) +
  theme_bw()

save_plot(
  plot_brier_heatmap,
  "05_brier_advantage_heatmap.png",
  width = 14,
  height = 8
)

plot_model_comparison <- MODEL_SUMMARY %>%
  filter(
    beta_true %in%
      unique(BETA_VALUES),
    n_processes ==
      max(N_VALUES),
    mean_length ==
      max(
        MEAN_LENGTH_VALUES
      )
  ) %>%
  ggplot(
    aes(
      x = model,
      y =
        departure_logloss_mean
    )
  ) +
  geom_col() +
  coord_flip() +
  facet_grid(
    baseline_design ~
      beta_true,
    labeller =
      label_both
  ) +
  labs(
    title =
      "Departure log loss under the largest simulated design",
    x = "Model",
    y =
      "Mean test log loss"
  ) +
  theme_bw()

save_plot(
  plot_model_comparison,
  "06_model_comparison_largest_design.png",
  width = 14,
  height = 9
)

plot_calibration <- MODEL_SUMMARY %>%
  filter(
    model %in% c(
      "Origin-specific homogeneous",
      "Stage-Spell SPAOM",
      "Spline spell-age model"
    ),
    n_processes ==
      max(N_VALUES),
    mean_length ==
      max(
        MEAN_LENGTH_VALUES
      )
  ) %>%
  ggplot(
    aes(
      x =
        calibration_intercept_mean,
      y =
        calibration_slope_mean,
      shape = model
    )
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 2
  ) +
  geom_hline(
    yintercept = 1,
    linetype = 2
  ) +
  geom_point(
    size = 3
  ) +
  facet_grid(
    baseline_design ~
      beta_true,
    labeller =
      label_both
  ) +
  labs(
    title =
      "Simulation calibration at the largest design",
    x =
      "Mean calibration intercept",
    y =
      "Mean calibration slope",
    shape =
      "Model"
  ) +
  theme_bw()

save_plot(
  plot_calibration,
  "07_calibration_summary.png",
  width = 13,
  height = 8
)

pdf(
  file.path(
    OUTPUT_DIR,
    "SS_SPAOM_Simulation_Graphs.pdf"
  ),
  width = 11,
  height = 8.5
)

for (plot_object in list(
  plot_beta_recovery,
  plot_beta_rmse,
  plot_beta_coverage,
  plot_advantage_heatmap,
  plot_brier_heatmap,
  plot_model_comparison,
  plot_calibration
)) {
  print(plot_object)
}

dev.off()


# -----------------------------------------------------------------------------
# 13. EXCEL WORKBOOK, INTERPRETATION, AND STATUS
# -----------------------------------------------------------------------------

writexl::write_xlsx(
  list(
    Scenario_Grid =
      SCENARIO_GRID,
    Calibration =
      CALIBRATION_TABLE,
    Parameter_Recovery =
      PARAMETER_SUMMARY,
    Eta_Recovery =
      ETA_SUMMARY,
    Model_Performance =
      MODEL_SUMMARY,
    Paired_Advantages =
      PAIRED_SUMMARY,
    Best_Model_Regions =
      BEST_MODEL_REGIONS,
    Errors =
      ERROR_RESULTS
  ),
  path = file.path(
    OUTPUT_DIR,
    "SS_SPAOM_Revised_Simulation_Results.xlsx"
  )
)

homogeneous_advantage <- PAIRED_SUMMARY %>%
  filter(
    comparator ==
      "Origin-specific homogeneous",
    metric ==
      "departure_logloss"
  )

moderate_strong <- homogeneous_advantage %>%
  filter(
    beta_true >= 0.52
  )

null_design <- homogeneous_advantage %>%
  filter(
    beta_true == 0
  )

interpretation_lines <- c(
  "REVISED SS-SPAOM MONTE CARLO STUDY",
  "=================================",
  "",
  paste0(
    "Simulation mode: ",
    SIMULATION_MODE
  ),
  paste0(
    "Replications per scenario: ",
    B_REPLICATIONS
  ),
  paste0(
    "Completed scenarios: ",
    length(
      checkpoint_files
    ),
    " / ",
    nrow(
      SCENARIO_GRID
    )
  ),
  "",
  paste0(
    "Mean SS-SPAOM log-loss advantage under beta >= 0.52: ",
    safe_mean(
      moderate_strong$
        mean_advantage
    )
  ),
  paste0(
    "Mean SS-SPAOM log-loss advantage under beta = 0: ",
    safe_mean(
      null_design$
        mean_advantage
    )
  ),
  "",
  "Interpretation standard:",
  "- beta = 0 should favour or tie the homogeneous model;",
  "- moderate or strong beta should increasingly favour SS-SPAOM;",
  "- recovery should improve with more processes and longer histories;",
  "- the spline may be marginally better, but SS-SPAOM should remain competitive;",
  "- unequal eta designs test separation of stage persistence from spell dynamics."
)

writeLines(
  interpretation_lines,
  file.path(
    OUTPUT_DIR,
    "SIMULATION_INTERPRETATION.txt"
  )
)

status_lines <- c(
  paste0(
    "Simulation mode: ",
    SIMULATION_MODE
  ),
  paste0(
    "Scenarios completed: ",
    length(
      checkpoint_files
    ),
    " / ",
    nrow(
      SCENARIO_GRID
    )
  ),
  paste0(
    "Replications per scenario: ",
    B_REPLICATIONS
  ),
  paste0(
    "Model-result rows: ",
    nrow(
      MODEL_RESULTS
    )
  ),
  paste0(
    "Parameter-result rows: ",
    nrow(
      PARAMETER_RESULTS
    )
  ),
  paste0(
    "Error rows: ",
    nrow(
      ERROR_RESULTS
    )
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

capture.output(
  sessionInfo(),
  file = file.path(
    OUTPUT_DIR,
    "R_session_info.txt"
  )
)

cat("\n============================================================\n")
cat("REVISED SS-SPAOM SIMULATION COMPLETED\n")
cat("============================================================\n")
cat("Results are saved in:\n")
cat(OUTPUT_DIR, "\n\n")
cat("Open these files first:\n")
cat(
  file.path(
    OUTPUT_DIR,
    "SIMULATION_INTERPRETATION.txt"
  ),
  "\n"
)
cat(
  file.path(
    OUTPUT_DIR,
    "SS_SPAOM_Revised_Simulation_Results.xlsx"
  ),
  "\n"
)
cat(
  file.path(
    OUTPUT_DIR,
    "SS_SPAOM_Simulation_Graphs.pdf"
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
