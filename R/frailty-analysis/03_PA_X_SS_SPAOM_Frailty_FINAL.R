# =============================================================================
# FRAILTY SENSITIVITY ANALYSIS FOR THE STAGE-SPELL SPAOM
# =============================================================================
#
# Frozen primary fixed-effects model:
#
#   log(mu_iar) = log(eta_a) - beta_s log(s_ir)
#
# Sensitivity extension:
#
#   u_i ~ N(0, sigma_u^2)
#
#   log(mu_iar)
#     = log(eta_a) + u_i - beta_s log(s_ir).
#
# This program asks whether the stage-spell effect remains after accounting
# for unobserved peace-process heterogeneity.  The fixed SS-SPAOM remains the
# primary model; the frailty model is explicitly secondary.
#
# In grouped process-level validation, held-out processes have no estimated
# random intercept.  Their predictions are therefore marginalised over the
# fitted normal frailty distribution using Gaussian quadrature.
#
# =============================================================================


# -----------------------------------------------------------------------------
# 0. USER SETTINGS
# -----------------------------------------------------------------------------

BASE_DIR <- "C:/Users/Hp/Downloads"

TRANSITION_FILE <- NA_character_

ANALYSIS_MODE <- "FINAL"  # "QUICK" or "FINAL"

OUTPUT_DIR <- file.path(
  BASE_DIR,
  paste0(
    "PA_X_SS_SPAOM_Frailty_",
    ANALYSIS_MODE
  )
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
  VARIANCE_BOOTSTRAP_B <- 49L
  PROCESS_BOOTSTRAP_B <- 99L
  CV_REPEATS <- 3L
  CV_FOLDS <- 5L
  PAIRED_BOOTSTRAP_B <- 500L
} else {
  VARIANCE_BOOTSTRAP_B <- 499L
  PROCESS_BOOTSTRAP_B <- 499L
  CV_REPEATS <- 10L
  CV_FOLDS <- 5L
  PAIRED_BOOTSTRAP_B <- 5000L
}

BASE_SEED <- 20260714L
PROBABILITY_FLOOR <- 1e-12
QUADRATURE_POINTS <- 25L
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
  "glmmTMB",
  "statmod"
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
  library(glmmTMB)
})


# -----------------------------------------------------------------------------
# 2. GENERAL HELPERS
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
    "frailty_errors.csv"
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


# -----------------------------------------------------------------------------
# 3. LOCATE AND PREPARE DATA
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
        log(
          spell_age
        )
    ) %>%
    select(
      -origin_character,
      -spell_start
    )
}

DATA_FILE <- locate_transition_file()

RAW_DATA <- read.csv(
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
      RAW_DATA$origin_stage
    ),
    as.character(
      RAW_DATA$destination_stage
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

DATA <- add_spell_age(
  RAW_DATA,
  STAGE_LEVELS
)

cat(
  "Transitions: ",
  nrow(DATA),
  "\nPeace processes: ",
  dplyr::n_distinct(
    DATA$PP
  ),
  "\n\n",
  sep = ""
)


# -----------------------------------------------------------------------------
# 4. MODEL FITTING
# -----------------------------------------------------------------------------

GLMM_CONTROL <- glmmTMBControl(
  optCtrl = list(
    iter.max = 2000,
    eval.max = 2000
  ),
  optimizer = optim,
  optArgs = list(
    method = "BFGS"
  )
)

fit_tmb_model <- function(
  formula_object,
  data
) {
  tryCatch(
    suppressWarnings(
      glmmTMB(
        formula_object,
        data = data,
        family =
          binomial(
            link =
              "cloglog"
          ),
        control =
          GLMM_CONTROL
      )
    ),
    error = function(e) NULL
  )
}

model_converged <- function(fit) {
  !is.null(fit) &&
    isTRUE(
      fit$sdr$pdHess
    ) &&
    is.finite(
      as.numeric(
        logLik(fit)
      )
    )
}

get_beta_spell <- function(fit) {
  if (is.null(fit)) {
    return(NA_real_)
  }

  coefficients <- fixef(fit)$cond

  if (
    !(
      "log_spell_age" %in%
        names(coefficients)
    )
  ) {
    return(NA_real_)
  }

  -unname(
    coefficients[
      "log_spell_age"
    ]
  )
}

get_beta_se <- function(fit) {
  if (is.null(fit)) {
    return(NA_real_)
  }

  covariance <- tryCatch(
    vcov(fit)$cond,
    error = function(e) NULL
  )

  if (
    is.null(covariance) ||
    !(
      "log_spell_age" %in%
        rownames(covariance)
    )
  ) {
    return(NA_real_)
  }

  sqrt(
    covariance[
      "log_spell_age",
      "log_spell_age"
    ]
  )
}

get_frailty_sd <- function(fit) {
  if (is.null(fit)) {
    return(0)
  }

  variance_component <- tryCatch(
    VarCorr(fit)$cond$PP,
    error = function(e) NULL
  )

  if (is.null(variance_component)) {
    return(0)
  }

  standard_deviation <- attr(
    variance_component,
    "stddev"
  )

  if (
    is.null(standard_deviation) ||
    length(standard_deviation) == 0L
  ) {
    return(0)
  }

  as.numeric(
    standard_deviation[1L]
  )
}

fit_all_models <- function(data) {
  list(
    fixed_homogeneous =
      fit_tmb_model(
        D ~ 0 +
          origin_stage,
        data
      ),
    fixed_stage_spell =
      fit_tmb_model(
        D ~ 0 +
          origin_stage +
          log_spell_age,
        data
      ),
    frailty_homogeneous =
      fit_tmb_model(
        D ~ 0 +
          origin_stage +
          (
            1 |
              PP
          ),
        data
      ),
    frailty_stage_spell =
      fit_tmb_model(
        D ~ 0 +
          origin_stage +
          log_spell_age +
          (
            1 |
              PP
          ),
        data
      ),
    frailty_spline_spell =
      fit_tmb_model(
        as.formula(
          paste0(
            "D ~ 0 + origin_stage + ",
            "splines::ns(log_spell_age, df = ",
            SPLINE_DF,
            ") + (1 | PP)"
          )
        ),
        data
      )
  )
}

cat(
  "Fitting fixed and frailty models...\n"
)

FULL_FITS <- fit_all_models(
  DATA
)

if (
  !model_converged(
    FULL_FITS$
      fixed_stage_spell
  )
) {
  stop(
    "The fixed Stage-Spell model did not converge."
  )
}

if (
  !model_converged(
    FULL_FITS$
      frailty_stage_spell
  )
) {
  warning(
    "The primary frailty Stage-Spell model has convergence concerns."
  )
}


# -----------------------------------------------------------------------------
# 5. FULL-DATA ESTIMATES AND MODEL COMPARISON
# -----------------------------------------------------------------------------

MODEL_NAMES <- c(
  fixed_homogeneous =
    "Fixed homogeneous",
  fixed_stage_spell =
    "Fixed Stage-Spell SPAOM",
  frailty_homogeneous =
    "Frailty homogeneous",
  frailty_stage_spell =
    "Frailty Stage-Spell SPAOM",
  frailty_spline_spell =
    "Frailty spline spell-age"
)

model_comparison_rows <- lapply(
  names(FULL_FITS),
  function(model_key) {
    fit <- FULL_FITS[[model_key]]

    loglikelihood <- if (
      is.null(fit)
    ) {
      NA_real_
    } else {
      as.numeric(
        logLik(fit)
      )
    }

    parameter_count <- if (
      is.null(fit)
    ) {
      NA_real_
    } else {
      attr(
        logLik(fit),
        "df"
      )
    }

    data.frame(
      model_key =
        model_key,
      model =
        MODEL_NAMES[
          model_key
        ],
      converged =
        model_converged(fit),
      loglik =
        loglikelihood,
      parameters =
        parameter_count,
      AIC = if (
        is.finite(loglikelihood)
      ) {
        -2 *
          loglikelihood +
          2 *
          parameter_count
      } else {
        NA_real_
      },
      BIC_cluster = if (
        is.finite(loglikelihood)
      ) {
        -2 *
          loglikelihood +
          parameter_count *
          log(
            dplyr::n_distinct(
              DATA$PP
            )
          )
      } else {
        NA_real_
      },
      beta_spell_hat =
        get_beta_spell(fit),
      beta_spell_se =
        get_beta_se(fit),
      frailty_sd =
        get_frailty_sd(fit),
      stringsAsFactors = FALSE
    )
  }
)

MODEL_COMPARISON <- bind_rows(
  model_comparison_rows
) %>%
  mutate(
    beta_ci_lower =
      beta_spell_hat -
      1.96 *
      beta_spell_se,
    beta_ci_upper =
      beta_spell_hat +
      1.96 *
      beta_spell_se,
    delta_AIC =
      AIC -
      min(
        AIC,
        na.rm = TRUE
      ),
    delta_BIC_cluster =
      BIC_cluster -
      min(
        BIC_cluster,
        na.rm = TRUE
      )
  ) %>%
  arrange(AIC)

FIXED_BETA <- get_beta_spell(
  FULL_FITS$
    fixed_stage_spell
)

FRAILTY_BETA <- get_beta_spell(
  FULL_FITS$
    frailty_stage_spell
)

FRAILTY_SD <- get_frailty_sd(
  FULL_FITS$
    frailty_stage_spell
)

ATTENUATION_RATIO <- FRAILTY_BETA /
  FIXED_BETA

save_csv(
  MODEL_COMPARISON,
  "01_frailty_model_comparison.csv"
)

PRIMARY_SENSITIVITY_ESTIMATES <- data.frame(
  fixed_beta_spell =
    FIXED_BETA,
  frailty_beta_spell =
    FRAILTY_BETA,
  frailty_standard_deviation =
    FRAILTY_SD,
  attenuation_ratio =
    ATTENUATION_RATIO,
  relative_intensity_s5_fixed =
    5^(
      -FIXED_BETA
    ),
  relative_intensity_s5_frailty =
    5^(
      -FRAILTY_BETA
    ),
  stringsAsFactors = FALSE
)

save_csv(
  PRIMARY_SENSITIVITY_ESTIMATES,
  "02_primary_frailty_sensitivity_estimates.csv"
)


# -----------------------------------------------------------------------------
# 6. MARGINAL PREDICTIONS FOR NEW PEACE PROCESSES
# -----------------------------------------------------------------------------

QUADRATURE <- statmod::gauss.quad.prob(
  QUADRATURE_POINTS,
  dist = "normal"
)

predict_fixed_response <- function(
  fit,
  new_data
) {
  prediction <- tryCatch(
    predict(
      fit,
      newdata = new_data,
      type = "response",
      allow.new.levels = TRUE
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

predict_frailty_marginal <- function(
  fit,
  new_data
) {
  fixed_linear_predictor <- tryCatch(
    predict(
      fit,
      newdata = new_data,
      type = "link",
      re.form = NA,
      allow.new.levels = TRUE
    ),
    error = function(e) {
      rep(
        NA_real_,
        nrow(new_data)
      )
    }
  )

  if (
    any(
      !is.finite(
        fixed_linear_predictor
      )
    )
  ) {
    return(
      rep(
        NA_real_,
        nrow(new_data)
      )
    )
  }

  random_sd <- get_frailty_sd(
    fit
  )

  if (
    !is.finite(random_sd) ||
    random_sd <= 1e-10
  ) {
    return(
      inverse_cloglog(
        fixed_linear_predictor
      )
    )
  }

  marginal_probability <- rep(
    0,
    length(
      fixed_linear_predictor
    )
  )

  for (
    quadrature_index in
    seq_along(
      QUADRATURE$nodes
    )
  ) {
    marginal_probability <-
      marginal_probability +
      QUADRATURE$weights[
        quadrature_index
      ] *
      inverse_cloglog(
        fixed_linear_predictor +
          random_sd *
          QUADRATURE$nodes[
            quadrature_index
          ]
      )
  }

  clip_probability(
    marginal_probability
  )
}


# -----------------------------------------------------------------------------
# 7. PARAMETRIC-BOOTSTRAP TEST OF THE FRAILTY VARIANCE
# -----------------------------------------------------------------------------

VARIANCE_BOOTSTRAP_FILE <- file.path(
  CHECKPOINT_DIR,
  paste0(
    "frailty_variance_bootstrap_",
    ANALYSIS_MODE,
    ".csv"
  )
)

run_variance_bootstrap <- function() {
  if (
    file.exists(
      VARIANCE_BOOTSTRAP_FILE
    )
  ) {
    results <- read.csv(
      VARIANCE_BOOTSTRAP_FILE,
      stringsAsFactors = FALSE
    )
  } else {
    results <- data.frame(
      replication_id =
        integer(0),
      LR =
        numeric(0),
      frailty_sd =
        numeric(0),
      converged =
        logical(0),
      stringsAsFactors = FALSE
    )
  }

  completed <-
    results$replication_id

  null_probability <-
    predict_fixed_response(
      FULL_FITS$
        fixed_stage_spell,
      DATA
    )

  observed_LR <- max(
    0,
    2 * (
      as.numeric(
        logLik(
          FULL_FITS$
            frailty_stage_spell
        )
      ) -
      as.numeric(
        logLik(
          FULL_FITS$
            fixed_stage_spell
        )
      )
    )
  )

  set.seed(
    BASE_SEED + 10000L
  )

  bootstrap_seeds <- sample.int(
    .Machine$integer.max,
    VARIANCE_BOOTSTRAP_B
  )

  for (
    replication_id in
    seq_len(
      VARIANCE_BOOTSTRAP_B
    )
  ) {
    if (
      replication_id %in%
        completed
    ) {
      next
    }

    cat(
      "Frailty null bootstrap ",
      replication_id,
      " of ",
      VARIANCE_BOOTSTRAP_B,
      "\n",
      sep = ""
    )

    set.seed(
      bootstrap_seeds[
        replication_id
      ]
    )

    bootstrap_data <- DATA

    bootstrap_data$D <- rbinom(
      nrow(bootstrap_data),
      size = 1L,
      prob =
        null_probability
    )

    null_fit <- fit_tmb_model(
      D ~ 0 +
        origin_stage +
        log_spell_age,
      bootstrap_data
    )

    frailty_fit <- fit_tmb_model(
      D ~ 0 +
        origin_stage +
        log_spell_age +
        (
          1 |
            PP
        ),
      bootstrap_data
    )

    converged <-
      model_converged(
        null_fit
      ) &&
      model_converged(
        frailty_fit
      )

    LR <- if (converged) {
      max(
        0,
        2 * (
          as.numeric(
            logLik(
              frailty_fit
            )
          ) -
          as.numeric(
            logLik(
              null_fit
            )
          )
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
        frailty_sd = if (
          converged
        ) {
          get_frailty_sd(
            frailty_fit
          )
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
      VARIANCE_BOOTSTRAP_FILE,
      row.names = FALSE
    )
  }

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
        VARIANCE_BOOTSTRAP_B,
      stringsAsFactors = FALSE
    )
  )
}

VARIANCE_BOOTSTRAP_RESULT <-
  run_variance_bootstrap()

save_csv(
  VARIANCE_BOOTSTRAP_RESULT$
    replications,
  "03_frailty_variance_bootstrap_replications.csv"
)

save_csv(
  VARIANCE_BOOTSTRAP_RESULT$
    summary,
  "04_frailty_variance_bootstrap_summary.csv"
)


# -----------------------------------------------------------------------------
# 8. PROCESS-BOOTSTRAP UNCERTAINTY FOR BETA AND SIGMA
# -----------------------------------------------------------------------------

PROCESS_BOOTSTRAP_FILE <- file.path(
  CHECKPOINT_DIR,
  paste0(
    "frailty_process_bootstrap_",
    ANALYSIS_MODE,
    ".csv"
  )
)

run_process_bootstrap <- function() {
  if (
    file.exists(
      PROCESS_BOOTSTRAP_FILE
    )
  ) {
    results <- read.csv(
      PROCESS_BOOTSTRAP_FILE,
      stringsAsFactors = FALSE
    )
  } else {
    results <- data.frame(
      replication_id =
        integer(0),
      beta_spell_hat =
        numeric(0),
      frailty_sd =
        numeric(0),
      converged =
        logical(0),
      stringsAsFactors = FALSE
    )
  }

  completed <-
    results$replication_id

  process_ids <- unique(
    as.character(
      DATA$PP
    )
  )

  set.seed(
    BASE_SEED + 20000L
  )

  bootstrap_seeds <- sample.int(
    .Machine$integer.max,
    PROCESS_BOOTSTRAP_B
  )

  for (
    replication_id in
    seq_len(
      PROCESS_BOOTSTRAP_B
    )
  ) {
    if (
      replication_id %in%
        completed
    ) {
      next
    }

    cat(
      "Frailty process bootstrap ",
      replication_id,
      " of ",
      PROCESS_BOOTSTRAP_B,
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
      size =
        length(process_ids),
      replace = TRUE
    )

    parts <- vector(
      "list",
      length(
        sampled_processes
      )
    )

    for (
      sample_index in
      seq_along(
        sampled_processes
      )
    ) {
      part <- DATA %>%
        filter(
          PP ==
            sampled_processes[
              sample_index
            ]
        )

      part$PP <- paste0(
        "bootstrap_",
        replication_id,
        "_",
        sample_index
      )

      parts[[sample_index]] <-
        part
    }

    bootstrap_data <- bind_rows(
      parts
    )

    fit <- fit_tmb_model(
      D ~ 0 +
        origin_stage +
        log_spell_age +
        (
          1 |
            PP
        ),
      bootstrap_data
    )

    converged <-
      model_converged(fit)

    results <- bind_rows(
      results,
      data.frame(
        replication_id =
          replication_id,
        beta_spell_hat = if (
          converged
        ) {
          get_beta_spell(fit)
        } else {
          NA_real_
        },
        frailty_sd = if (
          converged
        ) {
          get_frailty_sd(fit)
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
      PROCESS_BOOTSTRAP_FILE,
      row.names = FALSE
    )
  }

  results
}

PROCESS_BOOTSTRAP_RESULTS <-
  run_process_bootstrap()

summarise_bootstrap <- function(
  values,
  estimate,
  parameter_name
) {
  finite_values <- values[
    is.finite(values)
  ]

  interval <- if (
    length(finite_values) > 0L
  ) {
    quantile(
      finite_values,
      c(
        0.025,
        0.975
      )
    )
  } else {
    c(
      NA_real_,
      NA_real_
    )
  }

  data.frame(
    parameter =
      parameter_name,
    estimate =
      estimate,
    bootstrap_mean =
      safe_mean(
        finite_values
      ),
    bootstrap_standard_error =
      safe_sd(
        finite_values
      ),
    ci_lower =
      unname(
        interval[1L]
      ),
    ci_upper =
      unname(
        interval[2L]
      ),
    successful_replications =
      length(
        finite_values
      ),
    stringsAsFactors = FALSE
  )
}

PROCESS_BOOTSTRAP_SUMMARY <- bind_rows(
  summarise_bootstrap(
    PROCESS_BOOTSTRAP_RESULTS$
      beta_spell_hat,
    FRAILTY_BETA,
    "beta_spell"
  ),
  summarise_bootstrap(
    PROCESS_BOOTSTRAP_RESULTS$
      frailty_sd,
    FRAILTY_SD,
    "frailty_sd"
  )
)

save_csv(
  PROCESS_BOOTSTRAP_RESULTS,
  "05_frailty_process_bootstrap_replications.csv"
)

save_csv(
  PROCESS_BOOTSTRAP_SUMMARY,
  "06_frailty_process_bootstrap_summary.csv"
)


# -----------------------------------------------------------------------------
# 9. REPEATED GROUPED PROCESS-LEVEL CROSS-VALIDATION
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
      name =
        "transitions"
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
      nrow(
        process_sizes
      )
    )
  ) {
    candidates <- which(
      fold_load ==
        min(
          fold_load
        )
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

calibration_intercept_slope <- function(
  actual,
  predicted
) {
  fit <- tryCatch(
    suppressWarnings(
      glm(
        actual ~
          qlogis(
            clip_probability(
              predicted
            )
          ),
        family =
          binomial()
      )
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(
      c(
        intercept =
          NA_real_,
        slope =
          NA_real_
      )
    )
  }

  coefficients <- coef(fit)

  c(
    intercept =
      unname(
        coefficients[1L]
      ),
    slope =
      unname(
        coefficients[2L]
      )
  )
}

evaluate_departure_prediction <- function(
  data,
  probability
) {
  probability <- clip_probability(
    probability
  )

  calibration <-
    calibration_intercept_slope(
      data$D,
      probability
    )

  contributions <- data.frame(
    PP =
      as.character(
        data$PP
      ),
    transition_id =
      as.character(
        data$transition_id
      ),
    D =
      data$D,
    q =
      probability,
    departure_logloss =
      -(
        data$D *
          log(probability) +
        (
          1 -
            data$D
        ) *
          log(
            1 -
              probability
          )
      ),
    departure_brier =
      (
        data$D -
          probability
      )^2,
    stringsAsFactors = FALSE
  )

  list(
    metrics = c(
      departure_logloss =
        mean(
          contributions$
            departure_logloss
        ),
      departure_brier =
        mean(
          contributions$
            departure_brier
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
        )
    ),
    contributions =
      contributions
  )
}

CV_DIRECTORY <- file.path(
  CHECKPOINT_DIR,
  paste0(
    "CV_",
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
      DATA,
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
        "Frailty CV repeat ",
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

      training_data <- DATA %>%
        filter(
          !(PP %in%
              test_processes)
        )

      test_data <- DATA %>%
        filter(
          PP %in%
            test_processes
        )

      fits <- fit_all_models(
        training_data
      )

      prediction_list <- list(
        `Fixed homogeneous` =
          predict_fixed_response(
            fits$
              fixed_homogeneous,
            test_data
          ),
        `Fixed Stage-Spell SPAOM` =
          predict_fixed_response(
            fits$
              fixed_stage_spell,
            test_data
          ),
        `Frailty homogeneous` =
          predict_frailty_marginal(
            fits$
              frailty_homogeneous,
            test_data
          ),
        `Frailty Stage-Spell SPAOM` =
          predict_frailty_marginal(
            fits$
              frailty_stage_spell,
            test_data
          ),
        `Frailty spline spell-age` =
          predict_frailty_marginal(
            fits$
              frailty_spline_spell,
            test_data
          )
      )

      metric_rows <- list()
      contribution_rows <- list()
      result_index <- 1L

      for (
        model_name in
        names(
          prediction_list
        )
      ) {
        evaluation <-
          evaluate_departure_prediction(
            test_data,
            prediction_list[[model_name]]
          )

        metric_rows[[result_index]] <-
          data.frame(
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
            stringsAsFactors = FALSE
          )

        contribution <-
          evaluation$
            contributions

        contribution$repeat_id <-
          repeat_id

        contribution$fold_id <-
          fold_id

        contribution$model <-
          model_name

        contribution_rows[[result_index]] <-
          contribution

        result_index <-
          result_index + 1L
      }

      write.csv(
        bind_rows(
          metric_rows
        ),
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

  existing_metric_files <-
    metric_files[
      file.exists(
        metric_files
      )
    ]

  existing_contribution_files <-
    contribution_files[
      file.exists(
        contribution_files
      )
    ]

  if (
    length(
      existing_metric_files
    ) == 0L ||
    length(
      existing_contribution_files
    ) == 0L
  ) {
    stop(
      "No frailty cross-validation files were produced."
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
      mean(
        departure_logloss
      ),
    departure_brier =
      mean(
        departure_brier
      ),
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
      mean(
        departure_logloss
      ),
    departure_logloss_sd =
      safe_sd(
        departure_logloss
      ),
    departure_brier_mean =
      mean(
        departure_brier
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
    .groups = "drop"
  ) %>%
  arrange(
    departure_logloss_mean
  )

save_csv(
  CV_RESULT$fold_metrics,
  "07_frailty_cv_fold_metrics.csv"
)

save_csv(
  CV_REPEAT_METRICS,
  "08_frailty_cv_repeat_metrics.csv"
)

save_csv(
  CV_MODEL_SUMMARY,
  "09_frailty_cv_model_summary.csv"
)


# -----------------------------------------------------------------------------
# 10. PAIRED PROCESS-BOOTSTRAP ADVANTAGES
# -----------------------------------------------------------------------------

PRIMARY_MODEL_NAME <-
  "Frailty Stage-Spell SPAOM"

comparison_models <- setdiff(
  unique(
    CV_RESULT$
      contributions$model
  ),
  PRIMARY_MODEL_NAME
)

calculate_process_bootstrap_interval <- function(
  comparison_model,
  metric_name
) {
  primary <- CV_RESULT$
    contributions %>%
    filter(
      model ==
        PRIMARY_MODEL_NAME
    ) %>%
    select(
      repeat_id,
      transition_id,
      PP,
      primary_value =
        all_of(
          metric_name
        )
    )

  comparison <- CV_RESULT$
    contributions %>%
    filter(
      model ==
        comparison_model
    ) %>%
    select(
      repeat_id,
      transition_id,
      comparison_value =
        all_of(
          metric_name
        )
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

  observed <- sum(
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
      selected <- sample(
        seq_len(
          nrow(paired)
        ),
        size =
          nrow(paired),
        replace = TRUE
      )

      sum(
        paired$advantage_sum[
          selected
        ]
      ) /
        sum(
          paired$
            observation_count[
              selected
            ]
        )
    }
  )

  interval <- quantile(
    bootstrap_values,
    c(
      0.025,
      0.975
    ),
    na.rm = TRUE
  )

  c(
    advantage =
      observed,
    lower =
      unname(
        interval[1L]
      ),
    upper =
      unname(
        interval[2L]
      )
  )
}

paired_rows <- list()
paired_index <- 1L

for (
  comparison_model in
  comparison_models
) {
  for (
    metric_name in c(
      "departure_logloss",
      "departure_brier"
    )
  ) {
    interval <-
      calculate_process_bootstrap_interval(
        comparison_model,
        metric_name
      )

    repeat_primary <- CV_REPEAT_METRICS %>%
      filter(
        model ==
          PRIMARY_MODEL_NAME
      ) %>%
      select(
        repeat_id,
        primary_value =
          all_of(
            metric_name
          )
      )

    repeat_comparison <- CV_REPEAT_METRICS %>%
      filter(
        model ==
          comparison_model
      ) %>%
      select(
        repeat_id,
        comparison_value =
          all_of(
            metric_name
          )
      )

    repeat_difference <- inner_join(
      repeat_primary,
      repeat_comparison,
      by =
        "repeat_id"
    ) %>%
      mutate(
        difference =
          comparison_value -
          primary_value
      )

    paired_rows[[paired_index]] <-
      data.frame(
        comparison_model =
          comparison_model,
        metric =
          metric_name,
        mean_repeat_advantage =
          mean(
            repeat_difference$
              difference
          ),
        repeat_win_proportion =
          mean(
            repeat_difference$
              difference > 0
          ),
        process_bootstrap_advantage =
          unname(
            interval[
              "advantage"
            ]
          ),
        process_bootstrap_lower =
          unname(
            interval[
              "lower"
            ]
          ),
        process_bootstrap_upper =
          unname(
            interval[
              "upper"
            ]
          ),
        stringsAsFactors = FALSE
      )

    paired_index <-
      paired_index + 1L
  }
}

PAIRED_ADVANTAGES <- bind_rows(
  paired_rows
)

save_csv(
  PAIRED_ADVANTAGES,
  "10_frailty_paired_predictive_advantages.csv"
)


# -----------------------------------------------------------------------------
# 11. CALIBRATION TABLE
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
      transitions =
        n(),
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
    D =
      first(D),
    q =
      mean(q),
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
  "11_frailty_calibration_table.csv"
)


# -----------------------------------------------------------------------------
# 12. RANDOM-EFFECT ESTIMATES
# -----------------------------------------------------------------------------

random_effects <- tryCatch(
  ranef(
    FULL_FITS$
      frailty_stage_spell
  )$cond$PP,
  error = function(e) NULL
)

RANDOM_EFFECT_TABLE <- if (
  is.null(random_effects)
) {
  data.frame()
} else {
  data.frame(
    PP =
      rownames(
        random_effects
      ),
    random_intercept =
      as.numeric(
        random_effects[
          ,
          1L
        ]
      ),
    stringsAsFactors = FALSE
  )
}

save_csv(
  RANDOM_EFFECT_TABLE,
  "12_process_random_effects.csv"
)


# -----------------------------------------------------------------------------
# 13. AUTOMATIC SENSITIVITY INTERPRETATION
# -----------------------------------------------------------------------------

beta_bootstrap <- PROCESS_BOOTSTRAP_SUMMARY %>%
  filter(
    parameter ==
      "beta_spell"
  )

variance_p_value <- VARIANCE_BOOTSTRAP_RESULT$
  summary$
  bootstrap_p_value[1L]

fixed_comparison <- PAIRED_ADVANTAGES %>%
  filter(
    comparison_model ==
      "Fixed Stage-Spell SPAOM",
    metric ==
      "departure_logloss"
  )

beta_persists <-
  nrow(beta_bootstrap) > 0L &&
  beta_bootstrap$
    ci_lower[1L] > 0

substantial_attenuation <-
  is.finite(
    ATTENUATION_RATIO
  ) &&
  ATTENUATION_RATIO < 0.70

frailty_supported <-
  is.finite(
    variance_p_value
  ) &&
  variance_p_value < 0.05

frailty_predictively_better <-
  nrow(fixed_comparison) > 0L &&
  fixed_comparison$
    process_bootstrap_lower[1L] > 0

interpretation <- if (
  beta_persists &&
  frailty_supported &&
  frailty_predictively_better
) {
  "Frailty is supported and improves prediction, while the stage-spell effect remains positive. Process heterogeneity supplements rather than replaces stage entrenchment."
} else if (
  beta_persists &&
  frailty_supported
) {
  "Frailty is supported, but grouped prediction is not clearly improved. The stage-spell effect remains after accounting for process heterogeneity."
} else if (
  beta_persists
) {
  "The stage-spell effect remains positive after the frailty extension, but evidence for additional random-intercept heterogeneity is limited."
} else {
  "The stage-spell effect is substantially weakened after adding frailty; part of the apparent duration dependence may reflect unobserved process heterogeneity."
}

SENSITIVITY_VERDICT <- data.frame(
  criterion = c(
    "Frailty bootstrap p-value below 0.05",
    "Frailty beta bootstrap interval above zero",
    "Frailty model improves grouped log loss",
    "Frailty beta is less than 70% of fixed beta"
  ),
  satisfied = c(
    frailty_supported,
    beta_persists,
    frailty_predictively_better,
    substantial_attenuation
  ),
  stringsAsFactors = FALSE
)

save_csv(
  SENSITIVITY_VERDICT,
  "13_frailty_sensitivity_verdict.csv"
)

writeLines(
  c(
    "FRAILTY SENSITIVITY INTERPRETATION",
    "=================================",
    "",
    interpretation,
    "",
    paste0(
      "Fixed beta_spell: ",
      FIXED_BETA
    ),
    paste0(
      "Frailty beta_spell: ",
      FRAILTY_BETA
    ),
    paste0(
      "Attenuation ratio: ",
      ATTENUATION_RATIO
    ),
    paste0(
      "Frailty standard deviation: ",
      FRAILTY_SD
    ),
    paste0(
      "Frailty variance bootstrap p-value: ",
      variance_p_value
    ),
    "",
    "The fixed SS-SPAOM remains the frozen primary model. The random-intercept model is reported only as a sensitivity analysis."
  ),
  file.path(
    OUTPUT_DIR,
    "FRAILTY_INTERPRETATION.txt"
  )
)


# -----------------------------------------------------------------------------
# 14. GRAPHS
# -----------------------------------------------------------------------------

beta_plot_data <- MODEL_COMPARISON %>%
  filter(
    model %in% c(
      "Fixed Stage-Spell SPAOM",
      "Frailty Stage-Spell SPAOM"
    )
  )

plot_beta_comparison <- ggplot(
  beta_plot_data,
  aes(
    x = model,
    y =
      beta_spell_hat
  )
) +
  geom_point(
    size = 3
  ) +
  geom_errorbar(
    aes(
      ymin =
        beta_ci_lower,
      ymax =
        beta_ci_upper
    ),
    width = 0.15
  ) +
  coord_flip() +
  labs(
    title =
      "Stage-spell effect before and after process frailty",
    x = "Model",
    y = expression(
      hat(beta)[s]
    )
  ) +
  theme_bw()

save_plot(
  plot_beta_comparison,
  "01_fixed_vs_frailty_beta.png"
)

plot_random_effects <- if (
  nrow(
    RANDOM_EFFECT_TABLE
  ) > 0L
) {
  ggplot(
    RANDOM_EFFECT_TABLE,
    aes(
      x =
        random_intercept
    )
  ) +
    geom_histogram(
      bins = 30
    ) +
    geom_vline(
      xintercept = 0,
      linetype = 2
    ) +
    labs(
      title =
        "Estimated peace-process random intercepts",
      x =
        "Conditional random-intercept estimate",
      y =
        "Frequency"
    ) +
    theme_bw()
} else {
  ggplot() +
    labs(
      title =
        "Random-effect estimates were unavailable"
    ) +
    theme_void()
}

save_plot(
  plot_random_effects,
  "02_process_random_effect_distribution.png"
)

bootstrap_long <- PROCESS_BOOTSTRAP_RESULTS %>%
  filter(
    converged
  ) %>%
  select(
    beta_spell_hat,
    frailty_sd
  ) %>%
  pivot_longer(
    cols =
      everything(),
    names_to =
      "parameter",
    values_to =
      "estimate"
  )

plot_parameter_bootstrap <- ggplot(
  bootstrap_long,
  aes(
    x = estimate
  )
) +
  geom_histogram(
    bins = 30
  ) +
  facet_wrap(
    ~ parameter,
    scales = "free"
  ) +
  labs(
    title =
      "Process-bootstrap distributions for the frailty sensitivity model",
    x =
      "Bootstrap estimate",
    y =
      "Frequency"
  ) +
  theme_bw()

save_plot(
  plot_parameter_bootstrap,
  "03_frailty_parameter_bootstrap.png",
  width = 11,
  height = 6
)

plot_variance_bootstrap <- VARIANCE_BOOTSTRAP_RESULT$
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
      VARIANCE_BOOTSTRAP_RESULT$
        summary$
        observed_LR[1L],
    linetype = 2
  ) +
  labs(
    title =
      "Parametric-bootstrap test of the process-frailty variance",
    subtitle =
      "Dashed line marks the observed likelihood-ratio statistic",
    x =
      "Likelihood-ratio statistic",
    y =
      "Frequency"
  ) +
  theme_bw()

save_plot(
  plot_variance_bootstrap,
  "04_frailty_variance_bootstrap.png"
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
      "Grouped cross-validated departure log loss",
    x = "Model",
    y =
      "Departure log loss"
  ) +
  theme_bw()

save_plot(
  plot_cv_logloss,
  "05_frailty_cv_logloss.png"
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
      "Grouped cross-validated departure Brier score",
    x = "Model",
    y =
      "Departure Brier score"
  ) +
  theme_bw()

save_plot(
  plot_cv_brier,
  "06_frailty_cv_brier.png"
)

plot_calibration <- CALIBRATION_TABLE %>%
  filter(
    model %in% c(
      "Fixed Stage-Spell SPAOM",
      "Frailty Stage-Spell SPAOM",
      "Frailty spline spell-age"
    )
  ) %>%
  ggplot(
    aes(
      x =
        mean_predicted,
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
      size =
        transitions
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
      "Marginal cross-validated departure calibration",
    x =
      "Mean predicted probability",
    y =
      "Observed departure proportion",
    size =
      "Transitions"
  ) +
  theme_bw()

save_plot(
  plot_calibration,
  "07_frailty_calibration.png",
  width = 11,
  height = 8
)

advantage_plot_data <- PAIRED_ADVANTAGES %>%
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
      "Paired predictive advantages of the frailty SS-SPAOM",
    subtitle =
      "Positive values favour the frailty model",
    x =
      "Comparator and metric",
    y =
      "Predictive advantage"
  ) +
  theme_bw()

save_plot(
  plot_advantages,
  "08_frailty_paired_advantages.png",
  width = 11,
  height = 8
)

pdf(
  file.path(
    OUTPUT_DIR,
    "Frailty_Sensitivity_Graphs.pdf"
  ),
  width = 11,
  height = 8.5
)

for (plot_object in list(
  plot_beta_comparison,
  plot_random_effects,
  plot_parameter_bootstrap,
  plot_variance_bootstrap,
  plot_cv_logloss,
  plot_cv_brier,
  plot_calibration,
  plot_advantages
)) {
  print(plot_object)
}

dev.off()


# -----------------------------------------------------------------------------
# 15. WORKBOOK, LATEX BLOCK, STATUS, AND SESSION
# -----------------------------------------------------------------------------

writexl::write_xlsx(
  list(
    Model_Comparison =
      MODEL_COMPARISON,
    Primary_Estimates =
      PRIMARY_SENSITIVITY_ESTIMATES,
    Variance_Test =
      VARIANCE_BOOTSTRAP_RESULT$
        summary,
    Bootstrap_Summary =
      PROCESS_BOOTSTRAP_SUMMARY,
    CV_Model_Summary =
      CV_MODEL_SUMMARY,
    Paired_Advantages =
      PAIRED_ADVANTAGES,
    Calibration =
      CALIBRATION_TABLE,
    Random_Effects =
      RANDOM_EFFECT_TABLE,
    Sensitivity_Verdict =
      SENSITIVITY_VERDICT
  ),
  path = file.path(
    OUTPUT_DIR,
    "SS_SPAOM_Frailty_Sensitivity_Results.xlsx"
  )
)

beta_summary <- PROCESS_BOOTSTRAP_SUMMARY %>%
  filter(
    parameter ==
      "beta_spell"
  )

latex_lines <- c(
  "\\subsection{Process-Frailty Sensitivity Analysis}",
  "",
  "A secondary random-intercept extension was fitted:",
  "\\[",
  "\\log\\mu_{iar}=\\log\\eta_a+u_i-\\beta_s\\log s_{ir},",
  "\\qquad",
  "u_i\\sim\\mathcal N(0,\\sigma_u^2).",
  "\\]",
  "",
  paste0(
    "The fixed-effects estimate was \\(\\widehat{\\beta}_s=",
    sprintf(
      "%.4f",
      FIXED_BETA
    ),
    "\\), whereas the frailty estimate was \\(",
    sprintf(
      "%.4f",
      FRAILTY_BETA
    ),
    "\\)."
  ),
  paste0(
    "The estimated frailty standard deviation was \\(",
    sprintf(
      "%.4f",
      FRAILTY_SD
    ),
    "\\)."
  ),
  paste0(
    "The process-bootstrap interval for the frailty-adjusted spell effect was \\([",
    sprintf(
      "%.4f",
      beta_summary$ci_lower[1L]
    ),
    ",\\,",
    sprintf(
      "%.4f",
      beta_summary$ci_upper[1L]
    ),
    "]\\)."
  ),
  paste0(
    "The parametric-bootstrap likelihood-ratio test of \\(\\sigma_u^2=0\\) gave \\(p=",
    sprintf(
      "%.4f",
      variance_p_value
    ),
    "\\)."
  ),
  "",
  paste0(
    "\\textit{Sensitivity interpretation:} ",
    interpretation
  )
)

writeLines(
  latex_lines,
  file.path(
    OUTPUT_DIR,
    "Filled_Frailty_Sensitivity_Block.tex"
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
    nrow(DATA)
  ),
  paste0(
    "Peace processes: ",
    dplyr::n_distinct(
      DATA$PP
    )
  ),
  paste0(
    "Fixed beta: ",
    FIXED_BETA
  ),
  paste0(
    "Frailty beta: ",
    FRAILTY_BETA
  ),
  paste0(
    "Frailty SD: ",
    FRAILTY_SD
  ),
  paste0(
    "Variance bootstrap p-value: ",
    variance_p_value
  ),
  paste0(
    "Successful variance bootstraps: ",
    VARIANCE_BOOTSTRAP_RESULT$
      summary$
      successful_replications[1L],
    " / ",
    VARIANCE_BOOTSTRAP_RESULT$
      summary$
      requested_replications[1L]
  ),
  paste0(
    "Successful process bootstraps: ",
    sum(
      PROCESS_BOOTSTRAP_RESULTS$
        converged
    ),
    " / ",
    PROCESS_BOOTSTRAP_B
  ),
  paste0(
    "Completed grouped CV folds: ",
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
    "Interpretation: ",
    interpretation
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
cat("SS-SPAOM FRAILTY SENSITIVITY ANALYSIS COMPLETED\n")
cat("============================================================\n")
cat("Results are saved in:\n")
cat(OUTPUT_DIR, "\n\n")
cat("Open these files first:\n")
cat(
  file.path(
    OUTPUT_DIR,
    "FRAILTY_INTERPRETATION.txt"
  ),
  "\n"
)
cat(
  file.path(
    OUTPUT_DIR,
    "SS_SPAOM_Frailty_Sensitivity_Results.xlsx"
  ),
  "\n"
)
cat(
  file.path(
    OUTPUT_DIR,
    "Frailty_Sensitivity_Graphs.pdf"
  ),
  "\n"
)
cat(
  file.path(
    OUTPUT_DIR,
    "Filled_Frailty_Sensitivity_Block.tex"
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
