#' Run the HaloGen Shiny Application
#'
#' This function launches the main user interface for the HaloGen application,
#' allowing for Likelihood Ratio calculation based on DNA quantities.
#'
#' @export
#' @import shiny
#' @import dplyr
#' @import ggplot2
#' @import rstan
#' @import shinycssloaders
#' @import tidyr
#' @import scales
#' @import rmarkdown
#' @import gtools
#'
#' @return Does not return a value. This function is called for its side effect of launching the Shiny application.

run_halogen_app <- function() {

# Halogen Shiny App
#
# This script creates a fully functional Shiny app for the HaloGen framework.
# It requires a 'data' subfolder containing all the pre-computed model outputs.

# --- 0. Load Libraries ---
# This section now includes a check to automatically install missing packages.
if (!requireNamespace("shiny", quietly = TRUE)) { install.packages("shiny") }
if (!requireNamespace("dplyr", quietly = TRUE)) { install.packages("dplyr") }
if (!requireNamespace("ggplot2", quietly = TRUE)) { install.packages("ggplot2") }
if (!requireNamespace("rstan", quietly = TRUE)) { install.packages("rstan") }
if (!requireNamespace("shinycssloaders", quietly = TRUE)) { install.packages("shinycssloaders") }
if (!requireNamespace("tidyr", quietly = TRUE)) { install.packages("tidyr") }
if (!requireNamespace("scales", quietly = TRUE)) { install.packages("scales") }
if (!requireNamespace("rmarkdown", quietly = TRUE)) { install.packages("rmarkdown") }
if (!requireNamespace("gtools", quietly = TRUE)) { install.packages("gtools") }


library(shiny)
library(dplyr)
library(ggplot2)
library(rstan)
library(shinycssloaders)
library(tidyr)
library(scales)
library(rmarkdown)
library(gtools)


# --- 1. Helper Functions ---
# (These are the core calculation functions we developed)

get_mean_ts <- function(q_vector, post_direct, post_secondary) {
  q_vector <- as.numeric(q_vector)
  q_vector <- q_vector[!is.na(q_vector)]
  if (length(q_vector) == 0) return(list(t = NA, s = NA))

  if (is.null(post_direct) || is.null(post_secondary) || !is.data.frame(post_direct) || !is.data.frame(post_secondary) || nrow(post_direct) == 0 || nrow(post_secondary) == 0) {
    return(list(t = NA, s = NA))
  }

  mu_col <- if ("mu_lab" %in% names(post_direct)) "mu_lab" else "mu"
  log_sigma_col <- if ("log_sigma_lab" %in% names(post_direct)) "log_sigma_lab" else "log_sigma"
  k_col <- if ("k_lab" %in% names(post_direct)) "k_lab" else "k"

  n_total <- min(nrow(post_direct), nrow(post_secondary))
  if (!is.finite(n_total) || n_total == 0) return(list(t = NA, s = NA))

  get_joint_likelihoods <- function(samples, mu_col_name, log_sigma_col_name, k_col_name) {
    sapply(seq_len(nrow(samples)), function(i) {
      mu <- samples[[mu_col_name]][i]
      sigma_val <- exp(samples[[log_sigma_col_name]][i])
      k <- samples[[k_col_name]][i]

      eps <- 1e-9
      k_clamped <- max(eps, min(k, 1.0 - eps, na.rm = TRUE))
      sigma_adjusted <- if (is.na(sigma_val) || !is.finite(sigma_val) || sigma_val <= 0) eps else sigma_val

      per_stain_likelihoods <- numeric(length(q_vector))

      for (j in seq_along(q_vector)) {
        q_val <- q_vector[j]

        if (is.na(q_val)) {
          per_stain_likelihoods[j] <- 1e-300
        } else if (q_val == 0) {
          per_stain_likelihoods[j] <- k_clamped
        } else if (q_val > 0 && is.finite(q_val)) {
          density_val <- dlnorm(q_val, meanlog = mu, sdlog = sigma_adjusted)
          likelihood_pos <- (1 - k_clamped) * density_val
          per_stain_likelihoods[j] <- if (!is.finite(likelihood_pos) || likelihood_pos < 1e-300) 1e-300 else likelihood_pos
        } else {
          per_stain_likelihoods[j] <- 1e-300
        }
      }

      prod(per_stain_likelihoods, na.rm = TRUE)
    })
  }

  mean_t_full <- mean(get_joint_likelihoods(post_direct, mu_col, log_sigma_col, k_col), na.rm = TRUE)
  mean_s_full <- mean(get_joint_likelihoods(post_secondary, mu_col, log_sigma_col, k_col), na.rm = TRUE)

  return(list(t = max(mean_t_full, 1e-300), s = max(mean_s_full, 1e-300)))
}

compute_general_LRs <- function(known_t, known_s, unknown_tau, unknown_s, failrate, nSuspects) {
  known_t <- as.numeric(known_t)
  known_s <- as.numeric(known_s)
  unknown_tau <- as.numeric(unknown_tau)
  unknown_s <- as.numeric(unknown_s)

  nK <- length(known_t)
  nU <- length(unknown_tau)
  n <- nK + nU

  if (n == 0 && nSuspects > 0) {
    return(list(LR_suspects=numeric(0), LR_unknown=numeric(0), LR_none=Inf))
  } else if (n == 0 && nSuspects == 0) {
    return(list(LR_suspects=numeric(0), LR_unknown=numeric(0), LR_none=1))
  }

  totalHypotheses <- 2^n
  likelihoods <- numeric(totalHypotheses)
  subsetList <- vector("list", totalHypotheses)

  for(h_idx in 0:(totalHypotheses-1)){
    included <- as.integer(intToBits(h_idx))[1:n]
    subsetList[[h_idx+1]] <- included

    if(sum(included) > nSuspects){
      likelihoods[h_idx+1] <- NA
      next
    }

    L_known <- 1.0
    if(nK > 0){
      for(i_k in seq_len(nK)){
        term_val <- ifelse(included[i_k]==1, known_t[i_k], known_s[i_k])
        L_known <- L_known * term_val
      }
    }

    L_unknown <- 1.0
    if(nU > 0){
      for(j_u in seq_len(nU)){
        idx_loop <- nK + j_u
        term_val <- ifelse(included[idx_loop]==1, unknown_tau[j_u], unknown_s[j_u])
        L_unknown <- L_unknown * term_val
      }
    }

    L <- L_known * L_unknown
    x <- sum(included)

    if(x < nSuspects){
      safe_failrate <- max(min(failrate, 1.0 - 1e-9), 1e-9)
      power_term <- nSuspects - x
      L <- L * (safe_failrate^power_term)
    }

    if (is.na(L)) {
      likelihoods[h_idx + 1] <- 1e-300
    } else {
      likelihoods[h_idx + 1] <- L
    }
  }

  valid_indices <- !is.na(likelihoods)
  likelihoods_filtered <- likelihoods[valid_indices]
  subsetList_filtered <- subsetList[valid_indices]

  LR_suspects <- numeric(nK)
  if(nK > 0){
    for(i_s in seq_len(nK)){
      sapply_num_arg <- sapply(subsetList_filtered, function(sub) if(length(sub)>=i_s) sub[i_s]==1 else FALSE)
      sapply_den_arg <- sapply(subsetList_filtered, function(sub) if(length(sub)>=i_s) sub[i_s]==0 else FALSE)
      num_sum <- sum(likelihoods_filtered[sapply_num_arg], na.rm=T)
      den_sum <- sum(likelihoods_filtered[sapply_den_arg], na.rm=T)
      LR_suspects[i_s] <- ifelse(den_sum > 0, num_sum / den_sum, Inf)
    }
  }

  LR_unknown <- numeric(nU)
  if(nU > 0){
    for(j_u in seq_len(nU)){
      idx_loop <- nK + j_u
      sapply_num_arg <- sapply(subsetList_filtered, function(sub) if(length(sub)>=idx_loop) sub[idx_loop]==1 else FALSE)
      sapply_den_arg <- sapply(subsetList_filtered, function(sub) if(length(sub)>=idx_loop) sub[idx_loop]==0 else FALSE)
      num_sum <- sum(likelihoods_filtered[sapply_num_arg], na.rm=T)
      den_sum <- sum(likelihoods_filtered[sapply_den_arg], na.rm=T)
      LR_unknown[j_u] <- ifelse(den_sum > 0, num_sum / den_sum, Inf)
    }
  }

  none_index <- which(sapply(subsetList_filtered, sum) == 0)
  LR_none <- NA
  if(length(none_index) == 1){
    lik_none <- likelihoods_filtered[none_index]
    lik_others_sum <- sum(likelihoods_filtered[-none_index], na.rm=T)
    LR_none <- ifelse(lik_others_sum > 0, lik_none / lik_others_sum, Inf)
  } else if(length(none_index) == 0 && nSuspects > 0) {
    LR_none <- 0
  } else {
    LR_none <- NA
  }

  return(list(LR_suspects=LR_suspects, LR_unknown=LR_unknown, LR_none=LR_none))
}

get_F0_lab <- function(lab_id, params_df) {
  lab_info <- params_df[params_df$lab == lab_id, ]
  if (nrow(lab_info) == 0) return(NA)
  (lab_info$fail_count[1] + 0.5) / (lab_info$sample_size[1] + 1)
}
get_F0_group <- function(params_df) {
  nsamples <- sum(params_df$sample_size, na.rm = TRUE)
  m <- sum(params_df$sample_size * params_df$fail_rate, na.rm = TRUE)
  (m + 0.5) / (nsamples + 1)
}

# --- 2. User Interface (UI) ---
ui <- fluidPage(
  titlePanel("HaloGen LR Simulation App"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4("Instructions:"),
      p("Select a tab to either analyse a pre-loaded lab or upload your own data. Fill in the case details and click 'Calculate LR'."),

      tabsetPanel(
        id = "main_tabs",
        tabPanel("Analyse Pre-loaded Lab",
                 br(),
                 selectInput("lab_choice", "1. Select Laboratory:", choices = NULL),
                 hr(),
                 h4("Case Circumstances"),
                 numericInput("ns_offenders_preloaded", "Number of Assumed Offenders (Ns):", value = 1, min = 1),
                 numericInput("n_known_preloaded", "Number of Known Contributors (POIs):", value = 1, min = 0),
                 numericInput("n_unknown_preloaded", "Number of Unknown Contributors:", value = 0, min = 0),
                 numericInput("total_stains_preloaded", "Total Number of Stains Analysed:", value = 1, min = 1),
                 hr(),
                 h4("Evidence Quantities (ng)"),
                 uiOutput("known_inputs_preloaded"),
                 uiOutput("unknown_inputs_preloaded")
        ),
        tabPanel("Analyse New Lab Data",
                 br(),
                 p("Upload a CSV file with two columns: 'transfer_type' ('direct' or 'secondary') and 'quantity' (ng)."),
                 fileInput("new_lab_file", "1. Upload New Lab Data (.csv)", accept = ".csv"),
                 actionButton("fit_new_lab_model", "2. Fit New Lab Model", class = "btn-primary"),
                 hr(),
                 uiOutput("new_lab_ui")
        )
      ),
      hr(),
      actionButton("calculate_lr", "Calculate LR", icon = icon("calculator"), class = "btn-success btn-lg btn-block"),
      hr(),
      downloadButton("download_report", "Download Report", class = "btn-info btn-block")
    ),
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("Results Summary",
                 br(),
                 h3("Likelihood Ratio Results"),
                 withSpinner(tableOutput("lr_results_table")),
                 hr(),
                 h4("Contextual Plot"),
                 htmlOutput("contextual_plot_text"),
                 withSpinner(plotOutput("contextual_plot"))
        ),
        tabPanel("Model Parameters",
                 br(),
                 h3("Underlying Model Parameters"),
                 p("This table shows the mean parameter estimates for each model used in the calculation."),
                 withSpinner(tableOutput("parameter_table"))
        )
      )
    )
  )
)

# --- 3. Server Logic ---
server <- function(input, output, session) {

  # --- A. Dynamic UI and Data Loading ---

  observe({
    data_path <- system.file("app/data", package = "halogen")
      req(dir.exists(data_path))
    lab_ids <- list.dirs(data_path, full.names = FALSE, recursive = FALSE)
    if (require(gtools)) {
      sorted_lab_ids <- gtools::mixedsort(lab_ids)
    } else {
      sorted_lab_ids <- sort(lab_ids)
    }
    updateSelectInput(session, "lab_choice", choices = sorted_lab_ids)
  })

  rv <- reactiveValues(
    results_data = NULL,
    new_lab_posteriors = list(
      bayes_direct = NULL, bayes_secondary = NULL,
      vague_direct = NULL, vague_secondary = NULL,
      F0 = NULL
    )
  )

  observeEvent(input$main_tabs, {
    rv$results_data <- NULL
  })

  output$known_inputs_preloaded <- renderUI({
    req(input$n_known_preloaded > 0, input$total_stains_preloaded > 0)
    lapply(seq_len(input$n_known_preloaded), function(i) {
      tagList(
        h5(paste("Known Contributor", i)),
        lapply(seq_len(input$total_stains_preloaded), function(j) {
          numericInput(paste0("k_quant_preloaded_", i, "_", j), paste("Stain", j, "Quantity (ng):"), value = 0)
        })
      )
    })
  })

  output$unknown_inputs_preloaded <- renderUI({
    req(input$n_unknown_preloaded > 0, input$total_stains_preloaded > 0)
    lapply(seq_len(input$n_unknown_preloaded), function(i) {
      tagList(
        h5(paste("Unknown Contributor", i)),
        lapply(seq_len(input$total_stains_preloaded), function(j) {
          numericInput(paste0("u_quant_preloaded_", i, "_", j), paste("Stain", j, "Quantity (ng):"), value = 0)
        })
      )
    })
  })

  output$new_lab_ui <- renderUI({
    req(rv$new_lab_posteriors$bayes_direct)
    tagList(
      hr(),
      h4("3. Case Circumstances"),
      numericInput("ns_offenders_new", "Number of Assumed Offenders (Ns):", value = 1, min = 1),
      numericInput("n_known_new", "Number of Known Contributors (POIs):", value = 1, min = 0),
      numericInput("n_unknown_new", "Number of Unknown Contributors:", value = 0, min = 0),
      numericInput("total_stains_new", "Total Number of Stains Analysed:", value = 1, min = 1),
      hr(),
      h4("Evidence Quantities (ng)"),
      uiOutput("known_inputs_new"),
      uiOutput("unknown_inputs_new")
    )
  })

  output$known_inputs_new <- renderUI({
    req(input$n_known_new > 0, input$total_stains_new > 0)
    lapply(seq_len(input$n_known_new), function(i) {
      tagList(
        h5(paste("Known Contributor", i)),
        lapply(seq_len(input$total_stains_new), function(j) {
          numericInput(paste0("k_quant_new_", i, "_", j), paste("Stain", j, "Quantity (ng):"), value = 0)
        })
      )
    })
  })

  output$unknown_inputs_new <- renderUI({
    req(input$n_unknown_new > 0, input$total_stains_new > 0)
    lapply(seq_len(input$n_unknown_new), function(i) {
      tagList(
        h5(paste("Unknown Contributor", i)),
        lapply(seq_len(input$total_stains_new), function(j) {
          numericInput(paste0("u_quant_new_", i, "_", j), paste("Stain", j, "Quantity (ng):"), value = 0)
        })
      )
    })
  })

  # --- C. New Lab Model Fitting ---
  observeEvent(input$fit_new_lab_model, {
    req(input$new_lab_file)

data_path <- system.file("app/data", package = "halogen")
post_group <- readRDS(file.path(data_path, "post_group.rds"))
lab_model <- readRDS(file.path(data_path, "lab_model.rds"))
group_secondary_post <- readRDS(file.path(data_path, "group_secondary.rds"))

    if (is.null(post_group) || is.null(lab_model)) {
      showModal(modalDialog(
        title = "Error: Missing Core Model Files",
        "The 'post_group.rds' and 'lab_model.rds' files are required for fitting new lab data. Please ensure they are in the 'data' subfolder.",
        easyClose = TRUE
      ))
      return()
    }

    withProgress(message = 'Fitting new lab data...', value = 0, {

      setProgress(0.1, detail = "Reading data")
      new_data <- read.csv(input$new_lab_file$datapath)
      detection_limit_raw <- 0.001

      y_direct_raw <- new_data %>% filter(transfer_type == "direct") %>% pull(quantity)
      y_direct_stan <- y_direct_raw[y_direct_raw > detection_limit_raw & is.finite(y_direct_raw)]
      N_direct_stan <- length(y_direct_stan)

      y_secondary_raw <- new_data %>% filter(transfer_type == "secondary") %>% pull(quantity)
      y_secondary_stan <- y_secondary_raw[y_secondary_raw > detection_limit_raw & is.finite(y_secondary_raw)]
      N_secondary_stan <- length(y_secondary_stan)

      vague_priors <- list(mu0 = 0, tau_mu = 100, log_sigma0 = 0, tau_log_sigma = 10, mu_k = 0.5, phi_k = 1)

      setProgress(0.2, detail = "Fitting Bayes Direct model...")
      stan_data_lb_d <- list(N_lab = N_direct_stan, y_lab = y_direct_stan, mu0 = mean(post_group$mu0_dir), tau_mu = mean(post_group$tau_mu_dir), log_sigma0 = mean(post_group$log_sigma0_dir), tau_log_sigma = mean(post_group$tau_log_sigma_dir), mu_k = mean(post_group$mu_k_dir), phi_k = mean(post_group$phi_k_dir))
      fit_lb_d <- rstan::sampling(lab_model, data = stan_data_lb_d, iter = 2000, chains = 4, seed = 123)
      rv$new_lab_posteriors$bayes_direct <- as.data.frame(rstan::extract(fit_lb_d, pars = c("mu_lab", "log_sigma_lab", "k_lab")))

      setProgress(0.4, detail = "Fitting Vague Direct model...")
      stan_data_lv_d <- list(N_lab = N_direct_stan, y_lab = y_direct_stan, mu0 = vague_priors$mu0, tau_mu = vague_priors$tau_mu, log_sigma0 = vague_priors$log_sigma0, tau_log_sigma = vague_priors$tau_log_sigma, mu_k = vague_priors$mu_k, phi_k = vague_priors$phi_k)
      fit_lv_d <- rstan::sampling(lab_model, data = stan_data_lv_d, iter = 2000, chains = 4, seed = 123)
      rv$new_lab_posteriors$vague_direct <- as.data.frame(rstan::extract(fit_lv_d, pars = c("mu_lab", "log_sigma_lab", "k_lab")))

      if (N_secondary_stan > 0) {
        setProgress(0.6, detail = "Fitting Bayes Secondary model...")
        stan_data_lb_s <- list(N_lab = N_secondary_stan, y_lab = y_secondary_stan, mu0 = mean(post_group$mu0_sec), tau_mu = mean(post_group$tau_mu_sec), log_sigma0 = mean(post_group$log_sigma0_sec), tau_log_sigma = mean(post_group$tau_log_sigma_sec), mu_k = mean(post_group$mu_k_sec), phi_k = mean(post_group$phi_k_sec))
        fit_lb_s <- rstan::sampling(lab_model, data = stan_data_lb_s, iter = 2000, chains = 4, seed = 123)
        rv$new_lab_posteriors$bayes_secondary <- as.data.frame(rstan::extract(fit_lb_s, pars = c("mu_lab", "log_sigma_lab", "k_lab")))

        setProgress(0.8, detail = "Fitting Vague Secondary model...")
        stan_data_lv_s <- list(N_lab = N_secondary_stan, y_lab = y_secondary_stan, mu0 = vague_priors$mu0, tau_mu = vague_priors$tau_mu, log_sigma0 = vague_priors$log_sigma0, tau_log_sigma = vague_priors$tau_log_sigma, mu_k = vague_priors$mu_k, phi_k = vague_priors$phi_k)
        fit_lv_s <- rstan::sampling(lab_model, data = stan_data_lv_s, iter = 2000, chains = 4, seed = 123)
        rv$new_lab_posteriors$vague_secondary <- as.data.frame(rstan::extract(fit_lv_s, pars = c("mu_lab", "log_sigma_lab", "k_lab")))
      } else {
        rv$new_lab_posteriors$bayes_secondary <- group_secondary_post
        rv$new_lab_posteriors$vague_secondary <- group_secondary_post
      }

      fail_count <- sum(y_direct_raw <= detection_limit_raw)
      sample_size <- length(y_direct_raw)
      rv$new_lab_posteriors$F0 <- (fail_count + 0.5) / (sample_size + 1)

      setProgress(1, detail = "Done!")
    })

    showNotification("New lab model fitting complete.", type = "message")
  })

  # --- D. Main LR Calculation ---
  observeEvent(input$calculate_lr, {

    data_path <- system.file("app/data", package = "halogen")
     direct_params <- read.csv(file.path(data_path, "direct_params.csv"), stringsAsFactors = FALSE)
    group_direct_post <- readRDS(file.path(data_path, "group_direct.rds"))
  group_secondary_post <- readRDS(file.path(data_path, "group_secondary.rds"))

    if (input$main_tabs == "Analyse Pre-loaded Lab") {
      lab_id <- input$lab_choice
      req(lab_id)
      n_known <- input$n_known_preloaded
      n_unknown <- input$n_unknown_preloaded
      Ns <- input$ns_offenders_preloaded
      total_stains <- input$total_stains_preloaded

      known_quants <- if (n_known > 0) lapply(seq_len(n_known), function(i) {
        sapply(seq_len(total_stains), function(j) {
          input[[paste0("k_quant_preloaded_", i, "_", j)]]
        })
      }) else list()

      unknown_quants <- if (n_unknown > 0) lapply(seq_len(n_unknown), function(i) {
        sapply(seq_len(total_stains), function(j) {
          input[[paste0("u_quant_preloaded_", i, "_", j)]]
        })
      }) else list()

base_path <- file.path(data_path, lab_id)
post_bayes_d <- readRDS(file.path(base_path, paste0("post_lab_bayes_direct_", lab_id, ".rds")))
post_bayes_s <- readRDS(file.path(base_path, paste0("post_lab_bayes_secondary_", lab_id, ".rds")))
post_vague_d <- readRDS(file.path(base_path, paste0("post_lab_vague_direct_", lab_id, ".rds")))
post_vague_s <- readRDS(file.path(base_path, paste0("post_lab_vague_secondary_", lab_id, ".rds")))
      F0_lab <- get_F0_lab(lab_id, direct_params)

    } else { # New Lab Tab
      req(rv$new_lab_posteriors$bayes_direct)
      lab_id <- "New Lab"
      n_known <- input$n_known_new
      n_unknown <- input$n_unknown_new
      Ns <- input$ns_offenders_new
      total_stains <- input$total_stains_new

      known_quants <- if (n_known > 0) lapply(seq_len(n_known), function(i) {
        sapply(seq_len(total_stains), function(j) {
          val <- input[[paste0("k_quant_new_", i, "_", j)]]
          if (is.null(val)) 0 else val
        })
      }) else list()

      unknown_quants <- if (n_unknown > 0) lapply(seq_len(n_unknown), function(i) {
        sapply(seq_len(total_stains), function(j) {
          val <- input[[paste0("u_quant_new_", i, "_", j)]]
          if (is.null(val)) 0 else val
        })
      }) else list()

      post_bayes_d <- rv$new_lab_posteriors$bayes_direct
      post_bayes_s <- rv$new_lab_posteriors$bayes_secondary
      post_vague_d <- rv$new_lab_posteriors$vague_direct
      post_vague_s <- rv$new_lab_posteriors$vague_secondary
      F0_lab <- rv$new_lab_posteriors$F0
    }

    contributor_names <- c(if(n_known > 0) paste("Known", 1:n_known), if(n_unknown > 0) paste("Unknown", 1:n_unknown))

    validate(
      need((n_known + n_unknown) <= 14, "Error: Too many contributors. The app supports a maximum of 14 total contributors to ensure performance.")
    )

    withProgress(message = 'Calculating LRs...', value = 0, {

      ts_bayes_known <- lapply(known_quants, get_mean_ts, post_bayes_d, post_bayes_s)
      ts_bayes_unknown <- lapply(unknown_quants, get_mean_ts, post_bayes_d, post_bayes_s)
      lr_res_bayes <- compute_general_LRs(
        known_t = if (n_known > 0) unlist(lapply(ts_bayes_known, `[[`, "t"), use.names = FALSE) else numeric(0),
        known_s = if (n_known > 0) unlist(lapply(ts_bayes_known, `[[`, "s"), use.names = FALSE) else numeric(0),
        unknown_tau = if (n_unknown > 0) unlist(lapply(ts_bayes_unknown, `[[`, "t"), use.names = FALSE) else numeric(0),
        unknown_s = if (n_unknown > 0) unlist(lapply(ts_bayes_unknown, `[[`, "s"), use.names = FALSE) else numeric(0),
        failrate = F0_lab,
        nSuspects = Ns
      )

      ts_vague_known <- lapply(known_quants, get_mean_ts, post_vague_d, post_vague_s)
      ts_vague_unknown <- lapply(unknown_quants, get_mean_ts, post_vague_d, post_vague_s)
      lr_res_vague <- compute_general_LRs(
        known_t = if (n_known > 0) unlist(lapply(ts_vague_known, `[[`, "t"), use.names = FALSE) else numeric(0),
        known_s = if (n_known > 0) unlist(lapply(ts_vague_known, `[[`, "s"), use.names = FALSE) else numeric(0),
        unknown_tau = if (n_unknown > 0) unlist(lapply(ts_vague_unknown, `[[`, "t"), use.names = FALSE) else numeric(0),
        unknown_s = if (n_unknown > 0) unlist(lapply(ts_vague_unknown, `[[`, "s"), use.names = FALSE) else numeric(0),
        failrate = F0_lab,
        nSuspects = Ns
      )

      F0_group <- get_F0_group(direct_params)
      ts_group_known <- lapply(known_quants, get_mean_ts, group_direct_post, group_secondary_post)
      ts_group_unknown <- lapply(unknown_quants, get_mean_ts, group_direct_post, group_secondary_post)
      lr_res_group <- compute_general_LRs(
        known_t = if (n_known > 0) unlist(lapply(ts_group_known, `[[`, "t"), use.names = FALSE) else numeric(0),
        known_s = if (n_known > 0) unlist(lapply(ts_group_known, `[[`, "s"), use.names = FALSE) else numeric(0),
        unknown_tau = if (n_unknown > 0) unlist(lapply(ts_group_unknown, `[[`, "t"), use.names = FALSE) else numeric(0),
        unknown_s = if (n_unknown > 0) unlist(lapply(ts_group_unknown, `[[`, "s"), use.names = FALSE) else numeric(0),
        failrate = F0_group,
        nSuspects = Ns
      )
    })

    results_df <- data.frame(
      Contributor = contributor_names,
      `Bayes Lab` = log10(c(lr_res_bayes$LR_suspects, lr_res_bayes$LR_unknown)),
      `Standalone Lab` = log10(c(lr_res_vague$LR_suspects, lr_res_vague$LR_unknown)),
      `Group Model` = log10(c(lr_res_group$LR_suspects, lr_res_group$LR_unknown)),
      check.names = FALSE
    )

    lr_none_row <- data.frame(
      Contributor = "None",
      `Bayes Lab` = log10(lr_res_bayes$LR_none),
      `Standalone Lab` = log10(lr_res_vague$LR_none),
      `Group Model` = log10(lr_res_group$LR_none),
      check.names = FALSE
    )
    results_df <- bind_rows(results_df, lr_none_row)

    get_params <- function(post, mu_col = "mu_lab", log_sigma_col = "log_sigma_lab", k_col = "k_lab") {
      c(
        mu = mean(post[[mu_col]], na.rm = TRUE),
        sigma = exp(mean(post[[log_sigma_col]], na.rm = TRUE)),
        k = mean(post[[k_col]], na.rm = TRUE)
      )
    }

    params_bayes_d <- get_params(post_bayes_d)
    params_bayes_s <- get_params(post_bayes_s)
    params_vague_d <- get_params(post_vague_d)
    params_vague_s <- get_params(post_vague_s)
    params_group_d <- get_params(group_direct_post, "mu", "log_sigma", "k")
    params_group_s <- get_params(group_secondary_post, "mu", "log_sigma", "k")

    param_table <- data.frame(
      Parameter = rep(c("mu", "sigma", "k"), 2),
      Transfer = rep(c("Direct", "Secondary"), each = 3),
      `Bayes Lab` = c(params_bayes_d, params_bayes_s),
      `Standalone Lab` = c(params_vague_d, params_vague_s),
      `Group Model` = c(params_group_d, params_group_s),
      check.names = FALSE
    )

    f0_row <- data.frame(
      Parameter = "F0", Transfer = "-",
      `Bayes Lab` = F0_lab, `Standalone Lab` = F0_lab, `Group Model` = F0_group,
      check.names = FALSE
    )

    param_table <- bind_rows(param_table, f0_row)

    rv$results_data <- list(table = results_df, param_table = param_table, all_quants = c(known_quants, unknown_quants), contributor_names = contributor_names, lab = lab_id)
  })

  # --- E. Render Outputs ---
  output$lr_results_table <- renderTable({
    req(rv$results_data)
    rv$results_data$table
  }, digits = 2)

  output$parameter_table <- renderTable({
    req(rv$results_data)
    rv$results_data$param_table
  }, digits = 4)

  output$contextual_plot_text <- renderUI({
    HTML(
      "<p>The following plot illustrates the Mean Log LR curves for the selected laboratory across a range of DNA quantities from the first 'known' individual. Since it is a simple 2d plot, it is only conditional on one contributor varying, with all other contributors fixed. The vertical lines indicate the specific DNA quantities entered for the first 'known' individual in this case (i.e. for each stain), providing a visual context for the results in the table above.</p><p>If there is more than one stain, the values in the table will represent the joint LR which is the sum of logs of the individual stain LRs for the known individual. There is no informational context in the plot to help intepret results for other contributors (this is reserved for future iterations of the program).</p>"
    )
  })

  output$contextual_plot <- renderPlot({

    req(rv$results_data)

    res <- rv$results_data
    lab_id <- res$lab
    all_quants <- res$all_quants
    contributor_names <- res$contributor_names

    if (input$main_tabs == "Analyse Pre-loaded Lab") {
      Ns_for_plot <- input$ns_offenders_preloaded
      nK_val <- input$n_known_preloaded
    } else {
      Ns_for_plot <- input$ns_offenders_new
      nK_val <- input$n_known_new
    }
    # Define the path AGAIN, inside this block
    data_path <- system.file("app/data", package = "halogen")
    direct_params <- read.csv(file.path(data_path, "direct_params.csv"), stringsAsFactors = FALSE)
    group_direct_post <- readRDS(file.path(data_path, "group_direct.rds"))
    group_secondary_post <- readRDS(file.path(data_path, "group_secondary.rds"))

    quant_grid <- exp(seq(log(0.01), log(10), length.out = 40))

    plot_data_list <- lapply(quant_grid, function(q) {
      if (lab_id == "New Lab") {
        post_bayes_d <- rv$new_lab_posteriors$bayes_direct
        post_bayes_s <- rv$new_lab_posteriors$bayes_secondary
        post_vague_d <- rv$new_lab_posteriors$vague_direct
        post_vague_s <- rv$new_lab_posteriors$vague_secondary
        F0_lab <- rv$new_lab_posteriors$F0
      } else {
        base_path <- file.path(data_path, lab_id)
        post_bayes_d <- readRDS(file.path(base_path, paste0("post_lab_bayes_direct_", lab_id, ".rds")))
        post_bayes_s <- readRDS(file.path(base_path, paste0("post_lab_bayes_secondary_", lab_id, ".rds")))
        post_vague_d <- readRDS(file.path(base_path, paste0("post_lab_vague_direct_", lab_id, ".rds")))
        post_vague_s <- readRDS(file.path(base_path, paste0("post_lab_vague_secondary_", lab_id, ".rds")))
        F0_lab <- get_F0_lab(lab_id, direct_params)
      }
      F0_group <- get_F0_group(direct_params)

      known_quants_plot <- if(length(all_quants) > 0) all_quants else list(q)
      known_quants_plot[[1]] <- q

      ts_b <- lapply(known_quants_plot, get_mean_ts, post_bayes_d, post_bayes_s)
      lr_b <- compute_general_LRs(sapply(ts_b, "[[", "t"), sapply(ts_b, "[[", "s"), numeric(0), numeric(0), F0_lab, Ns_for_plot)$LR_suspects[1]

      ts_v <- lapply(known_quants_plot, get_mean_ts, post_vague_d, post_vague_s)
      lr_v <- compute_general_LRs(sapply(ts_v, "[[", "t"), sapply(ts_v, "[[", "s"), numeric(0), numeric(0), F0_lab, Ns_for_plot)$LR_suspects[1]

      ts_g <- lapply(known_quants_plot, get_mean_ts, group_direct_post, group_secondary_post)
      lr_g <- compute_general_LRs(sapply(ts_g, "[[", "t"), sapply(ts_g, "[[", "s"), numeric(0), numeric(0), F0_group, Ns_for_plot)$LR_suspects[1]

      data.frame(
        Quantity = q,
        Model = c("Bayes Lab", "Standalone Lab", "Group Model"),
        Mean_Log_LR = log10(c(lr_b, lr_v, lr_g))
      )
    })

    plot_data <- bind_rows(plot_data_list)

    vline_data_raw <- if (length(all_quants) > 0) {
      data.frame(
        Quantity    = unlist(all_quants),
        Contributor = rep(contributor_names, times = sapply(all_quants, length))
      )
    } else {
      data.frame(Quantity=numeric(0), Contributor=character(0))
    }

    known_contributor_names <- contributor_names[1:nK_val]
    vline_data_known <- subset(vline_data_raw,
                               Contributor %in% known_contributor_names &
                                 is.finite(Quantity) &
                                 Quantity > 0)

    model_colors <- c("Standalone Lab" = "red", "Bayes Lab" = "blue", "Group Model" = "black")
    model_linetypes <- c("Standalone Lab" = "solid", "Bayes Lab" = "solid", "Group Model" = "dashed")

    contributor_colors <- setNames(hcl(h = seq(15, 375, length = length(known_contributor_names) + 1)[1:length(known_contributor_names)], l = 65, c = 100), known_contributor_names)
    contributor_linetypes <- setNames(rep_len(c("dotted", "dotdash", "longdash", "twodash"), length.out = length(known_contributor_names)), known_contributor_names)

    ggplot(plot_data, aes(x = Quantity, y = Mean_Log_LR, color = Model)) +
      geom_line(aes(linetype = Model), linewidth = 1.2) +
      geom_vline(data = vline_data_known, aes(xintercept = Quantity, color = Contributor, linetype = Contributor), linewidth = 1, alpha = 0.8) +
      scale_x_log10(
        breaks = c(0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10),
        labels = c("0.01", "0.05", "0.1", "0.5", "1", "2", "5", "10")
      ) +
      scale_color_manual(
        name = "Trace",
        values = c(model_colors, contributor_colors)
      ) +
      scale_linetype_manual(
        name = "Trace",
        values = c(model_linetypes, contributor_linetypes)
      ) +
      labs(
        title = paste("LR Curves for", lab_id),
        x = "Quantity (ng) on Log10 Scale",
        y = "Mean Log Likelihood Ratio (log(LR))"
      ) +
      theme_bw(base_size = 16) +
      theme(legend.position = "bottom")
  })

  # --- F. Report Generation ---
  output$download_report <- downloadHandler(
    filename = function() {
      paste0("HaloGen_Report_", Sys.Date(), ".html")
    },
    content = function(file) {
      withProgress(message = 'Generating Report...', value = 0, {

        incProgress(0.1, detail = "Fetching Results...")

        res <- rv$results_data

        lab_id <- res$lab
        all_quants <- res$all_quants
        contributor_names <- res$contributor_names

        if (input$main_tabs == "Analyse Pre-loaded Lab") {
          Ns_val <- input$ns_offenders_preloaded
          nK_val <- input$n_known_preloaded
          nU_val <- input$n_unknown_preloaded
        } else {
          Ns_val <- input$ns_offenders_new
          nK_val <- input$n_known_new
          nU_val <- input$n_unknown_new
        }

        Ns_for_plot <- Ns_val

        incProgress(0.3, detail = "Preparing contextual plot data...")

        # *** FIX: Re-generate the plot object for the report ***
data_path <- system.file("app/data", package = "halogen")
direct_params <- read.csv(file.path(data_path, "direct_params.csv"), stringsAsFactors = FALSE)
group_direct_post <- readRDS(file.path(data_path, "group_direct.rds"))
group_secondary_post <- readRDS(file.path(data_path, "group_secondary.rds"))

        quant_grid <- exp(seq(log(0.01), log(10), length.out = 40))

        plot_data_list <- lapply(quant_grid, function(q) {
          if (lab_id == "New Lab") {
            post_bayes_d <- rv$new_lab_posteriors$bayes_direct
            post_bayes_s <- rv$new_lab_posteriors$bayes_secondary
            post_vague_d <- rv$new_lab_posteriors$vague_direct
            post_vague_s <- rv$new_lab_posteriors$vague_secondary
            F0_lab <- rv$new_lab_posteriors$F0
          } else {
            base_path <- file.path(data_path, lab_id)
            post_bayes_d <- readRDS(file.path(base_path, paste0("post_lab_bayes_direct_", lab_id, ".rds")))
            post_bayes_s <- readRDS(file.path(base_path, paste0("post_lab_bayes_secondary_", lab_id, ".rds")))
            post_vague_d <- readRDS(file.path(base_path, paste0("post_lab_vague_direct_", lab_id, ".rds")))
            post_vague_s <- readRDS(file.path(base_path, paste0("post_lab_vague_secondary_", lab_id, ".rds")))
            F0_lab <- get_F0_lab(lab_id, direct_params)
          }
          F0_group <- get_F0_group(direct_params)

          known_quants_plot <- if(length(all_quants) > 0) all_quants else list(q)
          known_quants_plot[[1]] <- q

          ts_b <- lapply(known_quants_plot, get_mean_ts, post_bayes_d, post_bayes_s)
          lr_b <- compute_general_LRs(sapply(ts_b, "[[", "t"), sapply(ts_b, "[[", "s"), numeric(0), numeric(0), F0_lab, Ns_for_plot)$LR_suspects[1]

          ts_v <- lapply(known_quants_plot, get_mean_ts, post_vague_d, post_vague_s)
          lr_v <- compute_general_LRs(sapply(ts_v, "[[", "t"), sapply(ts_v, "[[", "s"), numeric(0), numeric(0), F0_lab, Ns_for_plot)$LR_suspects[1]

          ts_g <- lapply(known_quants_plot, get_mean_ts, group_direct_post, group_secondary_post)
          lr_g <- compute_general_LRs(sapply(ts_g, "[[", "t"), sapply(ts_g, "[[", "s"), numeric(0), numeric(0), F0_group, Ns_for_plot)$LR_suspects[1]

          data.frame(
            Quantity = q,
            Model = c("Bayes Lab", "Standalone Lab", "Group Model"),
            Mean_Log_LR = log10(c(lr_b, lr_v, lr_g))
          )
        })

        plot_data <- bind_rows(plot_data_list)

        vline_data_raw <- if (length(all_quants) > 0) {
          data.frame(
            Quantity    = unlist(all_quants),
            Contributor = rep(contributor_names, times = sapply(all_quants, length))
          )
        } else {
          data.frame(Quantity=numeric(0), Contributor=character(0))
        }

        known_contributor_names <- contributor_names[1:nK_val]
        vline_data_known <- subset(vline_data_raw,
                                   Contributor %in% known_contributor_names &
                                     is.finite(Quantity) &
                                     Quantity > 0)

        model_colors <- c("Standalone Lab" = "red", "Bayes Lab" = "blue", "Group Model" = "black")
        model_linetypes <- c("Standalone Lab" = "solid", "Bayes Lab" = "solid", "Group Model" = "dashed")

        contributor_colors <- setNames(hcl(h = seq(15, 375, length = length(known_contributor_names) + 1)[1:length(known_contributor_names)], l = 65, c = 100), known_contributor_names)
        contributor_linetypes <- setNames(rep_len(c("dotted", "dotdash", "longdash", "twodash"), length.out = length(known_contributor_names)), known_contributor_names)

        contextual_plot <- ggplot(plot_data, aes(x = Quantity, y = Mean_Log_LR, color = Model)) +
          geom_line(aes(linetype = Model), linewidth = 1.2) +
          geom_vline(data = vline_data_known, aes(xintercept = Quantity, color = Contributor, linetype = Contributor), linewidth = 1, alpha = 0.8) +
          scale_x_log10(
            breaks = c(0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10),
            labels = c("0.01", "0.05", "0.1", "0.5", "1", "2", "5", "10")
          ) +
          scale_color_manual(
            name = "Trace",
            values = c(model_colors, contributor_colors)
          ) +
          scale_linetype_manual(
            name = "Trace",
            values = c(model_linetypes, contributor_linetypes)
          ) +
          labs(
            title = paste("LR Curves for", lab_id),
            x = "Quantity (ng) on Log10 Scale",
            y = "Mean Log Likelihood Ratio (log(LR))"
          ) +
          theme_bw(base_size = 16) +
          theme(legend.position = "bottom")

        incProgress(0.8, detail = "Copying Rmd template...")

        report_path <- system.file("app/report.Rmd", package = "halogen")
        temp_report <- file.path(tempdir(), "report.Rmd")
        file.copy(report_path, temp_report, overwrite = TRUE)

        incProgress(0.9, detail = "Rendering report...")

        rmarkdown::render(
          input = temp_report,
          output_file = file,
          params = list(
            lab_id            = res$lab,
            case_details      = paste("Ns =", Ns_val, "; Known =", nK_val, "; Unknown =", nU_val),
            results_table     = res$table,
            param_table       = res$param_table,
            context_plot      = contextual_plot,
            contributor_names = res$contributor_names,
            all_quants        = res$all_quants
          ),
          envir = new.env(parent = globalenv())
        )

        incProgress(1, detail = "Report generation complete.")
      })
    }
  )

} # This is the closing brace for the server function.

# --- 4. Run the App ---
shinyApp(ui = ui, server = server)

} # <-- Add this closing brace at the very end of the file.
