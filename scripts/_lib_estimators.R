# Shared estimator library of the R stages (30 estimators, 40 inference, 60 simulation).
#
# Purpose : one source for the design matrix, the learner library, the AIPW machine with cross-fitting,
#           the parametric rungs, seeds and I/O. The stages call these functions, so every reading and every
#           variant uses the same estimators.
# Thesis  : Sections 3.3 (covariates), 3.7 (estimators, cross-fitting, learners), 3.8 (influence-function SE)
# Usage   : sourced at the top of every R stage, not run on its own. A reading is one covariate specification:
#           "total" gives the total effect (the headline), "cond" the conditional contrast given the investor count
# Input   : config.R (constants), analysis.csv via read_analysis(), variants.csv via read_variants()
# Output  : none (stages write their own outputs; values go through <NN>_values.csv and _values.py merge)
# Seeds   : seed_for(name): h = (h * 31 + code point) mod (2^31 - 1), start SEED_BASE
suppressMessages({ library(glmnet); library(gbm); library(ranger) })

## ---- 1) Where am I: paths relative to this file, so every stage runs from any working directory
if (!exists("SCRIPT_DIR")) {
  .f <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)), fixed = TRUE)   # Rscript encodes spaces
  SCRIPT_DIR <- if (length(.f)) dirname(normalizePath(.f[1])) else getwd()
}
source(file.path(SCRIPT_DIR, "config.R"))
DATA_DIR <- dirname(SCRIPT_DIR)
OUT_DIR <- file.path(DATA_DIR, OUTPUT_DIR); OUT_DATA <- file.path(OUT_DIR, "data"); OUT_TABLES <- file.path(OUT_DIR, "tables")
for (.d in c(OUT_DATA, OUT_TABLES)) dir.create(.d, showWarnings = FALSE, recursive = TRUE)

seed_for <- function(name) {
  h <- SEED_BASE
  for (cp in utf8ToInt(name)) h <- (h * 31 + cp) %% 2147483647
  as.integer(h) }

## ---- 2) Inputs with loud checks: the analysis table must carry the dictionary columns
ANALYSIS_COLUMNS <- c("firm_id", "country", "cohort", "first_round_id", "first_round_date", "first_round_type", "first_round_code",
  "first_round_amount_usd", "amount_missing", "same_day_merged", "rounds_skipped", "n_kein_investor_first", "n_links_first",
  "n_public_first", "n_variant_first", "n_private_first", "has_links", "treated", "treated_broad", "treated_strict",
  "treated_solo", "treated_synd", "first_round_only_unresolved", "sector", "age_at_round", "age_missing", "pre_funding",
  "o1_36", "o1_24", "o1_48", "o1b_36", "o1_undisc_36", "o2_36", "o1_36_only_unresolved", "exit_acq_36", "status")

read_analysis <- function(path = file.path(OUT_DATA, "analysis.csv")) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  miss <- setdiff(ANALYSIS_COLUMNS, names(d))
  if (length(miss)) stop("analysis.csv lacks columns: ", paste(miss, collapse = ", "))
  d }

read_variants <- function(path = file.path(DATA_DIR, "variants.csv")) {
  v <- read.csv(path, stringsAsFactors = FALSE)
  stopifnot(all(c("name", "sample_filter", "treatment_col", "outcome_col", "status") %in% names(v)))
  v }

apply_filter <- function(d, expr_text) {
  keep <- eval(parse(text = expr_text), envir = d, enclos = environment())   # config constants are in scope here
  if (length(keep) == 1) keep <- rep(keep, nrow(d))
  d[!is.na(keep) & keep, , drop = FALSE] }

write_values_csv <- function(vals, labels, path) {
  stopifnot(length(vals) == length(labels), !is.null(names(vals)))
  fmt <- function(v) if (is.numeric(v)) sprintf("%.12g", v) else as.character(v)   # a list value would be splatted into columns
  write.csv(data.frame(key = names(vals), value = vapply(vals, fmt, ""), label = vapply(labels, as.character, ""), stringsAsFactors = FALSE),
            path, row.names = FALSE) }

## ---- 3) Design matrix per reading: age log1p with cohort-median imputation and missing flag,
##      pre-funding, country, sector and cohort dummies; reading "cond" adds the log1p investor count
build_X <- function(d, reading = "total") {
  stopifnot(reading %in% names(READINGS))
  age <- suppressWarnings(as.numeric(d$age_at_round))
  med <- ave(age, d$cohort, FUN = function(v) median(v, na.rm = TRUE))
  med[is.na(med)] <- median(age, na.rm = TRUE)   # fallback: cohort without any age
  X <- data.frame(
    log1p_age   = log1p(ifelse(is.na(age), med, age)),
    age_missing = as.integer(is.na(age)),
    pre_funding = d$pre_funding,
    log1p_ninv  = log1p(d$n_links_first),
    country = factor(d$country), sector = factor(d$sector),
    cohort  = factor(d$cohort))
  if (reading == "total") X$log1p_ninv <- NULL
  Xm <- model.matrix(~ ., data = X)[, -1, drop = FALSE]
  stopifnot(nrow(Xm) == nrow(d), all(is.finite(Xm)))          # model.matrix would drop rows with NA silently
  Xm }
clip <- function(p) pmin(pmax(p, PROPENSITY_CLIP[1]), PROPENSITY_CLIP[2])

## ---- 4) Learner library: fit(x, y) returns a prediction function. Parameters from config LEARNERS.
##      The logit learner and the parametric rungs use glmnet with one tiny lambda, practically an unpenalised logit.
mk_logit <- function(xtr, ytr) { f <- glmnet(xtr, ytr, family = "binomial", lambda = LEARNERS$logit$lambda)
  function(x) as.numeric(predict(f, x, type = "response")) }
mk_lasso <- function(xtr, ytr) { f <- cv.glmnet(xtr, ytr, family = "binomial", nfolds = LEARNERS$lasso$nfolds)
  function(x) as.numeric(predict(f, x, s = "lambda.min", type = "response")) }
mk_gbm <- function(xtr, ytr) { f <- suppressWarnings(gbm.fit(as.data.frame(xtr), ytr,
    distribution = "bernoulli", n.trees = LEARNERS$gbm$n_trees, interaction.depth = LEARNERS$gbm$interaction_depth,
    shrinkage = LEARNERS$gbm$shrinkage, verbose = FALSE))
  function(x) suppressWarnings(predict(f, as.data.frame(x), n.trees = LEARNERS$gbm$n_trees, type = "response")) }
mk_rf <- function(xtr, ytr) { f <- ranger(y = factor(ytr, levels = c(0, 1)),
    x = as.data.frame(xtr), probability = TRUE, num.trees = LEARNERS$rf$num_trees, seed = RF_SEED)
  function(x) predict(f, data = as.data.frame(x))$predictions[, "1"] }
RF_SEED <- SEED_BASE   # stages set RF_SEED <- seed_for(<unit>) before fitting, so every unit is reproducible on its own
LEARNERS_LIB <- list(logit = mk_logit, lasso = mk_lasso, rf = mk_rf, gbm = mk_gbm)

## ---- 5) AIPW machine: score phi, ATT and IF-SE; folds drawn separately per arm; cross-fitting with a fixed pair
# phi = (t - (1-t) e/(1-e)) (y - m0): treated carry their residual, controls are weighted by the propensity odds
aipw <- function(y, t, e, m0) {
  phi <- (t - (1 - t) * e / (1 - e)) * (y - m0)
  att <- mean(phi) / mean(t)
  psi <- (phi - t * att) / mean(t)
  c(att = att, se = sd(psi) / sqrt(length(y))) }

# id: firm ids of a bootstrap sample. Folds are then drawn per distinct firm and every copy of a firm gets the
# fold of its firm, so no firm is in the training and the test fold at once. Without id: one draw per row.
folds_strat <- function(t, K = K_FOLDS, id = NULL) {
  if (!is.null(id)) { first <- !duplicated(id); return(folds_strat(t[first], K)[match(id, id[first])]) }
  fold <- integer(length(t))
  fold[t == 1] <- sample(rep_len(1:K, sum(t == 1)))
  fold[t == 0] <- sample(rep_len(1:K, sum(t == 0)))
  fold }

crossfit_fixed <- function(Xm, t, y, spec, K = K_FOLDS, id = NULL) {
  mk_e <- LEARNERS_LIB[[spec[["e"]]]]; mk_m <- LEARNERS_LIB[[spec[["m0"]]]]
  fold <- folds_strat(t, K, id)
  e <- m0 <- rep(NA_real_, length(t))
  for (k in 1:K) {
    tr <- fold != k; te <- fold == k; ctr <- tr & t == 0
    e[te]  <- clip(mk_e(Xm[tr, , drop = FALSE], t[tr])(Xm[te, , drop = FALSE]))
    m0[te] <- clip(mk_m(Xm[ctr, , drop = FALSE], y[ctr])(Xm[te, , drop = FALSE]))
  }
  list(e = e, m0 = m0, fold = fold) }

## Learner selection (discrete super learner): every learner cross-fitted on the same folds,
## the smallest out-of-fold log-loss wins per nuisance. Returns the pair, the log-losses and the fitted columns.
select_learners <- function(Xm, t, y, K = K_FOLDS) {
  fold <- folds_strat(t, K)
  E <- M <- matrix(NA_real_, length(t), length(LEARNERS_LIB), dimnames = list(NULL, names(LEARNERS_LIB)))
  for (k in 1:K) {
    tr <- fold != k; te <- fold == k; ctr <- tr & t == 0
    for (L in names(LEARNERS_LIB)) {
      E[te, L] <- clip(LEARNERS_LIB[[L]](Xm[tr, , drop = FALSE], t[tr])(Xm[te, , drop = FALSE]))
      M[te, L] <- clip(LEARNERS_LIB[[L]](Xm[ctr, , drop = FALSE], y[ctr])(Xm[te, , drop = FALSE]))
    }
  }
  logloss <- function(p, obs) -mean(obs * log(p) + (1 - obs) * log(1 - p))
  ll_e <- sapply(colnames(E), function(L) logloss(E[, L], t))
  ll_m <- sapply(colnames(M), function(L) logloss(M[t == 0, L], y[t == 0]))
  list(spec = c(e = names(which.min(ll_e)), m0 = names(which.min(ll_m))), ll_e = ll_e, ll_m = ll_m, E = E, M = M, fold = fold) }

## ---- 6) Parametric rungs: naive, outcome regression, IPW (Hajek)
est_naive <- function(y, t, Xm) mean(y[t == 1]) - mean(y[t == 0])
est_or <- function(y, t, Xm) {
  f <- glmnet(Xm[t == 0, ], y[t == 0], family = "binomial", lambda = LEARNERS$logit$lambda)
  mean(y[t == 1]) - mean(as.numeric(predict(f, Xm[t == 1, ], type = "response"))) }
est_ipw <- function(y, t, Xm) {
  e <- clip(as.numeric(predict(glmnet(Xm, t, family = "binomial", lambda = LEARNERS$logit$lambda),
                               Xm, type = "response")))
  w <- e[t == 0] / (1 - e[t == 0])
  mean(y[t == 1]) - weighted.mean(y[t == 0], w) }

## ---- 7) Diagnostics on a fitted split: effective sample size of the controls and weighted balance
ess_controls <- function(e, t) { w <- e[t == 0] / (1 - e[t == 0]); sum(w)^2 / sum(w^2) }
weighted_smd <- function(v, t, w_c) {
  m1 <- mean(v[t == 1]); m0 <- weighted.mean(v[t == 0], w_c)
  s <- sqrt((var(v[t == 1]) + var(v[t == 0])) / 2); if (s == 0) 0 else (m1 - m0) / s }
