#!/usr/bin/env Rscript
# Stage 40: bootstrap and causal-forest cross-checks of the main AIPW estimate, both readings.
#
# Purpose : firm-level bootstrap of the complete main specification (design matrix, cross-fitting with the pair
#           chosen in stage 30, AIPW) for both readings; cluster bootstrap over country x cohort for the headline
#           reading "total"; grf causal forest as an independent implementation for both readings.
# Thesis  : Sections 3.8, 4.1; Appendix A.3
# Input   : outputs/data/analysis.csv, outputs/data/30_values.csv (learner pairs), config.R
# Output  : outputs/data/40_values.csv (-> _values.py merge), outputs/data/40_boot_draws.csv (one row per replicate,
#           resume file) with 40_boot_draws.csv.meta (fingerprint of the inputs), outputs/tables/40_bootstrap.csv
# Usage   : Rscript 40_inference.R, then python3 _values.py merge ../outputs/data/40_values.csv 40_inference.R
#           (run_all.sh does both)
# Seeds   : seed_for("40_boot_<mode>_<reading>_<rep>") with <rep> in four digits (0001), seed_for("40_grf_<reading>")
# Runtime : one to two hours on 8 cores (2 x B_BOOTSTRAP firm replicates, B_BOOTSTRAP cluster replicates); resumes after an interruption
# Checks  : every replicate with fewer than MIN_TREATED treated or MIN_CONTROLS controls is counted as failed, never
#           dropped silently; evaluation only when all replicates exist; cross-fitting folds inside a replicate are
#           drawn per firm, so the copies of a resampled firm never sit on both sides of a fold
.f <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)), fixed = TRUE)
SCRIPT_DIR <- if (length(.f)) dirname(normalizePath(.f[1])) else getwd()
source(file.path(SCRIPT_DIR, "_lib_estimators.R"))
suppressMessages({ library(parallel); library(grf) })
CORES <- max(1L, detectCores() - 2L); CHUNK <- 50L
Z <- qnorm(1 - (1 - CI_LEVEL) / 2)
DRAWS <- file.path(OUT_DATA, "40_boot_draws.csv"); META <- paste0(DRAWS, ".meta")   # fingerprint of the inputs the draws belong to

## ---- 1) Data, learner pairs from stage 30, cluster ids
d <- read_analysis(); prim <- d[d$has_links == 1, ]
v30 <- read.csv(file.path(OUT_DATA, "30_values.csv"), stringsAsFactors = FALSE)
spec_of <- lapply(setNames(names(READINGS), names(READINGS)), function(r)
  c(e = v30$value[v30$key == paste0("s30_learner_e_", r)], m0 = v30$value[v30$key == paste0("s30_learner_m0_", r)]))
stopifnot(all(sapply(spec_of, length) == 2))
.const <- tempfile(); writeLines(deparse(list(K_FOLDS, PROPENSITY_CLIP, LEARNERS, MIN_TREATED, MIN_CONTROLS)), .const)   # constants the replicates depend on
FINGERPRINT <- paste(unname(tools::md5sum(file.path(OUT_DATA, "analysis.csv"))), paste(unlist(spec_of), collapse = "/"), B_BOOTSTRAP, SEED_BASE, "folds-by-firm", unname(tools::md5sum(.const)))
if (file.exists(DRAWS)) { if (!file.exists(META) || readLines(META, n = 1) != FINGERPRINT) stop("40_boot_draws.csv belongs to other inputs (data, learner pairs, B, seeds, fold rule or constants): delete it to recompute") } else writeLines(FINGERPRINT, META)
cl_id <- interaction(prim$country, prim$cohort, drop = TRUE); n_cl <- nlevels(cl_id)

## ---- 2) Replicate plan: firm bootstrap for both readings, cluster bootstrap for the headline reading
plan <- rbind(expand.grid(mode = "firm", reading = names(READINGS), rep = seq_len(B_BOOTSTRAP), stringsAsFactors = FALSE),
              expand.grid(mode = "cluster", reading = "total", rep = seq_len(B_BOOTSTRAP), stringsAsFactors = FALSE))
plan$seed <- mapply(function(m, r, b) seed_for(sprintf("40_boot_%s_%s_%04d", m, r, b)), plan$mode, plan$reading, plan$rep)

one_rep <- function(mode, reading, rep, seed) {
  set.seed(seed)
  idx <- if (mode == "firm") sample.int(nrow(prim), replace = TRUE)
         else unlist(lapply(sample(levels(cl_id), n_cl, replace = TRUE), function(cc) which(cl_id == cc)), use.names = FALSE)
  dd <- prim[idx, ]
  att <- tryCatch({
    y <- dd$o1_36; t <- dd$treated
    if (sum(t) < MIN_TREATED || sum(t == 0) < MIN_CONTROLS) stop("degenerate world")
    RF_SEED <<- seed
    cf <- crossfit_fixed(build_X(dd, reading), t, y, spec_of[[reading]], id = dd$firm_id)   # copies of a firm share a fold
    aipw(y, t, cf$e, cf$m0)[["att"]] }, error = function(e) NA_real_)
  data.frame(mode = mode, reading = reading, rep = rep, seed = seed, att = att) }

## ---- 3) Run with checkpoint: finished replicates are in the draws file, only open ones are computed
done <- if (file.exists(DRAWS)) read.csv(DRAWS, stringsAsFactors = FALSE) else plan[0, ]
todo <- plan[!paste(plan$mode, plan$reading, plan$rep) %in% paste(done$mode, done$reading, done$rep), ]
cat(sprintf("bootstrap: %d of %d replicates open · %d cores · %d clusters\n", nrow(todo), nrow(plan), CORES, n_cl))
while (nrow(todo) > 0) {
  batch <- todo[seq_len(min(CHUNK, nrow(todo))), ]; t0 <- Sys.time()
  res <- do.call(rbind, mclapply(seq_len(nrow(batch)), function(i) one_rep(batch$mode[i], batch$reading[i], batch$rep[i], batch$seed[i]),
                                 mc.cores = CORES, mc.set.seed = FALSE))
  write.table(res, DRAWS, sep = ",", row.names = FALSE, append = file.exists(DRAWS), col.names = !file.exists(DRAWS))
  todo <- todo[-seq_len(nrow(batch)), ]
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("[%s] chunk done (%.1f min) · open %d · remaining ~%.0f min\n", format(Sys.time(), "%H:%M"), el, nrow(todo), el * nrow(todo) / nrow(batch)))
}

## ---- 4) Evaluation: bootstrap SE and percentile CI per mode and reading, failures counted
dr <- read.csv(DRAWS, stringsAsFactors = FALSE); stopifnot(nrow(dr) == nrow(plan))
vals <- c(); labels <- c(); put <- function(k, v, l) { vals[[k]] <<- v; labels[[k]] <<- l }
put("s40_boot_B", B_BOOTSTRAP, "bootstrap replicates per mode and reading"); put("s40_n_clusters", n_cl, "country x cohort clusters")
cs <- as.integer(table(cl_id)); put("s40_cluster_size_median", median(cs), "median firms per cluster"); put("s40_cluster_size_max", max(cs), "firms in the largest cluster"); put("s40_cluster_size_min", min(cs), "firms in the smallest cluster")
tab <- list()
for (m in unique(plan$mode)) for (r in unique(plan$reading[plan$mode == m])) {
  v <- dr$att[dr$mode == m & dr$reading == r]; ok <- v[!is.na(v)]
  st <- c(se = sd(ok), lo = unname(quantile(ok, (1 - CI_LEVEL) / 2)), hi = unname(quantile(ok, 1 - (1 - CI_LEVEL) / 2)), fail = sum(is.na(v)))
  put(sprintf("s40_boot_se_%s_%s", m, r), st[["se"]], sprintf("bootstrap SE, %s resampling, reading %s", m, r))
  put(sprintf("s40_boot_ci_lo_%s_%s", m, r), st[["lo"]], sprintf("lower percentile bootstrap CI, %s resampling, reading %s", m, r))
  put(sprintf("s40_boot_ci_hi_%s_%s", m, r), st[["hi"]], sprintf("upper percentile bootstrap CI, %s resampling, reading %s", m, r))
  put(sprintf("s40_boot_fail_%s_%s", m, r), st[["fail"]], sprintf("failed replicates, %s resampling, reading %s", m, r))
  tab[[length(tab) + 1]] <- data.frame(mode = m, reading = r, B = B_BOOTSTRAP, se = st[["se"]], ci_lo = st[["lo"]], ci_hi = st[["hi"]], fail = st[["fail"]])
  cat(sprintf("%-8s %-6s SE %.4f · CI [%.4f, %.4f] · %d failures\n", m, r, st[["se"]], st[["lo"]], st[["hi"]], st[["fail"]]))
}

## ---- 5) Causal forest (grf): own nuisance estimation, honest trees, doubly robust ATT on the treated
for (r in names(READINGS)) {
  Xm <- build_X(prim, r); y <- prim$o1_36; t <- prim$treated
  cf <- causal_forest(X = Xm, Y = y, W = t, num.trees = GRF_NUM_TREES, seed = seed_for(paste0("40_grf_", r)))
  a <- average_treatment_effect(cf, target.sample = "treated"); est <- a[["estimate"]]; se <- a[["std.err"]]
  put(paste0("s40_grf_att_", r), est, paste("causal-forest ATT on the treated, reading", r))
  put(paste0("s40_grf_se_", r), se, paste("causal-forest SE, reading", r))
  put(paste0("s40_grf_ci_lo_", r), est - Z * se, paste("lower causal-forest CI, reading", r)); put(paste0("s40_grf_ci_hi_", r), est + Z * se, paste("upper causal-forest CI, reading", r))
  put(paste0("s40_grf_what_max_treated_", r), max(cf$W.hat[t == 1]), paste("largest forest propensity among treated, reading", r))
  put(paste0("s40_grf_what_max_control_", r), max(cf$W.hat[t == 0]), paste("largest forest propensity among controls, reading", r))
  cat(sprintf("grf %-6s ATT %+.4f [%+.4f, %+.4f]\n", r, est, est - Z * se, est + Z * se))
}
write.csv(do.call(rbind, tab), file.path(OUT_TABLES, "40_bootstrap.csv"), row.names = FALSE)
write_values_csv(vals, labels, file.path(OUT_DATA, "40_values.csv"))
cat("wrote", length(vals), "keys -> 40_values.csv\n")
