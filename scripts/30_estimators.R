#!/usr/bin/env Rscript
# Stage 30: estimator ladder, learner choice, one-split AIPW with split band, robustness variants, diagnostics.
#
# Purpose : for both readings ("total": covariates without the investor count, the causal headline; "cond": with it,
#           the conditional contrast, fixed before the first estimate) the ladder naive -> outcome regression -> IPW
#           -> AIPW with cross-fitting, the learner pair chosen by cross-fitted log-loss, the one-split AIPW as the main number,
#           N_SPLITS further splits as the band, every row of variants.csv under the fixed pair, and the propensity
#           and balance diagnostics of the main split.
# Thesis  : Sections 3.7, 3.8 (estimator ladder, learner choice, split band), 4.1, 4.2; Appendix A.5
# Input   : outputs/data/analysis.csv, variants.csv (folder above scripts/), config.R
# Output  : outputs/data/30_values.csv (-> _values.py merge), outputs/tables/30_ladder.csv, 30_variants.csv,
#           30_splits_<reading>.csv, 30_firm_diagnostics_<reading>.csv, 30_balance_<reading>.csv
# Usage   : Rscript 30_estimators.R [<analysis.csv> <out_dir>], then python3 _values.py merge
#           ../outputs/data/30_values.csv 30_estimators.R (run_all.sh does both); without arguments the
#           defaults above are used, the two arguments run the stage on another table
# Seeds   : seed_for("30_boot_naive"), seed_for("30_boot_<rung>_<reading>"), seed_for("30_select_<reading>"), seed_for("30_split_<reading>_<s>"),
#           seed_for("30_variant_<name>_<reading>"); <s> is written with three digits (001)
# Runtime : about 6 minutes on 8 cores
# Checks  : column contract of analysis.csv, at least MIN_TREATED treated per variant, no missing outcome
.f <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)), fixed = TRUE)   # Rscript encodes spaces
SCRIPT_DIR <- if (length(.f)) dirname(normalizePath(.f[1])) else getwd()
source(file.path(SCRIPT_DIR, "_lib_estimators.R"))
suppressMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
ANALYSIS <- if (length(args) >= 1) args[1] else file.path(OUT_DATA, "analysis.csv")
if (length(args) >= 2) { OUT_DATA <- file.path(args[2], "data"); OUT_TABLES <- file.path(args[2], "tables")
  for (.d in c(OUT_DATA, OUT_TABLES)) dir.create(.d, showWarnings = FALSE, recursive = TRUE) }
CORES <- max(1L, detectCores() - 2L)
Z <- qnorm(1 - (1 - CI_LEVEL) / 2)

## ---- 1) Data: the primary sample is the analysis table restricted to firms with a first-round link
d <- read_analysis(ANALYSIS)
prim <- d[d$has_links == 1, ]
stopifnot(sum(prim$treated) >= MIN_TREATED, !anyNA(prim$o1_36))
variants <- read_variants()
vals <- c(); labels <- c()
put <- function(key, value, label) { vals[[key]] <<- value; labels[[key]] <<- label }

ladder_rows <- list(); variant_rows <- list(); spec_of <- list()
for (reading in names(READINGS)) {
  cat("\n==== reading", reading, "====\n")
  Xm <- build_X(prim, reading); y <- prim$o1_36; t <- prim$treated; n <- nrow(prim)

  ## ---- 2) Rungs 1 to 3 with a firm-level bootstrap (percentile CI), one seed per rung
  boot <- function(fn, name) {
    set.seed(seed_for(if (name == "naive") "30_boot_naive" else paste0("30_boot_", name, "_", reading)))   # the naive rung does not depend on the reading
    th <- fn(y, t, Xm)
    bs <- replicate(B_BOOTSTRAP, { i <- sample.int(n, n, replace = TRUE)
      tryCatch({ if (sum(t[i]) < MIN_TREATED) stop("degenerate world"); fn(y[i], t[i], Xm[i, , drop = FALSE]) }, error = function(e) NA) })
    c(att = th, lo = unname(quantile(bs, (1 - CI_LEVEL) / 2, na.rm = TRUE)), hi = unname(quantile(bs, 1 - (1 - CI_LEVEL) / 2, na.rm = TRUE)),
      se = sd(bs, na.rm = TRUE), fail = sum(is.na(bs))) }
  cat("bootstrap of rungs 1 to 3 ...\n")
  rungs <- list(naive = boot(est_naive, "naive"), or = boot(est_or, "or"), ipw = boot(est_ipw, "ipw"))
  for (r in names(rungs)) {
    lab <- c(naive = "naive", or = "outcome-regression", ipw = "IPW")[[r]]
    put(sprintf("s30_att_%s_%s", r, reading), rungs[[r]][["att"]], sprintf("ATT, %s rung, reading %s", lab, reading))
    put(sprintf("s30_ci_lo_%s_%s", r, reading), rungs[[r]][["lo"]], sprintf("lower bootstrap CI, %s rung, reading %s", lab, reading))
    put(sprintf("s30_ci_hi_%s_%s", r, reading), rungs[[r]][["hi"]], sprintf("upper bootstrap CI, %s rung, reading %s", lab, reading))
    put(sprintf("s30_boot_fail_%s_%s", r, reading), rungs[[r]][["fail"]], sprintf("failed bootstrap draws, %s rung, reading %s", lab, reading))
  }

  ## ---- 3) Learner choice on one split (discrete super learner), the same split gives the main AIPW estimate
  cat("learner selection ...\n")
  set.seed(seed_for(paste0("30_select_", reading))); RF_SEED <<- seed_for(paste0("30_select_", reading))
  sel <- select_learners(Xm, t, y); spec <- sel$spec; spec_of[[reading]] <- spec
  put(paste0("s30_learner_e_", reading), spec[["e"]], paste("learner chosen for the propensity, reading", reading))
  put(paste0("s30_learner_m0_", reading), spec[["m0"]], paste("learner chosen for the control outcome, reading", reading))
  for (L in names(LEARNERS_LIB)) {
    a <- aipw(y, t, sel$E[, L], sel$M[, L])
    put(sprintf("s30_att_aipw_%s_%s", L, reading), a[["att"]], sprintf("AIPW ATT with %s for both nuisances, reading %s", L, reading))
    put(sprintf("s30_ll_e_%s_%s", L, reading), sel$ll_e[[L]], sprintf("cross-fitted log-loss of %s for the propensity, reading %s", L, reading))
    put(sprintf("s30_ll_m0_%s_%s", L, reading), sel$ll_m[[L]], sprintf("cross-fitted log-loss of %s for the control outcome, reading %s", L, reading))
  }
  e_main <- sel$E[, spec[["e"]]]; m0_main <- sel$M[, spec[["m0"]]]
  main <- aipw(y, t, e_main, m0_main)
  w_c <- e_main[t == 0] / (1 - e_main[t == 0]); ess <- ess_controls(e_main, t)
  put(paste0("s30_att_main_", reading), main[["att"]], paste("ATT, one-split AIPW with the chosen pair, reading", reading))
  put(paste0("s30_se_main_", reading), main[["se"]], paste("influence-function SE of the main AIPW, reading", reading))
  put(paste0("s30_ci_lo_main_", reading), main[["att"]] - Z * main[["se"]], paste("lower IF-CI of the main AIPW, reading", reading))
  put(paste0("s30_ci_hi_main_", reading), main[["att"]] + Z * main[["se"]], paste("upper IF-CI of the main AIPW, reading", reading))
  put(paste0("s30_n_main_", reading), n, "firms in the primary sample"); put(paste0("s30_nt_main_", reading), sum(t), "treated firms in the primary sample")
  put(paste0("s30_ess_main_", reading), ess, paste("effective number of controls under ATT weights, reading", reading))
  ladder_rows[[reading]] <- data.frame(reading = reading,
    rung = c("naive", "outcome regression", "IPW", paste0("AIPW ", names(LEARNERS_LIB)), sprintf("AIPW main (e: %s, m0: %s)", spec[["e"]], spec[["m0"]])),
    att = c(sapply(rungs, `[[`, "att"), sapply(names(LEARNERS_LIB), function(L) aipw(y, t, sel$E[, L], sel$M[, L])[["att"]]), main[["att"]]),
    ci_lo = c(sapply(rungs, `[[`, "lo"), sapply(names(LEARNERS_LIB), function(L) { a <- aipw(y, t, sel$E[, L], sel$M[, L]); a[["att"]] - Z * a[["se"]] }), main[["att"]] - Z * main[["se"]]),
    ci_hi = c(sapply(rungs, `[[`, "hi"), sapply(names(LEARNERS_LIB), function(L) { a <- aipw(y, t, sel$E[, L], sel$M[, L]); a[["att"]] + Z * a[["se"]] }), main[["att"]] + Z * main[["se"]]),
    inference = c(rep("bootstrap percentile", 3), rep("influence function", length(LEARNERS_LIB) + 1)), stringsAsFactors = FALSE)

  ## ---- 4) Diagnostics of the main split: overlap, clipping, weighted balance, AUC, calibration, country weight share
  smd_u <- sapply(colnames(Xm), function(cn) { v <- Xm[, cn]; s <- sqrt((var(v[t == 1]) + var(v[t == 0])) / 2); if (s == 0) 0 else (mean(v[t == 1]) - mean(v[t == 0])) / s })
  smd_w <- sapply(colnames(Xm), function(cn) weighted_smd(Xm[, cn], t, w_c))
  auc <- { r <- rank(e_main); (sum(r[t == 1]) - sum(t) * (sum(t) + 1) / 2) / (sum(t) * sum(t == 0)) }   # Mann-Whitney
  cal <- unname(coef(suppressWarnings(glm(t ~ qlogis(e_main), family = binomial)))[2])
  w_country <- tapply(w_c, prim$country[t == 0], sum); w_top4 <- sum(sort(w_country, decreasing = TRUE)[1:4]) / sum(w_c)
  diag <- c(e_auc = auc, e_cal_slope = cal, e_max_control = max(e_main[t == 0]), e_max_treated = max(e_main[t == 1]),
            n_treated_e_gt_05 = sum(e_main[t == 1] > 0.5), n_clip_e_lo = sum(e_main <= PROPENSITY_CLIP[1] + 1e-12),
            n_clip_e_hi = sum(e_main >= PROPENSITY_CLIP[2] - 1e-12), max_abs_smd_weighted = max(abs(smd_w)),
            n_cols_smd_gt_01 = sum(abs(smd_w) > 0.1), n_design_cols = ncol(Xm), w_share_top4_countries = w_top4)
  diag_label <- c(e_auc = "AUC of the propensity score", e_cal_slope = "calibration slope of the propensity score",
                  e_max_control = "largest propensity among controls", e_max_treated = "largest propensity among treated",
                  n_treated_e_gt_05 = "treated firms with a propensity above 0.5",
                  n_clip_e_lo = "firms clipped at the lower propensity bound", n_clip_e_hi = "firms clipped at the upper propensity bound",
                  max_abs_smd_weighted = "largest absolute weighted SMD over the design columns",
                  n_cols_smd_gt_01 = "design columns with an absolute weighted SMD above 0.1", n_design_cols = "columns of the design matrix",
                  w_share_top4_countries = "share of the control weight in the four heaviest countries")
  for (k in names(diag)) put(sprintf("s30_%s_%s", k, reading), diag[[k]],
                             sprintf("%s%s, reading %s", diag_label[[k]], if (k == "n_design_cols") "" else ", main split", reading))
  w_att <- ifelse(t == 1, 1, e_main / (1 - e_main))
  write.csv(data.frame(firm_id = prim$firm_id, treated = t, y = y, e_main = e_main, m0_main = m0_main, w_att = w_att, fold = sel$fold),
            file.path(OUT_TABLES, sprintf("30_firm_diagnostics_%s.csv", reading)), row.names = FALSE)
  write.csv(data.frame(column = colnames(Xm), mean_treated = colMeans(Xm[t == 1, , drop = FALSE]), mean_control = colMeans(Xm[t == 0, , drop = FALSE]),
                       mean_control_weighted = apply(Xm[t == 0, , drop = FALSE], 2, weighted.mean, w = w_c), smd_unweighted = smd_u, smd_weighted = smd_w),
            file.path(OUT_TABLES, sprintf("30_balance_%s.csv", reading)), row.names = FALSE)

  ## ---- 5) Split band: N_SPLITS further cross-fitting splits under the fixed pair, one seed per split
  cat("split band ...\n")
  splits <- do.call(rbind, mclapply(seq_len(N_SPLITS), function(s) {
    set.seed(seed_for(sprintf("30_split_%s_%03d", reading, s))); RF_SEED <<- seed_for(sprintf("30_split_%s_%03d", reading, s))
    cf <- crossfit_fixed(Xm, t, y, spec); a <- aipw(y, t, cf$e, cf$m0)
    data.frame(split = s, att = a[["att"]], se = a[["se"]]) }, mc.cores = CORES, mc.set.seed = FALSE))
  write.csv(splits, file.path(OUT_TABLES, sprintf("30_splits_%s.csv", reading)), row.names = FALSE)
  med <- median(splits$att)
  put("s30_n_splits", N_SPLITS, "number of cross-fitting splits in the band")
  put(paste0("s30_att_split_median_main_", reading), med, paste("median ATT over the splits, reading", reading))
  put(paste0("s30_att_split_min_main_", reading), min(splits$att), paste("smallest ATT over the splits, reading", reading))
  put(paste0("s30_att_split_max_main_", reading), max(splits$att), paste("largest ATT over the splits, reading", reading))
  put(paste0("s30_att_split_sd_main_", reading), sd(splits$att), paste("SD of the ATT over the splits, reading", reading))
  se_med <- sqrt(median(splits$se^2 + (splits$att - med)^2))   # standard error of the median, split variation included
  put(paste0("s30_se_split_median_main_", reading), se_med, paste("median-adjusted SE over the splits (Chernozhukov et al. 2018), reading", reading))
  put(paste0("s30_share_splits_negative_main_", reading), mean(splits$att < 0), paste("share of splits with a negative ATT, reading", reading))
  put(paste0("s30_share_splits_ci_excl_zero_main_", reading), mean(splits$att + Z * splits$se < 0 | splits$att - Z * splits$se > 0), paste("share of splits whose IF-CI excludes zero, reading", reading))
  put(paste0("s30_ci_lo_split_median_main_", reading), med - Z * se_med, paste("lower CI of the median over the splits (median-adjusted SE), reading", reading))
  put(paste0("s30_ci_hi_split_median_main_", reading), med + Z * se_med, paste("upper CI of the median over the splits (median-adjusted SE), reading", reading))
  put(paste0("s30_share_splits_le_main_", reading), mean(splits$att <= main[["att"]]), paste("share of splits with an ATT at or below the main split, reading", reading))

  ## ---- 6) Variants: every row of variants.csv under the fixed pair; "main" reuses the main split above instead of fitting again
  cat("variants ...\n")
  for (i in seq_len(nrow(variants))) {
    v <- variants[i, ]
    if (v$name == "main") { r <- main; nn <- n; nt <- sum(t); es <- ess; nv <- mean(y[t == 1]) - mean(y[t == 0]) } else {
      dv <- apply_filter(d, v$sample_filter)
      if (reading == "cond" && any(dv$has_links == 0))   # no-link firms have no investor count: median of linked controls
        dv$n_links_first[dv$has_links == 0] <- median(prim$n_links_first[prim$treated == 0])
      tv <- dv[[v$treatment_col]]; yv <- dv[[v$outcome_col]]
      stopifnot(!anyNA(yv), !anyNA(tv), sum(tv) >= MIN_TREATED)
      set.seed(seed_for(sprintf("30_variant_%s_%s", v$name, reading))); RF_SEED <<- seed_for(sprintf("30_variant_%s_%s", v$name, reading))
      cf <- crossfit_fixed(build_X(dv, reading), tv, yv, spec); r <- aipw(yv, tv, cf$e, cf$m0)
      nn <- length(yv); nt <- sum(tv); es <- ess_controls(cf$e, tv); nv <- mean(yv[tv == 1]) - mean(yv[tv == 0])
    }
    cat(sprintf("%-16s %-5s N=%5d T=%3d naive=%+.1f AIPW=%+.1f [%+.1f, %+.1f]\n", v$name, reading, nn, nt, 100 * nv, 100 * r[["att"]],
                100 * (r[["att"]] - Z * r[["se"]]), 100 * (r[["att"]] + Z * r[["se"]])))
    if (v$name != "main") {                                   # the main keys and their labels are written in block 3 above
      put(sprintf("s30_att_%s_%s", v$name, reading), r[["att"]], sprintf("ATT, variant %s, reading %s", v$name, reading))
      put(sprintf("s30_se_%s_%s", v$name, reading), r[["se"]], sprintf("IF-SE, variant %s, reading %s", v$name, reading))
      put(sprintf("s30_ci_lo_%s_%s", v$name, reading), r[["att"]] - Z * r[["se"]], sprintf("lower IF-CI, variant %s, reading %s", v$name, reading))
      put(sprintf("s30_ci_hi_%s_%s", v$name, reading), r[["att"]] + Z * r[["se"]], sprintf("upper IF-CI, variant %s, reading %s", v$name, reading))
      put(sprintf("s30_n_%s_%s", v$name, reading), nn, sprintf("firms, variant %s, reading %s", v$name, reading))
      put(sprintf("s30_nt_%s_%s", v$name, reading), nt, sprintf("treated firms, variant %s, reading %s", v$name, reading))
      put(sprintf("s30_ess_%s_%s", v$name, reading), es, sprintf("effective controls, variant %s, reading %s", v$name, reading))
    }
    variant_rows[[length(variant_rows) + 1]] <- data.frame(name = v$name, reading = reading, description = v$description, status = v$status,
      n = nn, n_treated = nt, naive = nv, att = r[["att"]], se = r[["se"]],
      ci_lo = r[["att"]] - Z * r[["se"]], ci_hi = r[["att"]] + Z * r[["se"]], ess = es,
      key_att = sprintf("s30_att_%s_%s", v$name, reading), stringsAsFactors = FALSE)
  }
}

## ---- 7) Outputs: ladder and variant tables, values CSV for the merge
write.csv(do.call(rbind, ladder_rows), file.path(OUT_TABLES, "30_ladder.csv"), row.names = FALSE)
write.csv(do.call(rbind, variant_rows), file.path(OUT_TABLES, "30_variants.csv"), row.names = FALSE)
write_values_csv(vals, labels, file.path(OUT_DATA, "30_values.csv"))
cat("\nlearner pairs:", paste(names(spec_of), sapply(spec_of, function(s) paste(s, collapse = "/")), collapse = " · "), "\n")
cat("wrote", length(vals), "keys ->", file.path(OUT_DATA, "30_values.csv"), "\n")
