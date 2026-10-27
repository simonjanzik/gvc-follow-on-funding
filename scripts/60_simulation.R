#!/usr/bin/env Rscript
# Stage 60: simulation study (ADEMP) on the headline covariates, with misspecification and unobserved-confounder scenarios.
#
# Purpose : does the estimator ladder recover a known ATT on data like ours? Covariates are resampled rows of the
#           real design matrix (reading "total"), assignment and control outcome are logit fits on the real data with
#           intercepts recentred on the observed rates; the true effect is additive on the probability scale.
#           Scenarios: S0, S1, S2 (well specified, true ATT from SIM_SCENARIOS), M1 (an interaction and a squared age term
#           in assignment and outcome), M2 (a non-additive term from the candidate list, dampened to the asymptotic IPW bias
#           SIM_M2_IPW_TARGET or, if that is out of reach, to the propensity cap: key s60_m2_mode),
#           U scenarios rstar, shift5, shift10 (a binary unobserved confounder whose strength is calibrated so that
#           the X-adjusted bias makes the true ATT appear as the headline). Estimators: naive, OR, IPW, AIPW with the
#           learner pair of reading "total". Performance: bias and coverage, each with its Monte-Carlo SE, plus empirical
#           SE and RMSE.
# Thesis  : Sections 3.9, 4.3; Appendix A.4
# Input   : outputs/data/analysis.csv, outputs/data/30_values.csv (headline ATT for the U calibration, learner pair), config.R
# Output  : outputs/data/60_values.csv (-> _values.py merge), outputs/data/60_dgp_params.rds (deterministic cache with
#           data hash), outputs/data/60_sim_draws.csv (one row per replicate and estimator, resume file) with
#           60_sim_draws.csv.meta (fingerprint of the inputs), outputs/tables/60_simulation.csv, 60_simulation_u.csv,
#           60_dgp_params.csv (calibration of the U scenarios)
# Usage   : Rscript 60_simulation.R, then python3 _values.py merge ../outputs/data/60_values.csv 60_simulation.R
#           (run_all.sh does both)
# Seeds   : seed_for("60_<scenario>_<rep>") with <rep> in four digits (0001); the DGP calibration uses no random numbers
# Runtime : about one hour on 8 cores (SIM_REPS replicates x 8 scenarios); resumes after an interruption
# Checks  : cache and draws must carry the fingerprint of the current inputs (data, learner pair, headline, simulation
#           constants); a replicate with fewer than MIN_TREATED treated or MIN_CONTROLS controls is a failure, counted; every scenario must
#           hold SIM_REPS replicates before the evaluation; M2 stops if no candidate stays under the R2 cap
.f <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)), fixed = TRUE)
SCRIPT_DIR <- if (length(.f)) dirname(normalizePath(.f[1])) else getwd()
source(file.path(SCRIPT_DIR, "_lib_estimators.R"))
suppressMessages(library(parallel))
CORES <- max(1L, detectCores() - 2L); CHUNK <- 40L
Z <- qnorm(1 - (1 - CI_LEVEL) / 2)
DRAWS <- file.path(OUT_DATA, "60_sim_draws.csv"); PARAMS <- file.path(OUT_DATA, "60_dgp_params.rds"); META <- paste0(DRAWS, ".meta")

## ---- 1) Real data as template: primary sample, headline design matrix, learner pair and headline ATT from stage 30
d <- read_analysis(); prim <- d[d$has_links == 1, ]
n <- nrow(prim); Xm <- build_X(prim, "total"); y_real <- prim$o1_36; t_real <- prim$treated
v30 <- read.csv(file.path(OUT_DATA, "30_values.csv"), stringsAsFactors = FALSE)
spec <- c(e = v30$value[v30$key == "s30_learner_e_total"], m0 = v30$value[v30$key == "s30_learner_m0_total"])
HEADLINE <- as.numeric(v30$value[v30$key == "s30_att_main_total"]); stopifnot(length(spec) == 2, length(HEADLINE) == 1)
DATA_HASH <- unname(tools::md5sum(file.path(OUT_DATA, "analysis.csv")))
# fingerprint of everything the cache and the draws depend on: data, learner pair, headline, constants of the design
.const <- tempfile(); writeLines(deparse(list(SIM_REPS, SIM_SCENARIOS, SIM_M1_COEF, SIM_M2_CANDIDATES, SIM_M2_IPW_TARGET, SIM_M2_E_CAP, SIM_M2_R2_CAP,
  SIM_U_PREVALENCE, SIM_U_SHIFT_PP, SIM_P1_BOUNDS, LEARNERS, K_FOLDS, PROPENSITY_CLIP, MIN_TREATED, MIN_CONTROLS, SEED_BASE)), .const)
FINGERPRINT <- paste(DATA_HASH, paste(spec, collapse = "/"), sprintf("%.6f", HEADLINE), unname(tools::md5sum(.const)))
if (file.exists(DRAWS) || file.exists(PARAMS)) { if (!file.exists(META) || readLines(META, n = 1) != FINGERPRINT)
  stop("60_sim_draws.csv and 60_dgp_params.rds belong to other inputs (data, learner pair, headline or constants): delete both to recompute") } else writeLines(FINGERPRINT, META)

## ---- 2) DGP: deterministic calibration (cached with the data hash), misspecification terms, U scenarios
recenter <- function(lp, target) uniroot(function(c) mean(plogis(lp + c)) - target, c(-10, 10), tol = 1e-10)$root
calibrate <- function() {
  lp_e <- as.numeric(predict(glmnet(Xm, t_real, family = "binomial", lambda = LEARNERS$logit$lambda), Xm, type = "link"))
  lp_m <- as.numeric(predict(glmnet(Xm[t_real == 0, ], y_real[t_real == 0], family = "binomial", lambda = LEARNERS$logit$lambda), Xm, type = "link"))
  tgt <- c(e = mean(t_real), m = mean(y_real[t_real == 0]))
  age <- Xm[, "log1p_age"]; age_c <- age - mean(age); deu <- Xm[, "countryDEU"]; pre <- Xm[, "pre_funding"]
  # M1: interaction and curvature the logit rungs cannot represent, coefficients from config
  z1 <- deu * age; z2 <- age_c^2
  mis_e <- SIM_M1_COEF$e$deu_x_age * z1 + SIM_M1_COEF$e$age_sq * z2
  mis_m <- SIM_M1_COEF$m$deu_x_age * z1 + SIM_M1_COEF$m$age_sq * z2
  P <- list(lp_e_well = lp_e + recenter(lp_e, tgt[["e"]]), lp_m_well = lp_m + recenter(lp_m, tgt[["m"]]),
            lp_e_mis = lp_e + mis_e + recenter(lp_e + mis_e, tgt[["e"]]), lp_m_mis = lp_m + mis_m + recenter(lp_m + mis_m, tgt[["m"]]),
            targets = tgt, n = n, data_hash = DATA_HASH)
  # M2: the first candidate term under the R2 cap, dampened so that main-effect logit IPW is biased by SIM_M2_IPW_TARGET
  asym <- function(lpe, lpm) {
    e <- plogis(lpe); p0 <- plogis(lpm)
    fe <- suppressWarnings(glm.fit(cbind(1, Xm), e, family = quasibinomial())); ehat <- plogis(as.numeric(cbind(1, Xm) %*% fe$coefficients))
    fm <- suppressWarnings(glm.fit(cbind(1, Xm), p0, weights = 1 - e, family = quasibinomial())); mhat <- plogis(as.numeric(cbind(1, Xm) %*% fm$coefficients))
    EY0_t <- sum(e * p0) / sum(e)
    c(ipw = EY0_t - sum((1 - e) * ehat / (1 - ehat) * p0) / sum((1 - e) * ehat / (1 - ehat)), or = EY0_t - sum(e * mhat) / sum(e), max_e = max(e)) }
  q75 <- quantile(age, 0.75)
  cands <- list(pre_funding_x_age = pre * age_c, deu_x_age_high = deu * as.numeric(age >= q75))
  for (nm in SIM_M2_CANDIDATES) {
    z <- cands[[nm]]; r2 <- summary(lm(z ~ Xm))$r.squared
    world <- function(s) { lpe <- P$lp_e_well + s * z; lpe <- lpe + recenter(lpe, tgt[["e"]]); lpm <- P$lp_m_well - s * z; lpm <- lpm + recenter(lpm, tgt[["m"]]); list(lpe = lpe, lpm = lpm) }
    if (r2 >= SIM_M2_R2_CAP) { cat(sprintf("M2 candidate %s: R2 on main effects %.3f violates the cap\n", nm, r2)); next }
    s_hat <- tryCatch(uniroot(function(s) { w <- world(s); asym(w$lpe, w$lpm)[["ipw"]] - SIM_M2_IPW_TARGET }, c(0.05, 8), tol = 1e-6)$root, error = function(e) NA_real_)
    mode <- "target"
    if (is.na(s_hat) || asym(world(s_hat)$lpe, world(s_hat)$lpm)[["max_e"]] > SIM_M2_E_CAP) {
      # the target is out of reach within the propensity cap: dampen to the cap and report the bias that is reachable
      s_hat <- uniroot(function(s) { w <- world(s); asym(w$lpe, w$lpm)[["max_e"]] - SIM_M2_E_CAP }, c(0.01, 8), tol = 1e-6)$root; mode <- "cap" }
    w <- world(s_hat); a <- asym(w$lpe, w$lpm)
    cat(sprintf("M2 candidate %s: s = %.4f (%s), R2 on main effects %.3f, asymptotic IPW %+.4f, OR %+.4f, max e %.3f -> accepted\n", nm, s_hat, mode, r2, a[["ipw"]], a[["or"]], a[["max_e"]]))
    P$m2 <- list(term = nm, mode = mode, s = s_hat, r2 = r2, asym_ipw = a[["ipw"]], asym_or = a[["or"]], max_e = a[["max_e"]]); P$lp_e_mis2 <- w$lpe; P$lp_m_mis2 <- w$lpm; break
  }
  if (is.null(P$m2)) stop("no M2 candidate stays under the R2 cap: revise SIM_M2_CANDIDATES or SIM_M2_R2_CAP in config.py")
  # U scenarios: U ~ Bernoulli(p_U) independent of X, logit-additive with strength r in assignment and outcome,
  # marginal rates held fixed; r solved so that the X-adjusted population bias makes the true ATT appear as the headline
  lp_e <- P$lp_e_well; lp_m <- P$lp_m_well; pu <- SIM_U_PREVALENCE
  solve_c <- function(lp, g, target) uniroot(function(c) (1 - pu) * mean(plogis(lp + c)) + pu * mean(plogis(lp + g + c)) - target, c(-10, 10), tol = 1e-10)$root
  rr_of <- function(lp, g, c) mean(plogis(lp + g + c)) / mean(plogis(lp + c))
  solve_gamma <- function(lp, target, r, sign) {
    g <- uniroot(function(g) { c <- solve_c(lp, g, target); rr <- rr_of(lp, g, c); (if (sign > 0) rr else 1 / rr) - r }, sort(c(1e-6, 8) * sign), tol = 1e-10)$root
    c <- solve_c(lp, g, target); list(g = g, c = c, rr = rr_of(lp, g, c)) }
  xbias <- function(gA, cA, gY, cY) {
    e1 <- plogis(lp_e + gA + cA); e0 <- plogis(lp_e + cA); eb <- pu * e1 + (1 - pu) * e0
    p1 <- plogis(lp_m + gY + cY); p0 <- plogis(lp_m + cY)
    sum(eb * (pu * e1 / eb - pu * (1 - e1) / (1 - eb)) * (p1 - p0)) / sum(eb) }
  dir_y <- if (HEADLINE < 0) -1 else 1                    # U raises assignment and pushes the outcome towards the headline's sign
  u_scen <- function(r, delta) { A <- solve_gamma(lp_e, tgt[["e"]], r, +1); Y <- solve_gamma(lp_m, tgt[["m"]], r, dir_y)
    rr_y <- if (dir_y < 0) 1 / Y$rr else Y$rr                      # both associations as RR >= 1 for the bounding factor
    bf <- A$rr * rr_y / (A$rr + rr_y - 1)
    list(delta = delta, u = TRUE, r = r, gA = A$g, cA = A$c, gY = Y$g, cY = Y$c, rr_a = A$rr, rr_y = rr_y, xbias = xbias(A$g, A$c, Y$g, Y$c),
         e_equiv = bf + sqrt(bf * (bf - 1))) }
  solve_r <- function(delta) uniroot(function(r) u_scen(r, delta)$xbias - (HEADLINE - delta), c(1.02, 15), tol = 1e-8)$root
  P$u <- list(rstar = u_scen(solve_r(0), 0))
  for (pp in SIM_U_SHIFT_PP) P$u[[paste0("shift", pp)]] <- u_scen(solve_r(pp / 100), pp / 100)
  P }
if (file.exists(PARAMS)) { P <- readRDS(PARAMS); if (!identical(P$data_hash, DATA_HASH)) stop("60_dgp_params.rds does not match analysis.csv: delete the cache to recalibrate") ; cat("DGP parameters loaded from cache\n") } else {
  P <- calibrate(); saveRDS(P, PARAMS); cat("DGP calibrated and cached ->", PARAMS, "\n") }
cat(sprintf("M2 term %s (s = %.3f) · U: r* = %.3f, r(shift5) = %.3f, r(shift10) = %.3f\n", P$m2$term, P$m2$s, P$u$rstar$r, P$u$shift5$r, P$u$shift10$r))

## ---- 3) Scenario table and replicate plan
SCEN <- list(S0 = list(delta = SIM_SCENARIOS$S0, mis = "well"), S1 = list(delta = SIM_SCENARIOS$S1, mis = "well"), S2 = list(delta = SIM_SCENARIOS$S2, mis = "well"),
             M1 = list(delta = SIM_SCENARIOS$M1, mis = "mis"), M2 = list(delta = SIM_SCENARIOS$M2, mis = "mis2"))
for (nm in names(P$u)) SCEN[[nm]] <- c(P$u[[nm]], list(mis = "well"))
plan <- expand.grid(rep = seq_len(SIM_REPS), scen = names(SCEN), stringsAsFactors = FALSE)[, c("scen", "rep")]
plan$seed <- mapply(function(s, r) seed_for(sprintf("60_%s_%04d", s, r)), plan$scen, plan$rep)

one_rep <- function(scen, rep, seed) {
  s <- SCEN[[scen]]; set.seed(seed)
  idx <- sample.int(n, replace = TRUE); Xr <- Xm[idx, , drop = FALSE]
  U <- if (isTRUE(s$u)) rbinom(n, 1, SIM_U_PREVALENCE) else 0
  lpe <- switch(s$mis, well = P$lp_e_well, mis = P$lp_e_mis, mis2 = P$lp_e_mis2)[idx]
  lpm <- switch(s$mis, well = P$lp_m_well, mis = P$lp_m_mis, mis2 = P$lp_m_mis2)[idx]
  gA <- if (is.null(s$gA)) 0 else s$gA; cA <- if (is.null(s$cA)) 0 else s$cA; gY <- if (is.null(s$gY)) 0 else s$gY; cY <- if (is.null(s$cY)) 0 else s$cY
  e_t <- plogis(lpe + gA * U + cA); p0 <- plogis(lpm + gY * U + cY)
  A <- rbinom(n, 1, e_t); p1 <- pmin(pmax(p0 + s$delta, SIM_P1_BOUNDS[1]), SIM_P1_BOUNDS[2]); Y <- rbinom(n, 1, ifelse(A == 1, p1, p0))
  true_att <- mean((p1 - p0)[A == 1])
  out <- function(est, att, se = NA_real_) data.frame(scen = scen, rep = rep, seed = seed, true_att = true_att, est = est, att = att, se = se)
  tryCatch({
    if (sum(A) < MIN_TREATED || sum(A == 0) < MIN_CONTROLS) stop("degenerate world")
    n1 <- sum(A); n0 <- n - n1; p1h <- mean(Y[A == 1]); p0h <- mean(Y[A == 0])
    RF_SEED <<- seed
    cf <- crossfit_fixed(Xr, A, Y, spec); a <- aipw(Y, A, cf$e, cf$m0)
    rbind(out("naive", p1h - p0h, sqrt(p1h * (1 - p1h) / n1 + p0h * (1 - p0h) / n0)), out("or", est_or(Y, A, Xr)), out("ipw", est_ipw(Y, A, Xr)), out("aipw", a[["att"]], a[["se"]]))
  }, error = function(e) out("aipw", NA_real_)) }

## ---- 4) Run with checkpoint
done <- if (file.exists(DRAWS)) read.csv(DRAWS, stringsAsFactors = FALSE) else data.frame(scen = character(), rep = integer())
todo <- plan[!paste(plan$scen, plan$rep) %in% unique(paste(done$scen, done$rep)), ]
cat(sprintf("simulation: %d of %d replicates open · %d cores\n", nrow(todo), nrow(plan), CORES))
while (nrow(todo) > 0) {
  batch <- todo[seq_len(min(CHUNK, nrow(todo))), ]; t0 <- Sys.time()
  res <- do.call(rbind, mclapply(seq_len(nrow(batch)), function(i) one_rep(batch$scen[i], batch$rep[i], batch$seed[i]), mc.cores = CORES, mc.set.seed = FALSE))
  write.table(res, DRAWS, sep = ",", row.names = FALSE, append = file.exists(DRAWS), col.names = !file.exists(DRAWS))
  todo <- todo[-seq_len(nrow(batch)), ]
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("[%s] chunk done (%.1f min) · open %d · remaining ~%.0f min\n", format(Sys.time(), "%H:%M"), el, nrow(todo), el * nrow(todo) / nrow(batch)))
}

## ---- 5) Performance per scenario and estimator with Monte-Carlo SEs; keys and tables
dr <- read.csv(DRAWS, stringsAsFactors = FALSE)
stopifnot(nrow(unique(dr[, c("scen", "rep")])) == nrow(plan), !anyDuplicated(dr[, c("scen", "rep", "est")]))   # every scenario holds SIM_REPS replicates
perf <- do.call(rbind, lapply(names(SCEN), function(sc) do.call(rbind, lapply(c("naive", "or", "ipw", "aipw"), function(es) {
  v <- dr[dr$scen == sc & dr$est == es & !is.na(dr$att), ]; R <- nrow(v); err <- v$att - v$true_att
  cov <- if (all(is.na(v$se))) NA else mean(abs(err) <= Z * v$se)
  data.frame(scen = sc, true_att = SCEN[[sc]]$delta, true_att_realized = mean(v$true_att), est = es, R = R, bias = mean(err), bias_mcse = sd(err) / sqrt(R), empse = sd(v$att), rmse = sqrt(mean(err^2)),
             coverage = cov, cov_mcse = if (is.na(cov)) NA else sqrt(cov * (1 - cov) / R)) }))))
vals <- c(); labels <- c(); put <- function(k, v, l) { vals[[k]] <<- v; labels[[k]] <<- l }
EST_LABEL <- c(naive = "naive", or = "outcome regression", ipw = "IPW", aipw = "AIPW")
put("s60_R", SIM_REPS, "replicates per scenario"); put("s60_fail", sum(is.na(dr$att)), "failed replicates over all scenarios")
for (i in seq_len(nrow(perf))) { r <- perf[i, ]; pre <- if (isTRUE(SCEN[[r$scen]]$u)) "s60_u_" else "s60_"
  put(sprintf("%sbias_%s_%s", pre, r$est, r$scen), r$bias, sprintf("bias of %s in scenario %s", EST_LABEL[[r$est]], r$scen))
  put(sprintf("%smcse_bias_%s_%s", pre, r$est, r$scen), r$bias_mcse, sprintf("Monte-Carlo SE of the bias of %s in scenario %s", EST_LABEL[[r$est]], r$scen))
  if (r$est == "aipw") put(sprintf("%sempse_aipw_%s", pre, r$scen), r$empse, sprintf("empirical SE of AIPW in scenario %s", r$scen))
  if (!is.na(r$coverage)) { put(sprintf("%scov_%s_%s", pre, r$est, r$scen), r$coverage, sprintf("coverage of the %s CI in scenario %s", EST_LABEL[[r$est]], r$scen))
    put(sprintf("%smcse_cov_%s_%s", pre, r$est, r$scen), r$cov_mcse, sprintf("Monte-Carlo SE of the coverage of %s in scenario %s", EST_LABEL[[r$est]], r$scen)) } }
put("s60_m2_term", P$m2$term, "M2 term used"); put("s60_m2_mode", P$m2$mode, "M2 dampening reached the IPW target or stopped at the propensity cap"); put("s60_m2_s", P$m2$s, "scale factor s of the M2 term (+s in assignment, -s in outcome)"); put("s60_m2_r2", P$m2$r2, "R squared of the M2 term on the main effects")
put("s60_m2_asym_ipw", P$m2$asym_ipw, "asymptotic bias of main-effect logit IPW under M2"); put("s60_m2_asym_or", P$m2$asym_or, "asymptotic bias of main-effect logit OR under M2"); put("s60_m2_max_e", P$m2$max_e, "largest true propensity under M2")
put("s60_u_p_u", SIM_U_PREVALENCE, "prevalence of the unobserved confounder"); put("s60_u_rstar", P$u$rstar$r, "U strength (RR with treatment and outcome) that makes a true zero appear as the headline")
for (nm in names(P$u)) { u <- P$u[[nm]]
  put(paste0("s60_u_r_", nm), u$r, paste("U strength in scenario", nm)); put(paste0("s60_u_true_att_", nm), u$delta, paste("true ATT in scenario", nm))
  put(paste0("s60_u_xbias_", nm), u$xbias, paste("X-adjusted population bias in scenario", nm)); put(paste0("s60_u_e_equiv_", nm), u$e_equiv, paste("E-value equivalent of the U strength in scenario", nm)) }
# realised true ATT: the bounds SIM_P1_BOUNDS bind for firms with a very low control probability, the bias above is measured against it
for (sc in names(SCEN)) put(sprintf("%strue_att_realized_%s", if (isTRUE(SCEN[[sc]]$u)) "s60_u_" else "s60_", sc), perf$true_att_realized[perf$scen == sc & perf$est == "aipw"],
                            sprintf("realised true ATT in scenario %s (mean over replicates)", sc))
write.csv(perf[!sapply(perf$scen, function(s) isTRUE(SCEN[[s]]$u)), ], file.path(OUT_TABLES, "60_simulation.csv"), row.names = FALSE)
write.csv(perf[sapply(perf$scen, function(s) isTRUE(SCEN[[s]]$u)), ], file.path(OUT_TABLES, "60_simulation_u.csv"), row.names = FALSE)
write.csv(data.frame(scenario = names(P$u), r = sapply(P$u, `[[`, "r"), true_att = sapply(P$u, `[[`, "delta"), rr_treatment = sapply(P$u, `[[`, "rr_a"), rr_outcome = sapply(P$u, `[[`, "rr_y"),
                     x_adjusted_bias = sapply(P$u, `[[`, "xbias"), e_equivalent = sapply(P$u, `[[`, "e_equiv")), file.path(OUT_TABLES, "60_dgp_params.csv"), row.names = FALSE)
write_values_csv(vals, labels, file.path(OUT_DATA, "60_values.csv"))
print(perf[perf$est == "aipw", c("scen", "true_att", "bias", "bias_mcse", "coverage")])
cat("wrote", length(vals), "keys -> 60_values.csv\n")
