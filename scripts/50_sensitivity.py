#!/usr/bin/env python3
"""Stage 50: sensitivity to unobserved confounding and detectability, both readings.

Purpose : E-values (point estimate, CI bound nearer to zero, hypothetical true ATT of +5 and +10 pp), covariate
          benchmarks in the same convention, minimum detectable effect and power against hypothetical effects.
Thesis  : Sections 3.8 (E-values, covariate benchmarks, minimum detectable effect, power), 4.3 (results), 5.1
Input   : outputs/values.json (s20_o1_36_rate_treated, s30_att_main_*, s30_ci_*, s30_se_main_*, s30_nt_main_*,
          s30_ess_main_*), outputs/data/analysis.csv (benchmarks), config.py
Output  : outputs/tables/50_evalues.csv, 50_covariate_benchmarks.csv, 50_power.csv; keys s50_* in values.json
Usage   : python3 50_sensitivity.py   (no arguments; run_all.sh runs it as stage 50)
Seeds   : none
Checks  : stage 30 keys present; E-value convention: bias factor from the observed treated rate and the AIPW
          counterfactual, the sign of the ATT sets the direction, a CI that covers zero gives E = 1
"""
import math, statistics
from statistics import NormalDist
import config
from _lib import OUT_DATA, OUT_TABLES, read_csv, write_csv, require
import _values

ND = NormalDist(); Z = ND.inv_cdf(1 - (1 - config.CI_LEVEL) / 2)

def evalue(b):
    """bias factor B >= 1 -> smallest joint RR association of an unobserved confounder that could explain it"""
    return 1.0 if b is None or b <= 1 else b + math.sqrt(b * (b - 1))

def bias_factor(p1, p0_obs, true_att):
    """how much confounding turns a true ATT into the observed one: RR_obs / RR_true on the treated rate p1,
    inverted to >= 1 so that the direction follows from the sign"""
    p0_true = p1 - true_att
    if p0_true <= 0 or p0_obs <= 0: return None
    x = p0_true / p0_obs
    return max(x, 1 / x)

def main():
    vals = _values.load(); g = lambda k: vals[k]["value"]
    require("s30_att_main_total" in vals, "stage 30 keys missing: run 30_estimators.R and the merge first")
    p1 = g("s20_o1_36_rate_treated"); K = {}; ev_rows = []; pw_rows = []

    # ---- 1) E-values per reading: point, CI bound nearer to zero, hypothetical true ATT of +5 and +10 pp
    for r in config.READINGS:
        att, lo, hi, se = g(f"s30_att_main_{r}"), g(f"s30_ci_lo_main_{r}"), g(f"s30_ci_hi_main_{r}"), g(f"s30_se_main_{r}")
        p0 = p1 - att
        K[f"s50_rr_obs_{r}"] = (round(max(p1 / p0, p0 / p1), 6), f"observed risk ratio treated versus counterfactual (inverted to >= 1), reading {r}")
        e_point = evalue(bias_factor(p1, p0, 0.0))
        bound = hi if att < 0 else lo
        e_ci = 1.0 if lo <= 0 <= hi else evalue(bias_factor(p1, p1 - bound, 0.0))
        K[f"s50_evalue_point_{r}"] = (round(e_point, 6), f"E-value of the point estimate, reading {r}")
        K[f"s50_evalue_ci_{r}"] = (round(e_ci, 6), f"E-value of the CI bound nearer to zero (1 if the CI covers zero), reading {r}")
        ev_rows.append({"reading": r, "scenario": "point estimate versus zero", "evalue": round(e_point, 4)})
        ev_rows.append({"reading": r, "scenario": "CI bound versus zero", "evalue": round(e_ci, 4)})
        for pp in config.EVALUE_SHIFT_PP:
            e = evalue(bias_factor(p1, p0, pp / 100))
            K[f"s50_evalue_shift{pp}_{r}"] = (round(e, 6), f"E-value for a true ATT of +{pp} pp appearing as the observed one, reading {r}")
            ev_rows.append({"reading": r, "scenario": f"true ATT +{pp} pp", "evalue": round(e, 4)})

        # ---- 2) Detectability: SE decomposition, MDE, power against hypothetical effects (never observed power)
        n1, ess = g(f"s30_nt_main_{r}"), g(f"s30_ess_main_{r}")
        se_approx = math.sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / ess)
        K[f"s50_se_approx_ratio_{r}"] = (round(se_approx / se, 6), f"two-proportion SE approximation over the IF-SE, reading {r}")
        for tgt in config.POWER_TARGETS:
            K[f"s50_mde{int(tgt*100)}_{r}"] = (round((Z + ND.inv_cdf(tgt)) * se, 6), f"minimum detectable effect at {int(tgt*100)} percent power, as a share (0.06 = 6 pp), reading {r}")
        for pp in config.POWER_GRID_PP:
            lam = pp / 100 / se; pw = (1 - ND.cdf(Z - lam)) + ND.cdf(-Z - lam)
            K[f"s50_power_plus{pp}pp_{r}"] = (round(pw, 6), f"power against a true effect of {pp} pp, reading {r}")
            pw_rows.append({"reading": r, "true_effect_pp": pp, "power": round(pw, 4)})

    # ---- 3) Covariate benchmarks: unadjusted RR of each binary covariate with treatment and with the control outcome,
    #         bounding factor and E-equivalent (Ding and VanderWeele), compared with the E-values above
    d = [r for r in read_csv(OUT_DATA / "analysis.csv") if r["has_links"] == "1"]
    A = [int(r["treated"]) for r in d]; Y = [int(r["o1_36"]) for r in d]
    ages = [float(r["age_at_round"]) for r in d if r["age_at_round"] != ""]; med_age = statistics.median(ages)
    binaries = {"age_missing": [int(r["age_missing"]) for r in d], "pre_funding": [int(r["pre_funding"]) for r in d],
                "age_above_median": [int(r["age_at_round"] != "" and float(r["age_at_round"]) >= med_age) for r in d],
                "links_ge2": [int(int(r["n_links_first"]) >= 2) for r in d]}
    for kind in ("country", "sector", "cohort"):
        for lvl in sorted({r[kind] for r in d}):
            binaries[f"{kind}_{lvl}"] = [int(r[kind] == lvl) for r in d]
    rows = []
    for name, x in binaries.items():
        n1x, n0x = sum(x), len(x) - sum(x)
        if n1x < 20 or n0x < 20: continue                      # levels with fewer than 20 firms on either side are skipped
        pa1 = sum(a for a, xi in zip(A, x) if xi) / n1x; pa0 = sum(a for a, xi in zip(A, x) if not xi) / n0x
        c1 = [y for y, a, xi in zip(Y, A, x) if xi and not a]; c0 = [y for y, a, xi in zip(Y, A, x) if not xi and not a]
        if not c1 or not c0 or pa1 == 0 or pa0 == 0 or statistics.mean(c1) == 0 or statistics.mean(c0) == 0: continue   # zero rate: no ratio
        rr_a = max(pa1 / pa0, pa0 / pa1); rr_y = max(statistics.mean(c1) / statistics.mean(c0), statistics.mean(c0) / statistics.mean(c1))
        bf = rr_a * rr_y / (rr_a + rr_y - 1)
        rows.append({"covariate": name, "share": round(n1x / len(x), 4), "rr_treatment": round(rr_a, 4), "rr_outcome_controls": round(rr_y, 4), "bounding_factor": round(bf, 4), "e_equivalent": round(evalue(bf), 4)})
    rows.sort(key=lambda r: -r["e_equivalent"])
    write_csv(OUT_TABLES / "50_covariate_benchmarks.csv", rows)
    K["s50_bench_max_e_equiv"] = (rows[0]["e_equivalent"], f"largest E-equivalent among the observed covariates ({rows[0]['covariate']})")
    K["s50_bench_n_vars"] = (len(rows), "covariates in the benchmark table")
    for r in config.READINGS:
        K[f"s50_bench_n_vars_e_ge_evalue_point_{r}"] = (sum(1 for x in rows if x["e_equivalent"] >= K[f"s50_evalue_point_{r}"][0]), f"covariates whose E-equivalent reaches the point E-value, reading {r}")
        K[f"s50_bench_n_vars_e_ge_evalue_ci_{r}"] = (sum(1 for x in rows if x["e_equivalent"] >= K[f"s50_evalue_ci_{r}"][0]), f"covariates whose E-equivalent reaches the CI E-value, reading {r}")
    write_csv(OUT_TABLES / "50_evalues.csv", ev_rows); write_csv(OUT_TABLES / "50_power.csv", pw_rows)
    _values.update(K, "50_sensitivity.py")
    print(" · ".join(f"{r}: E point {K[f's50_evalue_point_{r}'][0]:.2f} CI {K[f's50_evalue_ci_{r}'][0]:.2f} MDE80 {100*K[f's50_mde80_{r}'][0]:.1f} pp" for r in config.READINGS))

if __name__ == "__main__":
    main()
