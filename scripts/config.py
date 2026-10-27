#!/usr/bin/env python3
"""Design constants of the pipeline. Single source of truth; config.R is generated from this file.

Purpose : every design constant that stages 10 to 60 use, grouped by topic. No estimation, no data access.
          Changing a constant changes the analysis: regenerate config.R and rerun from the first stage that reads it.
Thesis  : Sections 3.1 to 3.3 (sample, register, covariates), 3.7 to 3.9 (estimation, inference and sensitivity,
          simulation); Appendix A.3 (learner settings), A.4, A.5
Input   : none
Output  : config.R: the same constants in R syntax, do not edit by hand
Usage   : python3 config.py   (writes config.R; the Python stages import this file, the R stages read config.R)
"""
from pathlib import Path

# ---- 1) Sample: which firms and which rounds are in scope
COUNTRIES_EUROPE = sorted({"DEU", "FRA", "GBR", "ITA", "ESP", "NLD", "SWE", "DNK", "FIN", "NOR", "POL", "AUT",
                           "BEL", "CHE", "IRL", "PRT", "CZE", "HUN", "ROU", "GRC", "LUX", "EST", "LVA", "LTU",
                           "SVK", "SVN", "BGR", "HRV", "ISL"})   # EU without CYP and MLT, plus GBR, CHE, NOR, ISL: 29 codes
COHORT_FIRST, COHORT_LAST = 2005, 2012          # year of the first round
FIRST_ROUND_TYPES = ["seed", "angel"]           # qualify as first early-stage round unconditionally
FIRST_ROUND_VENTURE_CODES = ["A", ""]           # a venture round qualifies with code A or without code
EQUITY_TYPES = ["seed", "angel", "venture", "private_equity", "convertible_note", "equity_crowdfunding"]
EQUITY_TYPES_UNDISCLOSED = EQUITY_TYPES + ["undisclosed"]   # robustness column o1_undisc_36
WINDOW_DAYS = {"36": 1096, "24": 731, "48": 1461}           # months -> days
WINDOW_48_LAST_COHORT = 2011                                 # o1_48 is empty for later cohorts

# ---- 2) Register interface: what stage 10 reads from investors.csv
REGISTER_COLUMNS = ["permalink", "class", "manager", "name", "country", "entity", "confidence"]   # treatment, investor table, keys
REGISTER_FACT_COLUMNS = ["type", "investor_type", "public_capital", "instrument", "reason"]      # register keys and data flags
CLASSES_ALL = ["public", "public_variant", "kein_investor", "private"]
CLASS_PUBLIC_HEADLINE = ["public"]                           # treated
CLASS_PUBLIC_BROAD = ["public", "public_variant"]            # treated_broad
STRICT_MANAGER = "oeffentlich"                               # treated_strict = public with this manager value
CLASS_PRIVATE = ["private"]                                  # private link for O1 and for solo/syndicated
CLASS_IGNORED = ["kein_investor"]                            # link is dropped from every round
REGISTER_UNRESOLVED_VALUE = "unbekannt"                      # investor_type value counted as unresolved

# ---- 3) Covariates and filters
SECTOR_TOP_K = 25                                            # most frequent first categories, rest "other"
SECTOR_OTHER, SECTOR_UNKNOWN = "other", "unknown"
AGE_MIN_YEARS, AGE_MAX_YEARS = -0.5, 50.0                    # age at first round outside this range = missing
AGE_DECIMALS = 2                                             # age_at_round rounded in analysis.csv
CAP_FIRST_ROUND_USD = 10_000_000                             # variant cap_first_round
CAP_AGE_YEARS = 7.0

# ---- 4) Estimation and inference. A reading is one covariate specification: "total" gives the total effect (headline),
#         "cond" the conditional contrast given the investor count; key suffixes _total and _cond follow it
READINGS = {"total": "covariates without investor count (causal headline)",
            "cond": "covariates with log1p investor count (conditional contrast, fixed before the first estimate)"}
PROPENSITY_CLIP = [0.001, 0.99]                              # bounds of clip(): propensities and outcome-model predictions
K_FOLDS = 5                                                  # cross-fitting folds, drawn separately per arm
N_SPLITS = 100                                               # cross-fitting splits for the band
B_BOOTSTRAP = 999                                            # firm-level (both readings) and cluster (total)
CI_LEVEL = 0.95                                              # IF-CI = ATT +- z * SE, bootstrap CI = percentiles
MIN_TREATED = 30                                             # guard: a (sub)sample or replicate with fewer treated is invalid
MIN_CONTROLS = 100                                           # guard for bootstrap and simulation worlds
GRF_NUM_TREES = 2000                                         # causal forest cross-check
# learner choice: smallest cross-fitted log-loss, once per reading and nuisance, kept fixed for all variants
LEARNERS = {"logit": {"lambda": 1e-6},                       # glmnet binomial, practically unregularised
            "lasso": {"nfolds": 5},                          # cv.glmnet, lambda.min
            "rf": {"num_trees": 500},                        # ranger probability forest
            "gbm": {"n_trees": 300, "interaction_depth": 3, "shrinkage": 0.05}}
ROUND_DECIMALS = 6                                           # values.json rounding

# ---- 5) Sensitivity
EVALUE_SHIFT_PP = [5, 10]                                    # hypothetical true ATT of +5 and +10 pp (counterfactual shifted, treated rate observed)
POWER_TARGETS = [0.80, 0.90]                                 # MDE at these power levels
POWER_GRID_PP = [5, 8, 10, 15, 20]                           # power against hypothetical true effects

# ---- 6) Simulation: ADEMP on the headline covariates (reading "total")
SIM_REPS = 500
SIM_SCENARIOS = {"S0": 0.00, "S1": -0.06, "S2": 0.10, "M1": -0.06, "M2": -0.06}   # true ATT per scenario
SIM_M1_COEF = {"e": {"deu_x_age": 1.0, "age_sq": -0.6}, "m": {"deu_x_age": -0.8, "age_sq": 0.5}}
SIM_M2_CANDIDATES = ["pre_funding_x_age", "deu_x_age_high"]  # first term under the R2 cap is used (s60_m2_term); dampened to the IPW
#   target if reachable within the propensity cap, otherwise to the cap (s60_m2_mode, achieved bias s60_m2_asym_ipw)
#   pre_funding_x_age: pre_funding * centred log1p age; deu_x_age_high: country DEU * 1{log1p age >= upper quartile}
SIM_M2_IPW_TARGET = -0.06                                    # asymptotic IPW bias the M2 dampening is tuned to
SIM_M2_E_CAP, SIM_M2_R2_CAP = 0.95, 0.30                     # constraints of the M2 term (max propensity, R2 on main effects)
SIM_U_PREVALENCE = 0.30                                      # P(U = 1), U independent of X
SIM_U_SHIFT_PP = [5, 10]                                     # U scenarios shift5, shift10: true ATT +5 and +10 pp, U strength calibrated to land on the headline
SIM_P1_BOUNDS = [0.001, 0.999]                               # the treated outcome probability p0 + true ATT is kept inside these bounds,
#   so the realised true ATT can differ slightly from the nominal one (keys s60_true_att_realized_<scenario>, for the U scenarios s60_u_true_att_realized_<scenario>)

# ---- 7) Seeds and paths: one base, the unit name makes the seed (stage, variant, reading, replicate)
SEED_BASE = 2026                                             # seed_for(name) in _lib_estimators.R: fold of the name's code points
OUTPUT_DIR = "outputs"                                      # results folder next to scripts/

# ---- 8) Writer for config.R: keeps both languages on one source
def _r(v):
    if isinstance(v, bool): return "TRUE" if v else "FALSE"
    if isinstance(v, (int, float)): return repr(v)
    if isinstance(v, str): return '"' + v.replace('"', '\\"') + '"'
    if isinstance(v, (list, tuple)): return "c(" + ", ".join(_r(x) for x in v) + ")"
    if isinstance(v, dict): return "list(" + ", ".join(f"`{k}` = {_r(x)}" for k, x in v.items()) + ")"
    raise TypeError(type(v))

def write_r(path=None):
    path = Path(path or Path(__file__).with_name("config.R"))
    lines = ["# Generated by config.py, do not edit. Design constants, see config.py for comments."]
    for k, v in globals().items():
        if k.isupper() and not k.startswith("_"):
            lines.append(f"{k} <- {_r(v)}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path

if __name__ == "__main__":
    print("wrote", write_r())
