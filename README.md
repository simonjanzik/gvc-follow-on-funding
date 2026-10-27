# Certification or Crowding Out? A Doubly Robust Re-Examination of Governmental Venture Capital and Private Follow-on Funding in European Startups

**Type:** Bachelor's Thesis

**Author:** Simon Janzik

**1st Examiner:** Prof. Dr. Stefan Lessmann

**2nd Examiner:** Prof. Dr. Jan Mendling

Humboldt-Universität zu Berlin, School of Business and Economics, Chair of Information Systems

## Table of Content

- [Summary](#summary)
- [Working with the repo](#working-with-the-repo)
    - [Dependencies](#dependencies)
    - [Setup](#setup)
- [Reproducing results](#reproducing-results)
    - [Estimation (stages 10 to 30)](#estimation-stages-10-to-30)
    - [Inference, sensitivity, simulation (stages 40 to 60)](#inference-sensitivity-simulation-stages-40-to-60)
    - [Precomputed draws](#precomputed-draws)
- [Results](#results)
- [Investor register](#investor-register)
- [Project structure](#project-structure)
- [Data and license](#data-and-license)

## Summary

Governments invest public money in early funding rounds of startups to finance new technologies and for private capital to follow. Theory offers two predictions. Under certification, private investors read the public selection as a signal of quality. Under crowding out, public money replaces private money or keeps it away. Most European evidence for certification rests on data that are not open and on central assumptions argued rather than tested. This thesis asks whether a public investor in a firm's first early-stage round, in place of a private one, changes the probability of a private follow-on round within 36 months. The data are an open export of the startup database Crunchbase: 3,022 European firms with a first round between 2005 and 2012, 606 of them with a public investor. Because the export does not mark public investors, I classify every investor by a fixed rule from documented facts and release the register. I estimate the effect with a doubly robust estimator, check it in a simulation, and measure how strong a hidden cause would need to be to explain it away. The estimated effect of a public first investor on the probability of a private follow-on round is −5.0 percentage points, with a 95 percent confidence interval of [−9.3, −0.8]. The certification premium of the earlier evidence does not appear against a private first investor. Crowding out is not shown: the interval excludes zero only narrowly, and a hidden cause as strong as a typical recorded firm characteristic could move it to zero. The thesis states what open data can say about certification and crowding out, why they cannot say more, and what the records of the public funders would add.

**Keywords**: governmental venture capital, certification, crowding out, follow-on funding, doubly robust estimation, AIPW, E-value, Crunchbase

**Full text**: not published at the moment.

## Working with the repo

### Dependencies

- macOS or Linux (the R stages run in parallel with forked processes, `run_all.sh` needs bash).
- Python 3.9 or later, standard library only (the pipeline uses no third-party Python packages). `run_all.sh` calls `/usr/bin/python3`, the interpreter the results were computed with; set `PY` in that script to use another one.
- R 4.4 with the packages `glmnet`, `gbm`, `ranger`, `grf` and `Matrix`.
- `curl` and `shasum` for the data download.

The exact versions the results were computed with are in [`scripts/ENVIRONMENT.md`](scripts/ENVIRONMENT.md).

### Setup

1. Clone this repository.

2. Download the raw data (about 72 MB, see [Data and license](#data-and-license)). The script pins the source to one commit and verifies every file by checksum:
```bash
sh get_data.sh
```

3. Install the R packages:
```r
install.packages(c("glmnet", "gbm", "ranger", "grf", "Matrix"))
```

## Reproducing results

The analysis is a pipeline of six numbered stages in `scripts/`, run in order by one script. Every empirical number in the thesis is a key in `outputs/values.json` (value, label, writing script, compute date) or a row of a stage table in `outputs/tables/`; the design constants are in `scripts/config.py`. The scripts that draw the thesis figures are not part of this repository; they only read these files.

```bash
cd scripts
./run_all.sh              # everything in order
./run_all.sh --from 30    # from stage 30 on
```

| Stage | What it does | Thesis |
|---|---|---|
| `10_build_data.py` | sample, register lookup, treatment, outcomes, covariates → `outputs/data/analysis.csv` | 3.1 to 3.3, A.1, A.2 |
| `20_descriptives.py` | descriptive table by arm, balance, outcome rates per arm, firms with and without links, exits | 3.1, 3.3, 4.1, 4.2, A.1 |
| `30_estimators.R` | estimator ladder, learner choice, cross-fitted AIPW and split band, variants, diagnostics | 3.7, 3.8, 4.1, 4.2, A.5 |
| `40_inference.R` | firm and cluster bootstrap, causal forest | 3.8, 4.1 |
| `50_sensitivity.py` | E-values, covariate benchmarks, minimum detectable effect, power | 3.8, 4.3 |
| `60_simulation.R` | simulation study with misspecification and unobserved-confounder scenarios | 3.9, 4.3, A.4 |

Each stage runs without arguments from any working directory, reads only the raw data, the inputs and earlier outputs, and writes only its own outputs. Seeds are derived from stable names (`seed_for(name)` in `scripts/_lib_estimators.R`) and stored with the output, so a rerun on the same software versions reproduces every key.

### Estimation (stages 10 to 30)

Stage 10 builds the analysis table from the raw export and the investor register (`investors.csv`): firms with a European early-stage first round in the cohorts 2005 to 2012, treatment (public investor in the first round), the outcomes, the covariates and a data dictionary. Stage 20 writes the descriptive tables. Stage 30 fits the ladder: naive contrast, outcome regression, inverse probability weighting, and AIPW with cross-fitting, where the learner for each nuisance function (logit, lasso, random forest, gradient boosting) is chosen by cross-fitted log-loss; it also runs the robustness variants in `variants.csv` and the overlap and balance diagnostics.

### Inference, sensitivity, simulation (stages 40 to 60)

Stage 40 cross-checks the influence-function standard errors with a firm bootstrap and a cluster bootstrap and fits a causal forest as a cross-check that the thesis does not report. Stage 50 computes E-values, benchmarks them against observed covariates, and reports minimum detectable effect and power. Stage 60 is the simulation study: data-generating processes calibrated to the sample, with scenarios for misspecified nuisance models and an unobserved confounder, which tests whether the estimators recover a known truth.

### Precomputed draws

Stages 40 (bootstrap) and 60 (simulation) are by far the slowest. Their draws are shipped in `outputs/data/` (`40_boot_draws.csv`, `60_sim_draws.csv`) together with a fingerprint of their inputs (`*.meta`). On the same software versions, both stages resume from these files and finish in minutes; if data, learner pair or constants have changed, the stage stops and asks for the draw file (for stage 60 also `60_dgp_params.rds`) to be deleted. On another R or package version, delete these files and let both stages recompute.

## Results

Average treatment effect on the treated (ATT) of a public first investor on the probability of a private follow-on round within 36 months, in percentage points, with the 95 percent confidence interval (rungs 1 to 3: bootstrap percentile; AIPW: influence function). Full table: `outputs/tables/30_ladder.csv`; key names in `outputs/values.json`.

| Rung | Total contrast (`*_total`) | Conditional on the number of investors (`*_cond`) |
|---|---|---|
| Naive contrast (`s30_att_naive_*`) | −6.7 [−10.0, −2.7] | −6.7 [−10.0, −2.7] |
| Outcome regression (`s30_att_or_*`) | −4.8 [−9.0, −1.1] | −5.8 [−9.9, −1.5] |
| IPW (`s30_att_ipw_*`) | −4.3 [−8.5, −0.2] | −5.9 [−10.4, −1.2] |
| AIPW, cross-fitted, headline (`s30_att_main_*`) | −5.0 [−9.3, −0.8] | −7.4 [−12.1, −2.8] |

E-value of the headline estimate: 1.76 for the point estimate, 1.23 for the confidence bound (`s50_evalue_point_total`, `s50_evalue_ci_total`).

## Investor register

The Crunchbase export does not mark public investors. `investors.csv` therefore classifies every investor that appears in a round of a sample firm, one row per Crunchbase permalink; persons are private. For each organization the register records documented facts, and a fixed rule derives the class from them, never from the name or the outcomes (thesis Section 3.2, Figure 3.1):

| Investor type (`investor_type`) | Second fact | Class (`class`) |
|---|---|---|
| program or data artifact (`programm`, `artefakt`) | – | not an investor (`kein_investor`) |
| public body, or a vehicle it set up and controls (`stelle`, `vehikel_oeff`) | instrument equity or mixed (`equity`, `gemischt`) | `public` |
| | instrument mezzanine (`mezzanine`) | `public_variant` |
| | instrument grant, loan or unknown (`grant_darlehen`, `unbekannt`) | `kein_investor` |
| privately raised fund with public capital, university, or owner unknown (`vehikel_priv_oeffkap`, `uni`, `unbekannt`) | public majority of capital (`public_capital` = `mehrheit`) | `public_variant` |
| | no public majority | `private` |
| bank, corporate, VC, angel or other, also if state-owned (`bank`, `corporate`, `vc`, `angel`, `andere`) | – | `private` |

A row that the rule makes `public` or `public_variant` becomes `private` if its public owner is outside the 29 countries of the sample (`reason` = `nicht_europaeisch`).

**Treatment.** A firm is treated if its first round has at least one `public` investor (headline definition). The broad definition also counts `public_variant`; the strict cut keeps only `public` investors with a public manager (`manager` = `oeffentlich`). Links of `kein_investor` drop out of every round, and `private` investors count for the outcome. Counts per class: `outputs/tables/10_register_by_class.csv`; all `public` and `public_variant` investors with first-round links: `outputs/tables/10_public_investors.csv`.

**How it was built.** An AI assistant classified every organization, a second pass with sources covered contested rows, the result was checked against a list of public early-stage vehicles, and I decided each boundary case on its documented facts (thesis Appendix A.2). `checked_by` names the step that classified each row, `confidence` and `source` give the evidence.

**Columns.** Values are German as in the register; the rule above and this table translate them. Persons are private by rule (`investor_type` = `angel`, `reason` = `person`, `confidence` = `belegt`); their other fact columns and `source` are empty.

| Column | Content |
|---|---|
| `permalink`, `name`, `type`, `country` | Crunchbase identifier, name, `organization` or `person`, country code from the export; 56 names that the export double-encodes (for example `GestiÃ³n` for `Gestión`) are repaired |
| `entity` | legal entity behind the permalink, as researched; some rows add a short note in German, for example a later name or the managing company |
| `investor_type`, `public_capital`, `instrument` | the three facts; besides the values above, `public_capital` takes `anker` (anchor minority), `keins` (none) and `unbekannt` |
| `manager` | `oeffentlich` public, `privat` private, `unbekannt`; used only for the strict cut |
| `class`, `reason` | class from the rule, and one word for the reason: the decisive fact or, where the owner type settles the class, that type (for example `accelerator`, `lp_only`, `nicht_europaeisch`) |
| `confidence` | `belegt` primary source, `sekundaer` secondary source, `kenntnis` from knowledge, `unbekannt` |
| `checked_by` | the step whose reading the row keeps: `blind_web` first classification with web research, `blind_kenntnis` first classification from the assistant's knowledge (both without sight of an earlier name-based class), `verifikation_web` second pass over contested rows with sources, `konsolidierung` consolidation (one organization under several permalinks, rule updates, direct lookups), `auto_person` person, private by rule |
| `source` | URLs of the sources, separated by ` \| ` |

## Project structure

```bash
├── README.md
├── get_data.sh                   -- downloads and verifies the raw data into raw/cb-2015/
├── investors.csv                 -- investor register: one row per investor, class from documented facts
├── variants.csv                  -- robustness variants: filter, treatment and outcome column, status
├── raw/cb-2015/                  -- raw data, not in the repository (get_data.sh)
├── scripts
│   ├── run_all.sh                -- the six stages in order
│   ├── 10_build_data.py … 60_simulation.R
│   ├── config.py                 -- every design constant, single source; config.R is generated from it
│   ├── _lib.py, _lib_estimators.R -- shared helpers; the estimators the R stages share are in _lib_estimators.R
│   ├── _values.py                -- the single writer of outputs/values.json
│   ├── _env.py, ENVIRONMENT.md   -- software versions behind the outputs (run_all.sh rewrites it)
└── outputs
    ├── values.json               -- every reported number as a key
    ├── tables/                   -- stage tables (CSV)
    └── data/                     -- analysis table, data dictionary, follow-on rounds, per-stage values, bootstrap and
                                     simulation draws, simulation cache
```

## Data and license

**Raw data.** Crunchbase export of December 2015 (`companies`, `rounds`, `investments`, `acquisitions`), published at [github.com/notpeter/crunchbase-data](https://github.com/notpeter/crunchbase-data) and pinned here to commit `be8711f`. Crunchbase provides the dataset under the Creative Commons Attribution-NonCommercial license (CC BY-NC); the license text is downloaded with the data (`raw/cb-2015/LICENSE.txt`). The data are used for non-commercial research with attribution to Crunchbase and are not redistributed in this repository.

**Code and register.** The scripts and the investor register were written for this thesis and may be used for non-commercial research with attribution.
