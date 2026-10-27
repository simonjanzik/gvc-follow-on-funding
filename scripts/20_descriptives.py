#!/usr/bin/env python3
"""Stage 20: descriptives of the sample, no effects.

Purpose : rates per arm, the descriptive table by arm, unweighted balance, O1b decomposition, links-versus-no-links
          comparison, exits, amount fill rates, investor counts per arm, firms per cohort and arm, and the provenance
          markers of the export.
Thesis  : Sections 3.1 (no-link firms, country-code loss), 3.3 (Table 3.2, Figure 3.2, outcome rates), 4.1 (naive
          difference), 4.2 (O1b decomposition, Figure 4.3), 5.3 (no-link firms); Appendix A.1 (provenance markers)
Input   : outputs/values.json (s10_firms), outputs/data/analysis.csv, outputs/data/10_followons.csv, investors.csv
          (type and class of the investor), raw/cb-2015/{rounds,investments}.csv (provenance markers, country-code loss)
Output  : outputs/tables/20_table_by_arm.csv, 20_balance.csv, 20_links_vs_nolinks.csv, 20_o1b_decomposition.csv;
          keys s20_* in values.json
Usage   : python3 20_descriptives.py   (no arguments; run_all.sh runs it as stage 20)
Seeds   : none
Checks  : row count of analysis.csv equals s10_firms
"""
import bisect, math, statistics
from collections import Counter, defaultdict
import config
from _lib import OUT_DATA, OUT_TABLES, REGISTER, read_csv, write_csv, load_rows, parse_date, require
import _values

def load_register_maps():
    """permalink -> type and class of the investor"""
    reg = read_csv(REGISTER)
    return {r["permalink"]: r["type"] for r in reg}, {r["permalink"]: r["class"] for r in reg}

def rate(rows, key): return round(sum(int(r[key]) for r in rows) / len(rows), 6) if rows else ""
def smd(x1, x0):   # population variances; stage 30 uses sample variances, so its unweighted SMDs differ in the fourth digit
    s = math.sqrt((statistics.pvariance(x1) + statistics.pvariance(x0)) / 2) if len(x1) > 1 and len(x0) > 1 else 0
    return 0.0 if s == 0 else (statistics.mean(x1) - statistics.mean(x0)) / s

def main():
    vals = _values.load()
    d = read_csv(OUT_DATA / "analysis.csv")
    require(len(d) == vals["s10_firms"]["value"], "analysis.csv rows differ from s10_firms")
    for r in d:
        for k in ("treated", "has_links", "age_missing", "pre_funding", "o1_36", "o1_24", "o1b_36", "o1_undisc_36", "o2_36", "exit_acq_36", "amount_missing", "treated_synd", "treated_solo"):
            r[k] = int(r[k])
        r["n_links_first"] = int(r["n_links_first"]); r["cohort"] = int(r["cohort"])
    prim = [r for r in d if r["has_links"]]; nolink = [r for r in d if not r["has_links"]]
    tr = [r for r in prim if r["treated"]]; ct = [r for r in prim if not r["treated"]]
    K = {}

    # ---- 1) Rates per arm and sample shares
    for o in ("o1_36", "o1_24", "o1b_36", "o1_undisc_36", "o2_36"):
        K[f"s20_{o}_rate"] = (rate(prim, o), f"{o} rate in the primary sample")
        K[f"s20_{o}_rate_treated"] = (rate(tr, o), f"{o} rate among treated"); K[f"s20_{o}_rate_control"] = (rate(ct, o), f"{o} rate among controls")
    p48 = [r for r in prim if r["o1_48"] != ""]
    K["s20_o1_48_rate"] = (round(sum(int(r["o1_48"]) for r in p48) / len(p48), 6), f"o1_48 rate in the primary sample, cohorts up to {config.WINDOW_48_LAST_COHORT}")
    K["s20_o1_36_diff_naive"] = (round(rate(tr, "o1_36") - rate(ct, "o1_36"), 6), "naive difference of o1_36 rates, treated minus controls")
    K["s20_o1b_36_rate_by_cohort"] = ({str(c): rate([r for r in prim if r["cohort"] == c], "o1b_36") for c in sorted({r["cohort"] for r in prim})}, "o1b_36 rate per cohort (censoring check)")
    K["s20_pre_funding_share"] = (rate(prim, "pre_funding"), "share with a round before the first round")
    K["s20_age_missing_share"] = (rate(prim, "age_missing"), "share with missing age")
    K["s20_sector_top5"] = (dict(Counter(r["sector"] for r in prim).most_common(5)), "five most frequent sectors in the primary sample")
    K["s20_exit_acq_36_treated"] = (rate(tr, "exit_acq_36"), "acquired within 36 months, treated"); K["s20_exit_acq_36_control"] = (rate(ct, "exit_acq_36"), "acquired within 36 months, controls")
    cohorts = sorted({r["cohort"] for r in prim})   # firms per cohort and arm (primary sample)
    K["s20_firms_by_cohort_treated"] = ({str(c): sum(1 for r in tr if r["cohort"] == c) for c in cohorts}, "treated firms per cohort, primary sample")
    K["s20_firms_by_cohort_control"] = ({str(c): sum(1 for r in ct if r["cohort"] == c) for c in cohorts}, "control firms per cohort, primary sample")

    # ---- 2) Investor count per arm and control outcome by investor count
    K["s20_nlinks_mean_treated"] = (round(statistics.mean(r["n_links_first"] for r in tr), 6), "mean first-round links, treated")
    K["s20_nlinks_mean_control"] = (round(statistics.mean(r["n_links_first"] for r in ct), 6), "mean first-round links, controls")
    K["s20_nlinks_mean_treated_synd"] = (round(statistics.mean(r["n_links_first"] for r in tr if r["treated_synd"]), 6), "mean first-round links, syndicated treated")
    K["s20_nlinks_mean_treated_solo"] = (round(statistics.mean(r["n_links_first"] for r in tr if r["treated_solo"]), 6), "mean first-round links, treated without a private link")
    for name, lab, keep in (("1", "one link", lambda n: n == 1), ("2", "two links", lambda n: n == 2), ("3", "three links", lambda n: n == 3), ("4plus", "four or more links", lambda n: n >= 4)):
        cell = [r for r in ct if keep(r["n_links_first"])]   # outcome of controls by investor count, with the cell size n so that intervals can be computed
        K[f"s20_o1_36_rate_control_nlinks{name}"] = (rate(cell, "o1_36"), f"o1_36 rate of controls with {lab}")
        K[f"s20_n_control_nlinks{name}"] = (len(cell), f"controls with {lab}")

    # ---- 3) Descriptive table by arm and unweighted balance: the covariates of stage 30 with every dummy level and the
    #         log1p investor count, no weights
    def age_imputed(rows):
        med = {c: statistics.median([float(r["age_at_round"]) for r in rows if r["cohort"] == c and r["age_at_round"] != ""] or [0]) for c in {r["cohort"] for r in rows}}
        return [math.log1p(float(r["age_at_round"]) if r["age_at_round"] != "" else med[r["cohort"]]) for r in rows]
    cols = {"log1p_age": age_imputed(prim), "age_missing": [r["age_missing"] for r in prim], "pre_funding": [r["pre_funding"] for r in prim],
            "log1p_nlinks": [math.log1p(r["n_links_first"]) for r in prim]}
    for kind in ("country", "sector", "cohort"):
        for lvl in sorted({str(r[kind]) for r in prim}):
            cols[f"{kind}{lvl}"] = [int(str(r[kind]) == lvl) for r in prim]
    ti = [i for i, r in enumerate(prim) if r["treated"]]; ci = [i for i, r in enumerate(prim) if not r["treated"]]
    bal = [{"column": c, "mean_treated": round(statistics.mean(v[i] for i in ti), 6), "mean_control": round(statistics.mean(v[i] for i in ci), 6),
            "smd_unweighted": round(smd([v[i] for i in ti], [v[i] for i in ci]), 6)} for c, v in cols.items()]
    write_csv(OUT_TABLES / "20_balance.csv", bal)
    K["s20_balance_max_abs_smd"] = (round(max(abs(b["smd_unweighted"]) for b in bal), 6), "largest absolute unweighted SMD over the covariates (every dummy level and the log1p investor count)")
    def block(rows, label):
        ages = [float(r["age_at_round"]) for r in rows if r["age_at_round"] != ""]
        amts = [float(r["first_round_amount_usd"]) for r in rows if r["first_round_amount_usd"] != ""]
        return {"group": label, "n": len(rows), "o1_36_rate": rate(rows, "o1_36"), "o1b_36_rate": rate(rows, "o1b_36"), "mean_links": round(statistics.mean(r["n_links_first"] for r in rows), 3),
                "pre_funding_share": rate(rows, "pre_funding"), "age_median_years": round(statistics.median(ages), 2) if ages else "", "age_missing_share": rate(rows, "age_missing"),
                "amount_median_usd": round(statistics.median(amts)) if amts else "", "amount_missing_share": rate(rows, "amount_missing"),
                "share_seed": round(sum(r["first_round_type"] == "seed" for r in rows) / len(rows), 4), "share_angel": round(sum(r["first_round_type"] == "angel" for r in rows) / len(rows), 4),
                "share_DEU": round(sum(r["country"] == "DEU" for r in rows) / len(rows), 4), "share_GBR": round(sum(r["country"] == "GBR" for r in rows) / len(rows), 4), "share_FRA": round(sum(r["country"] == "FRA" for r in rows) / len(rows), 4),
                "cohort_mean": round(statistics.mean(r["cohort"] for r in rows), 2)}
    write_csv(OUT_TABLES / "20_table_by_arm.csv", [block(tr, "treated"), block(ct, "control"), block(prim, "primary_sample"), block(nolink, "no_link")])

    # ---- 4) Links versus no-link firms: the excluded firms differ, documented not hidden
    K["s20_nolink_share"] = (round(len(nolink) / len(d), 6), "share of firms without a first-round link")
    K["s20_nolink_o1b_gap"] = (round(rate(prim, "o1b_36") - rate(nolink, "o1b_36"), 6), "o1b_36 rate, linked minus no-link firms")
    K["s20_nolink_o1b_36_rate"] = (rate(nolink, "o1b_36"), "o1b_36 rate of no-link firms")
    K["s20_nolink_age_missing_share"] = (rate(nolink, "age_missing"), "share of no-link firms with missing age")
    K["s20_nolink_share_DEU"] = (round(sum(r["country"] == "DEU" for r in nolink) / len(nolink), 6), "share of no-link firms in Germany")
    amt = lambda rows: round(statistics.median([float(r["first_round_amount_usd"]) for r in rows if r["first_round_amount_usd"] != ""]))
    K["s20_nolink_amount_median_usd_nolink"] = (amt(nolink), "median first-round amount of no-link firms, USD")
    K["s20_nolink_amount_median_usd_control"] = (amt(ct), "median first-round amount of linked controls, USD")
    write_csv(OUT_TABLES / "20_links_vs_nolinks.csv", [block(prim, "linked"), block(nolink, "no_link")])

    # ---- 5) O1b decomposition: firms with an equity follow-on but no private link, by what the round looked like;
    #         "public only" means links of class public or public_variant only
    fo = defaultdict(list)
    for r in read_csv(OUT_DATA / "10_followons.csv"):
        if int(r["equity"]) and int(r["took_place"]) and int(r["days_since_first"]) <= config.WINDOW_DAYS["36"]: fo[r["firm_id"]].append(r)
    dec = []
    for arm, rows in (("treated", tr), ("control", ct)):
        only = [r for r in rows if r["o1b_36"] and not r["o1_36"]]
        nol = sum(1 for r in only if all(int(x["n_links"]) == 0 for x in fo[r["firm_id"]]))
        pub = len(only) - nol
        dec.append({"arm": arm, "n_arm": len(rows), "n_o1b_only": len(only), "share_o1b_only": round(len(only) / len(rows), 6), "follow_on_without_links": nol, "follow_on_public_only": pub})
        K[f"s20_o1bonly_n_arm_{arm}"] = (len(rows), f"firms in arm {arm}"); K[f"s20_o1bonly_n_o1bonly_{arm}"] = (len(only), f"firms with o1b_36 = 1 and o1_36 = 0, {arm}")
        K[f"s20_o1bonly_nolink_{arm}"] = (nol, f"of these: only follow-ons without links, {arm}"); K[f"s20_o1bonly_public_only_{arm}"] = (pub, f"of these: a follow-on with only public or variant links, {arm}")
    write_csv(OUT_TABLES / "20_o1b_decomposition.csv", dec)
    K["s20_o1bonly_diff_nolink_pp"] = (round(100 * (dec[0]["follow_on_without_links"] / dec[0]["n_arm"] - dec[1]["follow_on_without_links"] / dec[1]["n_arm"]), 6), "share with only link-less follow-ons, treated minus controls, pp")
    K["s20_o1bonly_diff_public_only_pp"] = (round(100 * (dec[0]["follow_on_public_only"] / dec[0]["n_arm"] - dec[1]["follow_on_public_only"] / dec[1]["n_arm"]), 6), "share with a follow-on with only public or variant links, treated minus controls, pp")

    # ---- 6) Follow-on round types, timing and amount fill
    allfo = read_csv(OUT_DATA / "10_followons.csv"); prim_ids = {r["firm_id"] for r in prim}; tr_ids = {r["firm_id"] for r in tr}
    K["s20_followon_types"] = (dict(Counter(r["round_type"] for r in allfo if r["firm_id"] in prim_ids).most_common(8)), "eight most frequent follow-on round types, primary sample")
    first_eq = {}
    for r in sorted(allfo, key=lambda r: int(r["days_since_first"])):
        if r["firm_id"] in prim_ids and int(r["equity"]) and int(r["took_place"]): first_eq.setdefault(r["firm_id"], int(r["days_since_first"]))
    K["s20_followon_median_days"] = (statistics.median(first_eq.values()), "median days from the first round to the first equity follow-on")
    fo36 = [r for r in allfo if r["firm_id"] in prim_ids and int(r["equity"]) and int(r["took_place"]) and int(r["days_since_first"]) <= config.WINDOW_DAYS["36"]]
    K["s20_amount_filled_share_followon_36"] = (round(1 - sum(int(r["amount_missing"]) for r in fo36) / len(fo36), 6), "share of equity follow-ons within 36 months with an amount")
    tr_with = {r["firm_id"] for r in fo36 if r["firm_id"] in tr_ids}; tr_miss = {r["firm_id"] for r in fo36 if r["firm_id"] in tr_ids and int(r["amount_missing"])}
    K["s20_amount_missing_share_followon_36_treated"] = (round(len(tr_miss) / len(tr_with), 6), "share of treated firms with an equity follow-on within 36 months whose amount is missing for at least one such round")

    # ---- 7) Provenance markers: where the export looks thin, per arm, with the outcome by marker.
    #         Links are those of the first round itself (first_round_id, without other rounds of the same day) after
    #         dropping kein_investor; a private link is one of class private.
    reg_type, reg_class = load_register_maps()
    kept = lambda ps: [p for p in ps if reg_class.get(p) not in config.CLASS_IGNORED]
    links_raw, inv_cc, inv_dates = defaultdict(list), {}, defaultdict(list)
    for i in load_rows("investments"):
        if i["investor_permalink"]:
            links_raw[i["funding_round_permalink"]].append(i["investor_permalink"]); inv_cc[(i["funding_round_permalink"], i["investor_permalink"])] = i["investor_country_code"]
            dd = parse_date(i["funded_at"])
            if dd: inv_dates[i["investor_permalink"]].append(dd)
    for v in inv_dates.values(): v.sort()
    for r in prim:
        rid = r["first_round_id"]; ps = kept(links_raw.get(rid, [])); d0 = parse_date(r["first_round_date"])
        r["_day1"] = int(r["first_round_date"][8:] == "01")             # dated on the first of a month
        r["_inv_cc_all_missing"] = int(bool(ps) and all(not inv_cc.get((rid, p)) for p in ps))
        r["_only_persons"] = int(bool(ps) and all(reg_type.get(p) == "person" for p in ps))
        r["_venture_nocode"] = int(r["first_round_type"] == "venture" and r["first_round_code"] == "")
        priv = [p for p in ps if reg_class.get(p) in config.CLASS_PRIVATE]   # experience = the investor's CB investments dated before the first round
        r["_coinv_prior_ge50"] = int(any(bisect.bisect_left(inv_dates.get(p, []), d0) >= 50 for p in priv))
    for m, lab in (("_day1", "first round dated on the first of a month"), ("_inv_cc_all_missing", "all investor country codes missing"), ("_only_persons", "only persons as investors"), ("age_missing", "age missing")):
        for arm, rows in (("treated", tr), ("control", ct)):
            K[f"s20_mark_{m.lstrip('_')}_{arm}"] = (rate(rows, m), f"share with marker: {lab}, {arm}")
    for arm, rows in (("treated", tr), ("control", ct)):
        K[f"s20_first_round_venture_nocode_share_{arm}"] = (rate(rows, "_venture_nocode"), f"share whose first round is a venture round without series code, {arm}")
        K[f"s20_coinv_prior_ge50_share_{arm}"] = (rate(rows, "_coinv_prior_ge50"), f"share with a private first-round investor that had 50 or more earlier investments in the export, {arm}")
    for m, name, lab in (("age_missing", "age_missing", "age missing"), ("_day1", "day1", "first round dated on the first of a month")):
        for v in (0, 1):
            K[f"s20_o1_36_rate_treated_{name}_{v}"] = (rate([r for r in tr if r[m] == v], "o1_36"), f"o1_36 rate of treated with {lab} = {v}")
    rew = sum(rate([r for r in tr if r["_day1"] == v], "o1_36") * sum(1 for r in ct if r["_day1"] == v) / len(ct) for v in (0, 1)) - rate(ct, "o1_36")
    K["s20_o1_36_diff_naive_reweighted_day1"] = (round(rew, 6), "naive difference after reweighting treated to the controls' share of first-of-month dates")
    # country-code loss (does the loss of the country code differ by treatment?): firms worldwide whose earliest early-stage round in the cohorts has a European investor
    #         (investor country code in the export) but whose own country code is empty; they never reach stage 10
    first_any = {}
    for r in load_rows("rounds"):
        dd = parse_date(r["funded_at"])
        if not (dd and config.COHORT_FIRST <= dd.year <= config.COHORT_LAST and (r["funding_round_type"] in config.FIRST_ROUND_TYPES or (r["funding_round_type"] == "venture" and r["funding_round_code"] in config.FIRST_ROUND_VENTURE_CODES))): continue
        key = (dd, r["funding_round_permalink"])
        if r["company_permalink"] not in first_any or key < first_any[r["company_permalink"]][0]: first_any[r["company_permalink"]] = (key, r)
    cc_firms, cc_treated, usa_pub = set(), set(), set()
    for fid, (_, r) in first_any.items():
        rid = r["funding_round_permalink"]; ps = kept(links_raw.get(rid, []))
        n_pub = sum(1 for p in ps if reg_class.get(p) in config.CLASS_PUBLIC_HEADLINE); n_eu = sum(1 for p in ps if inv_cc.get((rid, p)) in config.COUNTRIES_EUROPE)
        if r["company_country_code"] == "" and n_eu > 0:
            cc_firms.add(fid); cc_treated.add(fid) if n_pub > 0 else None
        if r["company_country_code"] == "USA" and n_pub > 0: usa_pub.add(fid)
    K["s20_ccmiss_firms_eu_inv"] = (len(cc_firms), "firms without country code whose earliest early-stage round in the cohorts has a European investor")
    K["s20_ccmiss_treated"] = (len(cc_treated), "of these: with a public investor in that earliest round (investors outside the register count as not public: a lower bound)")
    K["s20_usa_hq_with_eu_public_first_round"] = (len(usa_pub), "US-headquartered firms with a public investor in their earliest early-stage round in the cohorts (lower bound: investors outside the register count as not public)")
    _values.update(K, "20_descriptives.py")
    print(f"primary {len(prim)} · treated {len(tr)} · o1_36 treated {K['s20_o1_36_rate_treated'][0]} control {K['s20_o1_36_rate_control'][0]} · no-link {len(nolink)}")

if __name__ == "__main__":
    main()
