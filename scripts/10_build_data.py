#!/usr/bin/env python3
"""Stage 10: sample, register lookup, treatment, outcomes and covariates -> analysis.csv.

Purpose : one row per firm with its first European early-stage round, the treatment from the investor register,
          the follow-on outcomes and the covariates. Stages 30 to 60 read only analysis.csv; stage 20 also rereads
          the raw rounds and investments for the provenance markers and the country-code loss.
Thesis  : Sections 3.1 (sample, first round, links without an investor), 3.2 (treatment from the register, three
          definitions), 3.3 (outcomes, covariates); Appendix A.1 (sample rules in numbers), A.2 (public investors)
Input   : raw/cb-2015/{companies,rounds,investments,acquisitions}.csv, investors.csv, config.py
Output  : outputs/data/analysis.csv, outputs/data/analysis_dictionary.csv, outputs/data/10_followons.csv (one row per
          follow-on round of a sample firm, for stage 20), outputs/tables/10_funnel.csv,
          10_dropped_firms.csv, 10_public_investors.csv, 10_register_by_class.csv; keys s10_* in values.json
Usage   : python3 10_build_data.py   (no arguments; run_all.sh runs it as stage 10)
Seeds   : none (deterministic)
Checks  : register guard (every investor permalink in a relevant round is in the register, else stop with the list
          in outputs/tables/10_missing_permalinks.csv), dictionary = column registry.
"""
import datetime, hashlib
from collections import Counter, defaultdict
import config
from _lib import REGISTER, OUT_DATA, OUT_TABLES, parse_date, load_rows, read_csv, write_csv, require
import _values

W36, W24, W48 = config.WINDOW_DAYS["36"], config.WINDOW_DAYS["24"], config.WINDOW_DAYS["48"]
PUBLIC, BROAD = set(config.CLASS_PUBLIC_HEADLINE), set(config.CLASS_PUBLIC_BROAD)
PRIVATE, IGNORED = set(config.CLASS_PRIVATE), set(config.CLASS_IGNORED)

# ---- 0) Column registry: every column of analysis.csv is declared where it is built; the dictionary is written from it
DICT = []
def col(name, typ, definition, source):
    DICT.append({"column": name, "type": typ, "definition": definition, "source": source}); return name

# ---- 1) Register interface: what stage 10 knows about an investor permalink
def load_register(path=REGISTER):
    """permalink -> row (pipeline columns and fact columns). Rows without permalink are data artefacts and are skipped."""
    rows = [r for r in read_csv(path) if r["permalink"]]
    missing_cols = [c for c in config.REGISTER_COLUMNS + config.REGISTER_FACT_COLUMNS if c not in rows[0]]
    require(not missing_cols, f"register lacks columns {missing_cols}")
    for r in rows:
        require(r["class"] in config.CLASSES_ALL, f"register: unknown class {r['class']!r} for {r['permalink']}")
    reg = {r["permalink"]: r for r in rows}
    digest = hashlib.sha256(open(path, "rb").read()).hexdigest()[:12]
    return reg, len(rows), digest

def is_early(rtype, rcode):
    return rtype in config.FIRST_ROUND_TYPES or (rtype == "venture" and rcode in config.FIRST_ROUND_VENTURE_CODES)

def load_raw():
    companies = {c["permalink"]: c for c in load_rows("companies")}
    rounds_by_firm = defaultdict(list)
    for r in load_rows("rounds"):
        d = parse_date(r["funded_at"])
        if d:
            r["_date"] = d; rounds_by_firm[r["company_permalink"]].append(r)
    links_by_round = defaultdict(list)                      # links without permalink (no name either) do not count
    for i in load_rows("investments"):
        if i["investor_permalink"]:
            links_by_round[i["funding_round_permalink"]].append(i["investor_permalink"])
    acq_by_firm = defaultdict(list)
    for a in load_rows("acquisitions"):
        d = parse_date(a["acquired_at"])
        if d: acq_by_firm[a["company_permalink"]].append(d)
    return companies, rounds_by_firm, links_by_round, acq_by_firm

# ---- 2) Sample: scope rounds, the first round that took place, same-day rule
def build(reg, rounds_by_firm, links_by_round):
    class_of = lambda p: reg[p]["class"] if p in reg else None   # None: permalink not in the register
    missing = Counter()

    def kept_links(rp, dropped=None):
        """links of a round after dropping the ignored class; None if an unknown permalink leaves the round unclassified;
        the ignored links are appended to `dropped` when a list is given"""
        out = []
        for p in dict.fromkeys(links_by_round.get(rp, [])):      # unique permalinks, order kept
            c = class_of(p)
            if c is None: missing[p] += 1; return None
            if c not in IGNORED: out.append((p, c))
            elif dropped is not None: dropped.append(p)
        return out

    scope = defaultdict(list)
    for f, rs in rounds_by_firm.items():
        for r in rs:
            if (r["company_country_code"] in config.COUNTRIES_EUROPE and is_early(r["funding_round_type"], r["funding_round_code"])
                    and config.COHORT_FIRST <= r["_date"].year <= config.COHORT_LAST):
                scope[f].append(r)
    n_scope_rounds = sum(len(v) for v in scope.values())
    n_only_ignored = 0
    firms, dropped = [], []
    for f in sorted(scope):
        rs = sorted(scope[f], key=lambda r: (r["_date"], r["funding_round_permalink"]))
        happened = []                                          # (round, kept links) for rounds that took place
        for r in rs:
            kl = kept_links(r["funding_round_permalink"])
            if kl is None: continue                            # unknown permalink: the guard reports it below
            if links_by_round.get(r["funding_round_permalink"]) and not kl:
                n_only_ignored += 1; continue                  # only ignored links: this round did not take place
            happened.append((r, kl))
        if not happened:
            dropped.append({"firm_id": f, "name": rs[0]["company_name"], "country": rs[0]["company_country_code"],
                            "n_scope_rounds": len(rs), "first_scope_round_date": rs[0]["_date"].isoformat()})
            continue
        first, kl = happened[0]
        d0 = first["_date"]
        skipped = sum(1 for r in rs if r["_date"] < d0)     # scope rounds strictly before that did not take place; same-day ones are neither skipped nor merged
        shifted = int(first["funding_round_permalink"] != rs[0]["funding_round_permalink"])
        siblings = [(r, k) for r, k in happened[1:] if r["_date"] == d0]
        event = [(first, kl)] + siblings
        links = list(dict.fromkeys(p_c for _, k in event for p_c in k))
        dropped_links = []
        for r, _ in event: kept_links(r["funding_round_permalink"], dropped_links)
        dropped_links = list(dict.fromkeys(dropped_links))       # ignored links of the whole event, for the register keys
        amounts = [r["raised_amount_usd"] for r, _ in event]
        filled = [float(a) for a in amounts if a]
        firms.append({"firm": f, "first": first, "d0": d0, "links": links, "dropped": dropped_links, "same_day_merged": len(event),
                      "rounds_skipped": skipped, "shifted": shifted,
                      "amount": sum(filled) if filled else "", "amount_missing": int(len(filled) < len(amounts)),
                      "treated_single": int(any(c in PUBLIC for _, c in kl))})
    # guard: part 1 is counted in kept_links above; part 2 checks every round of a sample firm
    for f in firms:
        for r in rounds_by_firm[f["firm"]]:
            for p in links_by_round.get(r["funding_round_permalink"], []):
                if class_of(p) is None: missing[p] += 1
    return dict(firms=firms, dropped=dropped, scope=scope, n_scope_rounds=n_scope_rounds, n_only_ignored=n_only_ignored,
                missing=missing, kept_links=kept_links)

# ---- 3) Treatment, outcomes and covariates per firm
def assemble(b, reg, companies, rounds_by_firm, links_by_round, acq_by_firm):
    kept_links = b["kept_links"]
    unresolved = lambda p: p in reg and reg[p]["investor_type"] == config.REGISTER_UNRESOLVED_VALUE
    first_cat = {}
    for f in b["firms"]:
        c = companies.get(f["firm"])                        # companies.csv owns the categories; the round's copy only if the firm has no row
        cats = c["category_list"] if c is not None else f["first"].get("company_category_list", "")
        first_cat[f["firm"]] = cats.split("|")[0].strip() or config.SECTOR_UNKNOWN
    top = {s for s, _ in Counter(first_cat.values()).most_common(config.SECTOR_TOP_K)}

    rows, followon_rows = [], []
    for f in b["firms"]:
        fid, d0, links = f["firm"], f["d0"], f["links"]
        n_pub = sum(1 for _, c in links if c in PUBLIC); n_var = sum(1 for _, c in links if c == "public_variant")
        n_priv = sum(1 for _, c in links if c in PRIVATE)
        n_kein = len(f["dropped"])
        treated = int(n_pub > 0)
        strict = int(any(c in PUBLIC and reg.get(p, {}).get("manager") == config.STRICT_MANAGER for p, c in links))
        only_unres = int(bool(links) and all(unresolved(p) for p, _ in links))

        followons = [r for r in rounds_by_firm[fid] if r["_date"] > d0]
        for r in followons:                                    # table for stage 20: every later round with its link classes
            kl = kept_links(r["funding_round_permalink"]); raw = links_by_round.get(r["funding_round_permalink"], [])
            followon_rows.append({"firm_id": fid, "round_id": r["funding_round_permalink"], "days_since_first": (r["_date"] - d0).days,
                "round_type": r["funding_round_type"], "equity": int(r["funding_round_type"] in config.EQUITY_TYPES),
                "took_place": int(not (raw and not kl)), "n_links": len(kl or []), "n_private": sum(1 for _, c in (kl or []) if c in PRIVATE),
                "n_public": sum(1 for _, c in (kl or []) if c in PUBLIC), "n_variant": sum(1 for _, c in (kl or []) if c == "public_variant"),
                "amount_usd": r["raised_amount_usd"], "amount_missing": int(not r["raised_amount_usd"])})
        def outcome(window, kind, equity=config.EQUITY_TYPES, resolved_only=False):
            """kind o1: a follow-on with a private link · o1b: any equity follow-on, also without links ·
               o2: a follow-on with links and without public or variant link. A follow-on whose links are all
               ignored did not take place as an equity round."""
            for r in followons:
                if (r["_date"] - d0).days > window or r["funding_round_type"] not in equity: continue
                kl = kept_links(r["funding_round_permalink"])
                if links_by_round.get(r["funding_round_permalink"]) and not kl: continue
                if kind == "o1b": return 1
                if not kl: continue
                priv = [p for p, c in kl if c in PRIVATE and not (resolved_only and unresolved(p))]
                if kind == "o1" and priv: return 1
                if kind == "o2" and not any(c in BROAD for _, c in kl): return 1
            return 0
        o1_36 = outcome(W36, "o1")
        founded = parse_date(companies.get(fid, {}).get("founded_at", ""))
        age = (d0 - founded).days / 365.25 if founded else None
        if age is not None and not (config.AGE_MIN_YEARS <= age <= config.AGE_MAX_YEARS): age = None
        row = {
            col("firm_id", "text", "firm permalink", "companies"): fid,
            col("country", "text", "ISO3 of the headquarters", "rounds"): f["first"]["company_country_code"],
            col("cohort", "int", "year of the first round", "rounds"): d0.year,
            col("first_round_id", "text", "permalink of the first round (same-day event: the first in sort order)", "rounds"): f["first"]["funding_round_permalink"],
            col("first_round_date", "date", "date of the first round", "rounds"): d0.isoformat(),
            col("first_round_type", "text", "seed, angel or venture", "rounds"): f["first"]["funding_round_type"],
            col("first_round_code", "text", "series code as in the export (the A-or-empty restriction applies to venture rounds only)", "rounds"): f["first"]["funding_round_code"],
            col("first_round_amount_usd", "num", "amount in USD, summed over the filled parts of a same-day event, empty if no part has one", "rounds"): f["amount"],
            col("amount_missing", "0/1", "amount missing for the round or for at least one part of a same-day event", "rounds"): f["amount_missing"],
            col("same_day_merged", "int", "scope rounds merged on the first-round day (1 = none)", "rounds"): f["same_day_merged"],
            col("rounds_skipped", "int", "earlier scope rounds that did not take place (only ignored links)", "rounds, investments, register"): f["rounds_skipped"],
            col("n_kein_investor_first", "int", "links of the first-round event dropped as class kein_investor", "investments, register"): n_kein,
            col("n_links_first", "int", "links of the first-round event after dropping kein_investor", "investments, register"): len(links),
            col("n_public_first", "int", "links of class public in the first-round event", "investments, register"): n_pub,
            col("n_variant_first", "int", "links of class public_variant in the first-round event", "investments, register"): n_var,
            col("n_private_first", "int", "links of class private in the first-round event", "investments, register"): n_priv,
            col("has_links", "0/1", "n_links_first > 0 (primary sample)", "investments, register"): int(bool(links)),
            col("treated", "0/1", "at least one public link (headline definition)", "investments, register"): treated,
            col("treated_broad", "0/1", "at least one public or public_variant link", "investments, register"): int(n_pub + n_var > 0),
            col("treated_strict", "0/1", "at least one public link with a public manager", "investments, register"): strict,
            col("treated_solo", "0/1", "treated without a private link in the first round", "investments, register"): int(treated and n_priv == 0),
            col("treated_synd", "0/1", "treated with at least one private link", "investments, register"): int(treated and n_priv > 0),
            col("first_round_only_unresolved", "0/1", "all links of the first round are unresolved organisations (investor_type unbekannt)", "investments, register"): only_unres,
            col("sector", "text", f"first category, {config.SECTOR_TOP_K} most frequent kept (ties by first occurrence, an empty category is its own level), rest other", "companies"): first_cat[fid] if first_cat[fid] in top else config.SECTOR_OTHER,
            col("age_at_round", "num", f"years from founding to the first round, empty if missing or outside {config.AGE_MIN_YEARS} to {config.AGE_MAX_YEARS}", "companies"): round(age, config.AGE_DECIMALS) if age is not None else "",
            col("age_missing", "0/1", "age missing", "companies"): int(age is None),
            col("pre_funding", "0/1", "any round of any type before the first round, including earlier scope rounds with only kein_investor links", "rounds"): int(any(r["_date"] < d0 for r in rounds_by_firm[fid])),
            col("o1_36", "0/1", "equity follow-on round with at least one private link within 36 months (main outcome)", "rounds, investments, register"): o1_36,
            col("o1_24", "0/1", "as o1_36, within 24 months", "rounds, investments, register"): outcome(W24, "o1"),
            col("o1_48", "0/1", f"as o1_36, within 48 months, empty for cohorts after {config.WINDOW_48_LAST_COHORT}", "rounds, investments, register"): outcome(W48, "o1") if d0.year <= config.WINDOW_48_LAST_COHORT else "",
            col("o1b_36", "0/1", "any equity follow-on round within 36 months, also without links; a round whose links are all kein_investor does not count", "rounds, investments, register"): outcome(W36, "o1b"),
            col("o1_undisc_36", "0/1", "o1_36 with undisclosed rounds counted as equity", "rounds, investments, register"): outcome(W36, "o1", config.EQUITY_TYPES_UNDISCLOSED),
            col("o2_36", "0/1", "equity follow-on round with links, none of class public or public_variant, within 36 months", "rounds, investments, register"): outcome(W36, "o2"),
            col("o1_36_only_unresolved", "0/1", "o1_36 = 1 only through unresolved private links", "investments, register"): int(o1_36 and not outcome(W36, "o1", resolved_only=True)),
            col("exit_acq_36", "0/1", "acquired within 36 months after the first round", "acquisitions"): int(any(0 < (a - d0).days <= W36 for a in acq_by_firm.get(fid, []))),
            col("status", "text", "company status in the export", "companies"): companies.get(fid, {}).get("status", ""),
        }
        rows.append(row)
    return rows, followon_rows

# ---- 4) Register keys: counts over the register and over the first rounds of the sample
def register_keys(reg, first_links, dropped_links):
    orgs = [r for r in reg.values() if r["type"] == "organization"]
    links_of = Counter()                                       # permalink -> number of first-round links in this sample (kept and dropped)
    for links in first_links:
        for p, _ in links: links_of[p] += 1
    for links in dropped_links:
        for p in links: links_of[p] += 1
    def cnt(sel):
        s = [r for r in orgs if sel(r)]; return len(s), sum(links_of[r["permalink"]] for r in s)
    K = {}
    for c in config.CLASSES_ALL:
        n, l = cnt(lambda r, c=c: r["class"] == c)
        K[f"s10_reg_orgs_{c}"] = (n, f"organisations of class {c} in the register")
        K[f"s10_reg_links_first_{c}"] = (l, f"first-round links of organisations of class {c} in the sample" + (" (dropped links)" if c in IGNORED else ""))
    persons = [r for r in reg.values() if r["type"] == "person"]
    K["s10_reg_persons"] = (len(persons), "persons in the register (private by rule)")
    K["s10_reg_links_first_person"] = (sum(links_of[r["permalink"]] for r in persons), "first-round links of persons in the sample")
    n, l = cnt(lambda r: r["class"] == "public" and r["manager"] == config.STRICT_MANAGER)
    K["s10_reg_orgs_strict"] = (n, "public organisations with a public manager (strict cut)")
    K["s10_reg_links_first_strict"] = (l, "first-round links of the strict cut")
    for c in ("public", "public_variant"):
        n_all, l_all = cnt(lambda r, c=c: r["class"] == c); n_p, l_p = cnt(lambda r, c=c: r["class"] == c and r["confidence"] == "belegt")
        K[f"s10_reg_share_orgs_primary_source_{c}"] = (round(n_p / n_all, 6) if n_all else 0, f"share of {c} organisations with a primary source")
        K[f"s10_reg_share_links_first_primary_source_{c}"] = (round(l_p / l_all, 6) if l_all else 0, f"share of {c} first-round links with a primary source")
    n, l = cnt(lambda r: r["investor_type"] == config.REGISTER_UNRESOLVED_VALUE)
    K["s10_reg_orgs_unresolved"] = (n, "organisations without a resolved owner fact (private by rule)")
    K["s10_reg_links_first_unresolved"] = (l, "first-round links of unresolved organisations")
    n, l = cnt(lambda r: r["class"] == "public" and r["instrument"] == "gemischt")
    K["s10_reg_orgs_public_mixed_instrument"] = (n, "public organisations with a mixed instrument")
    K["s10_reg_links_first_public_mixed_instrument"] = (l, "first-round links of public organisations with a mixed instrument")
    n, l = cnt(lambda r: r["investor_type"] == "vehikel_priv_oeffkap" and r["public_capital"] == "anker")
    K["s10_reg_orgs_private_fund_anchor"] = (n, "privately raised funds with a public anchor investor (private)")
    K["s10_reg_links_first_private_fund_anchor"] = (l, "first-round links of privately raised funds with a public anchor")
    n, l = cnt(lambda r: r["class"] == "public" and r["investor_type"] == "vehikel_oeff" and r["public_capital"] != "mehrheit")
    K["s10_reg_orgs_public_vehicle_no_majority"] = (n, "publicly established vehicles without a documented public majority")
    K["s10_reg_links_first_public_vehicle_no_majority"] = (l, "first-round links of publicly established vehicles without a documented public majority")
    K["s10_reg_kein_investor_by_reason"] = (dict(Counter(r["reason"] for r in orgs if r["class"] == "kein_investor").most_common()), "organisations of class kein_investor by reason")
    table = [{"class": c, "entries": K[f"s10_reg_orgs_{c}"][0], "first_round_links": K[f"s10_reg_links_first_{c}"][0]} for c in config.CLASSES_ALL]
    table.append({"class": "person (private)", "entries": len(persons), "first_round_links": K["s10_reg_links_first_person"][0]})
    public = sorted((r for r in orgs if r["class"] in BROAD and links_of[r["permalink"]] > 0), key=lambda r: (-links_of[r["permalink"]], r["name"]))
    public_table = [{"permalink": r["permalink"], "name": r["name"], "entity": r["entity"], "country": r["country"], "class": r["class"], "manager": r["manager"],
                     "instrument": r["instrument"], "first_round_links": links_of[r["permalink"]], "confidence": r["confidence"]} for r in public]
    return K, table, public_table

def main():
    reg, n_reg, digest = load_register()
    companies, rounds_by_firm, links_by_round, acq_by_firm = load_raw()
    b = build(reg, rounds_by_firm, links_by_round)
    if b["missing"]:
        write_csv(OUT_TABLES / "10_missing_permalinks.csv", [{"permalink": p, "links": n} for p, n in b["missing"].most_common()])
        require(False, f"{len(b['missing'])} investor permalinks are not in the register -> classify them in "
                       f"investors.csv (list: outputs/tables/10_missing_permalinks.csv)")
    rows, followons = assemble(b, reg, companies, rounds_by_firm, links_by_round, acq_by_firm)
    write_csv(OUT_DATA / "10_followons.csv", followons)
    names = [d["column"] for d in DICT[:len(rows[0])]]
    require(names == list(rows[0].keys()), "column registry and row layout differ")
    write_csv(OUT_DATA / "analysis.csv", rows, names)
    seen = set(); dictionary = [d for d in DICT if not (d["column"] in seen or seen.add(d["column"]))]
    write_csv(OUT_DATA / "analysis_dictionary.csv", dictionary)
    write_csv(OUT_TABLES / "10_dropped_firms.csv", b["dropped"])

    prim = [r for r in rows if r["has_links"]]; tr = [r for r in prim if r["treated"]]
    funnel = [("scope_rounds", b["n_scope_rounds"], "early-stage rounds of European firms in the cohort years"),
              ("scope_rounds_only_kein_investor", b["n_only_ignored"], "scope rounds whose links are all of class kein_investor (did not take place)"),
              ("firms_scope", len(b["scope"]), "firms with at least one scope round"),
              ("firms_dropped_no_investor_round", len(b["dropped"]), "firms without a scope round that took place"),
              ("firms", len(rows), "firms in the sample"),
              ("firms_with_links", len(prim), "firms with at least one first-round link (primary sample)"),
              ("n_treated", len(tr), "treated firms in the primary sample (headline definition)")]
    write_csv(OUT_TABLES / "10_funnel.csv", [{"step": k, "n": v} for k, v, _ in funnel])
    K = {f"s10_{k}": (v, label) for k, v, label in funnel}
    K.update({
        "s10_firms_first_round_shifted": (sum(f["shifted"] for f in b["firms"]), "firms whose first round moved because an earlier scope round had only ignored links"),
        "s10_firms_same_day_merged": (sum(1 for f in b["firms"] if f["same_day_merged"] > 1), "firms with more than one scope round on the first-round day"),
        "s10_firms_same_day_treatment_changed": (sum(1 for f in b["firms"] if f["same_day_merged"] > 1 and f["treated_single"] != int(any(c in PUBLIC for _, c in f["links"]))), "firms whose treatment differs between the first round alone and the same-day event"),
        "s10_firms_by_cohort": ({str(k): v for k, v in sorted(Counter(r["cohort"] for r in rows).items())}, "firms per cohort"),
        "s10_firms_by_country_top5": (dict(Counter(r["country"] for r in rows).most_common(5)), "firms per country, five largest"),
        "s10_share_links": (round(len(prim) / len(rows), 6), "share of firms with at least one first-round link"),
        "s10_n_treated_broad": (sum(r["treated_broad"] for r in prim), "treated under the broad definition (primary sample)"),
        "s10_n_treated_strict": (sum(r["treated_strict"] for r in prim), "treated under the strict cut (primary sample)"),
        "s10_n_treated_solo": (sum(r["treated_solo"] for r in prim), "treated without a private co-investor"),
        "s10_n_treated_synd": (sum(r["treated_synd"] for r in prim), "treated with a private co-investor"),
        "s10_treated_share_of_linked": (round(len(tr) / len(prim), 6), "treated share of the primary sample"),
        "s10_treated_by_country_top5": (dict(Counter(r["country"] for r in tr).most_common(5)), "treated firms per country, five largest"),
        "s10_distinct_public_investors_first": (len({p for f in b["firms"] for p, c in f["links"] if c in PUBLIC}), "distinct public investor identifiers in first rounds"),
        "s10_firms_kein_investor_links_first": (sum(1 for r in rows if r["n_kein_investor_first"] > 0), "firms whose first round had at least one dropped kein_investor link"),
        "s10_firms_first_round_only_unresolved_links": (sum(r["first_round_only_unresolved"] for r in prim), "primary-sample firms whose first-round links are all unresolved organisations"),
        "s10_firms_o1_36_only_unresolved_links": (sum(r["o1_36_only_unresolved"] for r in prim), "primary-sample firms whose o1_36 rests only on unresolved private links"),
        "s10_reg_rows": (n_reg, "rows of the investor register"), "s10_reg_hash": (digest, "sha256 prefix of the register file"),
    })
    RK, table, public_table = register_keys(reg, [f["links"] for f in b["firms"]], [f["dropped"] for f in b["firms"]])
    export_end = max(r["_date"] for rs in rounds_by_firm.values() for r in rs)
    K["s10_export_last_round_date"] = (export_end.isoformat(), "date of the last round in the export")
    for w in ("36", "48"):
        late = [f for f in b["firms"] if (f["d0"] + datetime.timedelta(days=config.WINDOW_DAYS[w])) > export_end and (w == "36" or f["d0"].year <= config.WINDOW_48_LAST_COHORT)]
        K[f"s10_firms_window_censored_{w}"] = (len(late), f"firms whose {w}-month window ends after the last round in the export (outcome read as 0 if nothing is recorded)" + ("" if w == "36" else f", cohorts up to {config.WINDOW_48_LAST_COHORT}"))
        K[f"s10_firms_window_censored_{w}_linked"] = (sum(1 for f in late if f["links"]), f"of these: firms of the primary sample (with a link in the first round), {w}-month window")
    K.update(RK)
    write_csv(OUT_TABLES / "10_register_by_class.csv", table)
    write_csv(OUT_TABLES / "10_public_investors.csv", public_table)
    _values.update(K, "10_build_data.py")
    print("\n".join(f"{k}: {v}" for k, v, _ in funnel), f"\nanalysis.csv: {len(rows)} rows, {len(names)} columns; dictionary {len(dictionary)} entries")

if __name__ == "__main__":
    main()
