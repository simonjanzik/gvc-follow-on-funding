#!/usr/bin/env python3
"""values.json: the single writer and the merge for R stages.

Purpose : every reported number is a key with value, label, script and compute date. Python stages call
          update(); R stages write <NN>_values.csv (key,value,label) and
          run_all.sh merges them here.
Input   : outputs/values.json (created if missing), <NN>_values.csv for merge
Output  : outputs/values.json (sorted keys, atomic write)
Thesis  : Section 3.10, Appendix A.3
Usage   : python3 _values.py merge ../outputs/data/30_values.csv 30_estimators.R
          python3 _values.py show s30_att_main_total
"""
import csv, datetime, json, os, re, sys
import config
from _lib import OUT, DATA

VALUES = OUT / "values.json"

def load():
    return json.loads(VALUES.read_text(encoding="utf-8")) if VALUES.exists() else {}

def _round(v):
    if isinstance(v, bool) or not isinstance(v, float):
        return v
    return round(v, config.ROUND_DECIMALS)

def update(new, script):
    """new: {key: (value, label)}"""
    vals = load(); today = datetime.date.today().isoformat()
    for k, v in new.items():
        value, label = v
        if isinstance(value, dict):
            value = {kk: _round(vv) for kk, vv in value.items()}
        vals[k] = {"value": _round(value), "label": label, "script": script, "date": today}
    tmp = VALUES.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(vals, indent=1, ensure_ascii=False, sort_keys=True), encoding="utf-8")
    os.replace(tmp, VALUES)
    print(f"[values] {len(new)} keys -> {VALUES.relative_to(DATA)} ({len(vals)} total)")

def merge_csv(path, script):
    rows = list(csv.DictReader(open(path, newline="", encoding="utf-8")))
    new = {}
    for r in rows:
        raw = r["value"]
        try:
            value = int(raw) if re.fullmatch(r"-?\d+", raw) else float(raw)
        except ValueError:
            value = raw
        new[r["key"]] = (value, r.get("label", ""))
    update(new, script)

if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "merge":
        merge_csv(sys.argv[2], sys.argv[3])
    elif cmd == "show":
        print(json.dumps({k: load().get(k) for k in sys.argv[2:]}, indent=1, ensure_ascii=False))
    else:
        sys.exit("usage: _values.py merge <NN_values.csv> <script> | show <key> ...")
