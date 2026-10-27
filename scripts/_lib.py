#!/usr/bin/env python3
"""Shared helpers for the Python stages: paths, raw data, CSV I/O, stops. No analysis here.

Purpose : imported by 10_build_data.py, 20_descriptives.py, 50_sensitivity.py and _values.py.
Input   : none (paths are resolved from this file, so every stage runs from any working directory)
Output  : creates outputs/data and outputs/tables if missing; values.json is written only by _values.py
Usage   : imported by the stages, not run on its own
"""
import csv, datetime, os
from pathlib import Path
import config

SCRIPTS = Path(__file__).resolve().parent
DATA = SCRIPTS.parent                              # the folder above scripts/ (repository root)
RAW = DATA / "raw" / "cb-2015"
OUT = DATA / config.OUTPUT_DIR
OUT_DATA, OUT_TABLES = OUT / "data", OUT / "tables"
REGISTER = DATA / "investors.csv"
for _p in (OUT_DATA, OUT_TABLES):
    _p.mkdir(parents=True, exist_ok=True)

def parse_date(s):
    """ISO date or None; the raw export has empty and malformed dates."""
    try:
        return datetime.date.fromisoformat(s)
    except (ValueError, TypeError):
        return None

def load_rows(name):
    """All rows of a raw CB-2015 table as dicts (companies, rounds, investments, acquisitions)."""
    with open(RAW / f"{name}.csv", newline="", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f))

def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))

def write_csv(path, rows, fieldnames=None):
    """Atomic write (tmp + rename) so an interrupted run never leaves a half file."""
    path = Path(path); fieldnames = fieldnames or list(rows[0].keys())
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    os.replace(tmp, path)
    return path

def require(cond, msg):
    """Loud stop instead of a silent default."""
    if not cond:
        raise SystemExit(f"STOP: {msg}")
