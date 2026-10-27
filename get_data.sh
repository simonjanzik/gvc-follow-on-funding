#!/bin/sh
# get_data.sh: downloads the raw data the pipeline reads (raw/cb-2015/) and verifies every file.
# Source: Crunchbase export of December 2015, published at github.com/notpeter/crunchbase-data (last commit January 2016),
#         pinned to commit be8711fd4ee6b81110d243eb887ccc6bc8513601. License CC BY-NC 4.0 (raw/cb-2015/LICENSE.txt):
#         non-commercial use with attribution to Crunchbase. The files are not part of this repository.
# Usage:  sh get_data.sh      (from any directory; needs curl and shasum; about 72 MB)
set -eu
cd "$(dirname "$0")"
BASE=https://raw.githubusercontent.com/notpeter/crunchbase-data/be8711fd4ee6b81110d243eb887ccc6bc8513601
DIR=raw/cb-2015
mkdir -p "$DIR"
while read -r sum name src; do
  [ -f "$DIR/$name" ] || { echo "downloading $name"; curl -fsSL "$BASE/$src" -o "$DIR/$name.part" || { echo "download failed: $name"; exit 1; }; mv "$DIR/$name.part" "$DIR/$name"; }
  echo "$sum  $DIR/$name" | shasum -a 256 -c - >/dev/null || { echo "checksum mismatch: $name (delete $DIR/$name and run again)"; exit 1; }
  echo "ok  $name"
done <<EOF
5bef083fa9b47089ee084a5a1a24f7cda0d7bf6b873bfa7cd31db5d8f4f50330 acquisitions.csv acquisitions.csv
3d3a8f386f1cd64e65b6f7a4cb57b464d5ed448202c0df45d9f80645d70e2354 companies.csv companies.csv
9980cc2c62e22021b2c32b04ee71e5b76db8f27f47f0f7bd7f519f07789b3ca7 investments.csv investments.csv
89f575ed0d850e2890e903cde0bca8bdf66b028cac6fa298d534a252d6517315 rounds.csv rounds.csv
6be9af2cb7837b87f2947669a32b8ca2aa580b677d4ccaf68699944526d01d7d LICENSE.txt LICENSE
EOF
echo "raw data complete in $DIR"
