#!/usr/bin/env python3
import re
import sys

if len(sys.argv) != 4:
    raise SystemExit("usage: read_module_default.py NEXTFLOW_CONFIG KEY FALLBACK")
path, key, fallback = sys.argv[1:]
try:
    with open(path, encoding="utf-8") as handle:
        text = handle.read()
    match = re.search(rf"^\s*params\.{re.escape(key)}\s*=\s*(true|false)\s*$", text, re.I | re.M)
    value = match.group(1).lower() if match else fallback
except OSError:
    value = fallback
print(value)
