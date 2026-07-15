#!/usr/bin/env python3
"""Validate relative Markdown links in the documentation source."""

from __future__ import annotations

import re
import sys
from pathlib import Path
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs" / "content"
LINK = re.compile(r"(?<!!)\[[^\]]+\]\(([^)]+)\)")


def main() -> int:
    errors: list[str] = []
    for source in sorted(DOCS.rglob("*.md")):
        text = source.read_text(encoding="utf-8")
        for raw in LINK.findall(text):
            target = raw.strip().split(maxsplit=1)[0].strip("<>")
            if not target or target.startswith(("#", "http://", "https://", "mailto:")):
                continue
            path_text = unquote(target.split("#", 1)[0])
            if not path_text:
                continue
            candidate = (source.parent / path_text).resolve()
            if candidate.is_dir():
                candidate = candidate / "index.md"
            if not candidate.exists():
                errors.append(f"{source.relative_to(ROOT)} -> {target}")

    if errors:
        print("Broken internal Markdown links:", file=sys.stderr)
        print("\n".join(f"  {item}" for item in errors), file=sys.stderr)
        return 1
    print("Internal Markdown link checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
