#!/usr/bin/env python3
"""Add a Telescope-compatible locus attribute to a TE GTF once per run."""

import argparse
import re


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    transcript_re = re.compile(r'transcript_id "([^"]+)"')
    locus_re = re.compile(r'(?:^|;)\s*locus "')

    with open(args.input, encoding="utf-8") as src, open(args.output, "w", encoding="utf-8") as dst:
        for line in src:
            if line.startswith("#") or not line.strip() or locus_re.search(line):
                dst.write(line)
                continue
            match = transcript_re.search(line)
            if match:
                line = line.rstrip("\n") + f' locus "{match.group(1)}";\n'
            dst.write(line)


if __name__ == "__main__":
    main()
