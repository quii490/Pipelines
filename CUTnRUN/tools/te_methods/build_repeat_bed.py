#!/usr/bin/env python3
"""Build the BED4 RepeatMasker view required by the T3E parsers."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--te-saf", required=True)
    parser.add_argument("--te-anno", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    saf_path = Path(args.te_saf)
    anno_path = Path(args.te_anno)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with anno_path.open() as anno, saf_path.open() as saf, output_path.open("w") as output:
        header = anno.readline().rstrip("\n").split("\t")
        if header[:5] != ["GeneID", "repName", "Class", "Family", "milliDiv"]:
            raise ValueError(f"unexpected TE annotation header in {anno_path}")
        for line_number, saf_line in enumerate(saf, start=1):
            anno_line = anno.readline()
            if not anno_line:
                raise ValueError(f"annotation ended before SAF at line {line_number}")
            saf_fields = saf_line.rstrip("\n").split("\t")
            anno_fields = anno_line.rstrip("\n").split("\t")
            if len(saf_fields) != 5 or len(anno_fields) < 5:
                raise ValueError(f"malformed SAF/annotation row at line {line_number}")
            if saf_fields[0] != anno_fields[0]:
                raise ValueError(
                    f"SAF and TE annotation order differs at line {line_number}: "
                    f"{saf_fields[0]} != {anno_fields[0]}"
                )
            chrom, start, end, strand = saf_fields[1], saf_fields[2], saf_fields[3], saf_fields[4]
            rep_name = anno_fields[1] or saf_fields[0]
            # T3E's upstream parsers split every row into exactly four fields.
            # The strand remains available in the source SAF/annotation files.
            output.write(f"{chrom}\t{start}\t{end}\t{rep_name}\n")
            written += 1
        if anno.readline():
            raise ValueError("TE annotation has more rows than SAF")
    if written == 0:
        raise ValueError("no BED rows were written")
    print(f"rows={written}")


if __name__ == "__main__":
    main()
