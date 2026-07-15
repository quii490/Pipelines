#!/usr/bin/env python3
"""Convert transcript/exon records in a GTF file to Picard refFlat format."""

import argparse
import gzip
import os
import re


ATTR_RE = re.compile(r'([A-Za-z0-9_.:-]+)\s+"([^"]*)"')


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") else open(path, encoding="utf-8")


def parse_attrs(text):
    return dict(ATTR_RE.findall(text))


def format_record(transcript_id, record):
    exons = sorted(set(record["exons"]))
    if not exons or record.get("strand") not in {"+", "-"}:
        return None
    tx_start = min(start for start, _ in exons)
    tx_end = max(end for _, end in exons)
    if record["cds"]:
        cds_start = min(start for start, _ in record["cds"])
        cds_end = max(end for _, end in record["cds"])
    else:
        cds_start = tx_end
        cds_end = tx_end
    exon_starts = ",".join(str(start) for start, _ in exons) + ","
    exon_ends = ",".join(str(end) for _, end in exons) + ","
    return (
        record["gene_name"], transcript_id, record["chrom"], record["strand"],
        tx_start, tx_end, cds_start, cds_end, len(exons), exon_starts, exon_ends,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Gene annotation in GTF or GTF.GZ format")
    parser.add_argument("--output", required=True, help="Output Picard refFlat text file")
    args = parser.parse_args()

    output_tmp = f"{args.output}.tmp.{os.getpid()}"
    current_id = None
    current = None
    seen = set()
    row_count = 0

    try:
        with open_text(args.input) as handle, open(output_tmp, "w", encoding="utf-8") as out:
            for line_number, line in enumerate(handle, 1):
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9 or fields[2] not in {"transcript", "exon", "CDS"}:
                    continue
                attrs = parse_attrs(fields[8])
                transcript_id = attrs.get("transcript_id")
                if not transcript_id:
                    continue

                if transcript_id != current_id:
                    if current_id is not None:
                        row = format_record(current_id, current)
                        if row is not None:
                            out.write("\t".join(map(str, row)) + "\n")
                            row_count += 1
                        seen.add(current_id)
                    if transcript_id in seen:
                        raise ValueError(
                            f"Transcript {transcript_id} is not contiguous at {args.input}:{line_number}; "
                            "sort/group the GTF by transcript before conversion"
                        )
                    current_id = transcript_id
                    current = {"exons": [], "cds": []}

                try:
                    start0 = int(fields[3]) - 1
                    end = int(fields[4])
                except ValueError as exc:
                    raise ValueError(f"Invalid coordinates at {args.input}:{line_number}") from exc

                current["gene_name"] = attrs.get("gene_name") or attrs.get("gene_id") or transcript_id
                current["chrom"] = fields[0]
                current["strand"] = fields[6]
                if fields[2] == "exon":
                    current["exons"].append((start0, end))
                elif fields[2] == "CDS":
                    current["cds"].append((start0, end))

            if current_id is not None:
                row = format_record(current_id, current)
                if row is not None:
                    out.write("\t".join(map(str, row)) + "\n")
                    row_count += 1

        if row_count == 0:
            raise ValueError(f"No transcript exon records found in {args.input}")
        os.replace(output_tmp, args.output)
    except Exception:
        if os.path.exists(output_tmp):
            os.unlink(output_tmp)
        raise

    print(f"[gtf_to_refflat] transcripts={row_count} output={args.output}")


if __name__ == "__main__":
    main()
