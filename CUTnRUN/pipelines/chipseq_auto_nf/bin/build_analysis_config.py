#!/usr/bin/env python3
import argparse
import hashlib
import json
import logging
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger('build_analysis_config')


def r_string(value):
    return '"' + str(value).replace('\\', '\\\\').replace('"', '\\"') + '"'


def r_vec(items):
    items = [str(x) for x in items if str(x).strip()]
    if not items:
        return 'character(0)'
    return 'c(' + ', '.join(r_string(x) for x in items) + ')'


def r_named_vec(mapping):
    if not mapping:
        return 'character(0)'
    return 'c(' + ', '.join(
        f'{r_string(key)} = {r_string(value)}' for key, value in mapping.items()
    ) + ')'


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).resolve().open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def norm_bool(x):
    return str(x).strip().lower() in {'1', 'true', 'yes', 'y'}


def main() -> None:
    parser = argparse.ArgumentParser(description='Build auto config.R for ChIP/CUT&RUN downstream plotting')
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--gene-counts', required=True)
    parser.add_argument('--te-counts', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--results-dir', required=True)
    parser.add_argument('--repeat-anno', required=True)
    parser.add_argument('--gene-anno', required=True)
    parser.add_argument('--deny-list', required=True)
    parser.add_argument('--gene-saf', required=True)
    parser.add_argument('--te-saf', required=True)
    parser.add_argument('--te-classes', default='LTR,LINE,SINE')
    parser.add_argument('--min-total-normalized-counts', type=int, default=10)
    parser.add_argument('--n-labels-genes', type=int, default=20)
    parser.add_argument('--n-labels-te-family', type=int, default=20)
    parser.add_argument('--n-labels-te-repname', type=int, default=40)
    parser.add_argument('--max-te-plot-rows', type=int, default=200000)
    parser.add_argument('--max-te-scatter-points', type=int, default=100000)
    parser.add_argument('--max-te-facet-levels', type=int, default=100)
    parser.add_argument('--skip-te-deseq2', action='store_true', help='skip optional TE family DESeq2')
    parser.add_argument('--analysis-code', action='append', default=[])
    parser.add_argument('--fingerprint-param', action='append', default=[])
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, dtype=str).fillna('')
    if manifest.empty:
        raise ValueError('manifest is empty')

    species_set = set(manifest['species'].astype(str).str.strip())
    if len(species_set) != 1:
        raise ValueError(f'analysis currently expects one species per run, got: {sorted(species_set)}')
    species = next(iter(species_set))

    controls = manifest.loc[manifest['is_igg'].map(norm_bool), 'sample'].astype(str).str.strip().tolist()
    targets = manifest.loc[~manifest['is_igg'].map(norm_bool), 'sample'].astype(str).str.strip().tolist()
    known_samples = set(manifest['sample'].astype(str).str.strip())
    target_control_map = {}
    for _, row in manifest.loc[~manifest['is_igg'].map(norm_bool)].iterrows():
        target = str(row['sample']).strip()
        control = str(row.get('igg', '')).strip()
        if not control and len(controls) == 1:
            control = controls[0]
        if control and control not in known_samples:
            raise ValueError(f'target {target!r} refers to unknown control {control!r}')
        if not control and len(controls) > 1:
            raise ValueError(f'target {target!r} has no igg mapping but multiple controls exist')
        if control:
            target_control_map[target] = control
    te_classes = [x.strip() for x in args.te_classes.split(',') if x.strip()]

    fingerprint_files = {
        'manifest': args.manifest,
        'gene_counts': args.gene_counts,
        'te_counts': args.te_counts,
        'repeat_anno': args.repeat_anno,
        'gene_anno': args.gene_anno,
        'deny_list': args.deny_list,
        'gene_saf': args.gene_saf,
        'te_saf': args.te_saf,
    }
    for index, path in enumerate(args.analysis_code):
        fingerprint_files[f'analysis_code_{index}'] = path
    fingerprint_payload = {
        'files': {name: sha256(path) for name, path in fingerprint_files.items()},
        'params': sorted(args.fingerprint_param),
        'target_control_map': target_control_map,
        'te_classes': te_classes,
        'min_total_normalized_counts': args.min_total_normalized_counts,
        'labels': [args.n_labels_genes, args.n_labels_te_family, args.n_labels_te_repname],
        'max_te_plot_rows': args.max_te_plot_rows,
        'max_te_scatter_points': args.max_te_scatter_points,
        'max_te_facet_levels': args.max_te_facet_levels,
        'run_te_deseq2': not args.skip_te_deseq2,
    }
    analysis_fingerprint = hashlib.sha256(
        json.dumps(fingerprint_payload, sort_keys=True).encode()
    ).hexdigest()

    out = Path(args.output).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    results_dir = Path(args.results_dir).resolve()

    logger.info('writing auto config: %s', out)
    with open(out, 'w') as fh:
        fh.write('# Auto-generated by build_analysis_config.py\n')
        fh.write(f'repeat_annotations_file <- "{Path(args.repeat_anno).resolve()}"\n')
        fh.write(f'gene_annotations_file <- "{Path(args.gene_anno).resolve()}"\n')
        fh.write(f'deny_list_file <- "{Path(args.deny_list).resolve()}"\n')
        fh.write(f'gene_saf_file <- "{Path(args.gene_saf).resolve()}"\n')
        fh.write(f'te_saf_file <- "{Path(args.te_saf).resolve()}"\n')
        fh.write(f'gene_counts_file <- "{Path(args.gene_counts).resolve()}"\n')
        fh.write(f'te_counts_file <- "{Path(args.te_counts).resolve()}"\n')
        fh.write(f'output_dir <- "{results_dir}"\n')
        fh.write('intermediate_dir <- file.path(output_dir, "intermediate_data")\n')
        fh.write('figures_dir <- file.path(output_dir, "figures")\n')
        fh.write(f'MIN_TOTAL_NORMALIZED_COUNTS <- {args.min_total_normalized_counts}\n')
        fh.write(f'N_LABELS_GENES <- {args.n_labels_genes}\n')
        fh.write(f'N_LABELS_TE_FAMILY <- {args.n_labels_te_family}\n')
        fh.write(f'N_LABELS_TE_REPNAME <- {args.n_labels_te_repname}\n')
        fh.write(f'MAX_TE_PLOT_ROWS <- {args.max_te_plot_rows}\n')
        fh.write(f'MAX_TE_SCATTER_POINTS <- {args.max_te_scatter_points}\n')
        fh.write(f'MAX_TE_FACET_LEVELS <- {args.max_te_facet_levels}\n')
        fh.write(f'RUN_TE_DESEQ2 <- {str(not args.skip_te_deseq2).upper()}\n')
        fh.write(f'TE_CLASSES_OF_INTEREST <- {r_vec(te_classes)}\n')
        fh.write(f'CONTROL_SAMPLES <- {r_vec(controls)}\n')
        fh.write(f'TARGET_SAMPLES <- {r_vec(targets)}\n')
        fh.write(f'TARGET_CONTROL_MAP <- {r_named_vec(target_control_map)}\n')
        fh.write(f'ANALYSIS_FINGERPRINT <- "{analysis_fingerprint}"\n')
        fh.write('cache_files <- list(\n')
        fh.write('  fingerprint = file.path(intermediate_dir, "analysis_fingerprint.txt"),\n')
        fh.write('  gene_counts = file.path(intermediate_dir, "01_gene_counts.rds"),\n')
        fh.write('  te_counts   = file.path(intermediate_dir, "01_te_counts.rds"),\n')
        fh.write('  gene_plot_data = file.path(intermediate_dir, "02_gene_plot_data.rds"),\n')
        fh.write('  te_plot_data   = file.path(intermediate_dir, "02_te_plot_data.rds"),\n')
        fh.write('  plot_data    = file.path(intermediate_dir, "02_plot_data.rds")\n')
        fh.write(')\n')

    logger.info('controls: %s', ', '.join(controls) if controls else 'none')
    logger.info('targets : %s', ', '.join(targets) if targets else 'none')
    logger.info('target-control map: %s', target_control_map if target_control_map else 'none')
    logger.info('analysis fingerprint: %s', analysis_fingerprint)


if __name__ == '__main__':
    main()
