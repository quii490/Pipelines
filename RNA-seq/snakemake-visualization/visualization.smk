import csv
import os

configfile: "config.example.yaml"

OUTDIR = config.get("outdir", "rnaseq_snakemake_downstream")
SCRIPTS = config.get(
    "scripts_dir",
    os.path.abspath(os.path.join(workflow.basedir, "..", "tools")),
)
SPECIES = config.get("species", "hg38")
SAMPLE_TABLE = config["sample_table"]
CONTRAST_FILE = config["contrast_file"]
MATRICES = config.get("matrices", {})
R_CMD = config.get("rscript", "Rscript")

if not MATRICES:
    raise ValueError("config['matrices'] is empty")


def read_contrasts(path):
    out = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if "case" in row and "control" in row:
                out.append((row["case"], row["control"]))
            elif "group_col" in row and "case" in row and "control" in row:
                out.append((row["case"], row["control"]))
            else:
                raise ValueError("contrast_file must contain case,control or group_col,case,control")
    if not out:
        raise ValueError("contrast_file is empty")
    return out


CONTRASTS = read_contrasts(CONTRAST_FILE)


def tag(matrix_name, case, control):
    return f"{matrix_name}_{case}_vs_{control}"


def matrix_type(matrix_name):
    return MATRICES[matrix_name].get("type", "generic").lower()


def matrix_path(matrix_name):
    return MATRICES[matrix_name]["path"]


def matrix_extra(matrix_name):
    cfg = MATRICES[matrix_name]
    common = []
    if "matrix_format" in cfg:
        common += ["--matrix-format", cfg["matrix_format"]]
    return " ".join(common)


def annotation_args(matrix_name):
    cfg = MATRICES[matrix_name]
    typ = matrix_type(matrix_name)
    args = []
    if typ == "gene":
        tx = cfg.get("tx2gene_path", config.get("tx2gene_path", ""))
        if not tx:
            raise ValueError(f"matrix {matrix_name} type=gene requires tx2gene_path")
        args += ["--annotation-mode", "gene", "--tx2gene-path", tx]
    elif typ == "te":
        te = cfg.get("te_annotation_tsv", config.get("te_annotation_tsv", ""))
        if not te:
            raise ValueError(f"matrix {matrix_name} type=te requires te_annotation_tsv")
        args += ["--annotation-mode", "te", "--te-annotation-tsv", te]
        args += ["--te-label-level", cfg.get("te_label_level", config.get("te_label_level", "repName"))]
        args += ["--te-color-level", cfg.get("te_color_level", config.get("te_color_level", "repFamily"))]
    else:
        args += ["--annotation-mode", "generic"]
    return " ".join(args)


def visual_args(matrix_name):
    return annotation_args(matrix_name)


def pathway_args(matrix_name):
    cfg = MATRICES[matrix_name]
    tx = cfg.get("tx2gene_path", config.get("tx2gene_path", ""))
    if not tx:
        raise ValueError(f"matrix {matrix_name} pathway requires tx2gene_path")
    return " ".join([
        "--species", SPECIES,
        "--tx2gene-path", tx,
        "--run-go", str(config.get("run_go", True)).lower(),
        "--run-gsea", str(config.get("run_gsea", True)).lower(),
        "--disable-gseaplot2", str(config.get("disable_gseaplot2", False)).lower(),
    ])


def te_args(matrix_name):
    cfg = MATRICES[matrix_name]
    te = cfg.get("te_annotation_tsv", config.get("te_annotation_tsv", ""))
    if not te:
        raise ValueError(f"matrix {matrix_name} TE analysis requires te_annotation_tsv")
    return " ".join([
        "--te-annotation-tsv", te,
        "--te-label-level", cfg.get("te_label_level", config.get("te_label_level", "repName")),
        "--te-color-level", cfg.get("te_color_level", config.get("te_color_level", "repFamily")),
    ])


all_targets = []
for matrix_name in MATRICES:
    all_targets.append(f"{OUTDIR}/diff/{matrix_name}/{matrix_name}.diff.summary.txt")
    for case, control in CONTRASTS:
        t = tag(matrix_name, case, control)
        all_targets.append(f"{OUTDIR}/annotate/{matrix_name}/{t}/{t}.annotated_DE_matrix.csv")
        all_targets.append(f"{OUTDIR}/visuals/{matrix_name}/{t}/{t}.visuals.summary.txt")
        if matrix_type(matrix_name) == "gene" and config.get("run_pathway", True):
            all_targets.append(f"{OUTDIR}/pathway/{matrix_name}/{t}/{t}.pathway.summary.txt")
        if matrix_type(matrix_name) == "te" and config.get("run_te_analysis", True):
            all_targets.append(f"{OUTDIR}/te/{matrix_name}/{t}/{t}.te_analysis.summary.txt")


rule all:
    input:
        all_targets


rule diff_from_counts:
    input:
        matrix=lambda wc: matrix_path(wc.matrix)
    output:
        summary=f"{OUTDIR}/diff/{{matrix}}/{{matrix}}.diff.summary.txt"
    params:
        extra=lambda wc: matrix_extra(wc.matrix),
        threads_cfg=lambda wc: config.get("diff_threads", 1),
        padj=lambda wc: config.get("padj_cutoff", 0.05),
        lfc=lambda wc: config.get("lfc_cutoff", 0.58),
        base_mean=lambda wc: config.get("baseMean_min", 5)
    threads: lambda wc: int(config.get("diff_threads", 1))
    shell:
        """
        {R_CMD} {SCRIPTS}/run_diff_from_counts.R \
          --matrix {input.matrix} \
          --outdir {OUTDIR}/diff \
          --matrix-name {wildcards.matrix} \
          --sample-table {SAMPLE_TABLE} \
          --contrast-file {CONTRAST_FILE} \
          --threads {params.threads_cfg} \
          --padj-cutoff {params.padj} \
          --lfc-cutoff {params.lfc} \
          --baseMean-min {params.base_mean} \
          {params.extra}
        """


rule annotate_de:
    input:
        diff_summary=lambda wc: f"{OUTDIR}/diff/{wc.matrix}/{wc.matrix}.diff.summary.txt"
    output:
        annotated=f"{OUTDIR}/annotate/{{matrix}}/{{matrix}}_{{case}}_vs_{{control}}/{{matrix}}_{{case}}_vs_{{control}}.annotated_DE_matrix.csv"
    params:
        de=lambda wc: f"{OUTDIR}/diff/{wc.matrix}/{wc.matrix}_{wc.case}_vs_{wc.control}/{wc.matrix}_{wc.case}_vs_{wc.control}.DE_matrix.csv",
        args=lambda wc: annotation_args(wc.matrix),
        padj=lambda wc: config.get("padj_cutoff", 0.05),
        lfc=lambda wc: config.get("lfc_cutoff", 0.58),
        base_mean=lambda wc: config.get("baseMean_min", 5)
    shell:
        """
        {R_CMD} {SCRIPTS}/run_annotate_de.R \
          --de-matrix {params.de} \
          --outdir {OUTDIR}/annotate/{wildcards.matrix}/{wildcards.matrix}_{wildcards.case}_vs_{wildcards.control} \
          --prefix {wildcards.matrix}_{wildcards.case}_vs_{wildcards.control} \
          --padj-cutoff {params.padj} \
          --lfc-cutoff {params.lfc} \
          --baseMean-min {params.base_mean} \
          {params.args}
        """


rule de_visuals:
    input:
        diff_summary=lambda wc: f"{OUTDIR}/diff/{wc.matrix}/{wc.matrix}.diff.summary.txt"
    output:
        summary=f"{OUTDIR}/visuals/{{matrix}}/{{matrix}}_{{case}}_vs_{{control}}/{{matrix}}_{{case}}_vs_{{control}}.visuals.summary.txt"
    params:
        de=lambda wc: f"{OUTDIR}/diff/{wc.matrix}/{wc.matrix}_{wc.case}_vs_{wc.control}/{wc.matrix}_{wc.case}_vs_{wc.control}.DE_matrix.csv",
        args=lambda wc: visual_args(wc.matrix),
        padj=lambda wc: config.get("padj_cutoff", 0.05),
        lfc=lambda wc: config.get("lfc_cutoff", 0.58),
        base_mean=lambda wc: config.get("baseMean_min", 5),
        label_top_n=lambda wc: config.get("label_top_n", 40),
        volcano_orientation=lambda wc: config.get("volcano_orientation", "classic"),
        gray_nonsig=lambda wc: str(config.get("gray_nonsig", True)).lower()
    shell:
        """
        {R_CMD} {SCRIPTS}/run_de_matrix_visuals.R \
          --de-matrix {params.de} \
          --outdir {OUTDIR}/visuals/{wildcards.matrix}/{wildcards.matrix}_{wildcards.case}_vs_{wildcards.control} \
          --prefix {wildcards.matrix}_{wildcards.case}_vs_{wildcards.control} \
          --padj-cutoff {params.padj} \
          --lfc-cutoff {params.lfc} \
          --baseMean-min {params.base_mean} \
          --label-top-n {params.label_top_n} \
          --volcano-orientation {params.volcano_orientation} \
          --gray-nonsig {params.gray_nonsig} \
          {params.args}
        """


rule pathway_from_de:
    input:
        diff_summary=lambda wc: f"{OUTDIR}/diff/{wc.matrix}/{wc.matrix}.diff.summary.txt"
    output:
        summary=f"{OUTDIR}/pathway/{{matrix}}/{{matrix}}_{{case}}_vs_{{control}}/{{matrix}}_{{case}}_vs_{{control}}.pathway.summary.txt"
    params:
        de=lambda wc: f"{OUTDIR}/diff/{wc.matrix}/{wc.matrix}_{wc.case}_vs_{wc.control}/{wc.matrix}_{wc.case}_vs_{wc.control}.DE_matrix.csv",
        args=lambda wc: pathway_args(wc.matrix)
    shell:
        """
        {R_CMD} {SCRIPTS}/run_pathway_from_de.R \
          --de-matrix {params.de} \
          --outdir {OUTDIR}/pathway/{wildcards.matrix}/{wildcards.matrix}_{wildcards.case}_vs_{wildcards.control} \
          --prefix {wildcards.matrix}_{wildcards.case}_vs_{wildcards.control} \
          --case {wildcards.case} \
          --control {wildcards.control} \
          {params.args}
        """


rule te_analysis:
    input:
        diff_summary=lambda wc: f"{OUTDIR}/diff/{wc.matrix}/{wc.matrix}.diff.summary.txt"
    output:
        summary=f"{OUTDIR}/te/{{matrix}}/{{matrix}}_{{case}}_vs_{{control}}/{{matrix}}_{{case}}_vs_{{control}}.te_analysis.summary.txt"
    params:
        de=lambda wc: f"{OUTDIR}/diff/{wc.matrix}/{wc.matrix}_{wc.case}_vs_{wc.control}/{wc.matrix}_{wc.case}_vs_{wc.control}.DE_matrix.csv",
        args=lambda wc: te_args(wc.matrix),
        padj=lambda wc: config.get("padj_cutoff", 0.05),
        lfc=lambda wc: config.get("lfc_cutoff", 0.58),
        base_mean=lambda wc: config.get("baseMean_min", 5)
    shell:
        """
        {R_CMD} {SCRIPTS}/run_te_analysis.R \
          --de-matrix {params.de} \
          --outdir {OUTDIR}/te/{wildcards.matrix}/{wildcards.matrix}_{wildcards.case}_vs_{wildcards.control} \
          --prefix {wildcards.matrix}_{wildcards.case}_vs_{wildcards.control} \
          --padj-cutoff {params.padj} \
          --lfc-cutoff {params.lfc} \
          --baseMean-min {params.base_mean} \
          {params.args}
        """
