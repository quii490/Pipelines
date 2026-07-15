nextflow.enable.dsl=2

/*
 * CLI values are strings in recent Nextflow releases.  In Groovy, the
 * non-empty string "false" is truthy, so normalize every feature switch
 * before it controls a workflow branch.
 */
def enabled(value) {
    if (value instanceof Boolean) return value
    return value != null && value.toString().trim().toLowerCase() in ['true', '1', 'yes', 'y']
}

process CHECK_REFS {
    tag "${params.species}"
    publishDir "${params.outdir}/00_refcheck", mode: 'copy'

    output:
    path 'refcheck.ok'

    script:
    def required = [
        params.genome_fasta,
        params.gtf_genes,
        params.salmon_index,
        params.gtf_te,
        params.telocal_index,
        params.salmonte_dir
    ].findAll { it != null && it.toString().trim() != '' }

    if ((params.aligner ?: 'star') == 'hisat2' && enabled(params.run_star_fc)) {
        required << params.hisat2_index
    }
    if (((params.aligner ?: 'star') == 'star' && enabled(params.run_star_fc)) || enabled(params.run_tecount) || enabled(params.run_telocal) || enabled(params.run_tetranscripts)) {
        required << params.star_index
    }

    if (params.blacklist) {
        required << params.blacklist
    }

    if (enabled(params.run_rediscoverte) && params.species == 'hg38') {
        required << params.rediscoverte_salmon_index
        required << params.rediscoverte_rollup_dir
        required << params.rediscoverte_dir
    }
    if (enabled(params.run_rediscoverte_rollup)) {
        required << params.rediscoverte_rollup_conda_prefix
    }

    def checks = required.collect { "[[ -e '${it}' ]] || { echo 'Missing: ${it}' >&2; exit 1; }" }.join('\n')

    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/00_refcheck/logs"
    mkdir -p "\${LOG_DIR}"
    {
      echo "[CHECK_REFS] species=${params.species}"
      ${checks}
      echo OK > refcheck.ok
      echo "[CHECK_REFS] OK"
    } 2>&1 | tee "\${LOG_DIR}/refcheck.${params.species}.log"
    """
}

/********************
 * FASTP
 ********************/
process FASTP_PE {
    tag "$sample"
    publishDir "${params.outdir}/01_fastp", mode: 'copy'

    input:
    tuple val(sample), val(condition), val(replicate), path(r1), path(r2)

    output:
    tuple val(sample), val(condition), val(replicate), path("${sample}.R1.clean.fastq.gz"), path("${sample}.R2.clean.fastq.gz")

    script:
    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/01_fastp/logs"
    mkdir -p "\${LOG_DIR}"
    fastp \
      -i ${r1} \
      -I ${r2} \
      -o ${sample}.R1.clean.fastq.gz \
      -O ${sample}.R2.clean.fastq.gz \
      -h ${sample}.fastp.html \
      -j ${sample}.fastp.json \
      2>&1 | tee "\${LOG_DIR}/${sample}.fastp.log"
    """
}

process FASTP_SE {
    tag "$sample"
    publishDir "${params.outdir}/01_fastp", mode: 'copy'

    input:
    tuple val(sample), val(condition), val(replicate), path(r1)

    output:
    tuple val(sample), val(condition), val(replicate), path("${sample}.SE.clean.fastq.gz")

    script:
    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/01_fastp/logs"
    mkdir -p "\${LOG_DIR}"
    fastp \
      -i ${r1} \
      -o ${sample}.SE.clean.fastq.gz \
      -h ${sample}.fastp.html \
      -j ${sample}.fastp.json \
      2>&1 | tee "\${LOG_DIR}/${sample}.fastp.log"
    """
}

workflow PREP_INPUTS {
    take:
    pe_in
    se_in

    main:
    if (enabled(params.run_fastp)) {
        pe_reads = FASTP_PE(pe_in)
        se_reads = FASTP_SE(se_in)
    } else {
        pe_reads = pe_in
        se_reads = se_in
    }

    emit:
    pe_reads
    se_reads
}

/********************
 * FastQC (raw-read QC; fastp already reports post-filter quality)
 ********************/
process FASTQC_PE {
    tag "$sample"
    publishDir "${params.outdir}/01_fastqc", mode: 'copy'

    input:
    tuple val(sample), val(condition), val(replicate), path(r1), path(r2)

    output:
    path "${sample}_fastqc/*"

    script:
    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/01_fastqc/logs"
    mkdir -p "\${LOG_DIR}" "${sample}_fastqc"
    export LD_LIBRARY_PATH="/path/to/.conda/envs/rnaseq/lib:/opt/conda/lib:\${LD_LIBRARY_PATH:-}"
    if [[ -n "\${CONDA_PREFIX:-}" ]]; then
      export LD_LIBRARY_PATH="\${CONDA_PREFIX}/lib:\${LD_LIBRARY_PATH}"
    fi
    export JAVA_TOOL_OPTIONS="\${JAVA_TOOL_OPTIONS:-} -Djava.awt.headless=true"
    fastqc --threads ${task.cpus} --outdir "${sample}_fastqc" ${r1} ${r2} \
      2>&1 | tee "\${LOG_DIR}/${sample}.fastqc.log"
    html_files=("${sample}_fastqc/"*"_fastqc.html")
    zip_files=("${sample}_fastqc/"*"_fastqc.zip")
    [[ \${#html_files[@]} -eq 2 && \${#zip_files[@]} -eq 2 ]]
    for f in "\${html_files[@]}" "\${zip_files[@]}"; do test -s "\${f}"; done
    """
}

process FASTQC_SE {
    tag "$sample"
    publishDir "${params.outdir}/01_fastqc", mode: 'copy'

    input:
    tuple val(sample), val(condition), val(replicate), path(r1)

    output:
    path "${sample}_fastqc/*"

    script:
    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/01_fastqc/logs"
    mkdir -p "\${LOG_DIR}" "${sample}_fastqc"
    export LD_LIBRARY_PATH="/path/to/.conda/envs/rnaseq/lib:/opt/conda/lib:\${LD_LIBRARY_PATH:-}"
    if [[ -n "\${CONDA_PREFIX:-}" ]]; then
      export LD_LIBRARY_PATH="\${CONDA_PREFIX}/lib:\${LD_LIBRARY_PATH}"
    fi
    export JAVA_TOOL_OPTIONS="\${JAVA_TOOL_OPTIONS:-} -Djava.awt.headless=true"
    fastqc --threads ${task.cpus} --outdir "${sample}_fastqc" ${r1} \
      2>&1 | tee "\${LOG_DIR}/${sample}.fastqc.log"
    html_files=("${sample}_fastqc/"*"_fastqc.html")
    zip_files=("${sample}_fastqc/"*"_fastqc.zip")
    [[ \${#html_files[@]} -eq 1 && \${#zip_files[@]} -eq 1 ]]
    for f in "\${html_files[@]}" "\${zip_files[@]}"; do test -s "\${f}"; done
    """
}

/********************
 * STAR gene
 ********************/
process STAR_ALIGN_GENE_PE {
    tag "$sample"
    publishDir "${params.outdir}/02_gene_star", mode: 'copy'

    input:
    tuple val(sample), val(condition), val(replicate), path(r1), path(r2)
    val blacklist

    output:
    tuple val(sample), val('PE'), val(condition), val(replicate), path("${sample}.gene.bam"), path("${sample}.gene.bam.bai")

    script:
    def hasBlacklist = blacklist && blacklist.toString().trim()

    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/02_gene_star/logs"
    mkdir -p "\${LOG_DIR}"

    STAR \
      --runThreadN ${task.cpus} \
      --genomeDir ${params.star_index} \
      --readFilesIn ${r1} ${r2} \
      --readFilesCommand zcat \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode GeneCounts \
      --outFilterMultimapNmax ${params.star_gene_multimap_nmax} \
      --winAnchorMultimapNmax ${params.star_gene_win_anchor_multimap_nmax} \
      --outFilterMismatchNmax ${params.star_gene_mismatch_nmax} \
      --outFilterMismatchNoverReadLmax ${params.star_gene_mismatch_nover_read_lmax} \
      --outSAMattributes NH HI AS nM XS \
      --outSAMattrRGline ID:${sample} SM:${sample} LB:RNAseq_gene PL:ILLUMINA PU:${sample} \
      --outFileNamePrefix ${sample}.gene. \
      2>&1 | tee "\${LOG_DIR}/${sample}.gene.STAR.console.log"

    test -s ${sample}.gene.Aligned.sortedByCoord.out.bam

    if [[ -n "${hasBlacklist ? 'yes' : ''}" ]]; then
        bedtools intersect \
          -v \
          -abam ${sample}.gene.Aligned.sortedByCoord.out.bam \
          -b ${blacklist} \
          > ${sample}.gene.filtered.bam

        test -s ${sample}.gene.filtered.bam

        samtools sort -@ ${task.cpus} \
          -o ${sample}.gene.bam \
          ${sample}.gene.filtered.bam
    else
        samtools sort -@ ${task.cpus} \
          -o ${sample}.gene.bam \
          ${sample}.gene.Aligned.sortedByCoord.out.bam
    fi

    samtools index -@ ${task.cpus} ${sample}.gene.bam

    for f in \
      ${sample}.gene.Log.final.out \
      ${sample}.gene.Log.out \
      ${sample}.gene.Log.progress.out \
      ${sample}.gene.SJ.out.tab \
      ${sample}.gene.ReadsPerGene.out.tab; do
      [[ -f "\${f}" ]] && cp -f "\${f}" "\${LOG_DIR}/"
    done
    """
}

process STAR_ALIGN_GENE_SE {
    tag "$sample"
    publishDir "${params.outdir}/02_gene_star", mode: 'copy'

    input:
    tuple val(sample), val(condition), val(replicate), path(r1)
    val blacklist

    output:
    tuple val(sample), val('SE'), val(condition), val(replicate), path("${sample}.gene.bam"), path("${sample}.gene.bam.bai")

    script:
    def hasBlacklist = blacklist && blacklist.toString().trim()

    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/02_gene_star/logs"
    mkdir -p "\${LOG_DIR}"

    STAR \
      --runThreadN ${task.cpus} \
      --genomeDir ${params.star_index} \
      --readFilesIn ${r1} \
      --readFilesCommand zcat \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode GeneCounts \
      --outFilterMultimapNmax ${params.star_gene_multimap_nmax} \
      --winAnchorMultimapNmax ${params.star_gene_win_anchor_multimap_nmax} \
      --outFilterMismatchNmax ${params.star_gene_mismatch_nmax} \
      --outFilterMismatchNoverReadLmax ${params.star_gene_mismatch_nover_read_lmax} \
      --outSAMattributes NH HI AS nM XS \
      --outSAMattrRGline ID:${sample} SM:${sample} LB:RNAseq_gene PL:ILLUMINA PU:${sample} \
      --outFileNamePrefix ${sample}.gene. \
      2>&1 | tee "\${LOG_DIR}/${sample}.gene.STAR.console.log"

    test -s ${sample}.gene.Aligned.sortedByCoord.out.bam

    if [[ -n "${hasBlacklist ? 'yes' : ''}" ]]; then
        bedtools intersect \
          -v \
          -abam ${sample}.gene.Aligned.sortedByCoord.out.bam \
          -b ${blacklist} \
          > ${sample}.gene.filtered.bam

        test -s ${sample}.gene.filtered.bam

        samtools sort -@ ${task.cpus} \
          -o ${sample}.gene.bam \
          ${sample}.gene.filtered.bam
    else
        samtools sort -@ ${task.cpus} \
          -o ${sample}.gene.bam \
          ${sample}.gene.Aligned.sortedByCoord.out.bam
    fi

    samtools index -@ ${task.cpus} ${sample}.gene.bam

    for f in \
      ${sample}.gene.Log.final.out \
      ${sample}.gene.Log.out \
      ${sample}.gene.Log.progress.out \
      ${sample}.gene.SJ.out.tab \
      ${sample}.gene.ReadsPerGene.out.tab; do
      [[ -f "\${f}" ]] && cp -f "\${f}" "\${LOG_DIR}/"
    done
    """
}


process HISAT2_ALIGN_GENE_PE {
    tag "$sample"
    publishDir "${params.outdir}/02_gene_hisat2", mode: 'copy'

    input:
    tuple val(sample), val(condition), val(replicate), path(r1), path(r2)
    val blacklist

    output:
    tuple val(sample), val('PE'), val(condition), val(replicate), path("${sample}.gene.bam"), path("${sample}.gene.bam.bai")

    script:
    def hasBlacklist = blacklist && blacklist.toString().trim()

    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/02_gene_hisat2/logs"
    mkdir -p "\${LOG_DIR}"

    hisat2       -p ${task.cpus}       --dta       --rg-id ${sample}       --rg SM:${sample}       --rg LB:RNAseq_gene       --rg PL:ILLUMINA       --rg PU:${sample}       -x ${params.hisat2_index}       -1 ${r1}       -2 ${r2}       2> >(tee ${sample}.gene.hisat2.log "\${LOG_DIR}/${sample}.gene.hisat2.log" >&2)       | samtools view -@ ${task.cpus} -bS -       > ${sample}.gene.raw.bam

    test -s ${sample}.gene.raw.bam

    if [[ -n "${hasBlacklist ? 'yes' : ''}" ]]; then
        bedtools intersect           -v           -abam ${sample}.gene.raw.bam           -b ${blacklist}           > ${sample}.gene.filtered.bam

        test -s ${sample}.gene.filtered.bam

        samtools sort -@ ${task.cpus}           -o ${sample}.gene.bam           ${sample}.gene.filtered.bam
    else
        samtools sort -@ ${task.cpus}           -o ${sample}.gene.bam           ${sample}.gene.raw.bam
    fi

    samtools index -@ ${task.cpus} ${sample}.gene.bam
    echo "[HISAT2_ALIGN_GENE_PE] DONE sample=${sample}" | tee -a "\${LOG_DIR}/${sample}.gene.hisat2.log"
    """
}

process HISAT2_ALIGN_GENE_SE {
    tag "$sample"
    publishDir "${params.outdir}/02_gene_hisat2", mode: 'copy'

    input:
    tuple val(sample), val(condition), val(replicate), path(r1)
    val blacklist

    output:
    tuple val(sample), val('SE'), val(condition), val(replicate), path("${sample}.gene.bam"), path("${sample}.gene.bam.bai")

    script:
    def hasBlacklist = blacklist && blacklist.toString().trim()

    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/02_gene_hisat2/logs"
    mkdir -p "\${LOG_DIR}"

    hisat2       -p ${task.cpus}       --dta       --rg-id ${sample}       --rg SM:${sample}       --rg LB:RNAseq_gene       --rg PL:ILLUMINA       --rg PU:${sample}       -x ${params.hisat2_index}       -U ${r1}       2> >(tee ${sample}.gene.hisat2.log "\${LOG_DIR}/${sample}.gene.hisat2.log" >&2)       | samtools view -@ ${task.cpus} -bS -       > ${sample}.gene.raw.bam

    test -s ${sample}.gene.raw.bam

    if [[ -n "${hasBlacklist ? 'yes' : ''}" ]]; then
        bedtools intersect           -v           -abam ${sample}.gene.raw.bam           -b ${blacklist}           > ${sample}.gene.filtered.bam

        test -s ${sample}.gene.filtered.bam

        samtools sort -@ ${task.cpus}           -o ${sample}.gene.bam           ${sample}.gene.filtered.bam
    else
        samtools sort -@ ${task.cpus}           -o ${sample}.gene.bam           ${sample}.gene.raw.bam
    fi

    samtools index -@ ${task.cpus} ${sample}.gene.bam
    echo "[HISAT2_ALIGN_GENE_SE] DONE sample=${sample}" | tee -a "\${LOG_DIR}/${sample}.gene.hisat2.log"
    """
}

process MARKDUP_GENE_BAM {
    tag "$sample"
    publishDir "${params.outdir}/dedup", mode: 'copy'

    input:
    tuple val(sample), val(layout), val(condition), val(replicate), path(bam), path(bai)

    output:
    tuple val(sample), val(layout), val(condition), val(replicate),
          path("${sample}.gene.dedup.bam"),
          path("${sample}.gene.dedup.bam.bai")

    script:
    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/dedup/logs"
    mkdir -p "\${LOG_DIR}"
    echo "[MARKDUP_GENE_BAM] sample=${sample}" 2>&1 | tee "\${LOG_DIR}/${sample}.gene.markdup.log"

    bash ${projectDir}/scripts/prepare_picard_bam.sh \
      --bam ${bam} \
      --output ${sample}.gene.picard_input.bam \
      --sample ${sample} \
      --library RNAseq_gene \
      2>&1 | tee -a "\${LOG_DIR}/${sample}.gene.markdup.log"

    picard -Xmx${params.picard_markdup_java_heap} MarkDuplicates \
      I=${sample}.gene.picard_input.bam \
      O=${sample}.gene.dedup.bam \
      M=${sample}.gene.dedup.metrics.txt \
      REMOVE_DUPLICATES=false \
      ASSUME_SORTED=true \
      CREATE_INDEX=false \
      VALIDATION_STRINGENCY=LENIENT \
      2>&1 | tee -a "\${LOG_DIR}/${sample}.gene.markdup.log"

    rm -f ${sample}.gene.picard_input.bam
    samtools index -@ ${task.cpus} ${sample}.gene.dedup.bam ${sample}.gene.dedup.bam.bai

    mkdir -p ${params.outdir}/dedup
    cp ${sample}.gene.dedup.metrics.txt ${params.outdir}/dedup/${sample}.gene.dedup.metrics.txt
    """
}

process MARKDUP_TE_BAM {
    tag "$sample"
    publishDir "${params.outdir}/dedup", mode: 'copy'

    input:
    tuple val(sample), val(layout), val(condition), val(replicate), path(bam), path(bai)

    output:
    tuple val(sample), val(layout), val(condition), val(replicate),
          path("${sample}.te.dedup.bam"),
          path("${sample}.te.dedup.bam.bai")

    script:
    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/dedup/logs"
    mkdir -p "\${LOG_DIR}"
    echo "[MARKDUP_TE_BAM] sample=${sample}" 2>&1 | tee "\${LOG_DIR}/${sample}.te.markdup.log"

    bash ${projectDir}/scripts/prepare_picard_bam.sh \
      --bam ${bam} \
      --output ${sample}.te.picard_input.bam \
      --sample ${sample} \
      --library RNAseq_TE \
      2>&1 | tee -a "\${LOG_DIR}/${sample}.te.markdup.log"

    picard -Xmx${params.picard_markdup_java_heap} MarkDuplicates \
      I=${sample}.te.picard_input.bam \
      O=${sample}.te.dedup.bam \
      M=${sample}.te.dedup.metrics.txt \
      REMOVE_DUPLICATES=false \
      ASSUME_SORTED=true \
      CREATE_INDEX=false \
      VALIDATION_STRINGENCY=LENIENT \
      2>&1 | tee -a "\${LOG_DIR}/${sample}.te.markdup.log"

    rm -f ${sample}.te.picard_input.bam
    samtools index -@ ${task.cpus} ${sample}.te.dedup.bam ${sample}.te.dedup.bam.bai

    mkdir -p ${params.outdir}/dedup
    cp ${sample}.te.dedup.metrics.txt ${params.outdir}/dedup/${sample}.te.dedup.metrics.txt
    """
}

process MARKDUP_GENE_QC {
    tag "$sample"
    publishDir "${params.outdir}/dedup_qc", mode: 'copy'

    input:
    tuple val(sample), val(layout), val(condition), val(replicate), path(bam), path(bai)

    output:
    path "${sample}.gene.markdup.metrics.txt"

    script:
    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/dedup_qc/logs"
    mkdir -p "\${LOG_DIR}"
    echo "[MARKDUP_GENE_QC] sample=${sample} input=${bam}" 2>&1 | tee "\${LOG_DIR}/${sample}.gene.markdup_qc.log"

    bash ${projectDir}/scripts/prepare_picard_bam.sh \
      --bam ${bam} \
      --output ${sample}.gene.picard_input.bam \
      --sample ${sample} \
      --library RNAseq_gene \
      2>&1 | tee -a "\${LOG_DIR}/${sample}.gene.markdup_qc.log"

    picard -Xmx${params.picard_markdup_java_heap} MarkDuplicates \
      I=${sample}.gene.picard_input.bam \
      O=${sample}.gene.markdup_qc.tmp.bam \
      M=${sample}.gene.markdup.metrics.txt \
      REMOVE_DUPLICATES=false \
      ASSUME_SORTED=true \
      CREATE_INDEX=false \
      VALIDATION_STRINGENCY=LENIENT \
      2>&1 | tee -a "\${LOG_DIR}/${sample}.gene.markdup_qc.log"

    rm -f ${sample}.gene.markdup_qc.tmp.bam ${sample}.gene.picard_input.bam
    """
}

process MARKDUP_TE_QC {
    tag "$sample"
    publishDir "${params.outdir}/dedup_qc", mode: 'copy'

    input:
    tuple val(sample), val(layout), val(condition), val(replicate), path(bam), path(bai)

    output:
    path "${sample}.te.markdup.metrics.txt"

    script:
    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/dedup_qc/logs"
    mkdir -p "\${LOG_DIR}"
    echo "[MARKDUP_TE_QC] sample=${sample} input=${bam}" 2>&1 | tee "\${LOG_DIR}/${sample}.te.markdup_qc.log"

    bash ${projectDir}/scripts/prepare_picard_bam.sh \
      --bam ${bam} \
      --output ${sample}.te.picard_input.bam \
      --sample ${sample} \
      --library RNAseq_TE \
      2>&1 | tee -a "\${LOG_DIR}/${sample}.te.markdup_qc.log"

    picard -Xmx${params.picard_markdup_java_heap} MarkDuplicates \
      I=${sample}.te.picard_input.bam \
      O=${sample}.te.markdup_qc.tmp.bam \
      M=${sample}.te.markdup.metrics.txt \
      REMOVE_DUPLICATES=false \
      ASSUME_SORTED=true \
      CREATE_INDEX=false \
      VALIDATION_STRINGENCY=LENIENT \
      2>&1 | tee -a "\${LOG_DIR}/${sample}.te.markdup_qc.log"

    rm -f ${sample}.te.markdup_qc.tmp.bam ${sample}.te.picard_input.bam
    """
}

process PREP_REFFLAT {
    tag "gene_annotation"
    publishDir "${params.outdir}/rnaseq_metrics/reference", mode: 'copy'

    input:
    path gene_gtf

    output:
    path "generated.refFlat.txt"

    script:
    """
    set -euo pipefail
    python3 ${projectDir}/scripts/gtf_to_refflat.py \
      --input ${gene_gtf} \
      --output generated.refFlat.txt
    test -s generated.refFlat.txt
    """
}

process RNASEQ_METRICS {
    tag "$sample"
    publishDir "${params.outdir}/rnaseq_metrics", mode: 'copy'

    input:
    tuple val(sample), val(layout), val(condition), val(replicate), path(bam), path(bai)
    path ref_flat

    output:
    path "${sample}.rnaseq.metrics.txt"

    script:
    def metricStrand = params.rnaseq_metrics_strand ?: (
        params.strandedness == 'forward' ? 'FIRST_READ_TRANSCRIPTION_STRAND' :
        params.strandedness == 'reverse' ? 'SECOND_READ_TRANSCRIPTION_STRAND' :
        'NONE'
    )
    def rrnaArg = params.ribosomal_intervals ? "RIBOSOMAL_INTERVALS=${params.ribosomal_intervals}" : ''
    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/rnaseq_metrics/logs"
    mkdir -p "\${LOG_DIR}"
    echo "[RNASEQ_METRICS] sample=${sample} input=${bam} strand=${metricStrand}" 2>&1 | tee "\${LOG_DIR}/${sample}.rnaseq_metrics.log"

    picard -Xmx${params.rnaseq_metrics_java_heap} CollectRnaSeqMetrics \
      I=${bam} \
      O=${sample}.rnaseq.metrics.txt \
      REF_FLAT=${ref_flat} \
      STRAND_SPECIFICITY=${metricStrand} \
      ${rrnaArg} \
      VALIDATION_STRINGENCY=LENIENT \
      2>&1 | tee -a "\${LOG_DIR}/${sample}.rnaseq_metrics.log"

    test -s ${sample}.rnaseq.metrics.txt
    """
}

process FEATURECOUNTS_GENE {
    tag "$sample"
    publishDir "${params.outdir}/03_gene_featurecounts", mode: 'copy'

    input:
    tuple val(sample), val(layout), val(condition), val(replicate), path(bam), path(bai)

    output:
    path "${sample}.featureCounts.txt"

    script:
    def fcStrand = (
    params.strandedness == 'unstranded' ? '0' :
    params.strandedness == 'forward'    ? '1' :
    params.strandedness == 'reverse'    ? '2' :
    '0'
    )
    def fcPairArgs = (layout == 'PE') ? '-p --countReadPairs' : ''
    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/03_gene_featurecounts/logs"
    mkdir -p "\${LOG_DIR}"
    echo "[FEATURECOUNTS_GENE] sample=${sample} layout=${layout} strand=${fcStrand} pair_args='${fcPairArgs}' feature_type=exon gene_attr=gene_id" 2>&1 | tee "\${LOG_DIR}/${sample}.featureCounts.log"
    featureCounts \
        -T ${task.cpus} \
        ${fcPairArgs} \
        -s ${fcStrand} \
        -t exon \
        -g gene_id \
        -a ${params.gtf_genes} \
        -o ${sample}.featureCounts.txt \
        ${bam} \
        2>&1 | tee -a "\${LOG_DIR}/${sample}.featureCounts.log"
    [[ -f ${sample}.featureCounts.txt.summary ]] && cp -f ${sample}.featureCounts.txt.summary ${params.outdir}/03_gene_featurecounts/${sample}.featureCounts.txt.summary
    """
}

process STRINGTIE_ASSEMBLY {
    tag "$sample"
    publishDir "${params.outdir}/04_gene_stringtie", mode: 'copy'

    input:
    tuple val(sample), val(layout), val(condition), val(replicate), path(bam), path(bai)

    output:
    path "${sample}.stringtie.gtf"

    script:
    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/04_gene_stringtie/logs"
    mkdir -p "\${LOG_DIR}"
    stringtie ${bam} -p ${task.cpus} -G ${params.gtf_genes} -o ${sample}.stringtie.gtf \
      2>&1 | tee "\${LOG_DIR}/${sample}.stringtie.log"
    """
}

/********************
 * Salmon gene
 ********************/
process SALMON_GENE_PE {
    tag "$sample"
    publishDir "${params.outdir}/05_gene_salmon", mode: 'copy'

    input:
    tuple val(sample), val(condition), val(replicate), path(r1), path(r2)

    output:
    path "${sample}_salmon/quant.sf"

    script:
    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/05_gene_salmon/logs"
    mkdir -p "\${LOG_DIR}"
    salmon quant \
      -i ${params.salmon_index} \
      -l A \
      -1 ${r1} \
      -2 ${r2} \
      --validateMappings \
      --gcBias \
      --seqBias \
      -p ${task.cpus} \
      -o ${sample}_salmon \
      2>&1 | tee "\${LOG_DIR}/${sample}.salmon.log"
    if [[ -d ${sample}_salmon/logs ]]; then
      for f in ${sample}_salmon/logs/*; do
        [[ -f "\${f}" ]] && cp -f "\${f}" "\${LOG_DIR}/${sample}.salmon.\$(basename "\${f}")"
      done
    fi
    """
}

process SALMON_GENE_SE {
    tag "$sample"
    publishDir "${params.outdir}/05_gene_salmon", mode: 'copy'

    input:
    tuple val(sample), val(condition), val(replicate), path(r1)

    output:
    path "${sample}_salmon/quant.sf"

    script:
    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/05_gene_salmon/logs"
    mkdir -p "\${LOG_DIR}"
    salmon quant \
      -i ${params.salmon_index} \
      -l A \
      -r ${r1} \
      --validateMappings \
      --gcBias \
      --seqBias \
      -p ${task.cpus} \
      -o ${sample}_salmon \
      2>&1 | tee "\${LOG_DIR}/${sample}.salmon.log"
    if [[ -d ${sample}_salmon/logs ]]; then
      for f in ${sample}_salmon/logs/*; do
        [[ -f "\${f}" ]] && cp -f "\${f}" "\${LOG_DIR}/${sample}.salmon.\$(basename "\${f}")"
      done
    fi
    """
}

/********************
 * STAR TE
 ********************/
process STAR_ALIGN_TE_PE {
    tag "$sample"
    publishDir "${params.outdir}/06_te_star", mode: 'copy'

    input:
    tuple val(sample), val(condition), val(replicate), path(r1), path(r2)
    val blacklist

    output:
    tuple val(sample), val('PE'), val(condition), val(replicate), path("${sample}.te.bam"), path("${sample}.te.bam.bai")

    script:
    def hasBlacklist = blacklist && blacklist.toString().trim()

    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/06_te_star/logs"
    mkdir -p "\${LOG_DIR}"

    STAR \
      --runThreadN ${task.cpus} \
      --genomeDir ${params.star_index} \
      --readFilesIn ${r1} ${r2} \
      --readFilesCommand zcat \
      --outSAMtype BAM SortedByCoordinate \
      --outFilterMultimapNmax ${params.star_te_multimap_nmax} \
      --winAnchorMultimapNmax ${params.star_te_win_anchor_multimap_nmax} \
      --outSAMmultNmax ${params.star_te_out_sam_mult_nmax} \
      --outSAMprimaryFlag AllBestScore \
      --outFilterMismatchNmax ${params.star_te_mismatch_nmax} \
      --outFilterMismatchNoverReadLmax ${params.star_te_mismatch_nover_read_lmax} \
      --alignIntronMax ${params.star_te_align_intron_max} \
      --outSAMattributes NH HI AS nM \
      --outSAMattrRGline ID:${sample} SM:${sample} LB:RNAseq_TE PL:ILLUMINA PU:${sample} \
      --outFileNamePrefix ${sample}.te. \
      2>&1 | tee "\${LOG_DIR}/${sample}.te.STAR.console.log"

    test -s ${sample}.te.Aligned.sortedByCoord.out.bam

    samtools sort -@ ${task.cpus} \
          -o ${sample}.te.bam \
          ${sample}.te.Aligned.sortedByCoord.out.bam

    samtools index -@ ${task.cpus} ${sample}.te.bam

    for f in \
      ${sample}.te.Log.final.out \
      ${sample}.te.Log.out \
      ${sample}.te.Log.progress.out \
      ${sample}.te.SJ.out.tab; do
      [[ -f "\${f}" ]] && cp -f "\${f}" "\${LOG_DIR}/"
    done
    """
}

process STAR_ALIGN_TE_SE {
    tag "$sample"
    publishDir "${params.outdir}/06_te_star", mode: 'copy'

    input:
    tuple val(sample), val(condition), val(replicate), path(r1)
    val blacklist

    output:
    tuple val(sample), val('SE'), val(condition), val(replicate), path("${sample}.te.bam"), path("${sample}.te.bam.bai")

    script:
    def hasBlacklist = blacklist && blacklist.toString().trim()

    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/06_te_star/logs"
    mkdir -p "\${LOG_DIR}"

    STAR \
      --runThreadN ${task.cpus} \
      --genomeDir ${params.star_index} \
      --readFilesIn ${r1} \
      --readFilesCommand zcat \
      --outSAMtype BAM SortedByCoordinate \
      --outFilterMultimapNmax ${params.star_te_multimap_nmax} \
      --winAnchorMultimapNmax ${params.star_te_win_anchor_multimap_nmax} \
      --outSAMmultNmax ${params.star_te_out_sam_mult_nmax} \
      --outSAMprimaryFlag AllBestScore \
      --outFilterMismatchNmax ${params.star_te_mismatch_nmax} \
      --outFilterMismatchNoverReadLmax ${params.star_te_mismatch_nover_read_lmax} \
      --alignIntronMax ${params.star_te_align_intron_max} \
      --outSAMattributes NH HI AS nM \
      --outSAMattrRGline ID:${sample} SM:${sample} LB:RNAseq_TE PL:ILLUMINA PU:${sample} \
      --outFileNamePrefix ${sample}.te. \
      2>&1 | tee "\${LOG_DIR}/${sample}.te.STAR.console.log"

    test -s ${sample}.te.Aligned.sortedByCoord.out.bam

    samtools sort -@ ${task.cpus} \
          -o ${sample}.te.bam \
          ${sample}.te.Aligned.sortedByCoord.out.bam

    samtools index -@ ${task.cpus} ${sample}.te.bam

    for f in \
      ${sample}.te.Log.final.out \
      ${sample}.te.Log.out \
      ${sample}.te.Log.progress.out \
      ${sample}.te.SJ.out.tab; do
      [[ -f "\${f}" ]] && cp -f "\${f}" "\${LOG_DIR}/"
    done
    """
}

process TECOUNT_SAMPLE {
    tag "$sample"
    publishDir "${params.outdir}/07_tecount", mode: 'copy'
    conda "${projectDir}/envs/rediscoverte.yml"

    input:
    tuple val(sample), val(layout), val(condition), val(replicate), path(bam), path(bai)

    output:
    path "${sample}_tecount/*"

    script:
    def teToolStrandedness = params.strandedness == 'unstranded' ? 'no' : params.strandedness
    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/07_tecount/logs"
    mkdir -p "\${LOG_DIR}"
    mkdir -p ${sample}_tecount
    echo "[TECOUNT_SAMPLE] sample=${sample} strandedness=${params.strandedness} te_tool_strandedness=${teToolStrandedness}" 2>&1 | tee "\${LOG_DIR}/${sample}.TEcount.log"
    TEcount \
      -s ${teToolStrandedness} \
      -b ${bam} \
      -r ${params.rmsk} \
      --prefix ${sample}_tecount \
      --outdir ${sample}_tecount \
      2>&1 | tee -a "\${LOG_DIR}/${sample}.TEcount.log"
    """
}

process TE_QUANT {
    tag "TEtranscripts per-sample TEcount"
    publishDir "${params.outdir}/09_TEtranscripts",mode: 'copy'
    
    input:
    tuple val(sample), val(layout), val(condition), val(replicate), path(bam), path(bai)

    output:
    path "${sample}_tetranscripts/*" 

    script:
    def teToolStrandedness = params.strandedness == 'unstranded' ? 'no' : params.strandedness
    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/09_TEtranscripts/logs"
    mkdir -p "\${LOG_DIR}"
    mkdir -p ${sample}_tetranscripts
    echo "[TE_QUANT] sample=${sample} strandedness=${params.strandedness} te_tool_strandedness=${teToolStrandedness}" 2>&1 | tee "\${LOG_DIR}/${sample}.TEtranscripts.TEcount.log"
    TEcount -b ${bam} \
            --format BAM \
            --stranded ${teToolStrandedness} \
            --mode multi \
            --sortByPos \
            --GTF ${params.gtf_genes} \
            --TE ${params.gtf_te} \
            --project ${sample}_tetranscripts \
            --outdir ${sample}_tetranscripts \
            2>&1 | tee -a "\${LOG_DIR}/${sample}.TEtranscripts.TEcount.log"
    """
}

process TELOCAL_SAMPLE {
    tag "$sample"
    publishDir "${params.outdir}/08_telocal", mode: 'copy'

    input:
    tuple val(sample), val(layout), val(condition), val(replicate), path(bam), path(bai)

    output:
    path "${sample}_telocal.cntTable"

    script:
    def teToolStrandedness = params.strandedness == 'unstranded' ? 'no' : params.strandedness
    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/08_telocal/logs"
    mkdir -p "\${LOG_DIR}"
    echo "[TELOCAL_SAMPLE] sample=${sample} strandedness=${params.strandedness} te_tool_strandedness=${teToolStrandedness}" 2>&1 | tee "\${LOG_DIR}/${sample}.TElocal.log"
    TElocal \
      --sortByPos \
      -b ${bam} \
      --stranded ${teToolStrandedness} \
      --GTF ${params.gtf_genes} \
      --TE ${params.telocal_index} \
      --project ${sample}_telocal \
      2>&1 | tee -a "\${LOG_DIR}/${sample}.TElocal.log"
    """
}

/********************
 * Telescope locus-level TE quantification
 ********************/
process PREP_TELESCOPE_ANNOTATION {
    tag "telescope_annotation"

    input:
    path te_gtf

    output:
    path "telescope_annotation.gtf"

    script:
    """
    set -euo pipefail
    python3 ${projectDir}/scripts/make_telescope_gtf.py \
      --input ${te_gtf} \
      --output telescope_annotation.gtf
    test -s telescope_annotation.gtf
    """
}

process TELESCOPE_ASSIGN {
    tag "$sample"
    publishDir "${params.outdir}/12_telescope", mode: 'copy'

    input:
    tuple val(sample), val(layout), val(condition), val(replicate), path(bam), path(bai)
    path telescope_gtf

    output:
    path "${sample}_telescope_*"

    script:
    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/12_telescope/logs"
    mkdir -p "\${LOG_DIR}" telescope_work
    # Telescope 1.0.3 has a parallel-parser bug for coordinate-sorted BAMs:
    # reads left in a chromosome-region cache are emitted with code='cached'
    # and later cast to int.  Collate makes mates adjacent and lets Telescope
    # use its sequential reader, while the original coordinate-sorted BAM
    # remains the canonical TE BAM for all other modules.
    test -s ${bam}
    samtools quickcheck -v ${bam}
    samtools collate -@ ${task.cpus} -u -O ${bam} > ${sample}.telescope.collated.bam
    test -s ${sample}.telescope.collated.bam
    {
      echo "[TELESCOPE_ASSIGN] sample=${sample} input=${bam} telescope_input=${sample}.telescope.collated.bam ncpu=1"
      telescope assign \
        --attribute locus \
        --ncpu 1 \
        --outdir telescope_work \
        --exp_tag ${sample} \
        ${sample}.telescope.collated.bam telescope_annotation.gtf
    } 2>&1 | tee "\${LOG_DIR}/${sample}.telescope.log"

    shopt -s nullglob
    telescope_files=(telescope_work/${sample}_telescope_* telescope_work/${sample}-telescope_*)
    if (( \${#telescope_files[@]} == 0 )); then
      echo "[TELESCOPE_ASSIGN] ERROR: no Telescope output files found for ${sample}" >&2
      exit 1
    fi
    for f in "\${telescope_files[@]}"; do
      b=\$(basename "\${f}")
      b="\${b/${sample}-telescope_/${sample}_telescope_}"
      cp -f "\${f}" "\${b}"
    done
    test -s ${sample}_telescope_report.tsv
    """
}

/********************
 * REdiscoverTE
 ********************/
process REDISCOVERTE_QUANT_PE {
    tag "$sample"
    publishDir "${params.outdir}/10_rediscoverte", mode: 'copy'

    input:
    tuple val(sample), val(condition), val(replicate), path(r1), path(r2)

    output:
    path "${sample}_rediscoverte"

    script:
    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/10_rediscoverte/logs"
    mkdir -p "\${LOG_DIR}"
    salmon quant \
        --seqBias --gcBias \
        --index ${params.rediscoverte_salmon_index} \
        --libType A \
        --validateMappings \
        --threads ${task.cpus} \
        -1 ${r1} \
        -2 ${r2} \
        -o ${sample}_rediscoverte \
        2>&1 | tee "\${LOG_DIR}/${sample}.rediscoverte.salmon.log"
    if [[ -d ${sample}_rediscoverte/logs ]]; then
      for f in ${sample}_rediscoverte/logs/*; do
        [[ -f "\${f}" ]] && cp -f "\${f}" "\${LOG_DIR}/${sample}.rediscoverte.\$(basename "\${f}")"
      done
    fi
    """
}

process REDISCOVERTE_QUANT_SE {
    tag "$sample"
    publishDir "${params.outdir}/10_rediscoverte", mode: 'copy'

    input:
    tuple val(sample), val(condition), val(replicate), path(r1)

    output:
    path "${sample}_rediscoverte"

    script:
    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/10_rediscoverte/logs"
    mkdir -p "\${LOG_DIR}"
    salmon quant \
        --seqBias --gcBias \
        --index ${params.rediscoverte_salmon_index} \
        --libType A \
        --validateMappings \
        --threads ${task.cpus} \
        -r ${r1} \
        -o ${sample}_rediscoverte \
        2>&1 | tee "\${LOG_DIR}/${sample}.rediscoverte.salmon.log"
    if [[ -d ${sample}_rediscoverte/logs ]]; then
      for f in ${sample}_rediscoverte/logs/*; do
        [[ -f "\${f}" ]] && cp -f "\${f}" "\${LOG_DIR}/${sample}.rediscoverte.\$(basename "\${f}")"
      done
    fi
    """
}

process REDISCOVERTE_ROLLUP {
    publishDir "${params.outdir}/10_rediscoverte_rollup", mode: 'copy'
    conda "${params.rediscoverte_rollup_conda_prefix}"

    input:
    path quant_dirs

    output:
    path "rediscoverte_rollup/*"

    script:
    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/10_rediscoverte_rollup/logs"
    mkdir -p "\${LOG_DIR}"
    
    echo -e "sample\tquant_sf_path" > metadata.tsv
    for dir in *_rediscoverte; do
        sample=\${dir%_rediscoverte}
        test -s "\${dir}/quant.sf"
        echo -e "\${sample}\t\$(pwd)/\${dir}/quant.sf" >> metadata.tsv
    done

    Rscript ${params.rediscoverte_dir}/rollup.R \
        --metadata metadata.tsv \
        --datadir ${params.rediscoverte_rollup_dir} \
        --assembly ${params.species} \
        --nozero \
        --threads ${task.cpus} \
        --outdir rediscoverte_rollup \
        2>&1 | tee "\${LOG_DIR}/rediscoverte_rollup.log"
    """
}

/********************
 * SalmonTE
 ********************/
process SALMONTE_SAMPLE_PE {
    tag "$sample"
    publishDir "${params.outdir}/11_salmonte", mode: 'copy'

    input:
    tuple val(sample), val(condition), val(replicate), path(r1), path(r2)

    output:
    path "${sample}_SalmonTE_output"

    script:
    def refArg = params.salmonte_ref ? "--reference=${params.salmonte_ref}" : ''
    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/11_salmonte/logs"
    mkdir -p "\${LOG_DIR}"
    mkdir -p ${sample}_salmonte_input
    ln -sf "\$(realpath ${r1})" ${sample}_salmonte_input/${sample}_R1.fastq.gz
    ln -sf "\$(realpath ${r2})" ${sample}_salmonte_input/${sample}_R2.fastq.gz

    bash ${projectDir}/scripts/run_salmonte.sh \
      --salmonte-dir ${params.salmonte_dir} \
      --reference ${params.salmonte_ref} \
      --input-dir ${sample}_salmonte_input \
      --conda-prefix /path/to/.conda/envs/salmonte \
      --exprtype count \
      --outdir ${sample}_SalmonTE_output \
      2>&1 | tee "\${LOG_DIR}/${sample}.SalmonTE.log"
    """
}

process SALMONTE_SAMPLE_SE {
    tag "$sample"
    publishDir "${params.outdir}/11_salmonte", mode: 'copy'

    input:
    tuple val(sample), val(condition), val(replicate), path(r1)

    output:
    path "${sample}_SalmonTE_output"

    script:
    def refArg = params.salmonte_ref ? "--reference=${params.salmonte_ref}" : ''
    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/11_salmonte/logs"
    mkdir -p "\${LOG_DIR}"
    mkdir -p ${sample}_salmonte_input
    ln -sf "\$(realpath ${r1})" ${sample}_salmonte_input/${sample}.fastq.gz

    bash ${projectDir}/scripts/run_salmonte.sh \
      --salmonte-dir ${params.salmonte_dir} \
      --reference ${params.salmonte_ref} \
      --input-dir ${sample}_salmonte_input \
      --conda-prefix /path/to/.conda/envs/salmonte \
      --exprtype count \
      --outdir ${sample}_SalmonTE_output \
      2>&1 | tee "\${LOG_DIR}/${sample}.SalmonTE.log"
    """
}

process MULTIQC_REPORT {
    tag "multiqc"
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    val triggers

    output:
    path "multiqc_report.html"
    path "multiqc_data", optional: true

    script:
    def triggerCount = triggers instanceof Collection ? triggers.size() : 1
    """
    set -euo pipefail
    LOG_DIR="${params.outdir}/multiqc/logs"
    mkdir -p "\${LOG_DIR}"

    MULTIQC_CMD="${params.multiqc_cmd}"
    if [[ ! -x "\${MULTIQC_CMD}" ]]; then
      MULTIQC_CMD="\$(command -v multiqc || true)"
    fi

    if [[ -z "\${MULTIQC_CMD}" || ! -x "\${MULTIQC_CMD}" ]]; then
      {
        echo "[MULTIQC_REPORT] ERROR: multiqc not found in PATH"
        echo "[MULTIQC_REPORT] Install multiqc, set --multiqc_cmd with --extra, or keep --run-multiqc false."
      } 2>&1 | tee "\${LOG_DIR}/multiqc.log"
      exit 127
    fi

    echo "[MULTIQC_REPORT] trigger_count=${triggerCount} scan_dir=${params.outdir} multiqc_cmd=\${MULTIQC_CMD}" 2>&1 | tee "\${LOG_DIR}/multiqc.log"
    timeout --signal=TERM ${params.multiqc_time} "\${MULTIQC_CMD}" ${params.outdir} \
      --outdir . \
      --filename multiqc_report.html \
      --force \
      --dirs \
      --dirs-depth 3 \
      2>&1 | tee -a "\${LOG_DIR}/multiqc.log"
    """
}

workflow {
    // --- Samplesheet validation and channel setup ---
    if( !params.samplesheet ) exit 1, 'Please provide --samplesheet'

    /*
     * Read samplesheet
     * Required columns:
     * sample,layout,condition,replicate,r1,r2
     */
    Channel
    .fromPath(params.samplesheet)
    .splitCsv(header:true)
    .map { row ->
        def layout = row.layout.toString().trim().toUpperCase()
        assert layout in ['PE','SE'] : "Invalid layout for sample ${row.sample}: ${row.layout} (must be PE or SE)"
        def r1 = file(row.r1)
        def r2 = row.r2 && row.r2.toString().trim() ? file(row.r2) : null

        if( layout == 'PE' && !r2 )
            throw new IllegalArgumentException("Sample ${row.sample} is PE but r2 is missing in samplesheet")
        if( layout == 'SE' && r2 )
            log.warn("Sample ${row.sample} is SE but r2 is provided; it will be ignored")

        tuple(row.sample.toString(), layout, row.condition ?: 'NA', row.replicate ?: 'NA', r1, r2)
    }
    .set { ch_samples }

    // `ch_samples` is a queue channel.  Two independent filters would compete
    // for records, so route each row once before any downstream duplication.
    ch_samples
      .branch { sample, layout, condition, replicate, r1, r2 ->
        pe: layout == 'PE'
        se: layout == 'SE'
      }
      .set { reads_by_layout }

    ch_pe = reads_by_layout.pe
      .map { sample, layout, condition, replicate, r1, r2 -> tuple(sample, condition, replicate, r1, r2) }

    ch_se = reads_by_layout.se
      .map { sample, layout, condition, replicate, r1, r2 -> tuple(sample, condition, replicate, r1) }

    CHECK_REFS()
    multiqc_triggers = CHECK_REFS.out.map { "CHECK_REFS" }

    // Most analysis modules share the same cleaned-read channel.  FastQC is
    // deliberately run on raw reads, so fork only that branch when it is
    // requested.  Avoid declaring unconsumed multiMap outputs: they can
    // back-pressure a queue channel in small, module-specific reruns.
    def needs_prepared_reads = enabled(params.run_fastp) ||
                              enabled(params.run_star_fc) ||
                              enabled(params.run_salmon) ||
                              enabled(params.run_tecount) ||
                              enabled(params.run_telocal) ||
                              enabled(params.run_tetranscripts) ||
                              enabled(params.run_telescope) ||
                              enabled(params.run_rediscoverte) ||
                              enabled(params.run_salmonte)

    def prep_pe_input = ch_pe
    def prep_se_input = ch_se
    def prepared_pe = null
    def prepared_se = null

    if (needs_prepared_reads && enabled(params.run_fastqc)) {
        ch_pe
          .multiMap { sample, condition, replicate, r1, r2 ->
            prep: tuple(sample, condition, replicate, r1, r2)
            fastqc: tuple(sample, condition, replicate, r1, r2)
          }
          .set { pe_inputs }
        ch_se
          .multiMap { sample, condition, replicate, r1 ->
            prep: tuple(sample, condition, replicate, r1)
            fastqc: tuple(sample, condition, replicate, r1)
          }
          .set { se_inputs }
        prep_pe_input = pe_inputs.prep
        prep_se_input = se_inputs.prep
        fastqc_pe = FASTQC_PE(pe_inputs.fastqc)
        fastqc_se = FASTQC_SE(se_inputs.fastqc)
        multiqc_triggers = multiqc_triggers.mix(fastqc_pe.map { "FASTQC_PE" })
                                         .mix(fastqc_se.map { "FASTQC_SE" })
    } else if (enabled(params.run_fastqc)) {
        fastqc_pe = FASTQC_PE(ch_pe)
        fastqc_se = FASTQC_SE(ch_se)
        multiqc_triggers = multiqc_triggers.mix(fastqc_pe.map { "FASTQC_PE" })
                                         .mix(fastqc_se.map { "FASTQC_SE" })
    }

    if (needs_prepared_reads) {
        if (enabled(params.run_fastp)) {
            prepared = PREP_INPUTS(prep_pe_input, prep_se_input)
            prepared_pe = prepared.pe_reads
            prepared_se = prepared.se_reads
        } else {
            // Do not pass channels through an empty subworkflow.  Nextflow
            // 26 can leave that pass-through branch inactive in tool-only runs.
            prepared_pe = prep_pe_input
            prepared_se = prep_se_input
        }
    }
    def blacklistVal = params.blacklist ? params.blacklist.toString() : ''

    if (enabled(params.run_star_fc)) {
        if ((params.aligner ?: 'star') == 'hisat2') {
            gene_bams = HISAT2_ALIGN_GENE_PE(prepared_pe, blacklistVal)
                            .mix(HISAT2_ALIGN_GENE_SE(prepared_se, blacklistVal))
        } else {
            gene_bams = STAR_ALIGN_GENE_PE(prepared_pe, blacklistVal)
                            .mix(STAR_ALIGN_GENE_SE(prepared_se, blacklistVal))
        }
        if (enabled(params.run_markdup_qc)) {
            gene_markdup_qc = MARKDUP_GENE_QC(gene_bams)
            multiqc_triggers = multiqc_triggers.mix(gene_markdup_qc.map { "MARKDUP_GENE_QC" })
        }
        if (enabled(params.run_rnaseq_metrics)) {
            if (params.ref_flat) {
                ref_flat_for_metrics = Channel.value(file(params.ref_flat))
            } else {
                generated_ref_flat = PREP_REFFLAT(Channel.value(file(params.gtf_genes)))
                ref_flat_for_metrics = generated_ref_flat
            }
            rnaseq_metrics = RNASEQ_METRICS(gene_bams, ref_flat_for_metrics)
            multiqc_triggers = multiqc_triggers.mix(rnaseq_metrics.map { "RNASEQ_METRICS" })
        }
        gene_bams_for_count = enabled(params.run_dedup) ? MARKDUP_GENE_BAM(gene_bams) : gene_bams
        if (enabled(params.run_dedup)) {
            multiqc_triggers = multiqc_triggers.mix(gene_bams_for_count.map { "MARKDUP_GENE_BAM" })
        }
        gene_fc = FEATURECOUNTS_GENE(gene_bams_for_count)
        multiqc_triggers = multiqc_triggers.mix(gene_fc.map { "FEATURECOUNTS_GENE" })
        if (enabled(params.run_stringtie)) {
            stringtie_gtf = STRINGTIE_ASSEMBLY(gene_bams_for_count)
            multiqc_triggers = multiqc_triggers.mix(stringtie_gtf.map { "STRINGTIE_ASSEMBLY" })
        }
    }

    if (enabled(params.run_salmon)) {
        salmon_gene_pe = SALMON_GENE_PE(prepared_pe)
        salmon_gene_se = SALMON_GENE_SE(prepared_se)
        multiqc_triggers = multiqc_triggers.mix(salmon_gene_pe.map { "SALMON_GENE_PE" })
                                         .mix(salmon_gene_se.map { "SALMON_GENE_SE" })
    }

    if (enabled(params.run_tecount) || enabled(params.run_telocal) || enabled(params.run_tetranscripts) || enabled(params.run_telescope)) {
        te_bams = STAR_ALIGN_TE_PE(prepared_pe, blacklistVal)
                      .mix(STAR_ALIGN_TE_SE(prepared_se, blacklistVal))
        if (enabled(params.run_markdup_qc)) {
            te_markdup_qc = MARKDUP_TE_QC(te_bams)
            multiqc_triggers = multiqc_triggers.mix(te_markdup_qc.map { "MARKDUP_TE_QC" })
        }
        te_bams_for_count = enabled(params.run_dedup) ? MARKDUP_TE_BAM(te_bams) : te_bams
        if (enabled(params.run_dedup)) {
            multiqc_triggers = multiqc_triggers.mix(te_bams_for_count.map { "MARKDUP_TE_BAM" })
        }
        if (enabled(params.run_tecount)) {
            tecount_out = TECOUNT_SAMPLE(te_bams_for_count)
            multiqc_triggers = multiqc_triggers.mix(tecount_out.map { "TECOUNT_SAMPLE" })
        }
        if (enabled(params.run_tetranscripts)) {
            tetranscripts_out = TE_QUANT(te_bams_for_count)
            multiqc_triggers = multiqc_triggers.mix(tetranscripts_out.map { "TE_QUANT" })
        }
        if (enabled(params.run_telocal)) {
            telocal_out = TELOCAL_SAMPLE(te_bams_for_count)
            multiqc_triggers = multiqc_triggers.mix(telocal_out.map { "TELOCAL_SAMPLE" })
        }
        if (enabled(params.run_telescope)) {
            telescope_annotation = PREP_TELESCOPE_ANNOTATION(Channel.value(file(params.gtf_te)))
            telescope_out = TELESCOPE_ASSIGN(te_bams_for_count, telescope_annotation)
            multiqc_triggers = multiqc_triggers.mix(telescope_out.map { "TELESCOPE_ASSIGN" })
        }
    }

    if (enabled(params.run_rediscoverte)) {
        if (params.species == 'hg38') {
            quant_dirs = REDISCOVERTE_QUANT_PE(prepared_pe)
                           .mix(REDISCOVERTE_QUANT_SE(prepared_se))
            if (enabled(params.run_rediscoverte_rollup)) {
                rediscoverte_rollup = REDISCOVERTE_ROLLUP(quant_dirs.collect())
                multiqc_triggers = multiqc_triggers.mix(rediscoverte_rollup.map { "REDISCOVERTE_ROLLUP" })
            }
        } else {
            log.warn "REdiscoverTE is only enabled for hg38; current species=${params.species}, skipping."
        }
    }

    if (enabled(params.run_salmonte)) {
        salmonte_pe = SALMONTE_SAMPLE_PE(prepared_pe)
        salmonte_se = SALMONTE_SAMPLE_SE(prepared_se)
        multiqc_triggers = multiqc_triggers.mix(salmonte_pe.map { "SALMONTE_SAMPLE_PE" })
                                         .mix(salmonte_se.map { "SALMONTE_SAMPLE_SE" })
    }

    if (enabled(params.run_multiqc)) {
        MULTIQC_REPORT(multiqc_triggers.collect())
    }
}
