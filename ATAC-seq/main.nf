nextflow.enable.dsl=2

process CHECK_REFS {
    tag "${params.species}"
    publishDir "${params.outdir}/00_refcheck", mode: 'copy'
    output: path 'refcheck.ok'
    script:
    def required = [params.genome_fasta, params.bowtie2_index, params.chrom_sizes]
    if (params.remove_blacklist) required << params.blacklist
    if (params.run_downstream && params.gtf_genes) required << params.gtf_genes
    if (params.run_tss_enrich && params.tss_bed) required << params.tss_bed
    if (params.run_gene_body_profile && params.gene_body_bed) required << params.gene_body_bed
    if (params.run_te_heatmap && params.te_bed) required << params.te_bed
    if (params.run_te_relaxed_tracks && params.te_remove_blacklist) required << params.blacklist
    if (params.run_footprinting && params.motif_meme) required << params.motif_meme
    if (required.any { ref_path -> !ref_path }) {
        throw new IllegalArgumentException("Reference configuration is incomplete for species '${params.species}'")
    }
    def checks = required.collect { value ->
        def ref = value.toString().replace("'", "'\\\"'\\\"'")
        "[[ -e '${ref}' || -e '${ref}.1.bt2' || -e '${ref}.1.bt2l' ]] || { echo 'Missing: ${ref}' >&2; exit 1; }"
    }.join('\n')
    """
    set -euo pipefail
    ${checks}
    echo OK > refcheck.ok
    """
}

process FASTP_PE {
    tag "$sample"
    publishDir "${params.outdir}/01_fastp", mode: 'copy'
    input:
    path ref_ok
    tuple val(sample), val(layout), val(condition), val(replicate), path(r1), path(r2)
    output:
    tuple val(sample), val(layout), val(condition), val(replicate), path("${sample}.R1.clean.fastq.gz"), path("${sample}.R2.clean.fastq.gz")
    path "${sample}.fastp.html", optional: true
    path "${sample}.fastp.json", optional: true
    script:
    def r1Files = (r1 instanceof List) ? r1 : [r1]
    def r2Files = (r2 instanceof List) ? r2 : [r2]
    def r1Arg = r1Files.join(' ')
    def r2Arg = r2Files.join(' ')
    def prepR1 = (r1Files.size() == 1) ? "ln -sf ${r1Files[0]} ${sample}.R1.merged.fastq.gz" : "cat ${r1Arg} > ${sample}.R1.merged.fastq.gz"
    def prepR2 = (r2Files.size() == 1) ? "ln -sf ${r2Files[0]} ${sample}.R2.merged.fastq.gz" : "cat ${r2Arg} > ${sample}.R2.merged.fastq.gz"
    """
    set -euo pipefail
    ${prepR1}
    ${prepR2}
    timeout ${params.fastp_timeout} fastp --thread ${task.cpus} \
          -i ${sample}.R1.merged.fastq.gz -I ${sample}.R2.merged.fastq.gz \
          -o ${sample}.R1.clean.fastq.gz -O ${sample}.R2.clean.fastq.gz \
          -h ${sample}.fastp.html -j ${sample}.fastp.json
    """
}

process FASTP_SE {
    tag "$sample"
    publishDir "${params.outdir}/01_fastp", mode: 'copy'
    input:
    path ref_ok
    tuple val(sample), val(layout), val(condition), val(replicate), path(r1)
    output:
    tuple val(sample), val(layout), val(condition), val(replicate), path("${sample}.SE.clean.fastq.gz")
    path "${sample}.fastp.html", optional: true
    path "${sample}.fastp.json", optional: true
    script:
    def r1Files = (r1 instanceof List) ? r1 : [r1]
    def r1Arg = r1Files.join(' ')
    def prepR1 = (r1Files.size() == 1) ? "ln -sf ${r1Files[0]} ${sample}.SE.merged.fastq.gz" : "cat ${r1Arg} > ${sample}.SE.merged.fastq.gz"
    """
    set -euo pipefail
    ${prepR1}
    timeout ${params.fastp_timeout} fastp --thread ${task.cpus} \
          -i ${sample}.SE.merged.fastq.gz -o ${sample}.SE.clean.fastq.gz \
          -h ${sample}.fastp.html -j ${sample}.fastp.json
    """
}

process ALIGN_FILTER_PE {
    tag "$sample"
    publishDir "${params.outdir}/02_align", mode: 'copy', pattern: "${sample}.*"
    input:
    path ref_ok
    tuple val(sample), val(layout), val(condition), val(replicate), path(r1), path(r2)
    output:
    tuple val(sample), val(layout), val(condition), val(replicate), path("${sample}.clean.bam"), path("${sample}.clean.bam.bai"), path("${sample}.qc.tsv")
    path "${sample}.bowtie2.log", optional: true
    path "${sample}.markdup.metrics.txt", optional: true
    path "${sample}.library_complexity.tsv", optional: true
    tuple val(sample), val(layout), val(condition), val(replicate), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai")
    script:
    def removeDup = params.run_markdup ? 'true' : 'false'
    def removeMito = params.remove_mito ? 'true' : 'false'
    def removeBlacklist = params.remove_blacklist ? 'true' : 'false'
    def r1Arg = r1.join(',')
    def r2Arg = r2.join(',')
    """
    set -euo pipefail
    bowtie2 -p ${task.cpus} --very-sensitive --no-mixed --no-discordant -X 2000 \
      -x ${params.bowtie2_index} \
      --rg-id ${sample} --rg "SM:${sample}" --rg "LB:lib1" --rg "PL:ILLUMINA" \
      -1 ${r1Arg} -2 ${r2Arg} 2> ${sample}.bowtie2.log | \
      samtools view -@ ${task.cpus} -bS - > ${sample}.raw.bam
    samtools sort -@ ${task.cpus} -o ${sample}.sorted.bam ${sample}.raw.bam
    samtools index ${sample}.sorted.bam

    if [[ '${removeDup}' == 'true' ]]; then
      picard MarkDuplicates I=${sample}.sorted.bam O=${sample}.markdup.bam M=${sample}.markdup.metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true CREATE_INDEX=true
      mv ${sample}.markdup.bam ${sample}.work.bam
    else
      ln -s ${sample}.sorted.bam ${sample}.work.bam
    fi

    samtools view -@ ${task.cpus} -b -q ${params.mapq} -f 2 -F 1804 ${sample}.work.bam > ${sample}.mapq.bam
    samtools index -@ ${task.cpus} ${sample}.mapq.bam

    samtools view -q ${params.mapq} -f 66 -F 1804 ${sample}.sorted.bam | \
      awk 'BEGIN{OFS="\\t"} {
        if (\$9 > 0) {
          key=\$3":"\$4":"\$7":"\$8":"\$9
          c[key]++
          total++
        }
      } END {
        distinct=0; one=0; two=0
        for (k in c) {
          distinct++
          if (c[k] == 1) one++
          if (c[k] == 2) two++
        }
        nrf=(total>0 ? distinct/total : "NA")
        pbc1=(distinct>0 ? one/distinct : "NA")
        pbc2=(two>0 ? one/two : "NA")
        print "sample","layout","total_fragments_for_complexity","distinct_fragments","one_read_fragments","two_read_fragments","NRF","PBC1","PBC2"
        print "'${sample}'","PE",total+0,distinct+0,one+0,two+0,nrf,pbc1,pbc2
      }' > ${sample}.library_complexity.tsv

    if [[ '${removeMito}' == 'true' ]]; then
      mapfile -t keep_chroms < <(samtools idxstats ${sample}.mapq.bam | awk -v mt="${params.mito_chr}" '\$1 != mt && \$1 != "*" && \$1 != "" {print \$1}')
      printf '%s\n' "\${keep_chroms[@]}" > keep_chroms.txt
      if [[ \${#keep_chroms[@]} -gt 0 ]]; then
        samtools view -@ ${task.cpus} -b ${sample}.mapq.bam "\${keep_chroms[@]}" > ${sample}.nomito.bam
      else
        cp ${sample}.mapq.bam ${sample}.nomito.bam
      fi
    else
      cp ${sample}.mapq.bam ${sample}.nomito.bam
    fi

    if [[ '${removeBlacklist}' == 'true' ]]; then
      bedtools intersect -ubam -abam ${sample}.nomito.bam -b ${params.blacklist} | \
        samtools view - | cut -f1 | sort -u > ${sample}.blacklist.read_names.txt
      if [[ -s ${sample}.blacklist.read_names.txt ]]; then
        samtools view -h ${sample}.nomito.bam | \
          awk -v bad="${sample}.blacklist.read_names.txt" 'BEGIN{while((getline line < bad)>0) drop[line]=1} /^@/{print; next} !(\$1 in drop)' | \
          samtools view -@ ${task.cpus} -b - > ${sample}.clean.unsorted.bam
      else
        cp ${sample}.nomito.bam ${sample}.clean.unsorted.bam
      fi
    else
      cp ${sample}.nomito.bam ${sample}.clean.unsorted.bam
    fi

    samtools sort -@ ${task.cpus} -o ${sample}.clean.bam ${sample}.clean.unsorted.bam
    samtools index ${sample}.clean.bam

    total_raw=\$(samtools view -c ${sample}.raw.bam)
    total_clean=\$(samtools view -c ${sample}.clean.bam)
    mt_reads=\$(samtools idxstats ${sample}.mapq.bam | awk -v mito="${params.mito_chr}" '\$1==mito{s+=\$3} END{print s+0}')
    raw_records=\${total_raw}
    mapped_primary_records=\$(samtools view -c -F 260 ${sample}.raw.bam)
    mapq_records=\$(samtools view -c ${sample}.mapq.bam)
    mapq_fragments=\$(samtools view -c -f 64 ${sample}.mapq.bam)
    if samtools idxstats ${sample}.mapq.bam | cut -f1 | grep -Fxq "${params.mito_chr}"; then
      mitochondrial_fragments=\$(samtools view ${sample}.mapq.bam "${params.mito_chr}" | cut -f1 | sort -u | wc -l | awk '{print \$1}')
    else
      mitochondrial_fragments=0
    fi
    nuclear_fragments_before_blacklist=\$(samtools view -c -f 64 ${sample}.nomito.bam)
    final_nuclear_fragments=\$(samtools view -c -f 64 ${sample}.clean.bam)
    read nfr_fragments mono_fragments di_fragments tri_plus_fragments < <(
      samtools view -f 64 ${sample}.clean.bam | awk '{
        t=\$9; if (t<0) t=-t
        if (t > 0 && t < 100) nfr++
        else if (t >= 100 && t < 180) mono_free++
        else if (t >= 180 && t < 247) mono++
        else if (t >= 247 && t < 394) di++
        else if (t >= 394) tri++
      } END {print nfr+0, mono+0, di+0, tri+0}'
    )
    if [[ -f ${sample}.markdup.metrics.txt ]]; then
      duplicate_records=\$(awk 'BEGIN{v="NA"} /^READ_PAIR_DUPLICATES/ {getline; split(\$0,a,"\\t"); v=a[7]} END{print v}' ${sample}.markdup.metrics.txt)
    else
      duplicate_records=NA
    fi
    echo -e "sample\tlayout\ttotal_raw\ttotal_clean\tmt_reads\traw_records\tmapped_primary_records\tmapq_records\tmapq_fragments\tduplicate_records\tmitochondrial_records\tmitochondrial_fragments\tnuclear_fragments_before_blacklist\tfinal_nuclear_fragments\tnfr_fragments\tmono_fragments\tdi_fragments\ttri_plus_fragments" > ${sample}.qc.tsv
    echo -e "${sample}\tPE\t\${total_raw}\t\${total_clean}\t\${mt_reads}\t\${raw_records}\t\${mapped_primary_records}\t\${mapq_records}\t\${mapq_fragments}\t\${duplicate_records}\t\${mt_reads}\t\${mitochondrial_fragments}\t\${nuclear_fragments_before_blacklist}\t\${final_nuclear_fragments}\t\${nfr_fragments}\t\${mono_fragments}\t\${di_fragments}\t\${tri_plus_fragments}" >> ${sample}.qc.tsv
    """
}

process ALIGN_FILTER_SE {
    tag "$sample"
    publishDir "${params.outdir}/02_align", mode: 'copy', pattern: "${sample}.*"
    input:
    path ref_ok
    tuple val(sample), val(layout), val(condition), val(replicate), path(r1)
    output:
    tuple val(sample), val(layout), val(condition), val(replicate), path("${sample}.clean.bam"), path("${sample}.clean.bam.bai"), path("${sample}.qc.tsv")
    path "${sample}.bowtie2.log", optional: true
    path "${sample}.markdup.metrics.txt", optional: true
    path "${sample}.library_complexity.tsv", optional: true
    tuple val(sample), val(layout), val(condition), val(replicate), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai")
    script:
    def removeDup = params.run_markdup ? 'true' : 'false'
    def removeMito = params.remove_mito ? 'true' : 'false'
    def removeBlacklist = params.remove_blacklist ? 'true' : 'false'
    def r1Arg = r1.join(',')
    """
    set -euo pipefail
    bowtie2 -p ${task.cpus} --very-sensitive \
      -x ${params.bowtie2_index} \
      --rg-id ${sample} --rg "SM:${sample}" --rg "LB:lib1" --rg "PL:ILLUMINA" \
      -U ${r1Arg} 2> ${sample}.bowtie2.log | \
      samtools view -@ ${task.cpus} -bS - > ${sample}.raw.bam
    samtools sort -@ ${task.cpus} -o ${sample}.sorted.bam ${sample}.raw.bam
    samtools index ${sample}.sorted.bam

    if [[ '${removeDup}' == 'true' ]]; then
      picard MarkDuplicates I=${sample}.sorted.bam O=${sample}.markdup.bam M=${sample}.markdup.metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true CREATE_INDEX=true
      mv ${sample}.markdup.bam ${sample}.work.bam
    else
      ln -s ${sample}.sorted.bam ${sample}.work.bam
    fi

    samtools view -@ ${task.cpus} -b -q ${params.mapq} -F 1804 ${sample}.work.bam > ${sample}.mapq.bam
    samtools index -@ ${task.cpus} ${sample}.mapq.bam

    {
      samtools view -q ${params.mapq} -F 1820 ${sample}.sorted.bam | awk 'BEGIN{OFS="\\t"} {print \$3":"\$4":+"}'
      samtools view -q ${params.mapq} -f 16 -F 1804 ${sample}.sorted.bam | awk 'BEGIN{OFS="\\t"} {print \$3":"\$4":-"}'
    } | awk 'BEGIN{OFS="\\t"} {
        c[\$1]++
        total++
      } END {
        distinct=0; one=0; two=0
        for (k in c) {
          distinct++
          if (c[k] == 1) one++
          if (c[k] == 2) two++
        }
        nrf=(total>0 ? distinct/total : "NA")
        pbc1=(distinct>0 ? one/distinct : "NA")
        pbc2=(two>0 ? one/two : "NA")
        print "sample","layout","total_fragments_for_complexity","distinct_fragments","one_read_fragments","two_read_fragments","NRF","PBC1","PBC2"
        print "'${sample}'","SE",total+0,distinct+0,one+0,two+0,nrf,pbc1,pbc2
      }' > ${sample}.library_complexity.tsv

    if [[ '${removeMito}' == 'true' ]]; then
      mapfile -t keep_chroms < <(samtools idxstats ${sample}.mapq.bam | awk -v mt="${params.mito_chr}" '\$1 != mt && \$1 != "*" && \$1 != "" {print \$1}')
      printf '%s\n' "\${keep_chroms[@]}" > keep_chroms.txt
      if [[ \${#keep_chroms[@]} -gt 0 ]]; then
        samtools view -@ ${task.cpus} -b ${sample}.mapq.bam "\${keep_chroms[@]}" > ${sample}.nomito.bam
      else
        cp ${sample}.mapq.bam ${sample}.nomito.bam
      fi
    else
      cp ${sample}.mapq.bam ${sample}.nomito.bam
    fi

    if [[ '${removeBlacklist}' == 'true' ]]; then
      bedtools intersect -v -abam ${sample}.nomito.bam -b ${params.blacklist} > ${sample}.clean.unsorted.bam
    else
      cp ${sample}.nomito.bam ${sample}.clean.unsorted.bam
    fi

    samtools sort -@ ${task.cpus} -o ${sample}.clean.bam ${sample}.clean.unsorted.bam
    samtools index ${sample}.clean.bam

    total_raw=\$(samtools view -c ${sample}.raw.bam)
    total_clean=\$(samtools view -c ${sample}.clean.bam)
    mt_reads=\$(samtools idxstats ${sample}.mapq.bam | awk -v mito="${params.mito_chr}" '\$1==mito{s+=\$3} END{print s+0}')
    raw_records=\${total_raw}
    mapped_primary_records=\$(samtools view -c -F 260 ${sample}.raw.bam)
    mapq_records=\$(samtools view -c ${sample}.mapq.bam)
    if samtools idxstats ${sample}.mapq.bam | cut -f1 | grep -Fxq "${params.mito_chr}"; then
      mitochondrial_records=\$(samtools view -c ${sample}.mapq.bam "${params.mito_chr}")
    else
      mitochondrial_records=0
    fi
    nuclear_records_before_blacklist=\$(samtools view -c ${sample}.nomito.bam)
    final_nuclear_records=\$(samtools view -c ${sample}.clean.bam)
    if [[ -f ${sample}.markdup.metrics.txt ]]; then
      duplicate_records=\$(awk 'BEGIN{v="NA"} /^UNPAIRED_READ_DUPLICATES/ {getline; split(\$0,a,"\\t"); v=a[7]} END{print v}' ${sample}.markdup.metrics.txt)
    else
      duplicate_records=NA
    fi
    echo -e "sample\tlayout\ttotal_raw\ttotal_clean\tmt_reads\traw_records\tmapped_primary_records\tmapq_records\tmapq_fragments\tduplicate_records\tmitochondrial_records\tmitochondrial_fragments\tnuclear_fragments_before_blacklist\tfinal_nuclear_fragments\tnuclear_records_before_blacklist\tfinal_nuclear_records\tnfr_fragments\tmono_fragments\tdi_fragments\ttri_plus_fragments" > ${sample}.qc.tsv
    echo -e "${sample}\tSE\t\${total_raw}\t\${total_clean}\t\${mt_reads}\t\${raw_records}\t\${mapped_primary_records}\t\${mapq_records}\tNA\t\${duplicate_records}\t\${mitochondrial_records}\tNA\tNA\tNA\t\${nuclear_records_before_blacklist}\t\${final_nuclear_records}\tNA\tNA\tNA\tNA" >> ${sample}.qc.tsv
    """
}

process ATAC_QC {
    tag "$sample"
    publishDir "${params.outdir}/03_qc", mode: 'copy'
    input: tuple val(sample), val(layout), val(condition), val(replicate), path(bam), path(bai), path(qcts)
    output: tuple val(sample), path("${sample}.flagstat.txt"), path("${sample}.idxstats.txt"), path("${sample}.insert_size_metrics.txt"), path("${sample}.insert_size_hist.pdf"), path(qcts)
    script:
    def doInsert = (layout == 'PE') ? 'true' : 'false'
    """
    set -euo pipefail
    samtools flagstat ${bam} > ${sample}.flagstat.txt
    samtools idxstats ${bam} > ${sample}.idxstats.txt
    if [[ '${doInsert}' == 'true' ]]; then
      picard CollectInsertSizeMetrics I=${bam} O=${sample}.insert_size_metrics.txt H=${sample}.insert_size_hist.pdf M=0.5
    else
      echo -e "SAMPLE\tNOTE" > ${sample}.insert_size_metrics.txt
      echo -e "${sample}\tSE_library_insert_size_not_applicable" >> ${sample}.insert_size_metrics.txt
      : > ${sample}.insert_size_hist.pdf
    fi
    """
}

process BAMCOVERAGE {
    tag "$sample"
    publishDir "${params.outdir}/04_bw", mode: 'copy'
    input: tuple val(sample), val(layout), val(condition), val(replicate), path(bam), path(bai), path(qcts)
    output: path("${sample}.bw")
    script:
    def normalization = (params.track_normalization ?: 'RPGC').toString().toUpperCase()
    def useRpgc = normalization == 'RPGC' && params.effective_genome_size
    def normArgs = useRpgc ? "--normalizeUsing RPGC --effectiveGenomeSize ${params.effective_genome_size}" : "--normalizeUsing CPM"
    """
    set -euo pipefail
    bamCoverage --bam ${bam} -o ${sample}.bw --binSize 10 ${normArgs} -p ${task.cpus}
    """
}


process ATAC_TE_RELAXED_TRACKS {
    tag "$sample"
    publishDir "${params.outdir}", mode: 'copy'
    input:
    tuple val(sample), val(layout), val(condition), val(replicate), path(bam), path(bai)
    output:
    path "04_bw_te/${sample}.te.bw"
    path "02_align_te/${sample}.te.bam"
    path "02_align_te/${sample}.te.bam.bai"
    path "02_align_te/${sample}.te.clean_counts.tsv"
    path "02_align_te/${sample}.te.flagstat.txt"
    when:
    params.run_te_relaxed_tracks
    script:
    def pairedFilter = (layout == 'PE' && params.te_proper_pair_only) ? '-f 2' : ''
    def removeMito = params.te_remove_mito ? 'true' : 'false'
    def removeBlacklist = params.te_remove_blacklist ? 'true' : 'false'
    def runMarkdup = params.te_run_markdup ? 'true' : 'false'
    def normalization = (params.te_track_normalization ?: 'CPM').toString().toUpperCase()
    def rpgcArgs = normalization == 'RPGC' ? "--effectiveGenomeSize ${params.effective_genome_size}" : ''
    """
    set -euo pipefail
    mkdir -p 02_align_te 04_bw_te

    work_bam=${bam}
    if [[ '${runMarkdup}' == 'true' ]]; then
      picard MarkDuplicates I=${bam} O=${sample}.te.markdup.bam M=${sample}.te.markdup.metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true CREATE_INDEX=false
      work_bam=${sample}.te.markdup.bam
    fi

    samtools view -@ ${task.cpus} -b -q ${params.te_mapq} ${pairedFilter} -F ${params.te_exclude_flags} "\${work_bam}" > ${sample}.te.mapq.bam
    samtools index -@ ${task.cpus} ${sample}.te.mapq.bam

    if [[ '${removeMito}' == 'true' ]]; then
      mapfile -t keep_chroms < <(samtools idxstats ${sample}.te.mapq.bam | awk -v mt="${params.mito_chr}" '\$1 != mt && \$1 != "*" && \$1 != "" {print \$1}')
      if [[ \${#keep_chroms[@]} -gt 0 ]]; then
        samtools view -@ ${task.cpus} -b ${sample}.te.mapq.bam "\${keep_chroms[@]}" > ${sample}.te.nomito.bam
      else
        cp ${sample}.te.mapq.bam ${sample}.te.nomito.bam
      fi
    else
      cp ${sample}.te.mapq.bam ${sample}.te.nomito.bam
    fi

    if [[ '${removeBlacklist}' == 'true' ]]; then
      if [[ '${layout}' == 'PE' ]]; then
        bedtools intersect -ubam -abam ${sample}.te.nomito.bam -b ${params.blacklist} | \
          samtools view - | cut -f1 | sort -u > ${sample}.te.blacklist.read_names.txt
        if [[ -s ${sample}.te.blacklist.read_names.txt ]]; then
          samtools view -h ${sample}.te.nomito.bam | \
            awk -v bad="${sample}.te.blacklist.read_names.txt" 'BEGIN{while((getline line < bad)>0) drop[line]=1} /^@/{print; next} !(\$1 in drop)' | \
            samtools view -@ ${task.cpus} -b - > ${sample}.te.clean.unsorted.bam
        else
          cp ${sample}.te.nomito.bam ${sample}.te.clean.unsorted.bam
        fi
      else
        bedtools intersect -v -abam ${sample}.te.nomito.bam -b ${params.blacklist} > ${sample}.te.clean.unsorted.bam
      fi
    else
      cp ${sample}.te.nomito.bam ${sample}.te.clean.unsorted.bam
    fi

    samtools sort -@ ${task.cpus} -o 02_align_te/${sample}.te.bam ${sample}.te.clean.unsorted.bam
    samtools index -@ ${task.cpus} 02_align_te/${sample}.te.bam
    samtools flagstat 02_align_te/${sample}.te.bam > 02_align_te/${sample}.te.flagstat.txt

    bamCoverage -b 02_align_te/${sample}.te.bam -o 04_bw_te/${sample}.te.bw \
      --binSize ${params.te_bw_binsize} -p ${task.cpus} \
      --normalizeUsing ${normalization} ${rpgcArgs}

    {
      echo -e "step\talignments"
      echo -e "input_sorted\t\$(samtools view -c ${bam})"
      echo -e "post_markdup_step\t\$(samtools view -c \${work_bam})"
      echo -e "post_mapq_flag_filter\t\$(samtools view -c ${sample}.te.mapq.bam)"
      echo -e "post_mito_step\t\$(samtools view -c ${sample}.te.nomito.bam)"
      echo -e "final\t\$(samtools view -c 02_align_te/${sample}.te.bam)"
    } > 02_align_te/${sample}.te.clean_counts.tsv
    """
}

process CALL_PEAKS {
    tag "$sample"
    publishDir "${params.outdir}/05_peaks", mode: 'copy', pattern: "${sample}*"
    input: tuple val(sample), val(layout), val(condition), val(replicate), path(bam), path(bai), path(qcts)
    output: tuple val(sample), val(condition), val(replicate), path("${sample}.peaks.bed"), path("${sample}.summits.bed"), path("${sample}.frip.tsv")
    script:
    def doBroad = params.call_broad ? '--broad' : ''
    def macs_f = (layout == 'PE') ? 'BAMPE' : 'BAM'
    def macs_extra = (layout == 'SE') ? "--shift ${params.se_shift} --extsize ${params.se_extsize}" : ''
    """
    set -euo pipefail
    macs3 callpeak -t ${bam} -f ${macs_f} -g ${params.genome_size} -n ${sample} \
        --nomodel --keep-dup all --call-summits -q ${params.macs_qvalue} --outdir . ${doBroad} ${macs_extra}

    if [[ -f ${sample}_peaks.narrowPeak ]]; then
      cut -f1-3 ${sample}_peaks.narrowPeak > ${sample}.peaks.bed
    else
      cut -f1-3 ${sample}_peaks.broadPeak > ${sample}.peaks.bed
    fi

    if [[ -f ${sample}_summits.bed ]]; then
      cp ${sample}_summits.bed ${sample}.summits.bed
    else
      awk 'BEGIN{OFS="\\t"} {c=int((\$2+\$3)/2); s=c; e=c+1; if(s<0)s=0; print \$1,s,e}' ${sample}.peaks.bed > ${sample}.summits.bed
    fi

    if [[ '${layout}' == 'PE' ]]; then
      total=\$(samtools view -c -f 64 ${bam})
      inpeak=\$(bedtools intersect -u -abam ${bam} -b ${sample}.peaks.bed | samtools view - | cut -f1 | sort -u | wc -l | awk '{print \$1}')
      unit=fragments
    else
      total=\$(samtools view -c ${bam})
      inpeak=\$(bedtools intersect -u -abam ${bam} -b ${sample}.peaks.bed | samtools view -c)
      unit=reads
    fi
    frip=\$(awk -v tp=\$total -v ip=\$inpeak 'BEGIN {print (tp==0 ? 0 : ip/tp)}')
    echo -e "sample\tfrip_unit\ttotal_reads\tin_peak_reads\ttotal_fragments\tin_peak_fragments\tfrip" > ${sample}.frip.tsv
    if [[ "\$unit" == "fragments" ]]; then
      echo -e "${sample}\t\${unit}\t\${total}\t\${inpeak}\t\${total}\t\${inpeak}\t\${frip}" >> ${sample}.frip.tsv
    else
      echo -e "${sample}\t\${unit}\t\${total}\t\${inpeak}\tNA\tNA\t\${frip}" >> ${sample}.frip.tsv
    fi
    """
}

process BUILD_CONSENSUS {
    publishDir "${params.outdir}/06_consensus_peaks", mode: 'copy'
    input: path summit_files
    output: path 'consensus_peaks.bed'
    script:
    """
    set -euo pipefail
    cat *.summits.bed | cut -f1-3 | sort -k1,1 -k2,2n | \
      bedtools slop -i - -g ${params.chrom_sizes} -b ${params.consensus_half_width} | \
      bedtools sort -i - | bedtools merge -i - > consensus_peaks.bed
    """
}

process COUNT_CONSENSUS {
    publishDir "${params.outdir}/07_counts", mode: 'copy'
    input:
    path consensus_bed
    path bam_files
    output:
    path 'consensus_peak_counts.txt'
    path 'sample_metadata.csv'
    script:
    def sampleSheetPath = params.samplesheet.toString()
    def bamList = bam_files.join(' ')
    """
    set -euo pipefail
    awk 'BEGIN{OFS="\\t"; print "GeneID","Chr","Start","End","Strand"} {print "peak_"NR,\$1,\$2+1,\$3,"+"}' ${consensus_bed} > consensus_peaks.saf

    python3 - <<'PY' ${sampleSheetPath} > counts_config.txt
import csv, sys
pe, se = 0, 0
with open(sys.argv[1]) as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row['layout'].upper() == 'PE': pe += 1
        else: se += 1
print(f"{pe},{se}")
PY
    IFS=',' read -r pe_n se_n < counts_config.txt
    if [[ "\$pe_n" -gt 0 && "\$se_n" -gt 0 ]]; then
      echo "ERROR: PE and SE samples cannot be mixed in one consensus peak count matrix. Run them separately or split the samplesheet." >&2
      exit 1
    fi

    if [[ "\$pe_n" -gt 0 && "\$se_n" -eq 0 ]]; then
      featureCounts -a consensus_peaks.saf -F SAF -o consensus_peak_counts.txt -T ${task.cpus} -p --countReadPairs -B -C ${bamList}
    else
      featureCounts -a consensus_peaks.saf -F SAF -o consensus_peak_counts.txt -T ${task.cpus} ${bamList}
    fi

    python3 - <<'PY' ${sampleSheetPath}
import csv, sys
with open(sys.argv[1]) as rf, open('sample_metadata.csv', 'w') as wf:
    reader = csv.DictReader(rf)
    writer = csv.DictWriter(wf, fieldnames=['sample','condition','replicate'])
    writer.writeheader()
    for row in reader:
        writer.writerow({'sample': row['sample'], 'condition': row.get('condition','NA'), 'replicate': row.get('replicate','NA')})
PY
    """
}

process MERGE_QC_SUMMARY {
    publishDir "${params.outdir}/03_qc", mode: 'copy'
    input:
    path qc_tsvs
    path frip_tsvs
    path complexity_tsvs
    output:
    path 'atac_qc_summary.tsv'
    script:
    """
    set -euo pipefail
    echo "[MERGE_QC_SUMMARY] start" >&2
    python3 - <<'PY'
import pandas as pd
import glob
qc_files = sorted(glob.glob('*.qc.tsv'))
frip_files = sorted(glob.glob('*.frip.tsv'))
complexity_files = sorted(glob.glob('*.library_complexity.tsv'))
qc = pd.concat([pd.read_csv(x, sep='\t') for x in qc_files], ignore_index=True) if qc_files else pd.DataFrame()
frip = pd.concat([pd.read_csv(x, sep='\t') for x in frip_files], ignore_index=True) if frip_files else pd.DataFrame()
complexity = pd.concat([pd.read_csv(x, sep='\t') for x in complexity_files], ignore_index=True) if complexity_files else pd.DataFrame()
if not qc.empty and not frip.empty:
    out = qc.merge(frip, on='sample', how='outer')
elif not qc.empty:
    out = qc
else:
    out = frip
if not complexity.empty:
    out = out.merge(complexity, on=['sample', 'layout'], how='outer') if 'layout' in out.columns else out.merge(complexity, on='sample', how='outer')
out.to_csv('atac_qc_summary.tsv', sep='\t', index=False)
PY
    """
}

process QC_PLOTS {
    publishDir "${params.outdir}/03_qc", mode: 'copy'
    input:
    path qc_summary
    output:
    path 'FRiP_barplot.pdf'
    path 'FRiP_barplot.png', optional: true
    path 'Mito_fraction_barplot.pdf', optional: true
    path 'Mito_fraction_barplot.png', optional: true
    path 'Library_complexity_barplot.pdf', optional: true
    path 'Library_complexity_barplot.png', optional: true
    path 'Fragment_class_fraction_barplot.pdf', optional: true
    path 'Fragment_class_fraction_barplot.png', optional: true
    path 'qc_pass_fail.tsv'
    path 'atac_qc_summary_with_metrics.tsv'
    script:
    """
    set -euo pipefail
    Rscript ${projectDir}/scripts/plot_atac_qc.R --qc-summary ${qc_summary} --outdir .
    """
}

process TSS_ENRICH {
    publishDir "${params.outdir}/03_qc/tss_enrichment", mode: 'copy'
    input:
    path bw_files
    output:
    path 'TSS_enrichment_matrix.gz'
    path 'TSS_enrichment_profile.pdf'
    path 'TSS_enrichment_profile.png'
    path 'TSS_enrichment_heatmap.pdf'
    path 'TSS_enrichment_heatmap.png'
    when:
    params.run_tss_enrich && params.tss_bed
    script:
    def bwArg = bw_files.collect { bw -> bw.getName() }.join(' ')
    """
    set -euo pipefail
    echo "[TSS_ENRICH] start" >&2
    awk 'BEGIN{OFS="\t"} !/^#/ && NF>=3 { s=\$2; e=\$3; if (s<0) s=0; if (e<=s) e=s+1; print \$1,s,e,(NF>=4?\$4:"."),(NF>=5?\$5:"0"),(NF>=6?\$6:".") }' ${params.tss_bed} | sort -k1,1 -k2,2n > tss.valid.bed
    [[ -s tss.valid.bed ]] || { echo "[TSS_ENRICH] no valid regions after sanitization" >&2; exit 1; }
    computeMatrix reference-point \
      --referencePoint TSS \
      -R tss.valid.bed \
      -S ${bwArg} \
      -b 3000 -a 3000 \
      --skipZeros \
      -p ${task.cpus} \
      -o TSS_enrichment_matrix.gz
    plotProfile -m TSS_enrichment_matrix.gz -out TSS_enrichment_profile.pdf --perGroup
    plotProfile -m TSS_enrichment_matrix.gz -out TSS_enrichment_profile.png --perGroup
    plotHeatmap -m TSS_enrichment_matrix.gz -out TSS_enrichment_heatmap.pdf --sortRegions descend --whatToShow 'plot, heatmap and colorbar'
    plotHeatmap -m TSS_enrichment_matrix.gz -out TSS_enrichment_heatmap.png --sortRegions descend --whatToShow 'plot, heatmap and colorbar'
    """
}


process GENE_BODY_PROFILE {
    publishDir "${params.outdir}/03_qc/gene_body_profile", mode: 'copy'
    input:
    path bw_files
    output:
    path 'GeneBody_accessibility_matrix.gz'
    path 'GeneBody_accessibility_profile.pdf'
    path 'GeneBody_accessibility_profile.png'
    path 'GeneBody_accessibility_heatmap.pdf'
    path 'GeneBody_accessibility_heatmap.png'
    when:
    params.run_gene_body_profile && params.gene_body_bed
    script:
    def bwArg = bw_files.collect { bw -> bw.getName() }.join(' ')
    """
    set -euo pipefail
    echo "[GENE_BODY_PROFILE] start" >&2
    awk 'BEGIN{OFS="\t"} !/^#/ && NF>=3 { s=\$2; e=\$3; if (s<0) s=0; if (e<=s) next; print \$1,s,e,(NF>=4?\$4:"."),(NF>=5?\$5:"0"),(NF>=6?\$6:".") }' ${params.gene_body_bed} | sort -k1,1 -k2,2n > gene_body.valid.bed
    [[ -s gene_body.valid.bed ]] || { echo "[GENE_BODY_PROFILE] no valid regions after sanitization" >&2; exit 1; }
    computeMatrix scale-regions       -R gene_body.valid.bed       -S ${bwArg}       -b 3000 -a 3000       --regionBodyLength 5000       --skipZeros       -p ${task.cpus}       -o GeneBody_accessibility_matrix.gz
    plotProfile -m GeneBody_accessibility_matrix.gz -out GeneBody_accessibility_profile.pdf --perGroup
    plotProfile -m GeneBody_accessibility_matrix.gz -out GeneBody_accessibility_profile.png --perGroup
    plotHeatmap -m GeneBody_accessibility_matrix.gz -out GeneBody_accessibility_heatmap.pdf --sortRegions descend --whatToShow 'plot, heatmap and colorbar'
    plotHeatmap -m GeneBody_accessibility_matrix.gz -out GeneBody_accessibility_heatmap.png --sortRegions descend --whatToShow 'plot, heatmap and colorbar'
    """
}

process TE_HEATMAP {
    publishDir "${params.outdir}/03_qc/te_heatmap", mode: 'copy'
    input:
    path bw_files
    path consensus_bed
    output:
    path 'TE_accessibility_matrix.gz'
    path 'TE_accessibility_profile.pdf'
    path 'TE_accessibility_profile.png'
    path 'TE_accessibility_heatmap.pdf'
    path 'TE_accessibility_heatmap.png'
    path 'te.locus.bed', optional: true
    path 'te.filtered.bed', optional: true
    path 'te.sampled.bed', optional: true
    when:
    params.run_te_heatmap && params.te_bed
    script:
    """
    set -euo pipefail
    bash ${projectDir}/scripts/run_te_heatmap.sh \
      --bw-glob "*.bw" \
      --species ${params.species} \
      --te-bed ${params.te_bed} \
      --peak-bed ${consensus_bed} \
      --outdir . \
      --cores ${task.cpus} \
      --purpose global
    """
}

process RUN_DOWNSTREAM_PEAK {
    publishDir "${params.outdir}/08_downstream/peak_level", mode: 'copy'
    input:
    path count_file
    path sample_meta
    path peak_bed
    output:
    path 'downstream.ok'
    path 'contrast_beds', optional: true
    path 'qc', optional: true
    path 'figures', optional: true
    path 'contrasts', optional: true
    script:
    def contrastFile = params.contrast_file ? file(params.contrast_file) : null
    def contrastArg = (contrastFile && contrastFile.exists()) ? "--contrast-file ${contrastFile}" : ''
    def teArg = params.te_bed ? "--te-bed ${params.te_bed}" : ''
    def gtfArg = params.gtf_genes ? "--gtf ${params.gtf_genes}" : ''
    """
    set -euo pipefail
    Rscript ${projectDir}/atacseq-downstream/run_downstream_atac.R \
      --count-file ${count_file} \
      --sample-meta ${sample_meta} \
      --peak-bed ${peak_bed} \
      --outdir . \
      --species ${params.species} \
      ${contrastArg} \
      ${teArg} \
      ${gtfArg}
    echo OK > downstream.ok
    """
}

process BUILD_FIXED_BINS {
    publishDir "${params.outdir}/07_counts/bin_level", mode: 'copy'
    output:
    path "fixed_bins_${params.fixedbin_size}.bed"
    when:
    params.run_fixedbin
    script:
    """
    set -euo pipefail
    bedtools makewindows -g ${params.chrom_sizes} -w ${params.fixedbin_size} \
      | awk 'BEGIN{OFS="\\t"} \$1 ~ /^chr([0-9]+|X|Y|M)\$/ {print \$0}' > fixed_bins_${params.fixedbin_size}.bed
    [[ -s fixed_bins_${params.fixedbin_size}.bed ]] || { echo "No fixed bins generated from ${params.chrom_sizes}" >&2; exit 1; }
    """
}

process COUNT_FIXED_BINS {
    publishDir "${params.outdir}/07_counts/bin_level", mode: 'copy'
    input:
    path bins_bed
    path bam_files
    output:
    path 'fixedbin_counts.txt'
    path 'sample_metadata.csv'
    path "fixed_bins_${params.fixedbin_size}.bed"
    when:
    params.run_fixedbin
    script:
    def sampleSheetPath = params.samplesheet.toString()
    def bamList = bam_files.join(' ')
    """
    set -euo pipefail
    awk 'BEGIN{OFS="\\t"; print "GeneID","Chr","Start","End","Strand"} {print "bin_"NR,\$1,\$2+1,\$3,"+"}' ${bins_bed} > fixed_bins.saf

    python3 - <<'PY' ${sampleSheetPath} > counts_config.txt
import csv, sys
pe, se = 0, 0
with open(sys.argv[1]) as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row['layout'].upper() == 'PE': pe += 1
        else: se += 1
print(f"{pe},{se}")
PY
    IFS=',' read -r pe_n se_n < counts_config.txt
    if [[ "\$pe_n" -gt 0 && "\$se_n" -gt 0 ]]; then
      echo "ERROR: PE and SE samples cannot be mixed in one fixed-bin count matrix. Run them separately or split the samplesheet." >&2
      exit 1
    fi

    if [[ "\$pe_n" -gt 0 && "\$se_n" -eq 0 ]]; then
      featureCounts -a fixed_bins.saf -F SAF -o fixedbin_counts.txt -T ${task.cpus} -p --countReadPairs -B -C ${bamList}
    else
      featureCounts -a fixed_bins.saf -F SAF -o fixedbin_counts.txt -T ${task.cpus} ${bamList}
    fi

    python3 - <<'PY' ${sampleSheetPath}
import csv, sys
with open(sys.argv[1]) as rf, open('sample_metadata.csv', 'w') as wf:
    reader = csv.DictReader(rf)
    writer = csv.DictWriter(wf, fieldnames=['sample','condition','replicate'])
    writer.writeheader()
    for row in reader:
        writer.writerow({'sample': row['sample'], 'condition': row.get('condition','NA'), 'replicate': row.get('replicate','NA')})
PY
    """
}

process RUN_DOWNSTREAM_BIN {
    publishDir "${params.outdir}/08_downstream/bin_level", mode: 'copy'
    input:
    path count_file
    path sample_meta
    path bins_bed
    output:
    path 'downstream.ok'
    path 'contrast_beds', optional: true
    path 'qc', optional: true
    path 'figures', optional: true
    path 'contrasts', optional: true
    when:
    params.run_fixedbin && params.run_downstream
    script:
    def contrastFile = params.contrast_file ? file(params.contrast_file) : null
    def contrastArg = (contrastFile && contrastFile.exists()) ? "--contrast-file ${contrastFile}" : ''
    def teArg = params.te_bed ? "--te-bed ${params.te_bed}" : ''
    def gtfArg = params.gtf_genes ? "--gtf ${params.gtf_genes}" : ''
    """
    set -euo pipefail
    Rscript ${projectDir}/atacseq-downstream/run_downstream_atac.R \
      --count-file ${count_file} \
      --sample-meta ${sample_meta} \
      --peak-bed ${bins_bed} \
      --outdir . \
      --species ${params.species} \
      ${contrastArg} \
      ${teArg} \
      ${gtfArg}
    echo OK > downstream.ok
    """
}

process MOTIF_ENRICH {
    publishDir "${params.outdir}/09_motif", mode: 'copy'
    input:
    path contrast_beds_dir
    output:
    path 'motif_summary.tsv'
    path 'motif_results'
    when:
    params.run_motif && params.motif_genome
    script:
    """
    set -euo pipefail
    bash ${projectDir}/scripts/run_motif_homer.sh \
      --bed-dir contrast_beds \
      --genome ${params.motif_genome} \
      --outdir motif_results \
      --size 200 \
      --run-annotation true
    cp motif_results/motif_summary.tsv motif_summary.tsv
    """
}

process FOOTPRINTING {
    publishDir "${params.outdir}/10_footprinting", mode: 'copy'
    input:
    path bam_files
    path consensus_bed
    path sample_meta
    path contrast_beds_dir
    output:
    path 'tobias.done'
    when:
    params.run_footprinting && params.motif_meme
    script:
    def contrastFile = params.contrast_file ? file(params.contrast_file) : null
    def contrastArg = (contrastFile && contrastFile.exists()) ? "--contrast-file ${contrastFile}" : ''
    """
    set -euo pipefail
    bash ${projectDir}/scripts/run_tobias.sh \
      --bam-dir . \
      --peaks ${consensus_bed} \
      --sample-meta ${sample_meta} \
      --genome ${params.genome_fasta} \
      --motif-meme ${params.motif_meme} \
      --outdir . \
      --cores ${task.cpus} \
      ${contrastArg}
    echo OK > tobias.done
    """
}


process NUC_PHASING {
    publishDir "${params.outdir}/11_nuc_phasing", mode: "copy"
    input:
    path bam_files
    output:
    path "NRL_summary.csv"
    path "*_nuc_phasing.pdf", optional: true
    path "*_nuc_phasing.png", optional: true
    when:
    params.run_nuc_phasing
    script:
    def bam_list = bam_files.collect { bam -> bam.getName() }.join("\n")
    """
    set -euo pipefail
    if [[ -z "${bam_list}" ]]; then
      echo "sample,nrl,mono_peak,di_peak,tri_peak,n_fragments,note" > NRL_summary.csv
      echo "no_pe_samples,NA,NA,NA,NA,0,no paired-end BAMs were provided" >> NRL_summary.csv
      exit 0
    fi
    printf "%b\n" "${bam_list}" > bam_list.txt
    Rscript ${projectDir}/atacseq-downstream/nuc_phasing.R \
      --bam-list bam_list.txt \
      --outdir . \
      --cores ${task.cpus}
    """
}

workflow {
    if (!params.samplesheet) exit 1, 'Please provide --samplesheet'

    def supportedSpecies = ['hg38', 'mm10', 'mm39']
    if (!supportedSpecies.contains(params.species.toString())) {
        exit 1, "Unsupported --species '${params.species}'. Use one of: ${supportedSpecies.join(', ')}"
    }
    if (!params.genome_fasta || !params.bowtie2_index || !params.chrom_sizes) {
        exit 1, "Reference configuration is incomplete for species '${params.species}': genome_fasta, bowtie2_index and chrom_sizes are required"
    }
    if (params.run_downstream && !params.gtf_genes) {
        exit 1, "run_downstream=true requires gtf_genes in conf/species_refs.config or an explicit --gtf_genes override"
    }
    if (params.run_tss_enrich && !params.tss_bed) {
        exit 1, 'run_tss_enrich=true requires tss_bed'
    }
    if (params.run_gene_body_profile && !params.gene_body_bed) {
        exit 1, 'run_gene_body_profile=true requires gene_body_bed'
    }
    if (params.run_te_heatmap && !params.te_bed) {
        exit 1, 'run_te_heatmap=true requires te_bed'
    }
    if (params.run_motif && !params.motif_genome) {
        exit 1, 'run_motif=true requires motif_genome'
    }
    if (params.run_footprinting && !params.motif_meme) {
        exit 1, 'run_footprinting=true requires motif_meme'
    }
    if (params.run_downstream && !params.run_peak_calling) {
        exit 1, 'run_downstream=true requires run_peak_calling=true in this main workflow; use run_atac_downstream_only for matrix-only downstream.'
    }
    if (params.run_motif && !params.run_downstream) {
        exit 1, 'run_motif=true requires run_downstream=true because motif enrichment consumes contrast_beds from downstream output.'
    }
    if (params.run_footprinting && !params.run_downstream) {
        exit 1, 'run_footprinting=true requires run_downstream=true because footprinting currently consumes contrast outputs.'
    }

    def ch_samples = channel
        .fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map { row ->
            def layout = row.layout.toString().trim().toUpperCase()
            assert layout in ['PE','SE'] : "Invalid layout for sample ${row.sample}: ${row.layout}"
            def r1 = row.r1.toString().split(',').collect{ path -> file(path.trim()) }
            def r2 = row.r2 && row.r2.toString().trim() ? row.r2.toString().split(',').collect{ path -> file(path.trim()) } : []
            if( layout == 'PE' && !r2 ) throw new IllegalArgumentException("Sample ${row.sample} is PE but r2 missing")
            tuple(row.sample.toString(), layout, row.condition ?: 'NA', row.replicate ?: 'NA', r1, r2)
        }
    def ch_pe = ch_samples.filter { row -> row[1] == 'PE' }
    def ch_se = ch_samples.filter { row -> row[1] == 'SE' }.map { row -> tuple(row[0], row[1], row[2], row[3], row[4]) }

    CHECK_REFS()
    ref_ok = CHECK_REFS.out

    if (params.run_fastp) {
        FASTP_PE(ref_ok, ch_pe)
        FASTP_SE(ref_ok, ch_se)
        prep_pe = FASTP_PE.out[0]
        prep_se = FASTP_SE.out[0]
    } else {
        prep_pe = ch_pe
        prep_se = ch_se
    }
    ALIGN_FILTER_PE(ref_ok, prep_pe)
    ALIGN_FILTER_SE(ref_ok, prep_se)
    aligned_all = ALIGN_FILTER_PE.out[0].mix(ALIGN_FILTER_SE.out[0])
    ATAC_QC(aligned_all)

    if (params.run_nuc_phasing) {
        clean_bams_for_nrl = aligned_all.filter { sample, layout, condition, replicate, bam, bai, qcts -> layout == 'PE' }.map { sample, layout, condition, replicate, bam, bai, qcts -> bam }.collect()
        NUC_PHASING(clean_bams_for_nrl)
    }

    need_bw = params.run_bamcoverage || params.run_tss_enrich || params.run_gene_body_profile || (params.run_te_heatmap && !params.run_te_relaxed_tracks)
    if (need_bw) {
        BAMCOVERAGE(aligned_all)
        bw_files = BAMCOVERAGE.out.collect()
        if (params.run_tss_enrich && params.tss_bed) {
            TSS_ENRICH(bw_files)
        }
        if (params.run_gene_body_profile && params.gene_body_bed) {
            GENE_BODY_PROFILE(bw_files)
        }
    }

    if (params.run_te_relaxed_tracks) {
        sorted_bams_for_te = ALIGN_FILTER_PE.out[4].mix(ALIGN_FILTER_SE.out[4])
        ATAC_TE_RELAXED_TRACKS(sorted_bams_for_te)
        te_relaxed_bw_files = ATAC_TE_RELAXED_TRACKS.out[0].collect()
    }

    if (params.run_peak_calling) {
        peaks_all = CALL_PEAKS(aligned_all)

        qc_files = aligned_all.map { sample, layout, condition, replicate, bam, bai, qcts -> qcts }.collect()
        frip_files = peaks_all.map { sample, condition, replicate, peak, summits, frip -> frip }.collect()
        complexity_files = ALIGN_FILTER_PE.out[3].mix(ALIGN_FILTER_SE.out[3]).collect()
        MERGE_QC_SUMMARY(qc_files, frip_files, complexity_files)
        QC_PLOTS(MERGE_QC_SUMMARY.out)

        summit_files = peaks_all.map { sample, condition, replicate, peak, summits, frip -> summits }.collect()
        BUILD_CONSENSUS(summit_files)

        clean_bams = aligned_all.map { sample, layout, condition, replicate, bam, bai, qcts -> bam }.collect()
        COUNT_CONSENSUS(BUILD_CONSENSUS.out, clean_bams)

        if (params.run_fixedbin) {
            BUILD_FIXED_BINS()
            COUNT_FIXED_BINS(BUILD_FIXED_BINS.out, clean_bams)
        }

        if (params.run_downstream) {
            RUN_DOWNSTREAM_PEAK(COUNT_CONSENSUS.out[0], COUNT_CONSENSUS.out[1], BUILD_CONSENSUS.out)

            if (params.run_fixedbin) {
                RUN_DOWNSTREAM_BIN(COUNT_FIXED_BINS.out[0], COUNT_FIXED_BINS.out[1], COUNT_FIXED_BINS.out[2])
            }

            if (params.run_te_heatmap && params.te_bed) {
                if (params.run_te_relaxed_tracks) {
                    TE_HEATMAP(te_relaxed_bw_files, BUILD_CONSENSUS.out)
                } else if (need_bw) {
                    TE_HEATMAP(bw_files, BUILD_CONSENSUS.out)
                }
            }
            if (params.run_motif) {
                MOTIF_ENRICH(RUN_DOWNSTREAM_PEAK.out[1])
            }
            if (params.run_footprinting) {
                clean_bams2 = aligned_all.map { sample, layout, condition, replicate, bam, bai, qcts -> bam }.collect()
                FOOTPRINTING(clean_bams2, BUILD_CONSENSUS.out, COUNT_CONSENSUS.out[1], RUN_DOWNSTREAM_PEAK.out[1])
            }
        }
    }
}
