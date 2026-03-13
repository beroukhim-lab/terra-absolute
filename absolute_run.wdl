version 1.0

task absolute {
  input {
    File capseg_file
    File snp
    File indel
    Float skew = 0.99
    String sample_name

    Int mem_gb = 32
    Int cpu = 2
    Int disk_gb = 150
    String disk_type = "HDD"   # or "SSD"
  }

  command <<<
    set -euo pipefail
    set -x


    echo "=== INPUT FILES (as provided to task) ==="
    ls -lh "~{capseg_file}" "~{snp}" "~{indel}"
    echo ""

    # --- Clean seg: remove NA rows and strip chr from column 1 (keep header)
    grep -v "NA" "~{capseg_file}" > no_nan_segs.tsv
    awk 'BEGIN{FS=OFS="\t"} NR==1{print;next} {gsub(/^chr/, "", $1); print}' \
      no_nan_segs.tsv > reformat_seg.tsv


    # --- Use provided SNP/INDEL, only strip chr (Chromosome is column 2 in reduced MAF)
    awk 'BEGIN{FS=OFS="\t"} NR==1{print;next} {gsub(/^chr/, "", $2); print}' "~{snp}"   > reformat_snv.maf
    awk 'BEGIN{FS=OFS="\t"} NR==1{print;next} {gsub(/^chr/, "", $2); print}' "~{indel}" > reformat_indel.maf

    echo "=== PREVIEW: SEG (reformat_seg.tsv) ==="
    head -n 5 reformat_seg.tsv
    echo "Rows:"; wc -l reformat_seg.tsv
    echo ""

    echo "=== PREVIEW: SNV (reformat_snv.maf) ==="
    head -n 5 reformat_snv.maf
    echo "Rows:"; wc -l reformat_snv.maf
    echo ""

    echo "=== PREVIEW: INDEL (reformat_indel.maf) ==="
    head -n 5 reformat_indel.maf
    echo "Rows:"; wc -l reformat_indel.maf
    echo ""
    # --- HARD VALIDATION: ensure SNV/INDEL not empty (beyond header)
    seg_lines="$(wc -l < reformat_seg.tsv || true)"
    snv_lines="$(wc -l < reformat_snv.maf || true)"
    indel_lines="$(wc -l < reformat_indel.maf || true)"

    echo "seg_lines=${seg_lines}"
    echo "snv_lines=${snv_lines}"
    echo "indel_lines=${indel_lines}"

    test -n "${seg_lines}" || { echo "ERROR: failed to count lines in reformat_seg.tsv"; exit 1; }
    test -n "${snv_lines}" || { echo "ERROR: failed to count lines in reformat_snv.maf"; exit 1; }
    test -n "${indel_lines}" || { echo "ERROR: failed to count lines in reformat_indel.maf"; exit 1; }

    if [ "${seg_lines}" -lt 2 ]; then
      echo "ERROR: Seg file has no data rows after cleaning (only header)."
      exit 1
    fi

    if [ "${snv_lines}" -lt 2 ]; then
      echo "ERROR: SNV file has no data rows (only header). ABSOLUTE will fail."
      exit 1
    fi

    if [ "${indel_lines}" -lt 2 ]; then
      echo "WARNING: INDEL file has no data rows (only header). Continuing anyway."
    fi

    # --- Optional: install RColorBrewer if missing (won't stop run if install fails)
    echo "=== Checking R packages ==="
    Rscript -e 'if (!requireNamespace("RColorBrewer", quietly=TRUE)) { message("RColorBrewer missing; attempting install..."); try(install.packages("RColorBrewer", repos="https://cloud.r-project.org"), silent=TRUE) } else { message("RColorBrewer present") }' || true
    echo ""

    echo "=== RUNNING ABSOLUTE ==="
    Rscript /xchip/tcga/Tools/absolute/releases/v1.5/run/ABSOLUTE_cli_start.R \
      --seg_dat_fn reformat_seg.tsv \
      --maf_fn reformat_snv.maf \
      --indelmaf_fn reformat_indel.maf \
      --sample_name "~{sample_name}" \
      --results_dir . \
      --ssnv_skew "~{skew}" \
      --abs_lib_dir /xchip/tcga/Tools/absolute/releases/v1.5/

    echo ""
    echo "=== RUNNING ABSOLUTE EXTRACT ==="
    Rscript /xchip/tcga/Tools/absolute/releases/v1.5/run/ABSOLUTE_extract_cli_start.R \
      --solution_num 1 \
      --analyst_id force_called \
      --rdata_modes_fn "~{sample_name}.PP-modes.data.RData" \
      --sample_name "~{sample_name}" \
      --results_dir . \
      --abs_lib_dir /xchip/tcga/Tools/absolute/releases/v1.5/

    # Copy reviewed outputs if present (don't hard-fail if absent)
    cp -f reviewed/samples/~{sample_name}.ABSOLUTE.force_called.called.RData . 2>/dev/null || true
    cp -f reviewed/SEG_MAF/~{sample_name}_ABS_MAF.txt . 2>/dev/null || true
    cp -f reviewed/SEG_MAF/~{sample_name}.segtab.txt . 2>/dev/null || true

    echo "=== DONE ==="
    ls -lh
  >>>

  runtime {
    docker: "gcr.io/broad-getzlab-workflows/absolute_wolf:no_indel_filter_v6"
    cpu: cpu
    memory: "~{mem_gb}G"
    disks: "local-disk ~{disk_gb} ~{disk_type}"
  }

  output {
    File absolute_highres_plot = "~{sample_name}.ABSOLUTE_plot.pdf"
    File absolute_rdata = "~{sample_name}.PP-modes.data.RData"
    File? absolute_segdat_file = "~{sample_name}.ABSOLUTE.force_called.called.RData"
    File? absolute_annotated_maf_capture = "~{sample_name}_ABS_MAF.txt"
    File? absolute_seg_file = "~{sample_name}.segtab.txt"
  }
}

workflow absolute_wkflw {
  input {
    File capseg_file
    File snp
    File indel
    Float skew = 0.99
    String sample_name

    # workflow-level knobs users can set in Terra
    Int absolute_mem_gb = 32
    Int absolute_cpu = 2
    Int absolute_disk_gb = 150
    String absolute_disk_type = "HDD"
  }

  call absolute {
    input:
      capseg_file = capseg_file,
      snp = snp,
      indel = indel,
      skew = skew,
      sample_name = sample_name,
      mem_gb = absolute_mem_gb,
      cpu = absolute_cpu,
      disk_gb = absolute_disk_gb,
      disk_type = absolute_disk_type
  }

  output {
    File absolute_highres_plot = absolute.absolute_highres_plot
    File absolute_rdata = absolute.absolute_rdata
    File? absolute_segdat_file = absolute.absolute_segdat_file
    File? absolute_annotated_maf_capture = absolute.absolute_annotated_maf_capture
    File? absolute_seg_file = absolute.absolute_seg_file
  }
}