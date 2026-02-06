version 1.0

task absolute {
  input {
    File seg_file
    File snp
    File indel
    Float skew = 0.99
    String sample_name
  }

  command <<<
    set -euo pipefail

    echo "=== INPUT FILES (as provided to task) ==="
    ls -lh "~{seg_file}" "~{snp}" "~{indel}"
    echo ""

    # --- Clean seg: remove NA rows and strip chr from column 1 (keep header)
    grep -v "NA" "~{seg_file}" > no_nan_segs.tsv
    awk 'BEGIN{FS=OFS="\t"} NR==1{print;next} {gsub(/^chr/, "", $1); print}' \
      no_nan_segs.tsv > reformat_seg.tsv

    # --- Use provided SNP/INDEL, only strip chr (assumes Chromosome is column 5 like standard MAF)
    awk 'BEGIN{FS=OFS="\t"} NR==1{print;next} {gsub(/^chr/, "", $5); print}' "~{snp}"   > reformat_snv.maf
    awk 'BEGIN{FS=OFS="\t"} NR==1{print;next} {gsub(/^chr/, "", $5); print}' "~{indel}" > reformat_indel.maf

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
    snv_lines=$(wc -l < reformat_snv.maf)
    indel_lines=$(wc -l < reformat_indel.maf)
    seg_lines=$(wc -l < reformat_seg.tsv)

    if [ "${seg_lines}" -lt 2 ]; then
      echo "ERROR: Seg file has no data rows after cleaning (only header)."
      exit 1
    fi

    if [ "${snv_lines}" -lt 2 ]; then
      echo "ERROR: SNV file has no data rows (only header). ABSOLUTE will fail."
      exit 1
    fi

    # indels can be empty sometimes; warn but don’t fail
    if [ "${indel_lines}" -lt 2 ]; then
      echo "WARNING: INDEL file has no data rows (only header). Continuing anyway."
    fi

    # --- HARD VALIDATION: check chrom overlap between seg (col1) and snv (col5)
    cut -f1 reformat_seg.tsv | tail -n +2 | sort -u > seg.chroms.txt
    cut -f5 reformat_snv.maf | tail -n +2 | sort -u > snv.chroms.txt

    overlap_count=$(comm -12 seg.chroms.txt snv.chroms.txt | wc -l | tr -d ' ')
    echo "Chrom overlap count (SEG vs SNV): ${overlap_count}"
    if [ "${overlap_count}" -eq 0 ]; then
      echo "ERROR: No chromosome overlap between seg and SNV inputs."
      echo "This often means genome build mismatch (hg19 vs hg38) or wrong column assumption."
      echo "SEG chroms (first 20):"; head -n 20 seg.chroms.txt
      echo "SNV chroms (first 20):"; head -n 20 snv.chroms.txt
      exit 1
    fi
    echo ""

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
    memory: "7G"
    docker: "gcr.io/broad-getzlab-workflows/absolute_wolf:no_indel_filter_v6"
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
    File seg_file
    File snp
    File indel
    Float skew = 0.99
    String sample_name
  }

  call absolute {
    input:
      seg_file = seg_file,
      snp = snp,
      indel = indel,
      skew = skew,
      sample_name = sample_name
  }

  output {
    File absolute_highres_plot = absolute.absolute_highres_plot
    File absolute_rdata = absolute.absolute_rdata
    File? absolute_segdat_file = absolute.absolute_segdat_file
    File? absolute_annotated_maf_capture = absolute.absolute_annotated_maf_capture
    File? absolute_seg_file = absolute.absolute_seg_file
  }
}