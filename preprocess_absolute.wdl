version 1.0

# ------------------------
# Task: Runs your R script
# ------------------------
task maf_to_absolute_inputs {
  input {
    File maf
    String sample_id

    Int mem_gb = 8
    Int cpu = 1
    Int disk_gb = 20
  }

  command <<<
    set -euo pipefail
    # Create the R script at runtime
    cat << 'EOF' > maf_to_absolute_inputs.R
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript maf_to_absolute_inputs.R <input.maf> <sample_id>\n", call. = FALSE)
}

input_maf <- args[1]
sample_id <- args[2]

# Output file names
snp_out   <- paste0(sample_id, ".snp")
indel_out <- paste0(sample_id, ".indel")

# Converting input maf file to snp and indel files needed for ABSOLUTE

model_colnames <- c(
  "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build",
  "Chromosome", "Start_position", "End_position", "Strand",
  "Variant_Classification", "Variant_Type", "Reference_Allele",
  "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "dbSNP_RS",
  "dbSNP_Val_Status", "Tumor_Sample_Barcode",
  "Matched_Norm_Sample_Barcode", "Protein_Change",
  "UniProt_AApos", "t_alt_count", "t_ref_count"
)

# Read MAF
maf <- read.delim(input_maf, comment.char = "#", stringsAsFactors = FALSE, check.names = FALSE)

model_cols <- tolower(model_colnames)
maf_colnames <- tolower(colnames(maf))

retain_columns <- maf_colnames %in% model_cols
maf <- maf[, retain_columns, drop = FALSE]

# Fix potential column name variants
names(maf)[names(maf) == "Start_Position"] <- "Start_position"
names(maf)[names(maf) == "End_Position"]   <- "End_position"

# Strip 'chr' from chromosome column if present
if ("Chromosome" %in% names(maf)) {
  maf[["Chromosome"]] <- gsub("^chr", "", maf[["Chromosome"]])
}

# SNP-like events
snp_cols <- c("SNP", "DNP", "TNP", "MNP")
snp_maf <- maf[maf[["Variant_Type"]] %in% snp_cols, ]

write.table(
  snp_maf, snp_out,
  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

# INDEL events
indel_cols <- c("INS", "DEL")
indel_maf <- maf[maf[["Variant_Type"]] %in% indel_cols, ]

write.table(
  indel_maf, indel_out,
  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

message("Wrote SNP file:   ", snp_out)
message("Wrote INDEL file: ", indel_out)
EOF

    # Make it executable
    chmod +x maf_to_absolute_inputs.R

    # Run the R script
    Rscript maf_to_absolute_inputs.R ~{maf} ~{sample_id}
  >>>

  output {
    File snp   = "~{sample_id}.snp"
    File indel = "~{sample_id}.indel"
  }

  runtime {
    docker: "rocker/r-base:4.3.2"
    cpu: cpu
    memory: "~{mem_gb}G"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

task make_capseg {
  input {
    File processed_counts
    File segfile
    String participant_id

    Int cpu = 1
    Int mem_gb = 10
    Int disk_gb = 20
  }

  command <<<
    set -euo pipefail

    echo "Downloading capseg_conv.R from GitHub"
    wget -O capseg_conv.R \
      https://raw.githubusercontent.com/beroukhim-lab/terra-absolute/3390dbe6b7c629de136cfab54179c150c93b1d86/capseg_conv.R
    chmod +x capseg_conv.R

    echo "Detecting packages imported by capseg_conv.R and installing any that are missing"
    Rscript -e 'x <- readLines("capseg_conv.R"); \
      pkgs <- grep("^(library|require)\\(", x, value=TRUE); \
      pkgs <- gsub("^(library|require)\\(|\\).*", "", pkgs); \
      pkgs <- gsub(",.*$", "", pkgs); \
      pkgs <- gsub("[\"\\x27[:space:]]", "", pkgs); \
      pkgs <- pkgs[nzchar(pkgs)]; \
      pkgs <- unique(pkgs); \
      cat("Packages referenced:", paste(pkgs, collapse=", "), "\n"); \
      inst <- rownames(installed.packages()); \
      to_install <- setdiff(pkgs, inst); \
      if(length(to_install)) { \
        cat("Installing:", paste(to_install, collapse=", "), "\n"); \
        install.packages(to_install, repos="https://cloud.r-project.org"); \
      } else { \
        cat("No installs needed.\n"); \
      }'

    echo "Running capseg conversion"
    Rscript capseg_conv.R \
      --segfile ~{segfile} \
      --processed_cts ~{processed_counts} \
      --participant_id ~{participant_id}
  >>>

  output {
    File seg_file = "~{participant_id}.capseg.txt"
  }

  runtime {
    docker: "rocker/r-base:4.3.2"
    cpu: cpu
    memory: "~{mem_gb}G"
    disks: "local-disk ~{disk_gb} HDD"
  }
}



workflow preprocess_absolute_capseg {
  input {
    # --- ABSOLUTE inputs ---
    File maf
    String sample_id

    Int maf_to_abs_mem_gb = 8
    Int maf_to_abs_cpu = 1
    Int maf_to_abs_disk_gb = 20

    # --- CAPSEG conversion ---
    File processed_counts
    File segfile
    String participant_id

    Int capseg_cpu = 1
    Int capseg_mem_gb = 10
    Int capseg_disk_gb = 20
  }

  call maf_to_absolute_inputs {
    input:
      maf = maf,
      sample_id = sample_id,
      mem_gb = maf_to_abs_mem_gb,
      cpu = maf_to_abs_cpu,
      disk_gb = maf_to_abs_disk_gb
  }

  call make_capseg {
    input:
      processed_counts = processed_counts,
      segfile = segfile,
      participant_id = participant_id,
      cpu = capseg_cpu,
      mem_gb = capseg_mem_gb,
      disk_gb = capseg_disk_gb
  }

  output {
    File snp        = maf_to_absolute_inputs.snp
    File indel      = maf_to_absolute_inputs.indel
    File capseg_out = make_capseg.seg_file
  }
}
