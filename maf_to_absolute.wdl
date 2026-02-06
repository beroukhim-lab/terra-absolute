version 1.0

# ------------------------
# Task: Runs your R script
# ------------------------
task maf_to_absolute_inputs {
  input {
    File maf
    String sample_id
  }

  command <<<
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
    memory: "2G"
  }
}

# ------------------------
# Workflow: Calls the task
# ------------------------
workflow maf_to_absolute_workflow {
  input {
    File maf
    String sample_id
  }

  call maf_to_absolute_inputs {
    input:
      maf = maf,
      sample_id = sample_id
  }

  output {
    File snp   = maf_to_absolute_inputs.snp
    File indel = maf_to_absolute_inputs.indel
  }
}