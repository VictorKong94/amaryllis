# Run with command:
#   Rscript subset_fq_difference.R file1.fq.gz file2.fq.gz output_directory

# Pull relevant files' names from arguments
args = commandArgs(trailingOnly = T)
file1 = args[1]
file2 = args[2]
outpath = args[3]

# Determine the paths and names of the files to be written
if (!dir.exists(outpath)) dir.create(outpath, recursive = T)
if (!grepl("/$", outpath)) outpath = paste0(outpath, "/")


#################
# GENERAL SETUP #
#################

# Extract read information
fq1 = readLines(file1)
fq2 = readLines(file2)
print("Infiles read")

# Extract read IDs
fq1_ids = fq1[seq(1, length(fq1), by = 4)]
fq2_ids = fq2[seq(1, length(fq1), by = 4)]
print("Read IDs extracted")

# Sort reads by uniqueness to one file
fq1_only = setdiff(fq1_ids, fq2_ids)
fq2_only = setdiff(fq2_ids, fq1_ids)
both_fqs = intersect(fq1_ids, fq2_ids)

fq1_len = length(fq1_only)
fq2_len = length(fq2_only)
uniq_len = max(fq1_len, fq2_len)
both_len = length(both_fqs)

if (both_len > 20 * uniq_len) {
  both_fqs = sample(both_fqs, uniq_len)
} else if (both_len < 0.05 * uniq_len) {
  if (fq1_len > both_len) fq1_only = sample(fq1_only, both_len)
  if (fq2_len > both_len) fq2_only = sample(fq2_only, both_len)
} else if (both_len > 5 * uniq_len) {
  both_fqs = sample(both_fqs, uniq_len, replace = T)
} else if (both_len < 0.2 * uniq_len) {
  if (fq1_len > both_len) fq1_only = sample(fq1_only, both_len, replace = T)
  if (fq2_len > both_len) fq2_only = sample(fq2_only, both_len, replace = T)
}
print("Differential analysis complete")


####################################
# SORT READS INTO INDIVIDUAL FILES #
####################################

fq1_only_outfile = paste0(outpath, "fq1_only.fastq.gz")
fq2_only_outfile = paste0(outpath, "fq2_only.fastq.gz")
both_fqs_outfile = paste0(outpath, "both_fqs.fastq.gz")

sort_output = function(my_ids, all_ids, all_data) {
  if (is.na(my_ids)) {
    return(character(0))
  } else {
    return(sapply(my_ids, function(x)
      all_data[seq(4 * match(x, all_ids) - 3, length.out = 4)]))
  }
}

# Write subsets of data into separate files
con = gzfile(fq1_only_outfile, "wb")
writeLines(as.vector(sort_output(fq1_only, fq1_ids, fq1)), con)
close(con)
rm(fq1_only, fq1_ids, fq1)
print("Reads unique to first FASTQ written to file")

con = gzfile(fq2_only_outfile, "wb")
writeLines(as.vector(sort_output(fq2_only, fq2_ids, fq2)), con)
close(con)
rm(fq2_only)
print("Reads unique to second FASTQ written to file")

con = gzfile(both_fqs_outfile, "wb")
writeLines(as.vector(sort_output(both_fqs, fq2_ids, fq2)), con)
close(con)
rm(fq2_ids, fq2, both_fqs)
print("Reads belonging to both FASTQs written to file")
