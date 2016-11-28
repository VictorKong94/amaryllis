# Pull relevant files' names from arguments
args = commandArgs(trailingOnly = T)
method = args[1]
valid_methods = c("-b", "--both",
                  "-q", "--quality-assessment",
                  "-s", "--sort-by-uniqueness")
if (!(method %in% valid_methods)) {
  stop("Invalid method argument")
}
file1 = args[2]
file2 = args[3]
outpath = args[4]

if (dir.exists(outpath)) {
  stop("Output directory must not already exist")
} else {
  dir.create(outpath, recursive = T)
}

# Determine the paths and names of the files to be written
if (!grepl("/$", outpath)) {
  outpath = paste0(outpath, "/")
}

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

if (both_len > 10 * uniq_len) {
  both_fqs = sample(both_fqs, uniq_len)
} else if (both_len < 0.1 * uniq_len) {
  if (fq1_len > both_len) fq1_only = sample(fq1_only, both_len)
  if (fq2_len > both_len) fq2_only = sample(fq2_only, both_len)
}
print("Differential analysis complete")

if (method %in% c("-s", "--sort-options")) {
  
  ####################################
  # SORT READS INTO INDIVIDUAL FILES #
  ####################################
  
  fq1_only_outfile = paste0(outpath, "fq1_only.fastq.gz")
  fq2_only_outfile = paste0(outpath, "fq2_only.fastq.gz")
  both_fqs_outfile = paste0(outpath, "both_fqs.fastq.gz")
  
  sort_output = function(my_ids, all_ids, all_data) {
    sapply(my_ids, function(x) all_data[seq(match(x, all_ids), length.out = 4)])
  }
  
  # Write subsets of data into separate files
  con = gzfile(fq1_only_outfile, "wb")
  writeLines(as.vector(sort_output(fq1_only, fq1_ids, fq1)), con)
  close(con)
  rm(fq1_only, fq1_ids, fq1)
  print("Reads unique to first FASTQ written to file")
  
  con = gzfile(fq1_only_outfile, "wb")
  writeLines(as.vector(sort_output(fq2_only, fq2_ids, fq2)), con)
  close(con)
  rm(fq2_only)
  print("Reads unique to second FASTQ written to file")
  
  con = gzfile(fq1_only_outfile, "wb")
  writeLines(as.vector(sort_output(both_fqs, fq2_ids, fq2)), con)
  close(con)
  rm(fq2_ids, fq2, both_fqs)
  print("Reads belonging to both FASTQs written to file")
  
} else if (method %in% c("-q", "--quality-assessment")) {
  
  ######################
  # QUALITY ASSESSMENT #
  ######################
  
  quality_data_outfile = paste0(outpath, "quality_assessment.txt")
  quality_plot_outfile = paste0(outpath, "quality_assessment.pdf")
  
  quality = function(my_ids, all_ids, quality_scores) {
    # Define the the key for quality scores
    key = c("!", '"', "#", "$", "%", "&", "'", "(", ")", "*", "+", ",", "-",
            ".", "/", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", ":",
            ";", "<", "=", ">", "?", "@", "A", "B", "C", "D", "E", "F", "G",
            "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T",
            "U", "V", "W", "X", "Y", "Z", "[", "\\", "]", "^", "_", "`", "a",
            "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n",
            "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", "{",
            "|", "}", "~")
    all_quality = sapply(my_ids, function(x) quality_scores[match(x, all_ids)])
    quality_assessment = factor(unlist(strsplit(all_quality, "")), levels = key)
    return(quality_assessment)
  }
  
  # Extract read quality
  fq1_quality = fq1[seq(4, length(fq1), by = 4)]
  fq2_quality = fq2[seq(4, length(fq1), by = 4)]
  rm(fq1, fq2)
  print("Read quality extracted")
  
  # Determine quality of each class of reads
  qualities = list()
  
  qualities$"FASTQ 1 Only" = quality(fq1_only, fq1_ids, fq1_quality)
  rm(fq1_only, fq1_ids, fq1_quality)
  print("FASTQ 1 only read quality determined")
  
  qualities$"FASTQ 2 Only" = quality(fq2_only, fq2_ids, fq2_quality)
  rm(fq2_only)
  print("FASTQ 2 only read quality determined")
  
  qualities$"Both FASTQs" = quality(both_fqs, fq2_ids, fq2_quality)
  rm(fq2_ids, fq2_quality, both_fqs)
  print("Both FASTQs read quality determined")
  
  # Create a density plot of each set's qualities
  qualities = data.frame(
    "Score" = unlist(qualities, use.names = F),
    "Group" = rep(names(qualities), times = sapply(qualities, length)))
  qualities = qualities[!is.na(qualities$Score),]
  write.table(data.frame("Score" = as.numeric(qualities$Score),
                         "Group" = qualities$Group),
              file = quality_data_outfile, row.names = F)
  library("ggplot2", quietly = T)
  pdf(NULL)
  ggplot(qualities, aes(x = Score)) +
    labs(x = "Quality Score", y = "Density") +
    geom_density(aes(group = Group, colour = Group, fill = Group), alpha = 0.3)
  ggsave(quality_plot_outfile, width = 11, height = 8.5, units = "in")

} else if (method %in% c("-b", "--both")) {
  
  fq1_only_outfile = paste0(outpath, "fq1_only.fastq.gz")
  fq2_only_outfile = paste0(outpath, "fq2_only.fastq.gz")
  both_fqs_outfile = paste0(outpath, "both_fqs.fastq.gz")
  quality_data_outfile = paste0(outpath, "quality_assessment.txt")
  quality_plot_outfile = paste0(outpath, "quality_assessment.pdf")
  
  sort_output = function(my_ids, all_ids, all_data) {
    sapply(my_ids, function(x) all_data[seq(match(x, all_ids), length.out = 4)])
  }
  
  quality = function(my_ids, all_ids, quality_scores) {
    # Define the the key for quality scores
    key = c("!", '"', "#", "$", "%", "&", "'", "(", ")", "*", "+", ",", "-",
            ".", "/", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", ":",
            ";", "<", "=", ">", "?", "@", "A", "B", "C", "D", "E", "F", "G",
            "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T",
            "U", "V", "W", "X", "Y", "Z", "[", "\\", "]", "^", "_", "`", "a",
            "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n",
            "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", "{",
            "|", "}", "~")
    all_quality = sapply(my_ids, function(x) quality_scores[match(x, all_ids)])
    quality_assessment = factor(unlist(strsplit(all_quality, "")), levels = key)
    return(quality_assessment)
  }
  
  # Extract quality while writing to file
  qualities = list()
  
  fq1_quality = fq1[seq(4, length(fq1), by = 4)]
  qualities$"FASTQ 1 Only" = quality(fq1_only, fq1_ids, fq1_quality)
  rm(fq1_quality)
  con = gzfile(fq1_only_outfile, "wb")
  writeLines(as.vector(sort_output(fq1_only, fq1_ids, fq1)), con)
  close(con)
  rm(fq1_only, fq1_ids, fq1)
  print("Reads unique to first FASTQ written to file")
  
  fq2_quality = fq2[seq(4, length(fq1), by = 4)]
  qualities$"FASTQ 2 Only" = quality(fq2_only, fq2_ids, fq2_quality)
  con = gzfile(fq1_only_outfile, "wb")
  writeLines(as.vector(sort_output(fq2_only, fq2_ids, fq2)), con)
  close(con)
  rm(fq2_only)
  print("Reads unique to second FASTQ written to file")

  qualities$"Both FASTQs" = quality(both_fqs, fq2_ids, fq2_quality)
  rm(fq2_quality)
  con = gzfile(both_fqs_outfile, "wb")
  writeLines(as.vector(sort_output(both_fqs, fq2_ids, fq2)), con)
  close(con)
  rm(fq2_ids, fq2, both_fqs)
  print("Reads belonging to both FASTQs written to file")
  
  # Create a density plot of each set's qualities
  qualities = data.frame(
    "Score" = unlist(qualities, use.names = F),
    "Group" = rep(names(qualities), times = sapply(qualities, length)))
  qualities = qualities[!is.na(qualities$Score),]
  write.table(data.frame("Score" = as.numeric(qualities$Score),
                         "Group" = qualities$Group),
              file = quality_data_outfile, row.names = F)
  library("ggplot2", quietly = T)
  pdf(NULL)
  ggplot(qualities, aes(x = Score)) +
    labs(x = "Quality Score", y = "Density") +
    geom_density(aes(group = Group, colour = Group, fill = Group), alpha = 0.3)
  ggsave(quality_plot_outfile, width = 11, height = 8.5, units = "in")
  
}
