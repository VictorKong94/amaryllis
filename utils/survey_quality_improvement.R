# Run with command:
#   Rscript survey_quality_improvement.R \
#           <raw_fastq_file.fastq.gz> \
#           <processed_fastq_file.fastq.gz> \
#           <output_quality_comparison_graph.png>

# Set up utilities for quality assessment
library("ggplot2")
key = c("!", '"', "#", "$", "%", "&", "'", "(", ")", "*", "+", ",", "-", ".",
        "/", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", ":", ";", "<",
        "=", ">", "?", "@", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J",
        "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X",
        "Y", "Z", "[", "\\", "]", "^", "_", "`", "a", "b", "c", "d", "e", "f",
        "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t",
        "u", "v", "w", "x", "y", "z", "{", "|", "}", "~")
quality = function(read_quality) {
  scores = as.numeric(factor(unlist(strsplit(read_quality, "")), levels = key))
  return(mean(scores))
}

# Pull relevant files' names from arguments
args = commandArgs(trailingOnly = T)
raw_file = args[1]
processed_file = args[2]
outfile = args[3]

########################
# Extract quality data #
########################

# Establish connections to input files
raw_file = file(raw_file, "rt")
processed_file = file(processed_file, "rt")

# Extract read qualities and partition by uniqueness to file
dropped = character()
kept = character()
raw = readLines(raw_file, 4)
processed = readLines(processed_file, 4)
while (length(raw) != 0) {
  if (length(processed) == 0) {
    dropped = c(dropped, raw[4])
  } else if (raw[1] == processed[1]) {
    kept = c(kept, raw[4])
    processed = readLines(processed_file, 4)
  } else {
    dropped = c(dropped, raw[4])
  }
  raw = readLines(raw_file, 4)
}

# Randomly subset partitioned read qualities
# - Comment out set.seed for true randomness; leave uncommented for sanity check
set.seed(314159)
sample_length = min(length(dropped), length(kept), 5000)
dropped = sample(dropped, sample_length)
kept = sample(kept, sample_length)
quality_data = data.frame("Fate" = c(rep("Dropped", sample_length),
                                     rep("Kept", sample_length)),
                          "Quality_Score" = sapply(c(dropped, kept), quality))

################################################
# Create Density Plots of Quality Distribution #
################################################

pdf(NULL)
ggplot(quality_data, aes(x = Quality_Score)) +
  labs(x = "Quality Score", y = "Density") +
  geom_density(aes(group = Fate, colour = Fate, fill = Fate), alpha = 0.5) +
  ggtitle(paste("Quality Improvement in",
                rev(strsplit(outfile, "[\\./]")[[1]])[2]))
ggsave(outfile, width = 7.5, height = 6, units = "in")
