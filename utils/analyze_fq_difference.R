# Run with command:
#   Rscript analyze_fq_difference.R \
#           <file1.fastq.gz>:<group_1_name> \
#           <file2.fastq.gz>:<group_2_name> \
#           <output_directory>
options(stringsAsFactors = F)
library("ggplot2")
library("gridExtra")

quality = function(read_quality) {
  scores = as.numeric(factor(unlist(strsplit(read_quality, "")), levels = key))
  return(mean(scores))
}
GC_content = function(sequence) {
  1 - nchar(gsub("[GC]", "", sequence)) / nchar(sequence)
}
key = c("!", '"', "#", "$", "%", "&", "'", "(", ")", "*", "+", ",", "-", ".",
        "/", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", ":", ";", "<",
        "=", ">", "?", "@", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J",
        "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X",
        "Y", "Z", "[", "\\", "]", "^", "_", "`", "a", "b", "c", "d", "e", "f",
        "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t",
        "u", "v", "w", "x", "y", "z", "{", "|", "}", "~")

# Read in FASTQ data
args = strsplit(commandArgs(trailingOnly = T), ":")
fq1 = readLines(args[[1]][1])
group1 = gsub("_", " ", args[[1]][2])
fq2 = readLines(args[[2]][1])
group2 = gsub("_", " ", args[[2]][2])
outdir = args[[3]][1]

# Organize raw data
fq1 = data.frame("Group" = group1,
                 "sequence" = fq1[seq(2, length(fq1), by = 4)],
                 "quality" = fq1[seq(4, length(fq1), by = 4)])
fq2 = data.frame("Group" = group2,
                 "sequence" = fq2[seq(2, length(fq2), by = 4)],
                 "quality" = fq2[seq(4, length(fq2), by = 4)])
fq_data = rbind(fq1, fq2); rm(fq1, fq2)

# Compile useful quantitative data
fq_data$read_length = nchar(fq_data$sequence)
fq_data$GC_content = sapply(fq_data$sequence, GC_content)
fq_data$quality_score = sapply(fq_data$quality, quality)
fq_data[, c("sequence", "quality")] = NULL
rm(list = setdiff(ls(), c("fq_data", "group1", "group2", "outdir")))

# Create graphs for analysis
plots = list()
plots[[1]] = ggplot(fq_data, aes(x = read_length)) +
  labs(x = "Read Length", y = "Density") +
  geom_histogram(aes(group = Group, colour = Group, fill = Group), alpha = 0.5,
                 binwidth = 1, position = position_dodge(width = 0))
plots[[2]] = ggplot(fq_data, aes(x = GC_content)) +
  labs(x = "GC Content", y = "Density") +
  geom_density(aes(group = Group, colour = Group, fill = Group), alpha = 0.5)
plots[[3]] = ggplot(fq_data, aes(x = quality_score)) +
  labs(x = "Quality Score", y = "Density") +
  geom_density(aes(group = Group, colour = Group, fill = Group), alpha = 0.5)

# Save graphs
pdf(paste0(outdir, "/fq_difference.pdf"),
    height = 13, width = 7.5, onefile = T, paper = "us")
do.call("grid.arrange", c(plots, ncol = 1))
dev.off()

# Hypothesis Tests
fq1_mean = sapply(fq_data[fq_data$Group == group1, -1], mean)
fq2_mean = sapply(fq_data[fq_data$Group == group2, -1], mean)
diff_means = fq1_mean - fq2_mean
se = sapply(fq_data[, -1], sd) / sqrt(nrow(fq_data))
z = diff_means / se
p = 2 * pnorm(-abs(z))
results = data.frame("Read Length" = numeric(),
                     "GC Content" = numeric(),
                     "Quality Score" = numeric())
results[1,] = fq1_mean
results[2,] = fq2_mean
results[3,] = diff_means
results[4,] = se
results[5,] = z
results[6,] = p
row.names(results) = c(paste(group1, "mean"),
                       paste(group2, "mean"),
                       "Difference in means",
                       "Standard error",
                       "Test statistic",
                       "P-value")
write.table(results, file = paste0(outdir, "/fq_difference.txt"), sep = "\t")
