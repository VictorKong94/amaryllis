# Run with command:
#   Rscript analyze_fq_difference.R \
#           file1.fq.gz:group_1_name \
#           file2.fq.gz:group_2_name \
#           output_directory
options(stringsAsFactors = F)

quality = function(read_quality) {
  scores = as.numeric(factor(unlist(strsplit(read_quality, "")), levels = key))
  mean_score = mean(scores)
  mean_rank = mean(as.numeric(factor(scores, levels = sort(unique(scores)))))
  return(c(mean_score, mean_rank))
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
args = strsplit(commandArgs(), ":")
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
fq_data[, c("quality_score", "quality_rank")] =
  matrix(unlist(lapply(fq_data$quality, quality)), ncol = 2, byrow = T)
fq_data[, c("sequence", "quality")] = NULL
rm(list = setdiff(ls(), "fq_data"))

# Create graphs for analysis
plots = list()
plots[[1]] = ggplot(fq_data, aes(x = read_length)) +
  labs(x = "Read Length", y = "Density") +
  geom_histogram(aes(Group = Group, colour = Group, fill = Group), alpha = 0.5,
                 binwidth = 1, position = position_dodge(width = 0))
plots[[2]] = ggplot(fq_data, aes(x = GC_content)) +
  labs(x = "GC Content", y = "Density") +
  geom_density(aes(Group = Group, colour = Group, fill = Group), alpha = 0.5)
plots[[3]] = ggplot(fq_data, aes(x = quality_score)) +
  labs(x = "Quality Score", y = "Density") +
  geom_density(aes(Group = Group, colour = Group, fill = Group), alpha = 0.5)
plots[[4]] = ggplot(fq_data, aes(x = quality_rank)) +
  labs(x = "Quality Rank", y = "Density") +
  geom_density(aes(Group = Group, colour = Group, fill = Group), alpha = 0.5)

# Save graphs
library(ggplot2)
library(gridExtra)
pdf(paste0(outdir, "/fq_difference.pdf"),
    height = 7.5, width = 13, onefile = T, paper = "USr")
do.call("grid.arrange", c(plots, ncol = 2))
dev.off()

# Hypothesis Tests
fq1_mean = sapply(fq_data[fq_data$Group == group1, -1], mean)
fq2_mean = sapply(fq_data[fq_data$Group == group2, -1], mean)
diff_means = fq1_mean - both_mean
se = sapply(fq_data[, -1], sd) / sqrt(nrow(fq_data))
z = diff_means / se
p = 2 * pnorm(-abs(z))
write(c(paste(group1, "mean") = fq1_mean,
        paste(group2, "mean") = fq2_mean,
        "difference(mean)" = diff_means,
        "standard error" = se,
        "test statistic" = z,
        "p-value" = p),
      file = paste0(outdir, "/fq_difference.txt"))
