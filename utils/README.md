# Amaryllis: Utilities

### `subset_fq_difference.R`

Command: `Rscript subset_fq_difference.R <file1.fastq.gz> <file2.fastq.gz>
<output_directory>`  

- Used to partition reads by their uniqueness to one fastq file or the other.
- Creates three partitions: file1-unique, file2-unique, shared in both files.
- Will randomly sample one group if it is smaller than 5% the size of the other.

### `analyze_fq_difference.R`

Command: `Rscript analyze_fq_difference.R <file1.fastq.gz>:<group_1_name>
<file2.fastq.gz>:<group_2_name> <output_directory>`  

- Used to assess the difference between reads contained in two fastq files.
- Can be used with output from `subset_fq_difference.R`.
- Assesses differences in read length, GC content, and read quality.
- Creates a pdf containing three graphs and performs hypothesis tests.

### `survey_quality_improvement.R`

Command: `Rscript survey_quality_improvement.R <raw_fastq_file.fastq.gz>
<processed_fastq_file.fastq.gz> <output_quality_comparison_graph.png>`  

- Used in main pipeline.
- Randomly samples at most 5000 reads from each fastq file compared.
- Assesses difference in read quality between raw and processed fastq data.
- Creates a png file containing superimposed density plots.

### `truncate_fq.py`

Command: `python3 truncate_fq.py <original_file_to_truncate>
<new_truncated_file_to_create> <number_of_lines_to_keep>`  

- Used to remove all lines from a file after a specific line number.
- Compatible with gzipped files.
