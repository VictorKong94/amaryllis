# Locate principal working directory
BASE_DIR=

# Locate directory containing Trimmomatic
TRIM_BIN=

# Number of characters to read from filenames to determine corresponding sample
TRUNC_AT=

# For clipping:
# - Number of nucleotides to clip from 5' end
CLIP_LEN=8
# - Expected length of reads
EXP_LEN=

# For Trimmomatic:
TRIM_THREADS=25
LEADING=3
TRAILING=3
SLIDINGWINDOW=4:15
MINLEN=50

# Reference genome to which to map reads
REF_GENOME=

# For Bowtie2:
# - Available params:
#   - `--nofw` do not align forward (original) version of read (off)
#   - `--norc` do not align reverse-complement version of read (off)
BOWTIE_PARAMS=
BOWTIE_THREADS=5

# For dge_analysis.R:
# - Available methods:
#   - `--et` Fisher's exact test (suitable for experiments with a single factor)
#   - `--lrt` likelihood ratio test
#   - `--qlf` quasi-likelihood F-test
# - Jobs
#   - Format is "-j <baseline1>,<treatment1> <baseline2>,<treatment2> ..."
METHOD=
ANNOTATIONS=
SAMPLES=
JOBS=
