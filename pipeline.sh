#########
# SETUP #
#########

# Source in shell script containing parameters
. ./parameters.sh

# Locate and create directories used to store data
QA_DIR=$BASE_DIR/qa
RAW_DIR=$BASE_DIR/raw
GROUPED_DIR=$BASE_DIR/grouped
CLIP_DIR=$BASE_DIR/clipped
TRIM_DIR=$BASE_DIR/trimmed
BAM_DIR=$BASE_DIR/bam
COUNT_DIR=$BASE_DIR/counted
ANALYSIS_DIR=$BASE_DIR/analysis

# Locate directory containing adapter_clipper
ADAPTER_CLIPPER=adapter-clipper/clipper.py

# Locate directory containing Trimmomatic
TRIM_BIN=
TRIMMER=$TRIM_BIN/trimmomatic-0.36.jar
ADAPTERS=$TRIM_BIN/adapters/TruSeq3-PE-2.fa:2:30:10

# Locate directory containing read_counter
READ_COUNTER=read_counter/bin/simple_counts.pl

# Create join_by function, which concatenates multiple strings to one
function join_by { local IFS="$1"; shift; echo "$*"; }

# Locate directory containing dge
DGE=dge-analysis/dge.R

############
# PIPELINE #
############

# Stage 0: Create symbolic links back to raw files grouped by sample
for FILE in $(find $RAW_DIR -name "*.fastq*"); do
    FILENAME=${FILE/fq\/raw\//}
    SAMPLE=$(cut -c 1-$TRUNC_AT <<< $FILENAME)
    mkdir -p $GROUPED_DIR/$SAMPLE
    ln -s ../../../$FILE $GROUPED_DIR/$SAMPLE/$FILENAME
done
for SAMPLE in $(ls $GROUPED_DIR); do
    mkdir -p $QA_DIR/raw/$SAMPLE
    fastqc -o $QA_DIR/raw/$SAMPLE $(find $GROUPED_DIR/$SAMPLE -type l)
done

# Stage 1: Clip nucleotides off start of each read
for SAMPLE in $(ls -p $GROUPED_DIR); do
    python3 $ADAPTER_CLIPPER $GROUPED_DIR/$SAMPLE \
            $CLIP_LEN $EXP_LEN
    mkdir -p $QA_DIR/clipped/$SAMPLE
    fastqc -o $QA_DIR/clipped/$SAMPLE $(find $CLIP_DIR/$SAMPLE -name "*.gz")
done

# Stage 2: Use Trimmomatic to do quality trimming
for SAMPLE in $(ls $CLIP_DIR); do
    mkdir -p $TRIM_DIR/$SAMPLE $QA_DIR/trimmed_qc/$SAMPLE \
          $QA_DIR/trimmed_logs/$SAMPLE
    for FILE in $(ls $CLIP_DIR/$SAMPLE | grep -v QC); do
        java -jar $TRIMMER SE -phred33 -threads $TRIM_THREADS \
             $CLIP_DIR/$SAMPLE/$FILE $TRIM_DIR/$SAMPLE/$FILE \
             ILLUMINACLIP:$ADAPTERS \
             LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 \
             2>> $QA_DIR/trimmed_logs/$SAMPLE/${FILE/.fastq.gz/.log}
    done
    fastqc -o $QA_DIR/trimmed_qc/$SAMPLE $(find $TRIM_DIR/$SAMPLE -name "*.gz")
done

# Stage 3: Use Bowtie2 to map to cDNA reference
mkdir $BAM_DIR $QA_DIR/bam_logs # $QA_DIR/bam_qc
for SAMPLE in $(ls $TRIM_DIR); do
    FQ_LIST=$(join_by , $TRIM_DIR/$SAMPLE/*)
    bowtie2 $BOWTIE_PARAMS -p$BOWTIE_THREADS -x $REF_GENOME \
            -U $FQ_LIST \
            2> $QA_DIR/bam_logs/$SAMPLE.log \
            | samtools view -bS -o $BAM_DIR/$SAMPLE.bam
    samtools sort -@$BOWTIE_THREADS \
             $BAM_DIR/$SAMPLE.bam -o $BAM_DIR/$SAMPLE.sorted.bam
done
#fastqc -o $QA_DIR/bam_qc -f bam_mapped $(find $BAM_DIR -name "*.sorted.bam")


# Stage 4: Count reads per gene
mkdir -p $COUNT_DIR
perl $READ_COUNTER -o $COUNT_DIR $(find $BAM_DIR -name "*.sorted.bam")

# Stage 5: Use edgeR to perform differential expression analysis
mkdir -p $ANALYSIS_DIR
Rscript $DGE $METHOD $COUNT_DIR $ANNOTATIONS $SAMPLES $JOBS
