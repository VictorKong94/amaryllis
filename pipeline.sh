#########
# SETUP #
#########

# Source in base directory variable from plaintext file containing parameters
export $(grep 'BASE_DIR=' parameters.txt)

# Locate and establish directories used to store data
BASE_DIR=${BASE_DIR%/}
QA_DIR=$BASE_DIR/qa
RAW_DIR=$BASE_DIR/raw
GROUPED_DIR=$BASE_DIR/grouped
CLIP_DIR=$BASE_DIR/clipped
TRIM_DIR=$BASE_DIR/trimmed
BAM_DIR=$BASE_DIR/bam
COUNT_DIR=$BASE_DIR/counted
ANALYSIS_DIR=$BASE_DIR/analysis

# Source all other parameters necessary for pipeline
INDEX=1
PARAMS=$(grep -A1 '^SUBNAME' parameters.txt | grep -v '^SUBNAME')
while [[ $PARAMS ]]; do
  PARAMS=($PARAMS)
  SUBNAME[$INDEX]=${PARAMS[0]}
  SAMPLE[$INDEX]=${PARAMS[1]}
  CLIP_LEN[$INDEX]=${PARAMS[2]}
  if [[ ${#PARAMS[@]} -eq 11 ]]; then
    TRIM_THREADS[$INDEX]=${PARAMS[3]}
    LEADING[$INDEX]=${PARAMS[4]}
    TRAILING[$INDEX]=${PARAMS[5]}
    SLIDINGWINDOW[$INDEX]=${PARAMS[6]}
    MINLEN[$INDEX]=${PARAMS[7]}
    REF_GENOME[$INDEX]=${PARAMS[8]}
    BOWTIE_PARAMS[$INDEX]=${PARAMS[9]}
    BOWTIE_THREADS[$INDEX]=${PARAMS[10]}
  elif [[ ${#PARAMS[@]} -eq 12 ]]; then
    EXP_LEN[$INDEX]=${PARAMS[3]}
    TRIM_THREADS[$INDEX]=${PARAMS[4]}
    LEADING[$INDEX]=${PARAMS[5]}
    TRAILING[$INDEX]=${PARAMS[6]}
    SLIDINGWINDOW[$INDEX]=${PARAMS[7]}
    MINLEN[$INDEX]=${PARAMS[8]}
    REF_GENOME[$INDEX]=${PARAMS[9]}
    BOWTIE_PARAMS[$INDEX]=${PARAMS[10]}
    BOWTIE_THREADS[$INDEX]=${PARAMS[11]}
  else
    echo "Check format of pipeline parameters in parameters.txt"
    exit 3
  fi
  INDEX=$((INDEX+1))
  PARAMS=$(grep -A1 "$PARAMS" parameters.txt | grep -v "$PARAMS")
done

# Locate directory containing adapter_clipper
ADAPTER_CLIPPER=adapter-clipper/clipper.py

# Locate directory containing Trimmomatic
TRIM_BIN=/installs/Trimmomatic-0.36
TRIMMER=$TRIM_BIN/trimmomatic-0.36.jar
ADAPTERS=$TRIM_BIN/adapters/TruSeq3-PE-2.fa:2:30:10

# Locate directory containing read_counter
READ_COUNTER=read_counter/bin/simple_counts.pl

# Source in parameters necessary for DGE analysis
INDEX=1
PARAMS=$(grep -A1 '^EXP' parameters.txt | grep -v '^EXP')
while [[ $PARAMS ]]; do
  PARAMS=($PARAMS)
  if [[ ${#PARAMS[@]} -eq 5 ]]; then
    EXP[$INDEX]=${PARAMS[0]}
    METHOD[$INDEX]=${PARAMS[1]}
    ANNOTATIONS[$INDEX]=${PARAMS[2]}
    SAMPLES[$INDEX]=${PARAMS[3]}
    JOBS[$INDEX]=${PARAMS[4]}
  else
    echo "Check format of DGE analysis parameters in parameters.txt"
    exit 3
  fi
  INDEX=$((INDEX+1))
  PARAMS=$(grep -A1 "$PARAMS" parameters.txt | grep -v "$PARAMS")
done

# Locate directory containing dge_analysis
DGE=dge-analysis/dge.R

# Create join_by function, which concatenates multiple strings to one
function join_by { local IFS="$1"; shift; echo "$*"; }


############
# PIPELINE #
############

# Stage 0: Create symbolic links back to raw files grouped by sample
for INDEX in $(seq 1 ${#SAMPLE[@]}); do
  SAMPLE=${SAMPLE[$INDEX]}
  mkdir -p $GROUPED_DIR/$SAMPLE
  for FILE in $(find $RAW_DIR -name "${SUBNAME[$INDEX]}*.fastq*"); do
    ln -sf "../..${FILE/$BASE_DIR/}" $GROUPED_DIR/$SAMPLE/${FILE##*/}
  done
  mkdir -p $QA_DIR/raw/$SAMPLE
  fastqc -o $QA_DIR/raw/$SAMPLE $(find $GROUPED_DIR/$SAMPLE -type l)
done

# Stage 1: Clip nucleotides off start of each read
for INDEX in $(seq 1 ${#SAMPLE[@]}); do
  SAMPLE=${SAMPLE[$INDEX]}
  python3 $ADAPTER_CLIPPER $GROUPED_DIR/$SAMPLE \
          ${CLIP_LEN[$INDEX]} ${EXP_LEN[$INDEX]}
  mkdir -p $QA_DIR/clipped/$SAMPLE
  fastqc -o $QA_DIR/clipped/$SAMPLE $(find $CLIP_DIR/$SAMPLE -name "*.fastq*")
done

# Stage 2: Use Trimmomatic to do quality trimming
for INDEX in $(seq 1 ${#SAMPLE[@]}); do
  SAMPLE=${SAMPLE[$INDEX]}
  mkdir -p $TRIM_DIR/$SAMPLE $QA_DIR/trimmed_qc/$SAMPLE \
           $QA_DIR/trimmed_logs/$SAMPLE
  for FILE in $(ls $CLIP_DIR/$SAMPLE); do
    java -jar $TRIMMER SE -phred33 -threads ${TRIM_THREADS[$INDEX]} \
              $CLIP_DIR/$SAMPLE/$FILE $TRIM_DIR/$SAMPLE/$FILE \
              ILLUMINACLIP:$ADAPTERS \
              LEADING:${LEADING[$INDEX]} \
              TRAILING:${TRAILING[$INDEX]} \
              SLIDINGWINDOW:${SLIDINGWINDOW[$INDEX]} \
              MINLEN:${MINLEN[$INDEX]} \
              2>> $QA_DIR/trimmed_logs/$SAMPLE/${FILE/.fastq.gz/.log}
  done
  fastqc -o $QA_DIR/trimmed_qc/$SAMPLE $(find $TRIM_DIR/$SAMPLE -name "*.gz")
done

# Stage 3: Use Bowtie2 to map to cDNA reference
mkdir -p $BAM_DIR $QA_DIR/bam_logs
for INDEX in $(seq 1 ${#SAMPLE[@]}); do
  SAMPLE=${SAMPLE[$INDEX]}
  FQ_LIST=$(join_by , $TRIM_DIR/$SAMPLE/*)
  bowtie2 ${BOWTIE_PARAMS[$INDEX]} \
          -p${BOWTIE_THREADS[$INDEX]} \
          -x ${REF_GENOME[$INDEX]} \
          -U $FQ_LIST \
          2> $QA_DIR/bam_logs/$SAMPLE.log \
          | samtools view -bS -o $BAM_DIR/$SAMPLE.bam
  samtools sort -@$BOWTIE_THREADS[$INDEX] \
                $BAM_DIR/$SAMPLE.bam -o $BAM_DIR/$SAMPLE.sorted.bam
  rm $BAM_DIR/$SAMPLE.bam
done

# Stage 4: Count reads per gene
# sudo mkdir -p repro-archive
# sudo chown $USER:root repro-archive
mkdir -p $COUNT_DIR
perl $READ_COUNTER -o $COUNT_DIR $(find $BAM_DIR -name "*.sorted.bam")

# Stage 5: Use edgeR to perform differential expression analysis
for INDEX in $(seq 1 $NUM_ANALYSES)
  mkdir -p $ANALYSIS_DIR/${EXP[$INDEX]}
  Rscript $DGE $COUNT_DIR \
               ${EXP[$INDEX]} \
               ${METHOD[$INDEX]} \
               ${ANNOTATIONS[$INDEX]} \
               ${SAMPLES[$INDEX]} \
               -j ${JOBS[$INDEX]}
done

