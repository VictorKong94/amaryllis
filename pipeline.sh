#########
# SETUP #
#########

# Locate directory containing Trimmomatic
TRIM_BIN=/installs/Trimmomatic-0.36
TRIMMER=$TRIM_BIN/trimmomatic-0.36.jar
ADAPTERS=$TRIM_BIN/adapters/TruSeq3-PE-2.fa:2:30:10

# Locate survey_quality_improvement script
SURVEY_QUALITY_IMPROVEMENT=/data/amaryllis/utils/survey_quality_improvement.R

# Locate read_counter script
READ_COUNTER=/data/amaryllis/read_counter/bin/simple_counts.pl

# Locate dge_analysis script
DGE=/data/amaryllis/dge-analysis/dge.R


#######################################
# Nothing below here should be edited #
#######################################

# Run with command `sh pipeline.sh <parameters_file>`
PARAMETERS_FILE=$1

# Source in base directory variable from plaintext file containing parameters
export $(grep 'BASE_DIR=' $PARAMETERS_FILE)

# Locate and establish directories used to store data
BASE_DIR=${BASE_DIR%/}
LOG_DIR=$BASE_DIR/logs
QA_DIR=$BASE_DIR/qa
RAW_DIR=$BASE_DIR/raw
GROUPED_DIR=$BASE_DIR/grouped
TRIM_DIR=$BASE_DIR/trimmed
BAM_DIR=$BASE_DIR/bam
COUNT_DIR=$BASE_DIR/counted
ANALYSIS_DIR=$BASE_DIR/analysis

# Check format of parameters necessary for pipeline and source
INDEX=1
PARAMS=($(grep -A1 '^SUBNAME' $PARAMETERS_FILE | grep -v '^SUBNAME'))
if [[ ${#PARAMS[@]} -eq 11 ]]; then
  SUBNAME[1]=${PARAMS[0]}
  SAMPLE[1]=${PARAMS[1]}
  TRIM_THREADS[1]=${PARAMS[2]}
  HEADCROP[1]=${PARAMS[3]}
  LEADING[1]=${PARAMS[4]}
  TRAILING[1]=${PARAMS[5]}
  SLIDINGWINDOW[1]=${PARAMS[6]}
  MINLEN[1]=${PARAMS[7]}
  REF_GENOME[1]=${PARAMS[8]}
  BOWTIE_PARAMS[1]=${PARAMS[9]}
  BOWTIE_THREADS[1]=${PARAMS[10]}
else
  echo "Check format of pipeline parameters in $PARAMETERS_FILE" >&2
  exit 3
fi
INDEX=2
PARAMS=($(grep -A1 "$PARAMS" $PARAMETERS_FILE | grep -v "$PARAMS"))
# Assume the first set of parameters for all samples if only those are specified
if [[ ${#PARAMS[@]} -eq 2 ]]; then
  while [[ ${PARAMS[*]} ]]; do
    SUBNAME[$INDEX]=${PARAMS[0]}
    SAMPLE[$INDEX]=${PARAMS[1]}
    # Error if parameters are intermittently specified
    if [[ ${#PARAMS[@]} -ne 2 ]]; then
      echo "Check format of pipeline parameters in $PARAMETERS_FILE" >&2
      exit 3
    else
      TRIM_THREADS[$INDEX]=${TRIM_THREADS[1]}
      HEADCROP[$INDEX]=${HEADCROP[1]}
      LEADING[$INDEX]=${LEADING[1]}
      TRAILING[$INDEX]=${TRAILING[1]}
      SLIDINGWINDOW[$INDEX]=${SLIDINGWINDOW[1]}
      MINLEN[$INDEX]=${MINLEN[1]}
      REF_GENOME[$INDEX]=${REF_GENOME[1]}
      BOWTIE_PARAMS[$INDEX]=${BOWTIE_PARAMS[1]}
      BOWTIE_THREADS[$INDEX]=${BOWTIE_THREADS[1]}
    fi
    INDEX=$((INDEX+1))
    PARAMS=($(grep -A1 "$PARAMS" $PARAMETERS_FILE | grep -v "$PARAMS"))
  done
# Else source separately each set of parameters
else
  while [[ ${PARAMS[*]} ]]; do
    if [[ ${#PARAMS[@]} -eq 11 ]]; then
      SUBNAME[$INDEX]=${PARAMS[0]}
      SAMPLE[$INDEX]=${PARAMS[1]}
      TRIM_THREADS[$INDEX]=${PARAMS[2]}
      HEADCROP[$INDEX]=${PARAMS[3]}
      LEADING[$INDEX]=${PARAMS[4]}
      TRAILING[$INDEX]=${PARAMS[5]}
      SLIDINGWINDOW[$INDEX]=${PARAMS[6]}
      MINLEN[$INDEX]=${PARAMS[7]}
      REF_GENOME[$INDEX]=${PARAMS[8]}
      BOWTIE_PARAMS[$INDEX]=${PARAMS[9]}
      BOWTIE_THREADS[$INDEX]=${PARAMS[10]}
    else
      echo "Check format of pipeline parameters in $PARAMETERS_FILE" >&2
      exit 3
    fi
    INDEX=$((INDEX+1))
    PARAMS=($(grep -A1 "$PARAMS" $PARAMETERS_FILE | grep -v "$PARAMS"))
  done
fi

# Source in parameters necessary for DGE analysis
INDEX=1
PARAMS=($(grep -A1 '^EXP' $PARAMETERS_FILE | grep -v '^EXP'))
while [[ ${PARAMS[*]} ]]; do
  if [[ ${#PARAMS[@]} -eq 5 ]]; then
    EXP[$INDEX]=${PARAMS[0]}
    METHOD[$INDEX]=${PARAMS[1]}
    ANNOTATIONS[$INDEX]=${PARAMS[2]}
    SAMPLES[$INDEX]=${PARAMS[3]}
    JOBS[$INDEX]=${PARAMS[4]}
  else
    echo "Check format of DGE analysis parameters in $PARAMETERS_FILE" >&2
    exit 3
  fi
  INDEX=$((INDEX+1))
  PARAMS=($(grep -A1 "$PARAMS" $PARAMETERS_FILE | grep -v "$PARAMS"))
done

# Create join_by function, which concatenates multiple strings to one
function join_by { local IFS="$1"; shift; echo "$*"; }


############
# PIPELINE #
############

# Stage 0: Create symbolic links back to raw files grouped by sample
for INDEX in $(seq 1 ${#SAMPLE[@]}); do
  SAMPLE_I=${SAMPLE[$INDEX]}
  mkdir -p $GROUPED_DIR/$SAMPLE_I
  for FILE in $(find $RAW_DIR -name "${SUBNAME[$INDEX]}*.fastq*"); do
    ln -sf "../..${FILE/$BASE_DIR/}" $GROUPED_DIR/$SAMPLE_I/${FILE##*/}
  done
  mkdir -p $QA_DIR/raw/$SAMPLE_I
  fastqc -o $QA_DIR/raw/$SAMPLE_I $(find $GROUPED_DIR/$SAMPLE_I -type l) \
            &> /dev/null
done

# Stage 1: Use Trimmomatic to do quality trimming
for INDEX in $(seq 1 ${#SAMPLE[@]}); do
  SAMPLE_I=${SAMPLE[$INDEX]}
  mkdir -p $TRIM_DIR/$SAMPLE_I $LOG_DIR/trimmed/$SAMPLE_I \
           $QA_DIR/quality_improvement/$SAMPLE_I $QA_DIR/trimmed/$SAMPLE_I
  for FILE in $(ls $GROUPED_DIR/$SAMPLE_I); do
    LOG=$LOG_DIR/trimmed/$SAMPLE_I/${FILE/.fastq.gz/.log}
    printf '%s\n%s' '---' 'Java: Started: ' >> $LOG
    date >> $LOG
    printf '%s\n' 'Java: Version and Environment information:' >> $LOG
    java -version 2>> $LOG
    java -jar $TRIMMER SE -phred33 -threads ${TRIM_THREADS[$INDEX]} \
              $GROUPED_DIR/$SAMPLE_I/$FILE $TRIM_DIR/$SAMPLE_I/$FILE \
              HEADCROP:${HEADCROP[$INDEX]} \
              ILLUMINACLIP:$ADAPTERS \
              LEADING:${LEADING[$INDEX]} \
              TRAILING:${TRAILING[$INDEX]} \
              SLIDINGWINDOW:${SLIDINGWINDOW[$INDEX]} \
              MINLEN:${MINLEN[$INDEX]} \
              2>> $LOG
    printf '%s' 'Java: Finished: ' >> $LOG
    date >> $LOG
    Rscript $SURVEY_QUALITY_IMPROVEMENT \
            $GROUPED_DIR/$SAMPLE_I/$FILE $TRIM_DIR/$SAMPLE_I/$FILE \
            $QA_DIR/quality_improvement/$SAMPLE_I/${FILE/.fastq.gz/.png} \
            &> /dev/null
  done
  fastqc -o $QA_DIR/trimmed/$SAMPLE_I $(find $TRIM_DIR/$SAMPLE_I -type f) \
            &> /dev/null
done

# Stage 2: Use Bowtie2 to align reads to cDNA reference
mkdir -p $BAM_DIR $LOG_DIR/bam $QA_DIR/bam
for INDEX in $(seq 1 ${#SAMPLE[@]}); do
  SAMPLE_I=${SAMPLE[$INDEX]}
  FQ_LIST=$(join_by , $TRIM_DIR/$SAMPLE_I/*)
  LOG=$LOG_DIR/bam/$SAMPLE_I.log
  BT2_COMMAND="bowtie2 ${BOWTIE_PARAMS[$INDEX]} \
                       -p${BOWTIE_THREADS[$INDEX]} \
                       -x ${REF_GENOME[$INDEX]} \
                       -U $FQ_LIST \
                       2>> $LOG \
                       | samtools view -bS -o $BAM_DIR/$SAMPLE_I.bam"
  SAM_COMMAND="samtools sort -@$BOWTIE_THREADS[$INDEX] \
                             $BAM_DIR/$SAMPLE_I.bam \
                             -o $BAM_DIR/$SAMPLE_I.sorted.bam"
  # Align reads using Bowtie 2
  printf '%s\n%s' '---' 'Bowtie 2: Started: ' >> $LOG
  date >> $LOG
  printf '%s\n' 'Bowtie 2: Command line arguments:' >> $LOG
  printf '%s ' $BT2_COMMAND >> $LOG
  printf '\n%s\n' 'Bowtie 2: Version and Environment information:' >> $LOG
  bowtie2 --version >> $LOG
  eval "$BT2_COMMAND"
  printf '%s' 'Bowtie 2: Finished: ' >> $LOG
  date >> $LOG
  # Sort aligned reads using SAMtools
  printf '%s' 'SAMtools: Started: ' >> $LOG
  date >> $LOG
  printf '%s\n' 'SAMtools: Command line arguments:' >> $LOG
  printf '%s ' $SAM_COMMAND >> $LOG
  printf '\n%s\n' 'SAMtools: Version and Environment information:' >> $LOG
  samtools --version >> $LOG
  eval "$SAM_COMMAND"
  printf '%s' 'SAMtools: Finished: ' >> $LOG
  date >> $LOG
  rm $BAM_DIR/$SAMPLE_I.bam
done
fastqc -o $QA_DIR/bam -f bam_mapped $(find $BAM_DIR -type f) &> /dev/null

# Stage 3: Count reads per gene
# sudo mkdir -p repro-archive
# sudo chown $USER:root repro-archive
mkdir -p $COUNT_DIR
perl $READ_COUNTER -o $COUNT_DIR $(find $BAM_DIR -name "*.sorted.bam")

# Stage 4: Use edgeR to perform differential expression analysis
for INDEX in $(seq 1 ${#EXP[@]}); do
  mkdir -p $ANALYSIS_DIR/${EXP[$INDEX]} $LOG_DIR/analysis
  Rscript $DGE $COUNT_DIR \
               ${EXP[$INDEX]} \
               ${METHOD[$INDEX]} \
               ${ANNOTATIONS[$INDEX]} \
               ${SAMPLES[$INDEX]} \
               ${JOBS[$INDEX]} &> /dev/null
done
