# Create an index for the Mus_musculus.GRCm38.cdna reference genome
FILES=(Mus_musculus.GRCm38.cdna.all.1.bt2 \
       Mus_musculus.GRCm38.cdna.all.2.bt2 \
       Mus_musculus.GRCm38.cdna.all.3.bt2 \
       Mus_musculus.GRCm38.cdna.all.4.bt2 \
       Mus_musculus.GRCm38.cdna.all.rev.1.bt2 \
       Mus_musculus.GRCm38.cdna.all.rev.2.bt2)
cd reference
for FILE in $FILES; do
  if [[ ! -e $FILE ]]; then
    gunzip Mus_musculus.GRCm38.cdna.all.fa.gz
    bowtie2-build Mus_musculus.GRCm38.cdna.all.fa Mus_musculus.GRCm38.cdna.all
    gzip Mus_musculus.GRCm38.cdna.all.fa
    break
  fi
done
cd ..

# Confirm that checking of pipeline parameter checks exit properly
sh ../pipeline.sh parameters/parameters_error1.txt &> /dev/null
if [[ $? -eq 3 ]]; then
  sh ../pipeline.sh parameters/parameters_error2.txt &> /dev/null
  if [[ $? -eq 3 ]]; then
    echo "Pipeline parameters check errors correctly"
  else
    echo "Pipeline parameters check does not error as expected"
  fi
else
  echo "Pipeline parameters check does not error as expected"
fi

# Run full pipeline on datasets `fq` and `ftp-data`
sh ../pipeline.sh parameters/parameters_fq.txt
sh ../pipeline.sh parameters/parameters_ftp-data.txt
