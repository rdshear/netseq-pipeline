# make the STAR index
conda activate STAR

cd /n/groups/churchman/rds19

cd starRefFiles

rm -Rf *
  
gzcat ../../ucsc_golden_path/sacCer3.fa.gz > sacCer3.fa
   
STAR --runThreadN 3 --outFileNamePrefix genome/logs/Index. \
    --genomeSAindexNbases 10 \
        --runMode genomeGenerate \
        --genomeDir genome  \
        --genomeFastaFiles sacCer3.fa

conda activate cpa

picard CreateSequenceDictionary -R sacCer3.fa -O sacCer3.dict
