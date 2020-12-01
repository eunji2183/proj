
## trim_galore.sh
cat config | while read id
do
        fq1=./raw_fq/${id}_1.fastq.gz
        fq2=./raw_fq/${id}_2.fastq.gz
        echo "start trim_galore for ${id}" `date`
        trim_galore --nextera --gzip --cores 8 -o ./trim_fq  $fq1  $fq2 >> ./trim_fq/${id}_trim.log 2>&1
        echo "end trim_galore for ${id}" `date`
done
    
