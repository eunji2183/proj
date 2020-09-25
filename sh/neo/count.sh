## trim_galore.sh                                                               
cat config | while read id
do
        fq1=./1.raw_fq/${id}/${id}_1.fastq.gz
        fq2=./1.raw_fq/${id}/${id}_2.fastq.gz
        echo "start trim_galore for ${id}" `date`
        trim_galore  --paired --phred33 --length 30 --gzip --cores 8 -o ./2.cle\
an_fq  ${fq1}  ${fq2} >> ./2.clean_fq/${id}_trim.log 2>&1
        echo "end trim_galore for ${id}" `date`
done

## clean-qc.sh                                                                  
cat config | while read id
do
        fastqc --outdir ./3.qc/clean_qc/ --threads 16 ./2.clean_fq/${id}*.fq.gz\
 >> ./3.qc/clean_qc/${id}_fastqc.log 2>&1
done

multiqc  ./3.qc/clean_qc/*zip  -o ./3.qc/clean_qc/multiqc
rm ./3.qc/clean_qc/*zip
rm ./3.qc/clean_qc/*.html
rm ./3.qc/clean_qc/*.log

##hisat2_tran.sh                                                                
INDEX=/home/eunji/proj/0_sh/ref/grch38_tran/genome_tran

cat config | while read id
do
        echo "start hisat2_tran for ${id}" `date`
	fq1=./2.clean_fq/${id}_1_val_1.fq.gz
        fq2=./2.clean_fq/${id}_2_val_2.fq.gz
        time hisat2 -p 4 -x ${INDEX} -1 ${fq1} -2 ${fq2} -S ./4.align/hisat2_tr\
an/${id}.sam >>./4.align/hisat2_tran/${id}_hisat2.log 2>&1
        echo "end hisat2_tran for ${id}" `date`
done

##view.sh                                                                       
cat config | while read id
do
echo "start samtools view for ${id}" `date`
samtools view -Sb ./4.align/hisat2_tran/${id}.sam > ./4.align/hisat2_tran/${id}\
.bam
echo "end samtools view for ${id}" `date`
done

#sort_name.sh                                                                   
cat config | while read id
do
echo "start samtools sort name for ${id}" `date`
samtools sort -n -@ 10 ./4.align/hisat2_tran/${id}.bam -o ./4.align/hisat2_tran\
/${id}_nsort.bam
echo "end samtools sort name for ${id}" `date`
done

##featureCounts.sh                                                              
GTF=/home/eunji/proj/0_sh/ref/Homo_sapiens.GRCh38.96.gtf
cat config | while read id
do
 echo "start featureCounts for ${id}" `date`
 featureCounts -T 6 -p -t exon -g gene_id -a $GTF -o ./5.DEG/featurecounts/${id\
}.txt ./4.align/hisat2_tran/${id}_nsort.bam >>./5.DEG/featurecounts/${id}_featu\
recounts.log 2>&1
 echo "end featureCounts for ${id}" `date`
done



