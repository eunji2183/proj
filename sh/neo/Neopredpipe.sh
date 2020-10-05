## trim_galore.sh                                                               
cat config | while read id
do
        fq1=./1.raw_fq/${id}/${id}_1.fastq.gz
        fq2=./1.raw_fq/${id}/${id}_2.fastq.gz
        echo "start trim_galore for ${id}" `date`
        trim_galore  --paired --phred33 --length 30 --gzip --cores 8 -o ./2.clean_fq  ${fq1}  ${fq2} >> ./2.clean_fq/${id}_trim.log 2>&1
        echo "end trim_galore for ${id}" `date`
done

## clean-qc.sh                                                                  
cat config | while read id
do
        fastqc --outdir ./3.qc/clean_qc/ --threads 16 ./2.clean_fq/${id}*.fq.gz >> ./3.qc/clean_qc/${id}_fastqc.log 2>&1
done

multiqc  ./3.qc/clean_qc/*zip  -o ./3.qc/clean_qc/multiqc
rm ./3.qc/clean_qc/*zip
rm ./3.qc/clean_qc/*.html
rm ./3.qc/clean_qc/*.log

##hisat2_tran.sh                                                                
INDEX=/home/eunji/proj/0_sh/ref/rna/grch38_tran/genome_tran
#INDEX=/home/eunji/proj/0_sh/ref/rna/grcm38_tran/genome_tran

cat config | while read id
do
        echo "start hisat2_tran for ${id}" `date`
	fq1=./2.clean_fq/${id}_1_val_1.fq.gz
        fq2=./2.clean_fq/${id}_2_val_2.fq.gz
        time hisat2 -p 4 -x ${INDEX} -1 ${fq1} -2 ${fq2} -S ./4.align/hisat2_tran/${id}.sam >>./4.align/hisat2_tran/${id}_hisat2.log 2>&1
        echo "end hisat2_tran for ${id}" `date`
done

##view.sh                                                                       
cat config | while read id
do
echo "start samtools view for ${id}" `date`
samtools view -Sb ./4.align/hisat2_tran/${id}.sam > ./4.align/hisat2_tran/${id}.bam
echo "end samtools view for ${id}" `date`
done

# input bam file - HLA typing                                                   

ref=/home/eunji/proj/0_sh/ref/rna/hla_reference_rna.fasta

#bam2fq                                                                         
cat config | while read id
do
samtools bam2fq ./4.align/hisat2_tran/${id}.bam > ${id}.fastq
done

#mapping to hla reference                                                       
cat config | while read id
do
razers3 -i 95 -m 1 -dr 0 --thread-count 10 -o ${id}.bam ${ref} ${id}.fastq
done

#bam2fq                                                                         
cat config | while read id
do
samtools bam2fq ${id}.bam > ${id}.fastq
done

# separate paired fastq file (/1 /2)                                            
cat config | while read id
do
awk 'BEGIN{OFS="\n"} $0~/\/1$/{header = $0; getline seq; getline qheader; getline qseq; print header, seq, qheader, qseq}' ${id}.fastq > ${id}_1.fastq
awk 'BEGIN{OFS="\n"} $0~/\/2$/{header = $0; getline seq; getline qheader; getline qseq; print header, seq, qheader, qseq}' ${id}.fastq > ${id}_2.fastq
done

#HLA typing                                                                     
cat config | while read id
do
OptiTypePipeline.py -i ${id}_1.fastq ${id}_2.fastq --rna -v -o ./ -p ${id}.optitype.rna
done

python NeoPredPipe.py \
-I somatic.vcf \
-H hlatypes.txt \
-o ./ \
-n TestRun \
-c 1 2 -E 8 9 10
