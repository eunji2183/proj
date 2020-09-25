#RT-neoantigen prediction                                                       
#SRR12060755-HCC1143 , SRR12060763-HCC38, SRR12060810-HCC1187 SRR12060756-HCC1395                                                                              
#mkdir                                                                          
#mkdir -p 1.raw_fq 2.clean_fq 3.qc/{raw_qc,clean_qc} 4.align/{hisat2_tran}      


#cd ./0.sra                                                                     
## fasterq-dump                                                                 
#cat config | while read id                                                     
#do                                                                             
        #time fasterq-dump -3 -e 16 ${id}.sra -O ../1.raw_fq --outfile ${id}.fastq                                                                             
#done                                                                           



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
