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
awk 'BEGIN{OFS="\n"} $0~/\/1$/{header = $0; getline seq; getline qheader; getli\
ne qseq; print header, seq, qheader, qseq}' ${id}.fastq > ${id}_1.fastq
awk 'BEGIN{OFS="\n"} $0~/\/2$/{header = $0; getline seq; getline qheader; getli\
ne qseq; print header, seq, qheader, qseq}' ${id}.fastq > ${id}_2.fastq
done

#HLA typing                                                                     
cat config | while read id
do
OptiTypePipeline.py -i ${id}_1.fastq ${id}_2.fastq --rna -v -o ./ -p ${id}.opti\
type.rna
done
