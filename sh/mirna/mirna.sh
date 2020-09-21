#miRNA analysis                                                                
                                                                                

#fastx_trimmer                                                                 
                                                                                
cat config | while read id
do
fastq_quality_filter -v -q 20 -p 80 -Q33  -i ./1.fq/${id}.fastq -o ./8.mirna/1.trim/tmp/${id}.fastq
fastx_trimmer -v -f 1 -l 27 -m 15  -i ./8.mirna/1.trim/tmp/${id}.fastq  -Q33 -z -o ./8.mirna/1.trim/${id}_clean.fq.gz
done

#bowtie mapping     
index=/home/eunji/proj/0_sh/ref/mirna/hg38/hg38
cat config | while read id
do
bowtie2 -p 4 -x $index -1 ./8.mirna/1.trim/${id}_1_clean.fq.gz -2 ./8.mirna/1.trim/${id}_2_clean.fq.gz | samtools sort -@ 4 -o ${id}_genome.bam
done

#featureCounts 
gtf=/home/eunji/proj/0_sh/ref/mirna/hsa.gff3

featureCounts -T 4 -F gff  -M -t miRNA -g Name  -a $gtf -o all.counts.mature.txt   *genome* 1>counts.mature.log 2>&1
featureCounts -T 4 -F gff  -M -t miRNA_primary_transcript  -g Name  -a $gtf -o all.counts.hairpin.txt   *genome* 1>counts.hairpin.log 2>&1
