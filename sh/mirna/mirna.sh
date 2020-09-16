index=/home/eunji/proj/0_sh/ref/mirna/hg38/hg38
cat config2 | while read id
do
bowtie2 -p 4 -x $index -1 ./8.mirna/1.trim/${id}_1_clean.fq.gz -2 ./8.mirna/1.t\
rim/${id}_2_clean.fq.gz | samtools sort -@ 4 -o ${id}_genome.bam
done

gtf=/home/eunji/proj/0_sh/ref/mirna/hsa.gff3

featureCounts -T 4 -F gff  -M -t miRNA -g Name  -a $gtf -o all.counts.mature.tx\
t   *genome* 1>counts.mature.log 2>&1
featureCounts -T 4 -F gff  -M -t miRNA_primary_transcript  -g Name  -a $gtf -o \
all.counts.hairpin.txt   *genome* 1>counts.hairpin.log 2>&1
