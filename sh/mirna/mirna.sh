#miRNA analysis                                                                
                                                                                
wget ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz   ##　28645　reads
wget ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.zip   ##   35828 reads 
wget ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.zip
wget ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3 
wget ftp://mirbase.org/pub/mirbase/CURRENT/miFam.dat.zip

perl -alne '{if(/^>/){if(/Homo/){$tmp=1}else{$tmp=0}};next if $tmp!=1;s/U/T/g if !/>/;print }' hairpin.fa >hairpin.human.fa
perl -alne '{if(/^>/){if(/Homo/){$tmp=1}else{$tmp=0}};next if $tmp!=1;s/U/T/g if !/>/;print }' mature.fa >mature.human.fa

#index
bowtie2-build hairpin.human.fa hairpin_human
bowtie2-build mature.human.fa  mature_human

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
