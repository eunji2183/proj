#miRNA analysis                                                                 

#fastx_trimmer                                                                  

cat config | while read id
do
fastq_quality_filter -v -q 20 -p 80 -Q33  -i ./1.fq/${id}.fastq -o ./8.mirna/1.\
trim/tmp/${id}.fastq
fastx_trimmer -v -f 1 -l 27 -m 15  -i ./8.mirna/1.trim/tmp/${id}.fastq  -Q33 -z\
 -o ./8.mirna/1.trim/${id}_clean.fq.gz
done

mature=/home/eunji/proj/0_sh/ref/mirna/hsa-mature-bowtie-index
hairpin=/home/eunji/proj/0_sh/ref/mirna/hsa-hairpin-bowtie-index

cat config2 | while read id
do
bowtie -n 0 m1 --best --strata $mature -1 ./8.mirna/1.trim/${id}_1.clean.fq.gz \
-2 ./8.mirna/1.trim/${id}_2_clean.fq.gz -S ${id}_mature.sam
bowtie  -n 0 -m1 --best --strata $hairpin -1 ./8.mirna/1.trim/${id}_1.clean.fq.\
gz -2 ./8.mirna/1.trim/${id}_2_clean.fq.gz  -S ${id}_hairpin.sam
done

ls *.sam|while read id ;do (samtools sort -O bam -@ 5  -o $(basename ${id} ".sa\
\                                                                               
m").bam   ${id});done
rm *.sam

ls *.bam |xargs -i samtools index {}
ls *.bam|while read id ;do (samtools idxstats ${id} > ${id}.txt );done
# samtools view  matrue.bam |cut -f 3 |sort |uniq  -c 
