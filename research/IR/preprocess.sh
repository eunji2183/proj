cd /HDD8T/eunji/proj/IR/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/
STAR --runMode genomeGenerate --genomeDir /HDD8T/eunji/proj/IR/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/ --genomeFastaFiles /HDD8T/eunji/proj/IR/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa --sjdbGTFfile /HDD8T/eunji/proj/IR/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf --sjdbOverhang 75 --runThreadN 8

cd /HDD8T/eunji/proj/IR/STAR/

#STAR 2-pass 
cat config | while read id 
do
STAR \
--runThreadN 6 \
--genomeDir /HDD8T/eunji/proj/IR/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/ \
--twopassMode Basic \
--readFilesIn /HDD8T/eunji/proj/IR/STAR/trimmomatic_PE_genv33/${id}.fq.gz \
--readFilesCommand zcat \
--outSAMtype BAM Unsorted \
--chimOutType SeparateSAMold \
--quantMode GeneCounts \
--outFileNamePrefix /HDD8T/eunji/proj/IR/STAR/trimmomatic_PE_genv33/${id}
 done
#picard AddOrReplaceReadGroups
cat config | while read id 
do
picard AddOrReplaceReadGroups \
I=/HDD8T/eunji/proj/IR/STAR/trimmomatic_PE_genv33/${id}Aligned.out.bam \
O=/HDD8T/eunji/proj/IR/STAR/trimmomatic_PE_genv33/${id}_sorted.bam \
SO=coordinate RGID=${id} RGLB=mRNA RGPL=illumina RGPU=HiSeq2500 RGSM=${id}
done

#picard MarkDuplicates
cat config | while read id 
do
picard MarkDuplicates \
I=/HDD8T/eunji/proj/IR/STAR/trimmomatic_PE_genv33/${id}_sorted.bam \
O=/HDD8T/eunji/proj/IR/STAR/trimmomatic_PE_genv33/${id}_sorted_maked.bam \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT \
METRICS_FILE=/HDD8T/eunji/proj/IR/STAR/trimmomatic_PE_genv33/${id}.metrics
done

#gatk SplitNCigarReads AddOrReplaceReadGroups
cat config | while read id
do
gatk SplitNCigarReads \
-R /HDD8T/eunji/proj/IR/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa \
-I /HDD8T/eunji/proj/IR/STAR/trimmomatic_PE_genv33/${id}_sorted_maked.bam \
-O /HDD8T/eunji/proj/IR/STAR/trimmomatic_PE_genv33/${id}_sorted_maked_split.bam
done

cat config | while read id 
do
gatk AddOrReplaceReadGroups \
-I /HDD8T/eunji/proj/IR/STAR/trimmomatic_PE_genv33/${id}_sorted_maked_split.bam \
-O /HDD8T/eunji/proj/IR/STAR/trimmomatic_PE_genv33/${id}_split_add.bam \
-LB ${id} -PL ILLUMINA -PU ${id} -SM ${id}
done


#gatk BaseRecalibrator ApplyBQSR
GENOME=/HDD8T/eunji/proj/IR/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa
DBSNP=/HDD8T/eunji/proj/IR/ref/dbsnp_146.hg38.vcf.gz
kgSNP=/HDD8T/eunji/proj/IR/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz
kgINDEL=/HDD8T/eunji/proj/IR/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
cat config | while read id
do
gatk --java-options "-Xmx20G -Djava.io.tmpdir=./tmp" BaseRecalibrator \
 -I /HDD8T/eunji/proj/IR/STAR/trimmomatic_PE_genv33/${id}_sorted_maked_split.bam \
 -R $GENOME \
 -O /HDD8T/eunji/proj/IR/STAR/trimmomatic_PE_genv33/${id}_recal.table --known-sites $kgSNP --known-sites $kgINDEL

 gatk --java-options "-Xmx20G -Djava.io.tmpdir=./tmp" ApplyBQSR \
 -I /HDD8T/eunji/proj/IR/STAR/trimmomatic_PE_genv33/${id}_sorted_maked_split.bam \
 -R $GENOME \
 --output /HDD8T/eunji/proj/IR/STAR/trimmomatic_PE_genv33/${id}_recal.bam \
 -bqsr /HDD8T/eunji/proj/IR/STAR/trimmomatic_PE_genv33/${id}_recal.table
 done

GENOME=/HDD8T/eunji/proj/IR/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa
DBSNP=/HDD8T/eunji/proj/IR/ref/dbsnp_146.hg38.vcf.gz
cat config | while read id
do
gatk --java-options "-Xmx20G -Djava.io.tmpdir=./tmp" HaplotypeCaller \
-R $GENOME \
-I /HDD8T/eunji/proj/IR/STAR/trimmomatic_PE_genv33/${id}_recal.bam \
-O /HDD8T/eunji/proj/IR/STAR/trimmomatic_PE_genv33/${id}.vcf \
--dbsnp $DBSNP \
--dont-use-soft-clipped-bases \
-stand-call-conf 20.0
done

cat config | while read id 
do
gatk --java-options "-Xmx20G -Djava.io.tmpdir=./tmp" VariantFiltration \
-R $GENOME \
-V /HDD8T/eunji/proj/IR/STAR/trimmomatic_PE_genv33/${id}.vcf \
-window 35 \
-cluster 3 \
--filter-name FS -filter "FS > 30.0" \
--filter-name QD -filter "QD < 2.0" \
-O /HDD8T/eunji/proj/IR/STAR/trimmomatic_PE_genv33/${id}_filtered.vcf
done

