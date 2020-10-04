## WES
mkdir biosoft project data
cd project
mkdir -p 0.sra 1.raw_fq 2.clean_fq 3.qc/{raw_qc,clean_qc} 4.align/{qualimap,flagstat,stats} 5.gatk/gvcf 6.mutect 7.annotation/{vep,annovar,funcotator,snpeff} 8.cnv/{gatk,cnvkit,gistic,facet} 9.pyclone 10.signature

#download sra file 
cat SRR_Acc_List.txt | while read id
do
	prefetch ${id} -O  ./0.sra
done

cd ./0.sra
## fasterq-dump
cat config | while read id
do
	time fasterq-dump -3 -e 16 ${id}.sra -O ../1.raw_fq --outfile ${id}.fastq  >> ../1.raw_fq/${id}_sra2fq.log 2>&1
	time pigz -p 10 -f ../1.raw_fq/${id}_1.fastq >../1.raw_fq/${id}_1.fastq.gz
	time pigz -p 10 -f ../1.raw_fq/${id}_2.fastq >../1.raw_fq/${id}_2.fastq.gz
done 

# raw fastqc
cat config | while read id
do
	fastqc --outdir ./3.qc/raw_qc/ --threads 16 ./1.raw_fq/${id}*.fastq.gz >> ./3.qc/raw_qc/${id}_fastqc.log 2>&1 
done 

multiqc  ./3.qc/raw_qc/*zip  -o ./3.qc/raw_qc/multiqc

## trim_galore.sh
cat config | while read id
do
	fq1=./1.raw_fq/${id}_1.fastq.gz
	fq2=./1.raw_fq/${id}_2.fastq.gz
	trim_galore  --paired -q 28 --phred33 --length 30 --stringency 3 --gzip --cores 8 -o ./2.clean_fq  $fq1  $fq2 >> ./2.clean_fq/${id}_trim.log 2>&1
done

#trim fastqc
cat config | while read id
do
	fastqc --outdir ./3.qc/clean_qc/ --threads 16 ./2.clean_fq/${id}*.fq.gz >> ./3.qc/clean_qc/${id}_fastqc.log 2>&1 
done 

multiqc  ./3.qc/clean_qc/*zip  -o ./3.qc/clean_qc/multiqc

#cd ./data
#gunzip Homo_sapiens_assembly38.fasta.gz
#time bwa index -a bwtsw -p gatk_hg38 ~/wes_cancer/data/Homo_sapiens_assembly38.fasta
#cd ../project

## bwa.sh
INDEX=./data/gatk_hg38

cat config | while read id
do
	echo "start bwa for ${id}" `date`
	fq1=./2.clean_fq/${id}_1_val_1.fq.gz
	fq2=./2.clean_fq/${id}_2_val_2.fq.gz
	bwa mem -M -t 16 -R "@RG\tID:${id}\tSM:${id}\tLB:WXS\tPL:Illumina" ${INDEX} ${fq1} ${fq2} | samtools sort -@ 10 -m 1G  -o  ./4.align/${id}.bam -
	echo "end bwa for ${id}" `date`
done

#GATK - MarkDuplicates
GATK=./biosoft/gatk-4.1.4.1/gatk

cat config  | while read id
do
	BAM=./4.align/${id}.bam
	if [ ! -f ./5.gatk/ok.${id}_marked.status ]
	then
		echo "start MarkDuplicates for ${id}" `date`
		$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" MarkDuplicates \
		-I ${BAM} \
		--REMOVE_DUPLICATES=true \
		-O ./5.gatk/${id}_marked.bam \
		-M ./5.gatk/${id}.metrics \
		1>./5.gatk/${id}_log.mark 2>&1 
		
		if [ $? -eq 0 ]
		then
			touch ./5.gatk/ok.${id}_marked.status
		fi
		echo "end MarkDuplicates for ${id}" `date`
		samtools index -@ 16 -m 4G -b ./5.gatk/${id}_marked.bam ./5.gatk/${id}_marked.bai
	fi
done

#BQSR
GATK=./biosoft/gatk-4.1.4.1/gatk
snp=./data/dbsnp_146.hg38.vcf.gz
indel=./data/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
ref=./data/Homo_sapiens_assembly38.fasta
cat config  | while read id
do
	if [ ! -f ./5.gatk/${id}_bqsr.bam ]
	then
		echo "start BQSR for ${id}" `date`
		$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./"  BaseRecalibrator \
		-R $ref  \
		-I ./5.gatk/${id}_marked.bam  \
		--known-sites ${snp} \
		--known-sites ${indel} \
		-O ./5.gatk/${id}_recal.table \
		1>./5.gatk/${id}_log.recal 2>&1 
		
		$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./"  ApplyBQSR \
		-R $ref  \
		-I ./5.gatk/${id}_marked.bam  \
		-bqsr ./5.gatk/${id}_recal.table \
		-O ./5.gatk/${id}_bqsr.bam \
		1>./5.gatk/${id}_log.ApplyBQSR  2>&1 
		
		echo "end BQSR for ${id}" `date`
	fi
done

#Germline SNPs + Indels
GATK=./biosoft/gatk-4.1.4.1/gatk
snp==./data/dbsnp_146.hg38.vcf.gz
indel=./data/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
ref=./data/Homo_sapiens_assembly38.fasta
bed=./data/hg38.exon.bed

cat config  | while read id
do
	echo "start HC for ${id}" `date`
	$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" HaplotypeCaller -ERC GVCF \
	-R ${ref} \
	-I ./5.gatk/${id}_bqsr.bam \
	--dbsnp ${snp} \
	-L ${bed} \
	-O ./5.gatk/${id}_raw.vcf \
	1>./5.gatk/${id}_log.HC 2>&1
	echo "end HC for ${id}" `date`

done

cd ./5.gatk/gvcf
for chr in chr{1..22} chrX chrY chrM
do

time $GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" GenomicsDBImport \
-R ${ref} \
$(ls ./*raw.vcf | awk '{print "-V "$0" "}') \
-L ${chr} \
--genomicsdb-workspace-path gvcfs_${chr}.db

time $GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" GenotypeGVCFs \
-R ${ref} \
-V gendb://gvcfs_${chr}.db \
-O gvcfs_${chr}.vcf

done

$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" GatherVcfs \
$(for i in {1..22} X Y M;do echo "-I gvcfs_chr${i}.vcf" ;done) \
-O merge.vcf

## mutect.sh - somatic mutation 
GATK=./biosoft/gatk-4.1.4.1/gatk
ref=./data/Homo_sapiens_assembly38.fasta
bed=./data/hg38.exon.bed

cat config2 | while read id
do
	arr=(${id})
	sample=${arr[1]}
	T=./5.gatk/${arr[1]}_bqsr.bam
	N=./5.gatk/${arr[0]}_bqsr.bam
	echo "start Mutect2 for ${id}" `date`
	$GATK  --java-options "-Xmx20G -Djava.io.tmpdir=./"  Mutect2 -R ${ref} \
	-I ${T} -tumor  $(basename "$T" _bqsr.bam) \
	-I ${N} -normal $(basename "$N" _bqsr.bam) \
	-L ${bed}  \
	-O ./6.mutect/${sample}_mutect2.vcf

	$GATK  FilterMutectCalls \
  -R ${ref} \
	-V ./6.mutect/${sample}_mutect2.vcf \
	-O ./6.mutect/${sample}_somatic.vcf
	echo "end Mutect2 for ${id}" `date`

	cat ./6.mutect/${sample}_somatic.vcf | perl -alne '{if(/^#/){print}else{next unless $F[6] eq "PASS";next if $F[0] =~/_/;print } }' > ./6.mutect/${sample}_filter.vcf
done

#ANNOVAR - annotation 
cat config | while  read id
do
echo "start ANNOVAR for ${id} " `date`
~/biosoft/annovar/table_annovar.pl ./6.mutect/${id}_filter.vcf ~/biosoft/annovar/humandb/ \
-buildver hg38 \
-out ./7.annotation/annovar/${id} \
-remove \
-protocol refGene,knownGene,clinvar_20170905 \
-operation g,g,f \
-nastring . \
-vcfinput
echo "end ANNOVAR for ${id} " `date`
done
