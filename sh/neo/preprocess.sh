#WES & RNA both 

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
---------------------------------------------------------------------------------------------------------------------------------------------
#RNA

##hisat2_tran.sh                                                               \
                                                                                
INDEX=/home/eunji/proj/0_sh/ref/rna/grch38_tran/genome_tran

cat config3 | while read id
do
        echo "start hisat2_tran for ${id}" `date`
        fq1=./2.clean_fq/${id}-R_1_val_1.fq.gz
        fq2=./2.clean_fq/${id}-R_2_val_2.fq.gz
        time hisat2 -p 4 -x ${INDEX} -1 ${fq1} -2 ${fq2} | samtools view -bSh >\
 ./4.align/hisat2_tran/${id}.bam
        echo "end hisat2_tran for ${id}" `date`
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
---------------------------------------------------------------------------------------------------------------------------------------------------------
#WES 

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

##MarhDuplicates.sh 
GATK=/home/eunji/miniconda3/envs/gatk/share/gatk4-4.1.9.0-0/gatk-package-4.1.9.0-local.jar
cat config  | while read id
do
                echo "start MarkDuplicates for ${id}" `date`
      java -jar $GATK  MarkDuplicates -I ./${id}.bam -O ./5.gatk/${id}_marked.bam -M ./5.gatk/${id}.metrics.txt 1>./5.gatk/${id}_log.mark 2>&1
                echo "end MarkDuplicates for ${id}" `date`
done

##index_marked_bam.sh

cat config | while read id

 do
  echo "start samtools index mark duplicates bam ${id}" `date`
  time samtools index -@ 10 -b ./5.gatk/${id}_marked.bam ./5.gatk/${id}_marked.bai
  echo "end samtools index mark duplicates bam ${id}" `date`
done

#BQSR.sh
GATK=/home/eunji/miniconda3/envs/gatk/share/gatk4-4.1.9.0-0/gatk-package-4.1.9.0-local.jar
snp=/home/eunji/ref/dbsnp_146.hg38.vcf.gz
indel=/home/eunji/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
ref=/home/eunji/ref/Homo_sapiens_assembly38.fasta
cat config  | while read id
do
                echo "start BQSR for ${id}" `date`
      java -jar $GATK BaseRecalibrator \
                -R $ref  \
                -I ./${id}_marked.bam  \
                --known-sites ${snp} \
                --known-sites ${indel} \
                -O ./${id}_recal.table \
                1>./${id}_log.recal 2>&1
      java -jar $GATK ApplyBQSR \
                -R $ref  \
                -I ./${id}_marked.bam  \
                -bqsr ./${id}_recal.table \
                -O ./${id}_bqsr.bam \
                1>./${id}_log.ApplyBQSR  2>&1

                echo "end BQSR for ${id}" `date`
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

#strelka - call somatic (v2.9.10) 
ref=/home/eunji/proj/0_sh/ref/wgs/Homo_sapiens_assembly38.fasta
snp=/home/eunji/proj/0_sh/ref/wgs/dbsnp_146.hg38.vcf.gz

cat config2 | while read id
do
        arr=(${id})
        sample=${arr[1]}
        T=./5.gatk/BQSR/${arr[1]}_bqsr.bam
        N=./5.gatk/BQSR/${arr[0]}_bqsr.bam
        configureStrelkaSomaticWorkflow.py \
        --tumorBam ${T}  \
        --normalBam ${N}  \
        --referenceFasta ${ref} \
        --runDir 6.strelka \
        --disableEVS \
        --reportEVSFeatures \
        --config=/home/eunji/miniconda3/envs/strelka/share/strelka-2.9.10-0/bin/configureStrelkaSomaticWorkflow.py.ini \
        --snvScoringModelFile=/home/eunji/miniconda3/envs/strelka/share/strelka-2.9.10-0/share/config/somaticSNVScoringModels.json \
        --indelScoringModelFile=/home/eunji/miniconda3/envs/strelka/share/strelka-2.9.10-0/share/config/somaticIndelScoringModels.json \
        --outputCallableRegions
done

#ANNOVAR - annotation 
#cd ./biosoft
# wget 下载地址
#tar -zxvf annovar.latest.tar.gz
#cd annovar
#nohup ./annotate_variation.pl -downdb -webfrom annovar gnomad_genome --buildver hg38 humandb/ >down.log 2>&1 &  ##41G	humandb/

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


