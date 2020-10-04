#raw_qc
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

##clean_qc.sh
cat config | while read id
do
	fastqc --outdir ./3.qc/clean_qc/ --threads 16 ./2.clean_fq/${id}*.fq.gz >> ./3.qc/clean_qc/${id}_fastqc.log 2>&1 
done 

multiqc  ./3.qc/clean_qc/*zip  -o ./3.qc/clean_qc/multiqc

rm ./3.qc/clean_qc/*zip 
rm ./3.qc/clean_qc/*.html
rm ./3.qc/clean_qc/*.log

## bwa.sh
INDEX=../../0_sh/ref/wgs/Homo_sapiens_assembly38

cat config | while read id
do
	echo "start bwa for ${id}" `date`
	fq1=./2.clean_fq/${id}_R1_P.fq.gz
	fq2=./2.clean_fq/${id}_R2_P.fq.gz
	time bwa mem -M -Y -t 4 -R "@RG\tID:${id}\tSM:${id}\tLB:WGS\tPL:ILLUMINA" -o ./4.align/${id}.sam ${INDEX} ${fq1} ${fq2} >> ./4.align/${id}_bwa.log 2>&1
	echo "end bwa for ${id}" `date`
done

#sort.sh
cat config | while read id 
do
echo "start samtools sort for ${id}" `date`
samtools sort -@ 10 ./4.align/${id}.sam -o ./4.align/${id}_sort.bam 
echo "end samtools sort for ${id}" `date`
done 

#rm -rf ./4.align/${id}.sam 

##MarhDuplicates.sh 

cat config  | while read id
do
                echo "start MarkDuplicates for ${id}" `date`
      gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" MarkDuplicates -I ./4.align/${id}_sort.bam --REMOVE_DUPLICATES=true -O ./5.gatk/Mdup/${id}_marked.bam -M ./5.gatk/Mdup/${id}.metrics.txt 1>./5.gatk/Mdup/${id}_log.mark 2>&1
                echo "end MarkDuplicates for ${id}" `date`
done

##index_marked_bam.sh

cat config | while read id

 do
  echo "start samtools index mark duplicates bam ${id}" `date`
  time samtools index -@ 10 -b ./5.gatk/Mdup/${id}_marked.bam ./5.gatk/Mdup/${id}_marked.bai
  echo "end samtools index mark duplicates bam ${id}" `date`
 done

#rm -rf ./4.align/${id}_sort.bam

##BQSR.sh (base quality score recalibration) 
snp=../../0_sh/ref/wgs/dbsnp_146.hg38.vcf.gz
indel=../../0_sh/ref/wgs/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
ref=../../0_sh/ref/wgs/Homo_sapiens_assembly38.fasta
hc=../../0_sh/ref/wgs/1000G_phase1.snps.high_confidence.hg38.vcf.gz

cat config | while read id
do
echo "start BQSR for ${id}" `date`
gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" BaseRecalibrator -R ${ref} -I ./5.gatk/Mdup/${id}_marked.bam  --known-sites ${hc} --known-sites ${snp} --known-sites ${indel} -O ./5.gatk/BQSR/${id}_recal.table 1>./5.gatk/BQSR/${id}_log.recal 2>&1

gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" ApplyBQSR -R ${ref} -I ./5.gatk/Mdup/${id}_marked.bam -bqsr ./5.gatk/BQSR/${id}_recal.table -O ./5.gatk/BQSR/${id}_bqsr.bam 1>./5.gatk/BQSR/${id}_log.ApplyBQSR 2>&1

samtools index -@ 10 -b ./5.gatk/BQSR/${id}_bqsr.bam ./5.gatk/BQSR/${id}_bqsr.bai

echo "end BQSR for ${id}" `date`
done

#rm -rf ./5.gatk/BQSR/${id}_bqsr.bam 

##haplotypecaller.sh 
ref=../../0_sh/ref/wgs/Homo_sapiens_assembly38.fasta

cat config | while read id
do 
echo "start haplotypecaller for ${id}" `date`
gatk HaplotypeCaller --emit-ref-confidence GVCF -R ${ref} -I ./5.gatk/BQSR/${id}_bqsr.bam -O ./5.gatk/GVCF/${id}_HC.g.vcf.gz 1>./5.gatk/GVCF/${id}_log.HC 2>&1

gatk GenotypeGVCFs -R ${ref} -V ./5.gatk/GVCF/${id}_HC.g.vcf.gz -O ./5.gatk/VCF/${id}_HC.vcf.gz 1>./5.gatk/VCF/${id}_log.GG 2>&1

echo "end haplotypecaller for ${id}" `date`
done

##VQSR - SNP&INDEL 
ref=../../0_sh/ref/wgs/Homo_sapiens_assembly38.fasta
snp=../../0_sh/ref/wgs/dbsnp_146.hg38.vcf.gz
indel=../../0_sh/ref/wgs/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
hc=../../0_sh/ref/wgs/1000G_phase1.snps.high_confidence.hg38.vcf.gz
omini=../../0_sh/ref/wgs/1000G_omni2.5.hg38.vcf.gz







##mutect2.sh (somatic mutation)
##config2 - normal & tumor matching 
ref=../../0_sh/ref/wgs/Homo_sapiens_assembly38.fasta
list=../../0_sh/ref/wgs/hg38.interval_list

cat config2 | while read id
do
        arr=(${id})
        sample=${arr[1]}
        T=./5.gatk/BQSR/${arr[1]}_bqsr.bam
        N=./5.gatk/BQSR/${arr[0]}_bqsr.bam
        echo "start Mutect2 for ${id}" `date`
        gatk --java-options "-Xmx20G -Djava.io.tmpdir=./"  Mutect2 -R ${ref} \
        -I ${T} -tumor  $(basename "$T" _bqsr.bam) \
        -I ${N} -normal $(basename "$N" _bqsr.bam) \
        -L ${list}  \
        -O ./6.mutect/${sample}_mutect2.vcf

        gatk FilterMutectCalls \
  -R ${ref} \
        -V ./6.mutect/${sample}_mutect2.vcf \
        -O ./6.mutect/${sample}_somatic.vcf
        echo "end Mutect2 for ${id}" `date`

        cat ./6.mutect/${sample}_somatic.vcf | perl -alne '{if(/^#/){print}else{next unless $F[6] eq "PASS";next if $F[0] =~/_/;print } }' > ./6.mutect/${sample}_filter.vcf
done
