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


