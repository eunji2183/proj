
#sra_down.sh (sra-tools)
cat SRR_Acc_List.txt | while read id 
do
prefetch --max-size 100G ${id} -O ./sra 
done


#sra2fq.sh
cat SRR_Acc_List.txt | while read id
do
	time fasterq-dump -3 -e 16 ./sra/${id}/${id}.sra -O ./raw_fq --outfile ${id}.fastq  >> ./raw_fq/${id}_sra2fq.log 2>&1
	time pigz -p 10 -f ./raw_fq/${id}_1.fastq >./raw_fq/${id}_1.fastq.gz
	time pigz -p 10 -f ./raw_fq/${id}_2.fastq >./raw_fq/${id}_2.fastq.gz
done 

gunzip *.gz 
samtools fqidx *.fastq

barcode_splitter --bcfile ../barcode.txt ../raw_fq/ATM_UN2_2.fastq.gz --gzipin --gzipout --idxread 1 --suffix .fastq

## fastp trim 
cat config | while read id
do
        fastp -i ./raw_fq/${id}_1.fastq.gz -o ./trim_fq/${id}_1.fastq.gz -I ./raw_fq/${id}_2.fastq.gz -O ./trim_fq/${id}_2.fastq.gz --thread=16
done

##  trim_galore
cat config | while read id
do
        trim_galore --illumina -o ./trim_fq2 ./raw_fq/${id}_1.fastq.gz ./raw_fq/${id}_2.fastq.gz >> ./trim_fq2/${id}_trim.log 2>&1
done

#trimmomatic
cat config | while read id 
do
trimmomatic PE -phred33 ./raw_fq/${id}_1.fastq.gz ./raw_fq/${id}_2.fastq.gz \
./trimmomatic/${id}_R1_P.fq.gz ./trimmomatic/${id}_R1_U.fq.gz ./trimmomatic/${id}_R2_P.fq.gz ./trimmomatic/${id}_R2_U.fq.gz \
ILLUMINACLIP:/home/eunji/miniconda3/envs/IR/share/trimmomatic-0.39-1/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
