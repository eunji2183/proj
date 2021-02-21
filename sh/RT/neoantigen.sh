
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
