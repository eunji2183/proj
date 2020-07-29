### transcripts count (isoform) 
##hisat2 - stringtie 
#hisat2_genome.sh 
INDEX=../../0_sh/ref/grch38/genome

cat config | while read id
do
        echo "start hisat2_genome for ${id}" `date`
        fq1=./2.clean_fq/${id}_1_val_1.fq.gz
        fq2=./2.clean_fq/${id}_2_val_2.fq.gz
        time hisat2 -p 4 -x ${INDEX} -1 ${fq1} -2 ${fq2} -S ./4.align/hisat2_genome/${id}.sam >>./4.align/hisat2_genome/${id}_hisat2.log 2>&1
        echo "end hisat2_genome for ${id}" `date`
done

##hisat2_genome view 
cat config | while read id 
do
echo "start genome samtools view for ${id}" `date`
samtools view -Sb ./4.align/hisat2_genome/${id}.sam > ./4.align/hisat2_genome/${id}.bam
echo "end genome samtools view for ${id}" `date`
done

##hisat2_genome sort pos
cat config | while read id 
do
echo "start genome samtools sort pos for ${id}" `date`
samtools sort -@ 10 ./4.align/hisat2_genome/${id}.bam -o ./4.align/hisat2_genome/${id}_sort.bam
echo "end genome samtools sort pos for ${id}" `date`
done

##stringtie.sh : assemble transcripts for each sample > merge 
GTF=../../0_sh/ref/Homo_sapiens.GRCh38.96.gtf

cat config | while read id 
do 
 echo "start stringtie for ${id}" `date`
  stringtie -p 8 -G ${GTF} -o ./6.isoform/stringtie/${id}.gtf -l ./6.isoform/stringtie/${id} ./4.align/hisat2_genome/${id}_sort.bam >> ./6.isoform/stringtie/${id}_stringtie.log 2>&1 
 echo "end stringtie for ${id}" `date`
done 

ls -l ./6.isoform/stringtie/*.gtf | awk '{print $9}' > sample_assembly_gtf_list.txt 

stringtie --merge -p 8 -o ./6.isoform/stringtie/stringtie_merge.gtf -G ${GTF} sample_assembly_gtf_list.txt

##reassemble to stringtie_merge.gtf 
GTF=./6.isoform/stringtie/stringtie_merge.gtf 

cat config | while read id 
do 
 echo " start second stringtie for ${id}" `date`
 stringtie -e -B -p 4 ./4.align/hisat2_genome/${id}_sort.bam -G ${GTF} -o ./6.isoform/stringtie2/${id}/${id}.gtf -A ./6.isoform/stringtie2/${id}/${id}_gene_abun.txt 
 echo "end second stringtie for ${id}" `date`
done 

#salmon.sh 
##cd ../../0_sh/ref
##bash generateDecoyTranscriptome.sh -a Homo_sapiens.GRCh38.96.gtf -g GRCh38.primary_assembly.genome.fa -t Homo_sapiens.GRCh38.cdna.all.fa -o GRCh38_salmon_index
##salmon index -t gentrome.fa -d decoys.txt -i ./salmon_index 
salmon_index=../../0_sh/ref/salmon_index 
cat config | while read id 
do 
  echo "start salmon quant for ${id}" `date`
       fq1=./2.clean_fq/${id}_1_val_1.fq.gz
       fq2=./2.clean_fq/${id}_2_val_2.fq.gz
    salmon quant -i ${salmon_index} -l A -1 ${fq1} -2 ${fq2} -o ./6.isoform/salmon/${id}_quant >> ./6.isoform/salmon/${id}_salmon.log 2>&1 
  echo "end salmon quant for ${id}" `date` 
done 
