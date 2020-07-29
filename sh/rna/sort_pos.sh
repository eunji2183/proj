#sort_pos.sh
cat config | while read id 
do
echo "start samtools sort pos for ${id}" `date`
samtools sort -@ 10 ./4.align/hisat2_tran/${id}.bam -o ./4.align/hisat2_tran/${id}_sort.bam 
echo "end samtools sort pos for ${id}" `date`
done 
