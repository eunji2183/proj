##trinity.sh - isoform 
cat config | while read id
do 
cd ./6.isoform/trinity/${id}/
echo "start trinity for ${id}" `date`
fq1=./2.clean_fq/${id}_1_val_1.fq.gz
fq2=./2.clean_fq/${id}_2_val_2.fq.gz
Trinity --seqType fq --left ${fq1} --right ${fq2} --CPU 6 --max_memory 50G >> ./6.isoform/trinity/${id}/${id}_trinity.log 2>&1
echo "end trinity for ${id}" `date` 
done 
