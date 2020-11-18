##fusioncatcher.sh 
cat config | while read id 
do
echo "start fusioncatcher for ${id}" `date`
fusioncatcher.py -d /home/eunji/miniconda3/envs/python2/share/fusioncatcher-1.20/db/human_v98/ -i ./1.raw_fq/${id}/ -o ./7.fusion/${id}/ >> ./7.fusion/${id}_fusioncatcher.log 2>&1
echo "end fusioncatcher for ${id}" `date`
done 

#fusion만 뽑아내기 
grep "\*" summary_candidate_fusions.txt | awk '{printf ("%s\n",$2)}' > 15190-MLS.fusion.txt  
