perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20210501 humandb/

#annovar annotation 

cat config | while  read id
do
echo "start ANNOVAR for ${id} " `date`
~/tool/annovar/table_annovar.pl ./6.mutect/${id}_filter.vcf ~/tool/annovar/humandb/ \
-buildver hg38 \
-out ./7.annotation/annovar/${id} \
-remove \
-protocol refGene,knownGene,clinvar_20200316 \
-operation g,g,f \
-nastring . \
-vcfinput
echo "end ANNOVAR for ${id} " `date`
done

cat config | while read id
do 
	grep -v '^Chr' ./7.annotation/annovar/${id}.hg38_multianno.txt | cut -f 1-20 | awk -v T=${id} -v N=${id:0:5}_germline '{print $0"\t"T"\t"N}'  >./7.annotation/annovar/${id}.annovar.vcf 
done

head -1 ./7.annotation/annovar/${id}.hg38_multianno.txt| sed 's/Otherinfo/Tumor_Sample_Barcode\tMatched_Norm_Sample_Barcode/' >./7.annotation/annovar/header

cat ./7.annotation/annovar/header ./7.annotation/annovar/*annovar.vcf >./7.annotation/annovar/annovar_merge.vcf
