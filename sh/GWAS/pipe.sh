
# binomial phenotype 

# ls | cut -d"_" -f 1-2 > config
# cat config | while read id; do cp ../bfile/200624_1_HapMap_8840.bed ${id}.bed; done
# cat config | while read id; do cp ../bfile/200624_1_HapMap_8840.bim ${id}.bim; done

#cat.sh

### QC 
## step 1 : --geno (filter geno) --mind (filter SNP) cutoff (0.02) 
 cat config | while read id
 do
 mkdir ${id}

 plink --bfile ${id} --missing --no-sex --out ./${id}/${id}_2
 plink --bfile ${id} --geno 0.02 --mind 0.02 --no-sex --make-bed --out ./${id}/${id}_5

 awk '{ if($1 >=1 && $1 <= 22) print $2}' ./${id}/${id}_5.bim > ./${id}/${id}snp_1_22.txt
 plink --bfile ./${id}/${id}_5 --extract ./${id}/${id}snp_1_22.txt --allow-no-sex --make-bed --out ./${id}/${id}_6

## step2 : MAF < 0.05 
 plink --bfile ./${id}/${id}_6 --freq --allow-no-sex --out ./${id}/${id}MAF_check
 plink --bfile ./${id}/${id}_6 --maf 0.05 --allow-no-sex --make-bed --out ./${id}/${id}_7

## step 3 : HWE (hardy-weinberg) 
 plink --bfile ./${id}/${id}_7 --hardy --allow-no-sex --out ./${id}/${id}_7
 awk '{ if ($9 <0.00001) print $0 }' ./${id}/${id}_7.hwe > ./${id}/${id}zoomhwe.hwe
 plink --bfile ./${id}/${id}_7 --hwe 1e-6 --allow-no-sex --make-bed --out ./${id}/${id}_hwe_filter_step1
 plink --bfile ./${id}/${id}_hwe_filter_step1 --hwe 1e-10 --hwe-all --allow-no-sex --make-bed --out ./${id}/${id}_8

## step 4 : heterogeneity 
 plink --bfile ./${id}/${id}_8 --het --allow-no-sex --out R_check
 Rscript het.R R_check.het
 sed 's/"//g' fail-het-qc.txt |awk '{print $1,$2}' > ./${id}/het_fail_ind.txt
 plink --bfile ./${id}/${id}_8 --remove ./${id}/het_fail_ind.txt --allow-no-sex --make-bed --out ./${id}/${id}_9

## step 5 : pihat>0.2
 plink --bfile ./${id}/${id}_9 --genome --min 0.2 --allow-no-sex --out ./${id}/${id}pihat_min0.2
 awk '{if($8>0.9) print $0}' ./${id}/${id}pihat_min0.2.genome  > ./${id}/${id}zoom_pihat.genome
 plink --bfile ./${id}/${id}_9 --filter-founders --allow-no-sex --make-bed --out ./${id}/${id}_10
 done

#assoc.sh
##association 
cat config2 | while read id
do
plink --bfile ./${id}/${id}_10 --assoc --allow-no-sex --out ./${id}/${id}_assoc_results
plink --bfile ./${id}/${id}_10 -assoc --adjust --allow-no-sex --out ./${id}/${id}_adjusted_assoc_results

cd ${id}
Rscript ../PRSice.R --dir ../ --prsice ../PRSice_linux --base ${id}_adjusted_assoc_results.assoc --target ${id}_10 --thread 1 --stat OR --binary-target T

cd ..
done

-------------------------------------------------------------------------------------------------------------------------------------------------------------

#continuous 
# continuous phenotype 

# ls | cut -d"_" -f 1-2 > config
# cat config | while read id; do cp ../bfile/200624_1_HapMap_8840.bed ${id}.bed; done
# cat config | while read id; do cp ../bfile/200624_1_HapMap_8840.bim ${id}.bim; done

### QC 
## step 1 : --geno (filter geno) --mind (filter SNP) cutoff (0.02) 
 cat config | while read id
 do
 mkdir ${id}

 plink --bfile ${id} --missing --no-sex --out ./${id}/${id}_2
 plink --bfile ${id} --geno 0.02 --mind 0.02 --no-sex --make-bed --out ./${id}/${id}_5

 awk '{ if($1 >=1 && $1 <= 22) print $2}' ./${id}/${id}_5.bim > ./${id}/${id}snp_1_22.txt
 plink --bfile ./${id}/${id}_5 --extract ./${id}/${id}snp_1_22.txt --allow-no-sex --make-bed --out ./${id}/${id}_6

## step2 : MAF < 0.05 
 plink --bfile ./${id}/${id}_6 --freq --allow-no-sex --out ./${id}/${id}MAF_check
 plink --bfile ./${id}/${id}_6 --maf 0.05 --allow-no-sex --make-bed --out ./${id}/${id}_7

## step 3 : HWE (hardy-weinberg) 
 plink --bfile ./${id}/${id}_7 --hardy --allow-no-sex --out ./${id}/${id}_7
 awk '{ if ($9 <0.00001) print $0 }' ./${id}/${id}_7.hwe > ./${id}/${id}zoomhwe.hwe
 plink --bfile ./${id}/${id}_7 --hwe 1e-6 --allow-no-sex --make-bed --out ./${id}/${id}_hwe_filter_step1
 plink --bfile ./${id}/${id}_hwe_filter_step1 --hwe 1e-10 --hwe-all --allow-no-sex --make-bed --out ./${id}/${id}_8

## step 4 : heterogeneity 
 plink --bfile ./${id}/${id}_8 --het --allow-no-sex --out R_check
 Rscript het.R R_check.het
 sed 's/"//g' fail-het-qc.txt |awk '{print $1,$2}' > ./${id}/het_fail_ind.txt
 plink --bfile ./${id}/${id}_8 --remove ./${id}/het_fail_ind.txt --allow-no-sex --make-bed --out ./${id}/${id}_9

## step 5 : pihat>0.2
 plink --bfile ./${id}/${id}_9 --genome --min 0.2 --allow-no-sex --out ./${id}/${id}pihat_min0.2
 awk '{if($8>0.9) print $0}' ./${id}/${id}pihat_min0.2.genome  > ./${id}/${id}zoom_pihat.genome
 plink --bfile ./${id}/${id}_9 --filter-founders --allow-no-sex --make-bed --out ./${id}/${id}_10
 done

##association 
cat config2 | while read id
do
plink --bfile ./${id}/${id}_10 --linear --allow-no-sex --out ./${id}/${id}_assoc_results
plink --bfile ./${id}/${id}_10 --linear --adjust --allow-no-sex --out ./${id}/${id}_adjusted_assoc_results
done

cat config3 | while read id
do
cd ${id}
Rscript ../PRSice.R --dir ../ --prsice ../PRSice_linux --base ${id}_adjusted_assoc_results.assoc.linear --target ${id}_10 --thread 1 --stat BETA --binary-target F
cd ..
done

-----------------------------------------------------------------------------------------------------------------------------------------------------------------
#multinomial (trinculo) 
cat config | while read id
do
mkdir ${id}
./trinculo multinom --bfile ${id} --pheno ${id}.txt --phenoname ${id} --basepheno 0 --missingpheno -9 --out ./${id}/${id}
done



