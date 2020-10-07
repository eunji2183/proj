#rna raw file - trim-galore 
#dna file - wes analysis (ref-hg38)

#neopepsee.sh 
#neopepsee-3.0.0 안에 rna raw fastq file , vcf file , sh file 넣음 

java -jar ../neopepsee-3.0.0.jar \
-r ../database/ref/hg38/knownGene_hg38_v3.fasta \
-1 15190-MLS-R_1_val_1.fq \
-2 15190-MLS-R_2_val_2.fq \
-e 0 \
-v 15190-MLS-D.final.vcf \
-d . \
-n 9
