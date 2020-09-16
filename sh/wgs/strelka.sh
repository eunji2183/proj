#strelka - call somatic (v2.9.10) 
ref=/home/eunji/proj/0_sh/ref/wgs/Homo_sapiens_assembly38.fasta
snp=/home/eunji/proj/0_sh/ref/wgs/dbsnp_146.hg38.vcf.gz

cat config2 | while read id
do
        arr=(${id})
        sample=${arr[1]}
        T=./5.gatk/BQSR/${arr[1]}_bqsr.bam
        N=./5.gatk/BQSR/${arr[0]}_bqsr.bam
        configureStrelkaSomaticWorkflow.py \
        --tumorBam ${T}  \
        --normalBam ${N}  \
        --referenceFasta ${ref} \
        --runDir 6.strelka \
        --disableEVS \
        --reportEVSFeatures \
        --config=/home/eunji/miniconda3/envs/strelka/share/strelka-2.9.10-0/bin/configureStrelkaSomaticWorkflow.py.ini \
        --snvScoringModelFile=/home/eunji/miniconda3/envs/strelka/share/strelka-2.9.10-0/share/config/somaticSNVScoringModels.json \
        --indelScoringModelFile=/home/eunji/miniconda3/envs/strelka/share/strelka-2.9.10-0/share/config/somaticIndelScoringModels.json \
        --outputCallableRegions
done

