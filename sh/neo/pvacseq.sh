#
conda create -n pvac #(python3)
conda create -n pvac2 #(python2) 
conda install python=3
pip install pvactools
pip install pvactools --upgrade
#https://pvactools.readthedocs.io/en/latest/install.html

##annotation tool - VEP 
#https://m.blog.naver.com/PostView.nhn?blogId=glshh13&logNo=221565564066&proxyReferer=https:%2F%2Fwww.google.com%2F -> VEP install
conda create -n vep 
conda activate vep 
conda install -c bioconda ensembl-vep
mkdir ~/.vep
cd $HOME/.vep
wget ftp://ftp.ensembl.org/pub/release-101/variation/indexed_vep_cache/homo_sapiens_vep_101_GRCh38.tar.gz
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
git clone https://github.com/Ensembl/VEP_plugins.git

#vep.sh
cat config | while read id
do
vep --input_file ./6.mutect/${id}_filter.vcf --output_file ./7.annotation/vep/${id}_vep.vcf --format vcf --vcf 
--symbol --terms SO --tsl --hgvs --fasta /home/eunji/proj/0_sh/ref/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa 
--offline --cache $HOME/.vep --plugin Downstream --plugin Wildtype --dir_plugins /home/eunji/proj/0_sh/biosoft/VEP_plugins/ --pick --transcript_version
done

conda install vt 
pip install vcf-annotation-tools

#add_expression_vcf.sh (somatic)
cat config | while read id 
do
vcf-expression-annotator -s ${id} 
-o ../../8.pvacseq/vcf/${id}.tx.vcf 
--ignore-transcript-version 
${id}_vep.vcf 
/home/eunji/proj/RMLS/rna/6.isoform/kallisto/${id}/abundance.tsv kallisto transcript 
done 

#add_expression_vcf.sh (gvcf)
cat config | while read id
do
vcf-expression-annotator 
-o ./8.pvacseq/gvcf/${id}.tx.vcf 
--ignore-transcript-version 
./7.annotation/vep_gvcf/${id}_gvep.vcf 
/home/eunji/proj/RMLS/rna/6.isoform/kallisto/${id}/abundance.tsv kallisto transcript
done


conda create -n bam-read
conda activate bam-read 
conda install -c bioconda bam-readcount








