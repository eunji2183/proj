##Integrate-NEO (GRCh38.99 version 사용) 

#tophat mapping 
conda create -n tophat2
conda install tophat 
wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz (GRCh38.99.fa 로 고침)
bowtie2-build GRCh38.99.fa index
tophat -G Homo_sapiens.GRCh38.99.gtf --transcriptome-index=GRCh38.99.tr GRCh38.99 (dir GRCh38.99.tr 안에 Homo_sapiens.GRCh38.99 로 생김) 
#tophat.sh 
cat config | while read id
do
tophat -o /HDD2T/eunji/th/${id} --transcriptome-index=/home/eunji/ref/GRCh38.99.tr/Homo_sapiens.GRCh38.99 
/home/eunji/ref/GRCh38.99.tr/GRCh38.99 
/HDD2T/jeeh9/RMLS_RNA/trimmed_fq/Paired/${id}-R_R1_P.fq.gz 
/HDD2T/jeeh9/RMLS_RNA/trimmed_fq/Paired/${id}-R_R2_P.fq.gz
done

conda create -n integrate-neo
conda activate integrate-neo
conda install bwa python=2 gcc_linux-64 ucsc-gtftogenepred samtools bedtools matplotlib cmake 
wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz
gtfToGenePred -genePredExt -geneNameAsName2 Homo_sapiens.GRCh38.99.gtf Homo_sapiens.GRCh38.99.genePred
cut -f 1-10,12 Homo_sapiens.GRCh38.99.genePred > tmp.txt
echo -e "#GRCh38.ensGene.name\tGRCh38.ensGene.chrom\tGRCh38.ensGene.strand\tGRCh38.ensGene.txStart\tGRCh38.ensGene.txEnd\tGRCh38.ensGene.cdsStart\tGRCh38.ensGene.cdsEnd\tGRCh38.ensGene.exonCount\tGRCh38.ensGene.exonStarts\tGRCh38.ensGene.exonEnds\tGRCh38.ensemblToGeneName.value" > Homo_sapiens.GRCh38.99.tsv
cat tmp.txt >> Homo_sapiens.GRCh38.99.tsv 


git clone https://github.com/ChrisMaherLab/INTEGRATE-Vis.git
cd INTEGRATE-Vis.1.0.0
chmod +x install.sh
./install.sh -o /home/eunji/proj/0_sh/biosoft/integrate-vis
cd ~/proj/0_sh/biosoft/INTEGRATE-Vis/INTEGRATE-Vis.1.0.0/data/reference_genome/
wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{1..22}.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz
gunzip -c Homo_sapiens.GRCh38.dna.chromosome.* > GRCh38_r99.all.fa


cd INTEGRATE-Neo-V-1.2.1
chmod +x install.sh
./install.sh -o /home/eunji/proj/0_sh/biosoft/integrate-neo

wget https://sourceforge.net/projects/integrate-fusion/files/INTEGRATE.0.2.6.tar.gz
cd INTEGRATE_X_X
mkdir INTEGRATE-build
cd INTEGRATE-build
cmake ../Integrate/ -DCMAKE_BUILD_TYPE=release
make
cd ~/proj/0_sh/biosoft/INTEGRATE_0_2_6/INTERGRATE-build/bin
mkdir ./bwts 
./Integrate mkbwt /home/eunji/proj/0_sh/biosoft/INTEGRATE-Vis/INTEGRATE-Vis.1.0.0/data/reference_genome/GRCh38_r99.all.fa 

#integrate fusion Integrate fusion ref.fa annot.txt bwts accepted_hits.bam unmappeds.bam
./Integrate fusion /home/eunji/ref/GRCh38.99.tr/GRCh38_r99.all.fa \
/home/eunji/ref/GRCh38.99.tr/Homo_sapiens.GRCh38.99.tsv \
./bwts/ \
/HDD2T/eunji/th/${id}/accepted_hits.bam /HDD2T/eunji/th/${id}/unmappeds.bam
samtools index *.bam 

fusionBedpeAnnotator 
-r /home/eunji/ref/GRCh38.99.tr/GRCh38_r99.all.fa 
-g /home/eunji/ref/GRCh38.99.tr/Homo_sapiens.GRCh38.99.genePred 
-d /home/eunji/fusioncatcher/${id}.txt 
-i  /home/eunji/rna/70615-MLS/fusions.bedpe
-o ./fusions.annot.bedpe

#Hlaminer
bwa mem -a ../database_bwamem/HLA_ABC_CDS.fasta ${id}-R_1_val_1.fq.gz ${id}-R_2_val_2.fq.gz > ${id}.sam
perl ../bin/HLAminer.pl -a ${id}.sam -h ../database/HLA_ABC_CDS.fasta -s 500


#integrate-neo.py 
python integrate-neo.py 
-1 /home/eunji/tool/HLAminer-1.4/HLAminer_v1.4/70615-MLS/70615-MLS-R_1_val_1.fq.gz 
-2 /home/eunji/tool/HLAminer-1.4/HLAminer_v1.4/70615-MLS/70615-MLS-R_2_val_2.fq.gz 
-f /home/eunji/rna/70615-MLS/fusions.bedpe 
-r /home/eunji/ref/GRCh38.99.tr/GRCh38_r99.all.fa 
-g /home/eunji/ref/GRCh38.99.tr/Homo_sapiens.GRCh38.99.genePred -k

/home/eunji/tool/integrate/integrate-neo.py -t HLAminer_HPRA.csv -f fusions.bedpe -r ref.fa -g ref.genePred -k

pvacfuse run \
<example_data_dir>/fusions.bedpe.annot \
${id} \
HLA-A*02:01,HLA-B*35:01,DRB1*11:01 \
MHCflurry MHCnuggetsI MHCnuggetsII NNalign NetMHC PickPocket SMM SMMPMBEC SMMalign \
/home/eunji/proj/RMLS/rna/8.pvacfuse/integrate-neo/${id} \
-e 8,9,10

------------------------------------------------------------------------------------------------------------------------------------
#agfusion (GRCh38.95 version 사용) 
conda create -n agfusion 
conda activate agfusion 
conda install python=3.5
conda install pyensembl 
pyensembl install --species homo_sapiens --release 95
#pyensembl install --species mus_musculus --release 87
pip install agfusion
agfusion download -g hg38
#agfusion download -g mm10
conda install star-fusion 
#star-fusion index download(31GB)
wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play.tar.gz
#star-fusion.sh (RAM 32GB)
cat config | while read id
do
STAR-Fusion --genome_lib_dir /HDD2T/eunji/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir 
--left_fq /HDD2T/jeeh9/RMLS_RNA/trimmed_fq/Paired/${id}-R_R1_P.fq.gz 
--right_fq /HDD2T/jeeh9/RMLS_RNA/trimmed_fq/Paired/${id}-R_R2_P.fq.gz 
--output_dir /HDD2T/eunji/sf/${id}  --no_remove_dups
done

agfusion batch \
-f <star_fusion_tsv> \
-a starfusion \
-db agfusion.homo_sapiens.95.db \
- <output_directory> \
--middlestar \
--noncanonical

pvacfuse run \
<example_data_dir>/agfusion/ \
Test \
HLA-A*02:01,HLA-B*35:01,DRB1*11:01 \
MHCflurry MHCnuggetsI MHCnuggetsII NNalign NetMHC PickPocket SMM SMMPMBEC SMMalign \
<output_dir> \
-e 8,9,10

#fusioncatcher 
conda create -n fusioncatcher 
conda activate fusioncatcher 
conda install fusioncatcher
#fusioncatcher.sh (700GB)
cat config | while read id 
do
fusioncatcher.py -d /home/eunji/miniconda3/envs/fusioncatcher/share/fusioncatcher-1.20/db/human_v98/ 
-i /HDD2T/jeeh9/RMLS_RNA/trimmed_fq/Paired/ -o /HDD2T/eunji/fc/${id}/
done

