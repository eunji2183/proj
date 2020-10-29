#Integrate-NEO
conda create -n integrate-neo
conda activate integrate-neo
conda install bwa python gcc_linux-64 ucsc-gtftogenepred samtools bedtools matplotlib cmake 
wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz
gtfToGenePred -genePredExt -geneNameAsName2 Homo_sapiens.GRCh38.84.gtf Homo_sapiens.GRCh38.84.genePred
cut -f 1-10,12 Homo_sapiens.GRCh38.84.genePred > tmp.txt
echo -e "#GRCh38.ensGene.name\tGRCh38.ensGene.chrom\tGRCh38.ensGene.strand\tGRCh38.ensGene.txStart\tGRCh38.ensGene.txEnd\tGRCh38.ensGene.cdsStart\tGRCh38.ensGene.cdsEnd\tGRCh38.ensGene.exonCount\tGRCh38.ensGene.exonStarts\tGRCh38.ensGene.exonEnds\tGRCh38.ensemblToGeneName.value" > Homo_sapiens.GRCh38.84.tsv
cat tmp.txt >> Homo_sapiens.GRCh38.84.tsv 


git clone https://github.com/ChrisMaherLab/INTEGRATE-Vis.git
cd INTEGRATE-Vis.1.0.0
chmod +x install.sh
./install.sh -o /home/eunji/proj/0_sh/biosoft/integrate-vis
cd ~/proj/0_sh/biosoft/INTEGRATE-Vis/INTEGRATE-Vis.1.0.0/data/reference_genome/
wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{1..22}.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz
wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz
wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz
gunzip -c Homo_sapiens.GRCh38.dna.chromosome.* > GRCh38_r84.all.fa


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
./Integrate mkbwt /home/eunji/proj/0_sh/biosoft/INTEGRATE-Vis/INTEGRATE-Vis.1.0.0/data/reference_genome/GRCh38_r84.all.fa 

#integrate fusion 
./Integrate fusion /home/eunji/proj/0_sh/biosoft/INTEGRATE-Vis/INTEGRATE-Vis.1.0.0/data/reference_genome/GRCh38_r84.all.fa \
/home/eunji/proj/0_sh/biosoft/INTEGRATE-Vis/INTEGRATE-Vis.1.0.0/data/gene_model/Homo_sapiens.GRCh38.84.tsv \
./bwts/ \
/home/eunji/proj/RMLS/rna/4.align/hisat2_tran/${id}.sort.bam 

#tophat mapping 
wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
bowtie2-build Homo_sapiens.GRCh38.dna.primary_assembly.fa index
tophat -G Homo_sapiens.GRCh38.99.gtf --transcriptome-index=GRCh38.99.tr GRCh38.99
#tophat.sh 
cat config | while read id
do
tophat -o /HDD2T/eunji/th/${id} --transcriptome-index=/home/eunji/ref/GRCh38.99.tr/Homo_sapiens.GRCh38.99 
/home/eunji/ref/GRCh38.99.tr/GRCh38.99 
/HDD2T/jeeh9/RMLS_RNA/trimmed_fq/Paired/${id}-R_R1_P.fq.gz 
/HDD2T/jeeh9/RMLS_RNA/trimmed_fq/Paired/${id}-R_R2_P.fq.gz
done

#Integrate fusion 


------------------------------------------------------------------------------------------------------------------------------------
#agfusion 
conda create -n agfusion 
conda activate agfusion 
conda install pyensembl 
pyensembl install --species homo_sapiens --release 84
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

