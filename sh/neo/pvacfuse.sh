#Integrate-NEO
conda create -n integrate-neo
conda activate integrate-neo
conda install bwa python gcc_linux-64 ucsc-gtftogenepred samtools bedtools matplotlib cmake 
wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz
gtfToGenePred -genePredExt -geneNameAsName2 Homo_sapiens.GRCh38.84.gtf Homo_sapiens.GRCh38.84.genePred

git clone https://github.com/ChrisMaherLab/INTEGRATE-Vis.git
cd ~/proj/0_sh/biosoft/INTEGRATE-Vis/INTEGRATE-Vis.1.0.0/data/reference_genome/
wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{1..22}.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz
wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz
wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz

cd INTEGRATE-Neo-V-1.2.1
chmod +x install.sh
./install.sh -o /home/eunji/proj/0_sh/biosoft/



