
#python=2.7 annovar netMHCpan-4.1 netMHCIIpan-3.2 PeptideMatch
conda activate neopredpipe 
git clone https://github.com/MathOnco/NeoPredPipe.git
#Configure the 'usr_path.ini' file

#https://www.uniprot.org/news/2005/02/01/release > uniprot_sprot.fasta
#https://research.bioinformatics.udel.edu/peptidematch/commandlinetool.jsp , peptidematch
wget ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
java -jar /home/eunji/proj/0_sh/biosoft/PeptideMatchCMD_src_1.0/PeptideMatchCMD_1.0.jar -a index -d uniprot_sprot.fasta -i sprot_index




#https://github.com/MathOnco/NeoPredPipe
python NeoPredPipe.py \
-I somatic.vcf \
-H hlatypes.txt \
-o ./ \
-n TestRun \
-c 1 2 -E 8 9 10
