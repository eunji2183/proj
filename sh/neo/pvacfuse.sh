#Integrate-NEO
conda create -n integrate-neo
conda activate integrate-neo
conda install bwa python gcc ucsc-gtftogenepred


cd INTEGRATE-Neo-V-1.2.1
chmod +x install.sh
./install.sh -o /home/eunji/proj/0_sh/biosoft/
gtfToGenePred -genePredExt -geneNameAsName2 
