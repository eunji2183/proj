# mkdir                                                                         

mkdir -p biosoft project data
cd project
mkdir -p 1.raw_fq 2.clean_fq 3.qc/{raw_qc,clean_qc} 4.align/{hisat2_genome,hisa\
t2_tran} 5.DEG/{htseq-count,featurecounts} 6.isoform/{salmon,stringtie,stringti\
e2} 7.fusion
