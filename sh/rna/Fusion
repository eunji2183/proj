# STAR alignment 

#STAR-fusion 
 STAR --genomeDir ${star_index_dir} \                                                                                     
          --readFilesIn ${left_fq_filename} ${right_fq_filename} \                                                                      
          --outReadsUnmapped None \
          --twopassMode Basic \
          --readFilesCommand "gunzip -c" \
          --outSAMstrandField intronMotif \  # include for potential use with StringTie for assembly
          --outSAMunmapped Within 
          --chimSegmentMin 12 \  # ** essential to invoke chimeric read detection & reporting **
          --chimJunctionOverhangMin 8 \
          --chimOutJunctionFormat 1 \   # **essential** includes required metadata in Chimeric.junction.out file.
          --alignSJDBoverhangMin 10 \
          --alignMatesGapMax 100000 \   # avoid readthru fusions within 100k
          --alignIntronMax 100000 \
          --alignSJstitchMismatchNmax 5 -1 5 5 \   # settings improved certain chimera detections
          --outSAMattrRGline ID:GRPundef \
          --chimMultimapScoreRange 3 \
          --chimScoreJunctionNonGTAG -4 \
          --chimMultimapNmax 20 \
          --chimNonchimScoreDropMin 10 \
          --peOverlapNbasesMin 12 \
          --peOverlapMMp 0.1 \
          --alignInsertionFlush Right \
          --alignSplicedMateMapLminOverLmate 0 \
          --alignSplicedMateMapLmin 30
# out > Chimeric.out.junction
STAR-Fusion --genome_lib_dir /path/to/your/CTAT_resource_lib \
             -J Chimeric.out.junction \
             --output_dir star_fusion_outdir
             
#Arriba 
STAR \
    --runThreadN 8 \
    --genomeDir /path/to/STAR_index --genomeLoad NoSharedMemory \
    --readFilesIn read1.fastq.gz read2.fastq.gz --readFilesCommand zcat \
    --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
    --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \
    --chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 \
    --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 |
arriba \
    -x /dev/stdin \
    -o fusions.tsv -O fusions.discarded.tsv \
    -a /path/to/assembly.fa -g /path/to/annotation.gtf \
    -b /path/to/blacklist.tsv.gz -k /path/to/known_fusions.tsv.gz -t /path/to/known_fusions.tsv.gz -p /path/to/protein_domains.gff3

