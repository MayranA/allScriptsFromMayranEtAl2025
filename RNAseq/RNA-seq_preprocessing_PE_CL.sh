
# toolshed.g2.bx.psu.edu/repos/lparsons/cutadapt/cutadapt/4.0+galaxy1
# command_version:4.0
ln -f -s '48h_WT_rep1_R1.fastq.gz.fastqsanger.gz' '48h_WT_rep1_1.fq.gz' 
  ln -f -s '48h_WT_rep1_R2.fastq.gz.fastqsanger.gz' '48h_WT_rep1_2.fq.gz' 
    cutadapt  -j=${GALAXY_SLOTS:-4}      -a 'Please use: For R1: - For Nextera: CTGTCTCTTATACACATCTCCGAGCCCACGAGAC - For TrueSeq: GATCGGAAGAGCACACGTCTGAACTCCAGTCAC or AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'='CTGTCTCTTATACACATCTCCGAGCCCACGAGAC'         -A 'Please use: For R2: - For Nextera: CTGTCTCTTATACACATCTGACGCTGCCGACGA - For TruSeq: GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT or AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'='CTGTCTCTTATACACATCTGACGCTGCCGACGAX'      --output='out1.fq.gz' --paired-output='out2.fq.gz'  --error-rate=0.1 --times=1 --overlap=3   --action=trim      --minimum-length=15 --pair-filter=any   --quality-cutoff=30     '48h_WT_rep1_1.fq.gz' '48h_WT_rep1_2.fq.gz'  > report.txt
# toolshed.g2.bx.psu.edu/repos/iuc/rgrnastar/rna_star/2.7.8a+galaxy1
# command_version:
STAR  --runThreadN ${GALAXY_SLOTS:-4} --genomeLoad NoSharedMemory --genomeDir '/data/galaxy/galaxy/var/tool-data/rnastar/2.7.4a/mm10_UCSC/mm10_UCSC/dataset_163252_files' --sjdbOverhang 99 --sjdbGTFfile 'mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gtf'   --readFilesIn 'Cutadapt on data 92 and data 91.fastqsanger.gz' 'Cutadapt on data 92 and data 91.fastqsanger.gz'   --readFilesCommand zcat   --outSAMtype BAM SortedByCoordinate  --twopassMode None ''  --quantMode GeneCounts   --outSAMattrIHstart 1 --outSAMattributes NH HI AS nM  --outSAMprimaryFlag OneBestScore  --outSAMmapqUnique 255    --outSAMunmapped None  --outFilterType BySJout --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.3 --outFilterMismatchNoverReadLmax 0.04 --outFilterScoreMin 0 --outFilterScoreMinOverLread 0.66 --outFilterMatchNmin 0 --outFilterMatchNminOverLread 0.66 --outSAMmultNmax -1 --outSAMtlen 1   --seedSearchStartLmax 50 --seedSearchStartLmaxOverLread 1.0 --seedSearchLmax 0 --seedMultimapNmax 10000 --seedPerReadNmax 1000 --seedPerWindowNmax 50 --seedNoneLociPerWindow 10  --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJstitchMismatchNmax 0 -1 0 0 --alignSJDBoverhangMin 1 --alignSplicedMateMapLmin 0 --alignSplicedMateMapLminOverLmate 0.66 --alignWindowsPerReadNmax 10000 --alignTranscriptsPerWindowNmax 100 --alignTranscriptsPerReadNmax 10000 --alignEndsType Local --peOverlapNbasesMin 0 --peOverlapMMp 0.01  --limitOutSJoneRead 1000 --limitOutSJcollapsed 1000000 --limitSjdbInsertNsj 1000000  --outBAMsortingThreadN ${GALAXY_SLOTS:-4} --outBAMsortingBinsN 50 --winAnchorMultimapNmax 50 --limitBAMsortRAM $((${GALAXY_MEMORY_MB:-0}*1000000))  
 samtools view -b -o 'RNA STAR on data 147, data 338, and data 337: mapped.bam.bam' Aligned.sortedByCoord.out.bam
# toolshed.g2.bx.psu.edu/repos/devteam/bamtools_filter/bamFilter/2.5.1+galaxy0
# command_version:
cp '/data/galaxy/galaxy/jobs/000/118/118588/configs/tmpau9z6gkd' 'Filter on data 429: JSON filter rules.txt' 
 ln -s 'RNA STAR on data 147, data 338, and data 337: mapped.bam.bam' localbam.bam 
 ln -s '/data/galaxy/data/_metadata_files/016/metadata_16757.dat' localbam.bam.bai 
 cat '/data/galaxy/galaxy/jobs/000/118/118588/configs/tmpau9z6gkd' 
 bamtools filter -script '/data/galaxy/galaxy/jobs/000/118/118588/configs/tmpau9z6gkd' -in localbam.bam -out 'Filter on data 429: Filtered BAM.bam'
# toolshed.g2.bx.psu.edu/repos/lldelisle/revertr2orientationinbam/revertR2orientationInBam/0.0.2
# command_version:
set -o pipefail; bash /data/galaxy/galaxy/var/shed_tools/toolshed.g2.bx.psu.edu/repos/lldelisle/revertr2orientationinbam/21ddefab2e4f/revertr2orientationinbam/revertR2orientationInBam.sh 'Filter on data 429: Filtered BAM.bam' 'Filter on data 429: Filtered BAM with R2 orientation reversed.bam'
# toolshed.g2.bx.psu.edu/repos/iuc/bedtools/bedtools_genomecoveragebed/2.30.0
# command_version:bedtools v2.30.0
bedtools genomecov  -ibam 'Filter on data 429: Filtered BAM with R2 orientation reversed.bam'   -split -strand +  -bg  -scale 0.0575374     > 'bedtools Genome Coverage on data 623.bedgraph'
# wig_to_bigWig
# command_version:
grep -v "^track" 'bedtools Genome Coverage on data 623.bedgraph' | wigToBigWig stdin '/data/galaxy/galaxy/var/tool-data/mm10_UCSC/len/mm10_UCSC.len' 'negative strand coverage.bigwig' -clip 2>&1 || echo "Error running wigToBigWig." >&2

# toolshed.g2.bx.psu.edu/repos/devteam/cufflinks/cufflinks/2.2.1.3
# command_version:cufflinks v2.2.1
cufflinks -q --no-update-check 'RNA STAR on data 147, data 338, and data 337: mapped.bam.bam' --num-threads "${GALAXY_SLOTS:-4}" --no-effective-length-correction -G 'mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gtf' -b '/data/galaxy/galaxy/var/tool-data/mm10_UCSC/sam_indexes/mm10_UCSC/mm10_UCSC.fa' -u --library-type fr-firststrand --mask-file 'chrM_mm10.gtf.gtf' --max-bundle-length 10000000 --max-bundle-frags 1000000 2> stderr 
 python '/data/galaxy/galaxy/var/shed_tools/toolshed.g2.bx.psu.edu/repos/devteam/cufflinks/d080005cffe1/cufflinks/mass.py' stderr 'None' "transcripts.gtf"
