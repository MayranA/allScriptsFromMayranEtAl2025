# toolshed.g2.bx.psu.edu/repos/lparsons/cutadapt/cutadapt/4.0+galaxy1
# command_version:4.0
ln -f -s 'SRR6224362.fastqsanger.gz' 'SRR6224362.fq.gz' 
 cutadapt -j=${GALAXY_SLOTS:-4} -a 'Please use: For R1: - For Nextera: CTGTCTCTTATACACATCTCCGAGCCCACGAGAC - For TrueSeq: GATCGGAAGAGCACACGTCTGAACTCCAGTCAC or AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'='CTGTCTCTTATACACATCTCCGAGCCCACGAGACX' --output='out1.fq.gz' --action=trim --minimum-length=15 --quality-cutoff=30 'SRR6224362.fq.gz' > report.txt

# toolshed.g2.bx.psu.edu/repos/iuc/rgrnastar/rna_star/2.7.8a+galaxy1
# command_version:
STAR --runThreadN ${GALAXY_SLOTS:-4} --genomeDir '/data/galaxy/galaxy/var/tool-data/rnastar/2.7.4a/mm10_UCSC/mm10_UCSC/dataset_163252_files' --sjdbOverhang 99 --sjdbGTFfile 'mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gtf' --readFilesIn 'Cutadapt on data 55: Read 1 Output.fastqsanger.gz' --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate '' --quantMode GeneCounts --outSAMattributes NH HI AS nM --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJstitchMismatchNmax 0 -1 0 0 --alignSJDBoverhangMin 1 --outBAMsortingThreadN ${GALAXY_SLOTS:-4} --limitBAMsortRAM $((${GALAXY_MEMORY_MB:-0}*1000000)) 
 samtools view -b -o 'RNA STAR on data 7 and data 68: mapped.bam.bam' Aligned.sortedByCoord.out.bam

# toolshed.g2.bx.psu.edu/repos/devteam/bamtools_filter/bamFilter/2.5.2+galaxy1
# command_version:
cp '/data/galaxy/galaxy/jobs/000/118/118228/configs/tmpq19hsxpr' 'Filter BAM on data 99: JSON filter rules.txt' 
 ln -s 'RNA STAR on data 7 and data 68: mapped.bam.bam' localbam.bam 
 ln -s '/data/galaxy/data/_metadata_files/016/metadata_16674.dat' localbam.bam.bai 
 cat '/data/galaxy/galaxy/jobs/000/118/118228/configs/tmpq19hsxpr' 
 bamtools filter -script '/data/galaxy/galaxy/jobs/000/118/118228/configs/tmpq19hsxpr' -in localbam.bam -out 'Filter BAM on data 99: Filtered BAM.bam'

# toolshed.g2.bx.psu.edu/repos/iuc/bedtools/bedtools_genomecoveragebed/2.30.0
# command_version:bedtools v2.30.0
bedtools genomecov  -ibam 'Filter BAM on data 99: Filtered BAM.bam'   -split -strand -  -bg  -scale 0.0320191     > 'bedtools Genome Coverage on data 143.bedgraph'

# wig_to_bigWig
# command_version:
grep -v "^track" 'bedtools Genome Coverage on data 143.bedgraph' | wigToBigWig stdin '/data/galaxy/galaxy/var/tool-data/mm10_UCSC/len/mm10_UCSC.len' 'negative strand coverage.bigwig' -clip 2>&1 || echo "Error running wigToBigWig." >&2

# toolshed.g2.bx.psu.edu/repos/devteam/cufflinks/cufflinks/2.2.1.3
# command_version:cufflinks v2.2.1
cufflinks -q --no-update-check 'RNA STAR on data 7 and data 68: mapped.bam.bam' --num-threads "${GALAXY_SLOTS:-4}" --no-effective-length-correction -G 'mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gtf' -b '/data/galaxy/galaxy/var/tool-data/mm10_UCSC/sam_indexes/mm10_UCSC/mm10_UCSC.fa' -u --library-type fr-unstranded --mask-file 'chrM_mm10.gtf.gtf' --max-bundle-length 10000000 --max-bundle-frags 1000000 2> stderr 
 python '/data/galaxy/galaxy/var/shed_tools/toolshed.g2.bx.psu.edu/repos/devteam/cufflinks/d080005cffe1/cufflinks/mass.py' stderr 'None' "transcripts.gtf"
