# toolshed.g2.bx.psu.edu/repos/lparsons/cutadapt/cutadapt/4.4+galaxy0
# command_version:4.4
ln -f -s '48h_WT_rep1_R1.fastq.gz.fastqsanger.gz' '48h_WT_rep1_1.fq.gz' 
  ln -f -s '48h_WT_rep1_R2.fastq.gz.fastqsanger.gz' '48h_WT_rep1_2.fq.gz' 
    cutadapt  -j=${GALAXY_SLOTS:-4}      -a 'Nextera R1'='CTGTCTCTTATACACATCTCCGAGCCCACGAGAC'         -A 'Nextera R2'='CTGTCTCTTATACACATCTGACGCTGCCGACGA'      --output='out1.fq.gz' --paired-output='out2.fq.gz'  --error-rate=0.1 --times=1 --overlap=3   --action=trim      --minimum-length=15 --pair-filter=any   --quality-cutoff=30     '48h_WT_rep1_1.fq.gz' '48h_WT_rep1_2.fq.gz'  > report.txt

# toolshed.g2.bx.psu.edu/repos/devteam/bowtie2/bowtie2/2.5.0+galaxy0
# command_version:/data/galaxy/galaxy/var/dependencies/_conda/envs/mulled-v1-536425fba3670dfc88459e23d8cc068013e0b6d61ae3e2101de6bda1da4a4c6a/bin/bowtie2-align-s version 2.5.0
64-bit
Built on fv-az123-980
Tue Nov  1 03:44:13 UTC 2022
Compiler: gcc version 10.4.0 (conda-forge gcc 10.4.0-19) 
Options: -O3 -msse2 -funroll-loops -g3 -fvisibility-inlines-hidden -std=c++17 -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /data/galaxy/galaxy/var/dependencies/_conda/envs/mulled-v1-536425fba3670dfc88459e23d8cc068013e0b6d61ae3e2101de6bda1da4a4c6a/include -fdebug-prefix-map=/opt/conda/conda-bld/bowtie2_1667273633358/work=/usr/local/src/conda/bowtie2-2.5.0 -fdebug-prefix-map=/data/galaxy/galaxy/var/dependencies/_conda/envs/mulled-v1-536425fba3670dfc88459e23d8cc068013e0b6d61ae3e2101de6bda1da4a4c6a=/usr/local/src/conda-prefix -std=c++11 -DPOPCNT_CAPABILITY -DNO_SPINLOCK -DWITH_QUEUELOCK=1 -DWITH_ZSTD
Sizeof {int, long, long long, void*, size_t, off_t}: {4, 8, 8, 8, 8, 8}
set -o | grep -q pipefail 
 set -o pipefail;   ln -s 'Cutadapt on data 19 and data 18.fastqsanger.gz' input_f.fastq.gz 
  ln -s 'Cutadapt on data 19 and data 18.fastqsanger.gz' input_r.fastq.gz 
    bowtie2  -p ${GALAXY_SLOTS:-4}  -x '/data/galaxy/galaxy/var/tool-data/mm10_UCSC/bowtie2_index/mm10_UCSC/mm10_UCSC'   -1 'input_f.fastq.gz' -2 'input_r.fastq.gz'              --very-sensitive   2> 'mapping stats.txt'  | samtools sort --no-PG -@${GALAXY_SLOTS:-2} -T "${TMPDIR:-.}" -O bam -o 'bowtie2 output (BAM).bam'

# toolshed.g2.bx.psu.edu/repos/devteam/bamtools_filter/bamFilter/2.5.2+galaxy1
# command_version:
cp '/data/galaxy/galaxy/jobs/000/123/123980/configs/tmp18cax_3d' 'Filter BAM on data 417: JSON filter rules.txt' 
 ln -s 'bowtie2 output (BAM).bam' localbam.bam 
 ln -s '/data/galaxy/data/_metadata_files/017/metadata_17094.dat' localbam.bam.bai 
 cat '/data/galaxy/galaxy/jobs/000/123/123980/configs/tmp18cax_3d' 
 bamtools filter -script '/data/galaxy/galaxy/jobs/000/123/123980/configs/tmp18cax_3d' -in localbam.bam -out 'filtered BAM.bam'

# toolshed.g2.bx.psu.edu/repos/devteam/picard/picard_MarkDuplicates/2.18.2.4
# command_version:
_JAVA_OPTIONS=${_JAVA_OPTIONS:-"-Xmx2048m -Xms256m -Djava.io.tmpdir=${TMPDIR:-${_GALAXY_JOB_TMPDIR}}"} 
 export _JAVA_OPTIONS 
   ln -f -s 'filtered BAM.bam' '48h_WT_rep1' 
  picard MarkDuplicates  INPUT='48h_WT_rep1' OUTPUT='BAM filtered rmDup.bam'  METRICS_FILE='MarkDuplicates metrics.txt'  REMOVE_DUPLICATES='true' ASSUME_SORTED='true'  DUPLICATE_SCORING_STRATEGY='SUM_OF_BASE_QUALITIES'  OPTICAL_DUPLICATE_PIXEL_DISTANCE='100'   VALIDATION_STRINGENCY='LENIENT' TAGGING_POLICY=All QUIET=true VERBOSITY=ERROR

# toolshed.g2.bx.psu.edu/repos/iuc/bedtools/bedtools_bamtobed/2.30.0+galaxy2
# command_version:bedtools v2.30.0
ln -s 'BAM filtered rmDup.bam' ./input.bam 
  bedtools bamtobed    -i ./input.bam > 'BED filtered rmDup.bed'

# toolshed.g2.bx.psu.edu/repos/iuc/macs2/macs2_callpeak/2.2.9.1+galaxy0
# command_version:macs2 2.2.9.1
export PYTHON_EGG_CACHE=`pwd` 
   (macs2 callpeak   -t 'BED filtered rmDup.bed'  --name 48h_WT_rep1    --format BED   --gsize '1870000000'           --call-summits  --keep-dup 'all'  --d-min 20 --buffer-size 100000  --bdg  --qvalue '0.05'  --nomodel --extsize '200' --shift '-100'  2>&1 > macs2_stderr) 
 cp 48h_WT_rep1_peaks.xls 'MACS2 peaks xls.tabular'   
 ( count=`ls -1 48h_WT_rep1* 2>/dev/null | wc -l`; if [ $count != 0 ]; then mkdir '/data/galaxy/galaxy/jobs/000/124/124036/working/dataset_530937_files' 
 cp -r 48h_WT_rep1* '/data/galaxy/galaxy/jobs/000/124/124036/working/dataset_530937_files' 
 python '/data/galaxy/galaxy/var/shed_tools/toolshed.g2.bx.psu.edu/repos/iuc/macs2/86e2413cf3f8/macs2/dir2html.py' '/data/galaxy/galaxy/jobs/000/124/124036/working/dataset_530937_files' macs2_stderr > 'MACS2 callpeak on data 488 (html report).html'; fi; ) 
 exit_code_for_galaxy=$? 
 cat macs2_stderr 2>&1 
 (exit $exit_code_for_galaxy)

# wig_to_bigWig
# command_version:
grep -v "^track" 'MACS2 treatment coverage.bedgraph' | wigToBigWig stdin '/data/galaxy/galaxy/var/tool-data/mm10_UCSC/len/mm10_UCSC.len' 'Coverage from MACS2 (bigwig).bigwig' -clip 2>&1 || echo "Error running wigToBigWig." >&2

# toolshed.g2.bx.psu.edu/repos/bgruening/deeptools_bigwig_average/deeptools_bigwig_average/3.5.4+galaxy0
# command_version:bigwigAverage 3.5.4
bigwigAverage --numberOfProcessors "${GALAXY_SLOTS:-4}" --bigwigs Coverage from MACS2 (bigwig).bigwig --outFileName 'bigwig normalized per million reads in peaks.bigwig' --outFileFormat 'bigwig'     --scaleFactors '0.13428370386683353' --binSize 10
