wd=/scratch/ldelisle/Alex_plots/getData/
pathRNA=/scratch/ldelisle/Alex_plots/toGEO/RNAseq
dirPathWithDependencies="/home/ldelisle/softwares/"
chrsToRemove="chrM"
condaEnvName=RNAseq_R_202302
filePathForGTF=/scratch/ldelisle/Alex_plots/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz

# Check everything is set correctly:
if [ ! -z ${condaEnvName} ]; then
    # Conda environment:
    # This line is to adapt the conda to the shell
    source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
    # We check if the conda environment exists
    exists=$(conda info --envs | awk -v ce=${condaEnvName} '$1==ce{print}' | wc -l)
    # It if does not exists an error is raised
    if [ $exists -ne 1 ]; then
    echo "conda environment ${condaEnvName} does not exists. Create it before."
    exit 1
    fi
    # Activate the conda environment
    conda activate ${condaEnvName}
fi

# Check all softwares are present and write version to stdout:
R --version
if [ $? -ne 0 ]
then
  echo "R is not installed but required. Please install it for example in the conda environment."
  exit 1
fi
R --version > ${wd}/versions.txt
# Check  are installed:
Rscript -e "library(colorspace);library(DESeq2);library(ggplot2);library(pheatmap);library(RColorBrewer);library(rtracklayer)"
if [ $? -ne 0 ]
then
  echo "Some R packages are missing check rtracklayer, colorspace, DESeq2, ggplot2, pheatmap and RColorBrewer are installed."
  exit 1
fi

# Check the 2 github:
if [ ! -e ${dirPathWithDependencies}rnaseq_rscripts/ ]; then
  echo "${dirPathWithDependencies}rnaseq_rscripts/ does not exists please clone https://github.com/lldelisle/rnaseq_rscripts"
  exit 1
fi
cd ${dirPathWithDependencies}rnaseq_rscripts/
echo "Version of rnaseq_rscripts" >> ${wd}/versions.txt
git rev-parse HEAD >> ${wd}/versions.txt

if [ ! -e ${dirPathWithDependencies}toolBoxForMutantAndWTGenomes/ ]; then
  echo "${dirPathWithDependencies}toolBoxForMutantAndWTGenomes/ does not exists please clone https://github.com/lldelisle/toolBoxForMutantAndWTGenomes"
  exit 1
fi
cd ${dirPathWithDependencies}toolBoxForMutantAndWTGenomes/
echo "Version of toolBoxForMutantAndWTGenomes" >> ${wd}/versions.txt
git rev-parse HEAD >> ${wd}/versions.txt


## START
cd "${wd}"

### Single gastruloids

filePathForSamplesPlan=samples.plan.single.gastruloids.txt
mkdir -p single_gastruloids/
dirWithFPKMs=$(ls -d $wd/91*/*gene_expression)
echo -e "sample\tcufflinks_file\thtseq_count_file" > $filePathForSamplesPlan
for f in $dirWithFPKMs/*; do
  count_file=$(ls $wd/91*/*HTS*/$(basename $f))
  echo -e "$(basename $f .tabular)\t${f}\t${count_file}" >> $filePathForSamplesPlan
done

configFile="configFileRNAseq_step1.R"

echo "
### Required for all steps ###
RNAseqFunctionPath <- \"${dirPathWithDependencies}/rnaseq_rscripts/RNAseqFunctions.R\"
# This file should be a tabulated file with at least one column called
# 'sample'. Optionnaly, the paths to the counts tables and FPKM tables can be
# provided under the column called: htseq_count_file and cufflinks_file.
samplesPlan <- \"${filePathForSamplesPlan}\"

#### STEP 1 - MERGE TABLES ####
# If the merged tables are not already generated:
outputFolderForStep1 <- \"${wd}/single_gastruloids/\"
# Needed for DESeq2: Do you want to merge counts? T=yes F or commented=no
mergeCounts <- T
# Optional: subset the count table Do you want to remove some genes from the
# count table
subsetCounts <- F
# If the table with counts have already been generated and you just want to
# remove some genes.
# initialTableWithCount<-'${dirPathWithDependencies}/rnaseq_rscripts/example/mergedTablesE/AllHTSeqCounts.txt'
# If you provide the initialTableWithCount you need to provide the name of the
# column with the ensembl id.
# geneIDColInInitialTable<-'Ens_ID'
# List of genes id to remove (one per line with no header).
genesToRmFromCounts <- \"$PWD/genes${chrsToRemove//,/_}.txt\"
# Optional:
mergeFPKM <- T
# By default cufflinks split the transcripts which do not overlap in different
# locus and so different lines, put T if you want to sum the FPKM for non
# overlapping transcripts (put F if not).
oneLinePerEnsemblID <- T
# Optional: subset the FPKM table Do you want to remove some genes from the
# FPKM table
subsetFPKM <- F
chrToRemove <- c(\"${chrsToRemove//,/\",\"}\")
# Anouk method: Genes that have the less variable rank should have the same
# expression.
normFPKMWithAnoukMethod <- F
# If the table with FPKM have already been generated and you just want to
# normalize it.
# initialTableWithFPKM<-'${dirPathWithDependencies}/rnaseq_rscripts/example/mergedTablesE/AllCufflinks_Simplified.txt'
# Usually, it is recommanded to remove mitochondrial genes before doing the
# normalization. In some cases, it can also be useful to remove the sex
# chromosomes (put c('chrX','chrY','chrM')). If you do not want to remove any
# gene put NA or comment the line.
chrToRemoveBeforeNormWithAnoukMethod <- c(\"${chrsToRemove//,/\",\"}\")
# Default is 1000, you can change here.
nbOfGenesWithAnoukMethod <- 1000
# If you want to keep the genes used in the normalization from Anouk, they will
# be written in a file.
keepGenesUsedForNorm <- F
" > ${configFile}

if [ ! -e genes${chrsToRemove//,/_}.txt ]; then
  Rscript ${dirPathWithDependencies}toolBoxForMutantAndWTGenomes/scripts/getGeneListFromChrAndGTF.R $filePathForGTF ${chrsToRemove} ./
  mv genesIn${chrsToRemove}from* genes${chrsToRemove//,/_}.txt
fi
Rscript ${dirPathWithDependencies}rnaseq_rscripts/step1-generateTables.R $configFile

### Mutants

filePathForSamplesPlan=samples.plan.mutants.txt
mkdir -p mutants/
dirWithFPKMs=$(dirname $(ls $wd/df*/*gene_expression/* | grep WT | head -n 1))
echo -e "sample\tcufflinks_file\thtseq_count_file" > $filePathForSamplesPlan
for f in $dirWithFPKMs/*; do
  count_file=$(ls $wd/df*/*HTS*/$(basename $f))
  echo -e "$(basename $f .tabular)\t${f}\t${count_file}" >> $filePathForSamplesPlan
done

configFile="configFileRNAseq_step1.R"

echo "
### Required for all steps ###
RNAseqFunctionPath <- \"${dirPathWithDependencies}/rnaseq_rscripts/RNAseqFunctions.R\"
# This file should be a tabulated file with at least one column called
# 'sample'. Optionnaly, the paths to the counts tables and FPKM tables can be
# provided under the column called: htseq_count_file and cufflinks_file.
samplesPlan <- \"${filePathForSamplesPlan}\"

#### STEP 1 - MERGE TABLES ####
# If the merged tables are not already generated:
outputFolderForStep1 <- \"${wd}/mutants/\"
# Needed for DESeq2: Do you want to merge counts? T=yes F or commented=no
mergeCounts <- T
# Optional: subset the count table Do you want to remove some genes from the
# count table
subsetCounts <- F
# If the table with counts have already been generated and you just want to
# remove some genes.
# initialTableWithCount<-'${dirPathWithDependencies}/rnaseq_rscripts/example/mergedTablesE/AllHTSeqCounts.txt'
# If you provide the initialTableWithCount you need to provide the name of the
# column with the ensembl id.
# geneIDColInInitialTable<-'Ens_ID'
# List of genes id to remove (one per line with no header).
genesToRmFromCounts <- \"$PWD/genes${chrsToRemove//,/_}.txt\"
# Optional:
mergeFPKM <- T
# By default cufflinks split the transcripts which do not overlap in different
# locus and so different lines, put T if you want to sum the FPKM for non
# overlapping transcripts (put F if not).
oneLinePerEnsemblID <- T
# Optional: subset the FPKM table Do you want to remove some genes from the
# FPKM table
subsetFPKM <- F
chrToRemove <- c(\"${chrsToRemove//,/\",\"}\")
# Anouk method: Genes that have the less variable rank should have the same
# expression.
normFPKMWithAnoukMethod <- F
# If the table with FPKM have already been generated and you just want to
# normalize it.
# initialTableWithFPKM<-'${dirPathWithDependencies}/rnaseq_rscripts/example/mergedTablesE/AllCufflinks_Simplified.txt'
# Usually, it is recommanded to remove mitochondrial genes before doing the
# normalization. In some cases, it can also be useful to remove the sex
# chromosomes (put c('chrX','chrY','chrM')). If you do not want to remove any
# gene put NA or comment the line.
chrToRemoveBeforeNormWithAnoukMethod <- c(\"${chrsToRemove//,/\",\"}\")
# Default is 1000, you can change here.
nbOfGenesWithAnoukMethod <- 1000
# If you want to keep the genes used in the normalization from Anouk, they will
# be written in a file.
keepGenesUsedForNorm <- F
" > ${configFile}

if [ ! -e genes${chrsToRemove//,/_}.txt ]; then
  Rscript ${dirPathWithDependencies}toolBoxForMutantAndWTGenomes/scripts/getGeneListFromChrAndGTF.R $filePathForGTF ${chrsToRemove} ./
  mv genesIn${chrsToRemove}from* genes${chrsToRemove//,/_}.txt
fi
Rscript ${dirPathWithDependencies}rnaseq_rscripts/step1-generateTables.R $configFile


### Time course

filePathForSamplesPlan=samples.plan.time.course.txt
mkdir -p time_course/
dirWithFPKMs=$(dirname $(ls $wd/df*/*gene_expression/* | grep 48h | head -n 1))
echo -e "sample\tcufflinks_file\thtseq_count_file" > $filePathForSamplesPlan
for f in $dirWithFPKMs/*; do
  count_file=$(ls $wd/df*/*HTS*/$(basename $f))
  echo -e "$(basename $f .tabular)\t${f}\t${count_file}" >> $filePathForSamplesPlan
done

configFile="configFileRNAseq_step1.R"

echo "
### Required for all steps ###
RNAseqFunctionPath <- \"${dirPathWithDependencies}/rnaseq_rscripts/RNAseqFunctions.R\"
# This file should be a tabulated file with at least one column called
# 'sample'. Optionnaly, the paths to the counts tables and FPKM tables can be
# provided under the column called: htseq_count_file and cufflinks_file.
samplesPlan <- \"${filePathForSamplesPlan}\"

#### STEP 1 - MERGE TABLES ####
# If the merged tables are not already generated:
outputFolderForStep1 <- \"${wd}/time_course/\"
# Needed for DESeq2: Do you want to merge counts? T=yes F or commented=no
mergeCounts <- T
# Optional: subset the count table Do you want to remove some genes from the
# count table
subsetCounts <- F
# If the table with counts have already been generated and you just want to
# remove some genes.
# initialTableWithCount<-'${dirPathWithDependencies}/rnaseq_rscripts/example/mergedTablesE/AllHTSeqCounts.txt'
# If you provide the initialTableWithCount you need to provide the name of the
# column with the ensembl id.
# geneIDColInInitialTable<-'Ens_ID'
# List of genes id to remove (one per line with no header).
genesToRmFromCounts <- \"$PWD/genes${chrsToRemove//,/_}.txt\"
# Optional:
mergeFPKM <- T
# By default cufflinks split the transcripts which do not overlap in different
# locus and so different lines, put T if you want to sum the FPKM for non
# overlapping transcripts (put F if not).
oneLinePerEnsemblID <- T
# Optional: subset the FPKM table Do you want to remove some genes from the
# FPKM table
subsetFPKM <- F
chrToRemove <- c(\"${chrsToRemove//,/\",\"}\")
# Anouk method: Genes that have the less variable rank should have the same
# expression.
normFPKMWithAnoukMethod <- F
# If the table with FPKM have already been generated and you just want to
# normalize it.
# initialTableWithFPKM<-'${dirPathWithDependencies}/rnaseq_rscripts/example/mergedTablesE/AllCufflinks_Simplified.txt'
# Usually, it is recommanded to remove mitochondrial genes before doing the
# normalization. In some cases, it can also be useful to remove the sex
# chromosomes (put c('chrX','chrY','chrM')). If you do not want to remove any
# gene put NA or comment the line.
chrToRemoveBeforeNormWithAnoukMethod <- c(\"${chrsToRemove//,/\",\"}\")
# Default is 1000, you can change here.
nbOfGenesWithAnoukMethod <- 1000
# If you want to keep the genes used in the normalization from Anouk, they will
# be written in a file.
keepGenesUsedForNorm <- F
" > ${configFile}

if [ ! -e genes${chrsToRemove//,/_}.txt ]; then
  Rscript ${dirPathWithDependencies}toolBoxForMutantAndWTGenomes/scripts/getGeneListFromChrAndGTF.R $filePathForGTF ${chrsToRemove} ./
  mv genesIn${chrsToRemove}from* genes${chrsToRemove//,/_}.txt
fi
Rscript ${dirPathWithDependencies}rnaseq_rscripts/step1-generateTables.R $configFile

# gzip all tables
gzip */All*.txt
