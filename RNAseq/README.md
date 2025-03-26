# RNA-seq

## Experiments

Two experiments of RNA-seq were used in the manuscript:

- Reanalysis of GSE106226 (GSM2832454 to GSM2832463), single-read datasets. Fastqs were retrived from SRA using fasterq-dump version 3.0.5.

- Analysis of paired-end datasets generated in this study

## Fastq to counts/FPKM/coverage

For both, the first part of analysis, fastq to counts/coverage/FPKM values were computed using a local [galaxy](https://doi.org/10.1093/nar/gkac247) server.

The workflows used have been exported [here for single-read](./Galaxy-Workflow-RNAseq_SingleRead.ga) and [here for paired-end](./Galaxy-Workflow-RNAseq_PE.ga). They have been run with the following inputs:
- SR fastq input (for single-read): from SRA
- forward_adapter (for single-read): `CTGTCTCTTATACACATCTCCGAGCCCACGAGAC`
- strandness (for single-read): `unstranded`
- PE fastq input (for paired-end): from this study will be available on SRA/GEO
- forward_adapter (for paired-end): `CTGTCTCTTATACACATCTCCGAGCCCACGAGAC`
- reverse_adapter (for paired-end): `CTGTCTCTTATACACATCTGACGCTGCCGACGA`
- strandness (for paired-end): `stranded - reverse`
- reference genome (for both workflows): `mm10` (from UCSC)
- gtf (for both workflows): https://doi.org/10.5281/zenodo.7510406
- gtf with regions to exclude from FPKM normalization  (for both workflows): [chrM.gtf](./chrM.gtf)

The command lines used have been written [here](./RNA-seq_preprocessing_CL.sh).

## Combine counts and FPKM

Combination of individual counts / FPKM have been computing using the R scripts available on the [rnaseq_rscript repository](https://github.com/lldelisle/rnaseq_rscripts/). The detailed commands are available [here](./merge_tables.sh) and the merged tables are [here](./mergedTables/).

## Launch DESeq2

All differential gene expression analysis was computed with DESeq2. For mutants, the R script is available [here](./DESeq2_mutant.R) and the results are [here](./../output.files/RNAseq/DESeq2_pairwise/). For time-course, the R script is available [here](./DESeq2_time-course.R) and the results are [here](./../output.files/RNAseq/DESeq2_pairwise_time-course/).

## Average of replicates' coverage

For datasets generated in this study, normalized coverage of replicates were averaged using the galaxy workflow available [here](./Galaxy-Workflow-Average_Bigwig_between_replicates.ga) with `bin_size` to 5.

The command lines are [here](./Average_Bigwig.sh)

## Display coverage

In order to generate plots for figures. We used [pyGenomeTracks](https://github.com/deeptools/pyGenomeTracks). All command lines are available in [this bash script](./pgt.sh).
