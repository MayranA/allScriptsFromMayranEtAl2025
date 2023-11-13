# ATAC-seq

## Experiments

Two experiments of ATAC-seq were used in the manuscript:

- Analysis of time-course, duplicates of the following time-points: 48h, 72h, 96h, 120h

- Analysis of Snai1 mutants of the following time-points and genotypes:

    - 72h: KO clone 1 and 2 and WT
    - 96h: KO clone 1 and 2 and WT
    - 120h: KO clone 1 and WT

## By sample analysis: Fastq to peaks/coverage

The first part of analysis, fastq to peaks/coverage values were computed using a local [galaxy](https://doi.org/10.1093/nar/gkac247) server.

The workflow used have been exported [here](./ATACseq.ga). It has been run with the following inputs:

- PE fastq input: will be available on SRA/GEO
- reference_genome: `mm10`
- effective_genome_size: `1870000000`
- bin_size: `10`

The command lines used have been written [here](./ATACseq_CL.sh).


## Average of replicates' coverage

For the analysis of time-course, normalized coverage (per million reads in peaks) of replicates were averaged using the galaxy workflow available [here](../RNAseq/Galaxy-Workflow-Average_Bigwig_between_replicates.ga) with `bin_size` to 10.

The command lines are [here](./average_Bigwig_CL.sh)

## Identification of peaks with differential accessibility

### Quantification on shared peaks

#### Identification of consensus peaks

For each time point of the time-course, consensus peaks were called using [this galaxy workflow](./Get_Confident_Peaks_From_ATAC_or_CUTandRUN_duplicates.ga) with the BED files from the first workflow. The corresponding command lines are available [here](./Consensus_CL.sh) The strategy is described [here](https://github.com/iwc-workflows/consensus-peaks/blob/v0.1/README.md).

#### Merge and quantification

The consensus peaks were concatenated and merged to generate a single list of non-overlapping intervals. This list has been used to quantify with multiBAMsummary from deeptools. The Galaxy workflow is available [here](./ATAC_correlation_Time_Course.ga). The corresponding command lines are [here](./ATAC_correlation_CL.sh). The final output table is [here](./counts_on_peaks.txt.gz).

### Differential accessibility

DESeq2 was used to determine differential accessibility on consecutive timepoints. The R script is [here](./DESeq2_ATAC.R). The outputs are [here](./DESeq2_pairwise/). Among outputs, the list of regions significantly more accessible at 96h compared to 72h sorted by normalized coverage at 96h and the list of regions significantly less accessible at 96h compared to 72h sorted by normalized coverage at 72h.

## Heatmap

The two list of regions mentioned above were intersected with the consensus peaks of 48h and 120h respectively in order to split them into two categories. Then deeptools computeMatrix was used to get coverage on WT and Snai1 mutant coverage in a 5kb window centered on peaks and finally a heatmap was plotted.
The galaxy workflow used is [here](./Compute_heatmap.ga) and the corresponding command lines [here](./heatmap_CL.sh). Outputs are [here](../output.files/ATACseq/).

## Display coverage

In order to generate plots for figures. We used [pyGenomeTracks](https://github.com/deeptools/pyGenomeTracks). All command lines are available in [this bash script](./pgt.sh) and the outputs are [here](../output.files/ATACseq/pgt_outputs/).
