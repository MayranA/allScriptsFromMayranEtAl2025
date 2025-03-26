# Single-cell RNA-seq analysis FASTQ to matrices

The first steps of single-cell RNA-seq analysis was performed using a local [galaxy](https://doi.org/10.1093/nar/gkac247) server.

There are two types of experiments:

- Regular 10X (22 experiments)
- CellPlex 10X (7 experiments)

## FASTQ to matrices

The workflows used have been exported [here for regular 10X](./scRNA-seq_preprocessing_10X_v3_Bundle.ga) and [here for cellPlex](./scRNA-seq_preprocessing_10X_cellPlex_UPDUB.ga). They have been run with the following inputs:

- fastqPE collection GEX (for both workflows): fastqs are available on GEO/SRA.
- reference genome (for both workflows): `mm10` (from UCSC)
- gtf (for both workflows): available on [zenodo](https://doi.org/10.5281/zenodo.10079673).
- cellranger_barcodes_3M-february-2018.txt (for both workflows): downloaded from [zenodo](https://zenodo.org/record/3457880/files/3M-february-2018.txt.gz).
- fastqPE collection CMO (only CellPlex): fastqs are available on GEO/SRA
- cmo_10X_seq.txt (only CellPlex): available [here](./cmo_10X_seq.txt)
- sample name and CMO sequence collection (only CellPlex): see [CMO_samples](./CMO_samples) directory
- Number of expected cells (used by CITE-seq-Count) (only CellPlex): `24000`

The command lines used have been written [here](./scRNA-seq_preprocessing_CL.sh).

## Velocyto analysis

The workflow used have been exported [here](./Velocyto_on10X.ga). It has been run with the following inputs:

- BAM files from STAR solo (see above)
- Filtered barcodes from DropletUtils (see above)
- gtf same as above

The command lines used have been written [here](./velocyto_CL.sh).
