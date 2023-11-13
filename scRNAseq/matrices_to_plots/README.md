# Single-cell RNA-seq analysis matrices to plots

## R configuration

The scripts have been launched on a RStudio server.

The session information is available [here](./sessionInfo.txt).

## Pipeline

All samples have been registered into a [csv file](./metadata.sample.galaxy.csv).

Once matrices have been generated in galaxy they were downloaded and organized like this:

```bash
.
├── CMO
│   └── Directory1
│       ├── barcodes.tsv
│       ├── genes.tsv
│       └── matrix.mtx
└── GEX
    ├── Directory1
    │   ├── barcodes.tsv
    │   ├── genes.tsv
    │   └── matrix.mtx
    └── Directory2
        ├── barcodes.tsv
        ├── genes.tsv
        └── matrix.mtx
```

where Directory1 and Directory2 would be in the `Directory` column of the csv file.

[Step1](./Step1.Seurat.Demultiplexing.Analysis.R) has been run on all available samples. This R script generates a RDS for each sample. If needed it demultiplexes CellPlex experiments.

[Step2](./Step2.Seurat.Analysis.and.Merging.R) has been run twice. Once to generate a merged Seurat object with all wild-type samples (`samples.used = 1:22`, `desiredName.for.RDS = combined.WT`). Then to generate a merged Seurat object with mutant and wild-type samples (`samples.used = 2:32`, `desiredName.for.RDS = combined.Snai1`).

[Step3 for WT](./Step3.WT.qmd) contains all code used to generate figures relative to the wild-type samples.

[Step3 for Snai1](./Step3.Snai1.qmd) contains all code used to generate figures relative to the mutant vs wild-type samples.

All common functions to both Step3 have been collected into a [single file](./scRNAseqFunctions.R).
