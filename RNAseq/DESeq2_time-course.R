# Load dependencies:
if (!"devtools" %in% installed.packages()) {
    install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)

suppressPackageStartupMessages(safelyLoadAPackageInCRANorBioconductor("DESeq2"))
suppressPackageStartupMessages(safelyLoadAPackageInCRANorBioconductor("rtracklayer"))


# Functions:
# Wrapper for DESeq2
deseqAnaWithCovariates <- function(count.table, factorForAna, covariates,
                                   pathOutput, samplesPlan, LRT = FALSE, pvalT = 0.05,
                                   lfcT = 1.5, writeRLOG = FALSE, gene_id = "gene_id", ...) {
    # Checking the conditions
    if (!(factorForAna %in% colnames(samplesPlan))) {
        stop("The factor is is not part of the column names.")
    }
    if (!is.null(covariates) && !(all(unlist(covariates) %in% colnames(samplesPlan)))) {
        stop("Not all covariates are part of the column names.")
    }
    if (length(levels(samplesPlan[, factorForAna])) == 1) {
        stop("The factor you chose have only 1 value. The analysis is not possible.")
    }
    if (length(levels(samplesPlan[, factorForAna])) > 2 && !LRT) {
        print("The factor you chose have more than 2 values. LRT will be applied.")
        LRT <- TRUE
    }
    # Lauching DESeq2
    dds <- DESeqDataSetFromMatrix(
        countData = count.table[, match(rownames(samplesPlan), colnames(count.table))],
        colData = samplesPlan,
        design = as.formula(paste0(
            "~",
            paste(
                c(
                    unlist(covariates),
                    factorForAna
                ),
                collapse = " + "
            )
        ))
    )
    print("Design is:")
    print(design(dds))
    print("Genes that are never expressed are removed")
    dds <- dds[rowSums(counts(dds)) > 1, ]
    if (LRT) {
        # Here I am really not sure about the reduced
        reduced.formula <- as.formula("~1")
        if (!is.null(covariates)) {
            reduced.formula <- as.formula(paste0("~", paste(unlist(covariates), collapse = " + ")))
        }
        dds <- DESeq(dds, minReplicatesForReplace = Inf, test = "LRT", reduced = reduced.formula)
    } else {
        dds <- DESeq(dds, minReplicatesForReplace = Inf)
    }
    res <- results(dds, ...)
    resOrdered <- res[order(res$padj), ]
    # Subsetting the annotation file
    ann <- subset(
        count.table,
        select = intersect(
            colnames(count.table),
            c(gene_id, "gene_id", "gene_short_name", "locus")
        )
    )
    rownames(ann) <- ann[, gene_id]
    resToExport <- data.frame(
        ann[rownames(resOrdered), ],
        counts(dds, normalized = TRUE)[rownames(resOrdered), ],
        resOrdered,
        check.names = FALSE
    )
    if (ncol(ann) == 1) {
        colnames(resToExport)[1] <- colnames(ann)
    }
    write.table(
        resToExport,
        file = paste0(pathOutput, "DESeq2Results.txt"),
        sep = "\t", row.names = FALSE, quote = FALSE
    )
    dfDiffExp <- subset(resToExport, resToExport$padj < pvalT & abs(resToExport$log2FoldChange) > lfcT)
    write.table(
        dfDiffExp,
        file = paste0(pathOutput, "DESeq2significant.txt"),
        sep = "\t", row.names = FALSE, quote = FALSE
    )
    rld <- rlog(dds)
    rlogdata <- assay(rld)
    if (writeRLOG) {
        resToExport2 <- data.frame(
            ann[rownames(resOrdered), ],
            rlogdata[rownames(resOrdered), ],
            check.names = FALSE
        )
        if (ncol(ann) == 1) {
            colnames(resToExport2)[1] <- colnames(ann)
        }
        write.table(
            resToExport2,
            file = paste0(pathOutput, "rlog.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE
        )
    }
    return(invisible(dfDiffExp))
}

# Fixed variables:
path <- "RNAseq/"
pathForDESeq2 <- file.path("output.files", "RNAseq", "DESeq2_pairwise_time-course")
tableWithCounts <- file.path(pathForDESeq2, "..", "mergedTables", "time_course", "AllHTSeqCounts.txt.gz")
gtf.file <- file.path(path, "mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz")
gtf.url <- "https://zenodo.org/record/7510406/files/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz?download=1"
samplesPlan <- file.path(path, "samples.plan.time.course.txt")
log2FC.threshold <- 1
factorToStudy <- "Time"

# Prepare inputs
samples.plan.df <- read.delim(samplesPlan, check.names = FALSE)
rownames(samples.plan.df) <- samples.plan.df$sample
samples.plan.df$Time <- sapply(strsplit(samples.plan.df$sample, "_"), head, n = 1)
samples.plan.df$Time <- factor(samples.plan.df$Time, levels = paste0(sort(as.numeric(gsub("h", "", unique(samples.plan.df$Time)))), "h"))
samples.plan.df$Replicate <- sapply(strsplit(samples.plan.df$sample, "_"), tail, n = 1)

if (!file.exists(gtf.file)) {
    download.file(gtf.url, gtf.file)
}

if (!dir.exists(pathForDESeq2)) {
    dir.create(pathForDESeq2, recursive = TRUE)
}

count.table <- read.delim(tableWithCounts, check.names = FALSE)
colnames(count.table)[1] <- "gene_id"

big.table.fn.long <- "summary_long.txt"

if (!file.exists(file.path(pathForDESeq2, "summary.txt"))) {
    # Prepare a big table with the results of all DESeq2
    gtf <- readGFF(gtf.file)
    big.annot <- unique(gtf[, c("gene_id", "gene_name", "seqid", "gene_biotype")])
    colnames(big.annot) <- c("gene_id", "gene_short_name", "chr", "gene_biotype")
    ##### !  WARNING #### I am using only protein coding genes
    ### And excluding chrM ###

    big.annot <- subset(big.annot, gene_biotype == "protein_coding" & chr != "chrM")
    count.table <- merge(count.table, big.annot)
    rownames(count.table) <- count.table[, "gene_id"]
    rm(gtf)
    big.annot2 <- NULL
    # I will loop over the values of factorToStudy:
    for (increment in 1:2) {
        for (i.ref in seq_along(levels(samples.plan.df[, factorToStudy]))) {
            if (i.ref + increment > length(levels(samples.plan.df[, factorToStudy]))) {
                next
            }
            ref.value <- levels(samples.plan.df[, factorToStudy])[i.ref]
            my.value <- levels(samples.plan.df[, factorToStudy])[i.ref + increment]
            new.samples.plan <- subset(
                samples.plan.df,
                samples.plan.df[, factorToStudy] %in% c(ref.value, my.value)
            )
            # Drop levels for factorToStudy
            new.samples.plan[, factorToStudy] <- factor(new.samples.plan[, factorToStudy],
                levels = intersect(
                    levels(new.samples.plan[, factorToStudy]),
                    unique(new.samples.plan[, factorToStudy])
                )
            )
            base.filename <- paste0(factorToStudy, "_", my.value, "vs", ref.value, "_")

            if (
                length(
                    intersect(
                        new.samples.plan$Batch[new.samples.plan[, factorToStudy] == ref.value],
                        new.samples.plan$Batch[new.samples.plan[, factorToStudy] == my.value]
                    )
                ) == 0
            ) {
                covariates <- NULL
                print("Could not use covariates")
            } else {
                covariates <- list("Batch")
                base.filename <- paste0(base.filename, "_covarBatch_")
            }
            print(base.filename)
            # Run or read DESeq2 results with Wald test with Batch as covariate
            if (!file.exists(file.path(pathForDESeq2, paste0(base.filename, "DESeq2significant.txt")))) {
                print(new.samples.plan)
                deseqAnaWithCovariates(count.table, factorToStudy, covariates,
                    file.path(pathForDESeq2, base.filename),
                    new.samples.plan,
                    LRT = FALSE,
                    lfcT = log2FC.threshold,
                    writeRLOG = FALSE,
                    gene_id = "gene_id"
                )
                # theta = c(0.15, 0.99))
            } else {
                print("Exists")
            }
            # Add results to the dataframe
            all.res <- read.delim(file.path(pathForDESeq2, paste0(base.filename, "DESeq2Results.txt")))
            rownames(all.res) <- all.res[, "gene_id"]
            # Add results to the dataframe
            big.annot[, paste0(base.filename, "l2fc")] <-
                all.res$log2FoldChange[match(big.annot[, "gene_id"], all.res[, "gene_id"])]
            big.annot[, paste0(base.filename, "padj")] <- all.res$padj[match(big.annot[, "gene_id"], all.res[, "gene_id"])]
            big.annot[, paste0(base.filename, "signif")] <-
                with(
                    all.res[match(big.annot[, "gene_id"], all.res[, "gene_id"]), ],
                    !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > log2FC.threshold
                )
            all.res.fmt <- subset(
                all.res,
                select = intersect(
                    c("gene_id", "gene_short_name", "baseMean", "log2FoldChange", "padj"),
                    colnames(all.res)
                )
            )
            all.res.fmt$factor <- factorToStudy
            all.res.fmt$value <- my.value
            all.res.fmt$ref.value <- ref.value
            big.annot2 <- rbind(big.annot2, all.res.fmt)
        }
    }
    write.table(
        big.annot,
        file.path(pathForDESeq2, "summary.txt"),
        sep = "\t", row.names = FALSE, quote = FALSE
    )
    write.table(
        big.annot2,
        file.path(pathForDESeq2, big.table.fn.long),
        sep = "\t", row.names = FALSE, quote = FALSE
    )
}
