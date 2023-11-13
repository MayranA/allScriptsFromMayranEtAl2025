# Load dependencies:
if (!"devtools" %in% installed.packages()) {
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)

suppressPackageStartupMessages(safelyLoadAPackageInCRANorBioconductor("DESeq2"))
suppressPackageStartupMessages(safelyLoadAPackageInCRANorBioconductor("rtracklayer"))
suppressPackageStartupMessages(safelyLoadAPackageInCRANorBioconductor("ggplot2"))
suppressPackageStartupMessages(safelyLoadAPackageInCRANorBioconductor("ggrepel"))


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
  dds <- DESeqDataSetFromMatrix(countData = count.table[, match(rownames(samplesPlan), colnames(count.table))],
                                colData = samplesPlan,
                                design = as.formula(paste0("~",
                                                           paste(c(unlist(covariates),
                                                                   factorForAna),
                                                                 collapse = " + ")
                                )))
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
    select = intersect(colnames(count.table),
                       c(gene_id, "gene_id", "gene_short_name", "locus"))
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
pathForDESeq2 <- file.path("output.files", "RNAseq", "DESeq2_pairwise")
tableWithCounts <- file.path(pathForDESeq2, "..", "mergedTables", "mutants", "AllHTSeqCounts.txt.gz")
gtf.file <- file.path(path, "mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz")
gtf.url <- "https://zenodo.org/record/7510406/files/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz?download=1"
samplesPlan <- file.path(path, "samples.plan.mutants.txt")
log2FC.threshold <- 1


# list(factorToStudy = list(loopingVariable = list(subsetting)))
all.analyses <- list(
  "genotype" = list(
    "genotype" = list("gene" = "Cdh2"),
    "genotype" = list("gene" = "Snai1")
  )
)

# Prepare inputs
samples.plan.df <- read.delim(samplesPlan, check.names = FALSE)
rownames(samples.plan.df) <- samples.plan.df$sample

samples.plan.df$stage_genotype <- sapply(strsplit(samples.plan.df$sample, "_rep"), head, n = 1)
samples.plan.df$genotype <- sapply(strsplit(samples.plan.df$stage_genotype, "[0-9]h_"), tail, n = 1)
samples.plan.df$genotype <- factor(samples.plan.df$genotype,
                                   levels = c(
                                     grep("KO", unique(samples.plan.df$genotype), invert = TRUE, value = TRUE),
                                     grep("KO", unique(samples.plan.df$genotype), value = TRUE)
                                   )
)
samples.plan.df$gene <- "Snai1"
samples.plan.df$gene[grep("Cdh2", samples.plan.df$genotype)] <- "Cdh2"

if (!file.exists(gtf.file)) {
  download.file(gtf.url, gtf.file)
}

count.table <- read.delim(tableWithCounts, check.names = FALSE)
colnames(count.table)[1] <- "gene_id"

big.table.fn.long <- "summary_long.txt"

if (!file.exists(file.path(pathForDESeq2, "summary.txt"))) {
  # Prepare a big table with the results of all DESeq2
  gtf <- readGFF(gtf.file)
  # ! Only use gene_name
  big.annot <- unique(gtf[, c("gene_id", "gene_name", "seqid", "gene_biotype")])
  colnames(big.annot) <- c("gene_id", "gene_short_name", "chr", "gene_biotype")
  ##### !  WARNING #### I am using only protein coding genes
  ### And excluding chrM ###
  
  big.annot <- subset(big.annot, gene_biotype == "protein_coding" & chr != "chrM")
  count.table <- merge(count.table, big.annot)
  rownames(count.table) <- count.table[, "gene_id"]
  rm(gtf)
  big.annot2 <- NULL
  if (!dir.exists(pathForDESeq2)) {
    dir.create(pathForDESeq2, recursive = TRUE)
  }
  # I will loop over the list all.analyses:
  # First is the factorToStudy
  for (factorToStudy in names(all.analyses)) {
    print("FACTOR TO STUDY")
    print(factorToStudy)
    # Second indicates the looping variable
    for (i in 1:length(all.analyses[[factorToStudy]])) {
      loopingVariable <- names(all.analyses[[factorToStudy]][i])
      print("LOOPING VARIABLE")
      print(loopingVariable)
      # Third indicates the subsetting
      pre.samples.plan <- samples.plan.df
      subsetting.list <- all.analyses[[factorToStudy]][i][[loopingVariable]]
      subsetting.name <- ""
      for (rn in names(subsetting.list)) {
        print("SUBSET")
        print(rn)
        pre.samples.plan <- pre.samples.plan[pre.samples.plan[, rn] %in% subsetting.list[[rn]], ]
        subsetting.name <- paste0(subsetting.name, paste(subsetting.list[[rn]], collapse = "or"), "_")
      }
      looping.values <- unique(pre.samples.plan[, loopingVariable])
      ref.value <- NULL
      if (factorToStudy == loopingVariable) {
        ref.value <- intersect(levels(pre.samples.plan[, loopingVariable]), pre.samples.plan[, loopingVariable])[1]
        looping.values <- setdiff(looping.values, ref.value)
      }
      for (my.value in looping.values) {
        new.samples.plan <- pre.samples.plan[pre.samples.plan[, loopingVariable] %in% c(ref.value, my.value), ]
        # Drop levels for factorToStudy
        new.samples.plan[, factorToStudy] <- factor(new.samples.plan[, factorToStudy],
                                                    levels = intersect(levels(new.samples.plan[, factorToStudy]),
                                                                       unique(new.samples.plan[, factorToStudy])))
        base.filename <- paste0(factorToStudy, "_", subsetting.name, my.value)
        if (factorToStudy == loopingVariable) {
          base.filename <- paste0(base.filename, "vs", ref.value, "_")
        } else {
          base.filename <- paste0(base.filename, "_")
        }
        print(base.filename)
        # Run or read DESeq2 results with Wald test threshold of FC at 1.5
        if (!file.exists(file.path(pathForDESeq2, paste0(base.filename, "DESeq2significant.txt")))) {
          print(new.samples.plan)
          deseqAnaWithCovariates(count.table, factorToStudy, NULL,
                                 file.path(pathForDESeq2, base.filename),
                                 new.samples.plan,
                                 LRT = FALSE,
                                 lfcT = log2FC.threshold,
                                 writeRLOG = FALSE,
                                 gene_id = "gene_id")
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
          with(all.res[match(big.annot[, "gene_id"], all.res[, "gene_id"]), ],
               !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > log2FC.threshold)
        tail(big.annot)
        all.res.fmt <- subset(
          all.res,
          select = intersect(c("gene_id", "gene_short_name", "baseMean", "log2FoldChange", "padj"),
                             colnames(all.res))
        )
        all.res.fmt$factor <- factorToStudy
        all.res.fmt$subsetting <- subsetting.name
        all.res.fmt$value <- my.value
        all.res.fmt$ref.value <- ref.value
        big.annot2 <- rbind(big.annot2, all.res.fmt)
        tail(big.annot2)
      }
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
} else {
  big.annot <- read.delim(file.path(pathForDESeq2, "summary.txt"), check.names = FALSE)
  big.annot2 <- read.delim(file.path(pathForDESeq2, big.table.fn.long))
}

## Volcano plots
fig.path <- "output.files/RNAseq/mutants/"
dir.create(fig.path, showWarnings = FALSE, recursive = TRUE)
figs.width <- 5
figs.heights <- list("Cdh2Het" = 19 * figs.width / 14, "WT" = 34 * figs.width / 17)
# Default is 0.5
# Increase if labels overlap
max.time <- 10
genes.to.highlight.if.signif <- c("Cdh2", "Cdh1",
                                  "Snai1", "Pou5f1", "Dppa5a", "Foxc1",
                                  "Tbx6", "Pcdh19", "Foxc2")

conds.summary <- unique(big.annot2[, c("ref.value", "value")])
for (i in 1:nrow(conds.summary)) {
  cur.ref.value <- conds.summary$ref.value[i]
  cur.value <- conds.summary$value[i]
  title <- paste0(cur.value, "VS", cur.ref.value)
  df <- subset(big.annot2, ref.value == cur.ref.value & value == cur.value)
  df$signif <- !is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) > log2FC.threshold
  df$diffexpressed <- "Not significant"
  df$diffexpressed[df$signif & df$log2FoldChange > 0] <- "Upregulated"
  df$diffexpressed[df$signif & df$log2FoldChange < 0] <- "Downregulated"
  df$diffexpressed <- factor(df$diffexpressed, levels = c("Downregulated", "Not significant", "Upregulated"))
  df$delabel <- ""
  df$delabel[df$signif & df$gene_short_name %in% genes.to.highlight.if.signif] <-
    df$gene_short_name[df$signif & df$gene_short_name %in% genes.to.highlight.if.signif]
  # Reorder:
  df <- rbind(subset(df, !signif), subset(df, signif))
  # From https://biostatsquid.com/volcano-plots-r-tutorial/
  g <- ggplot(data = df, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = delabel)) +
    geom_vline(xintercept = c(-log2FC.threshold, log2FC.threshold), col = "gray", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = "dashed") +
    geom_point(size = 2) +
    scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00")) + # to set the colours of our variable
    scale_y_continuous(expand = c(0, NA)) +
    labs(
      color = "", # color legend_title
      x = expression("log"[2] * "FC"),
      y = expression("-log"[10] * "adjusted p-value")
    ) +
    geom_text_repel(
      seed = 1,
      max.time = max.time,
      max.overlaps = Inf
    ) +
    theme_classic() +
    theme(legend.position = "bottom")
  ggsave(file.path(fig.path, paste0(title, ".pdf")),
         width = figs.width, height = figs.heights[[cur.ref.value]]
  )
}

## Make a short selection of genes to compare with single cell analysis:
big.annot$AverageSnai1log2FC <- with(big.annot,
                                     (genotype_Snai1_Snai1KO_Clone1vsWT_l2fc +
                                        genotype_Snai1_Snai1KO_Clone2vsWT_l2fc) /2)
big.annot$ConcordentSnai1log2FCSign <- with(big.annot,
                                            sign(genotype_Snai1_Snai1KO_Clone1vsWT_l2fc) ==
                                              sign(genotype_Snai1_Snai1KO_Clone2vsWT_l2fc))
big.annot$MaxpValueSnai1 <- with(big.annot,
                                 pmax(genotype_Snai1_Snai1KO_Clone1vsWT_padj,
                                      genotype_Snai1_Snai1KO_Clone2vsWT_padj))
simplifiedTable <- big.annot[, c("gene_short_name", "AverageSnai1log2FC", "ConcordentSnai1log2FCSign", "MaxpValueSnai1",
                                 paste0("genotype_Snai1_Snai1KO_Clone", c(1, 2), "vsWT_l2fc"))]
colnames(simplifiedTable)[1] <- "gene"
colnames(simplifiedTable)[ncol(simplifiedTable) - (1:0)] <- paste0("Clone", c(1, 2))

filteredSimplifiedTable <- subset(simplifiedTable,
                                  !is.na(AverageSnai1log2FC) &
                                    !is.na(MaxpValueSnai1) &
                                    abs(AverageSnai1log2FC) > 2 &
                                    ConcordentSnai1log2FCSign &
                                    MaxpValueSnai1 < 0.00005)

write.table(
  filteredSimplifiedTable,
  file = file.path(pathForDESeq2, "filteredSimplified.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)


writeLines(capture.output(sessionInfo()), file.path(path, "sessionInfo.txt"))
