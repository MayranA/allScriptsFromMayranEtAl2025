# Load dependencies:
if (!"devtools" %in% installed.packages()) {
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)

suppressPackageStartupMessages(safelyLoadAPackageInCRANorBioconductor("DESeq2"))
suppressPackageStartupMessages(safelyLoadAPackageInCRANorBioconductor("ggplot2"))

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
path <- "ATACseq/"
tableWithCounts <- file.path(path, "counts_on_peaks.txt.gz")
pathForDESeq2 <- file.path(path, "DESeq2_pairwise")
log2FC.threshold <- 1


# list(factorToStudy = list(loopingVariable = list(subsetting)))
all.analyses <- list(
  "time" = list(
    "time" = list("time" = c("48h", "72h")),
    "time" = list("time" = c("72h", "96h")),
    "time" = list("time" = c("96h", "120h"))
  )
)

# Prepare inputs
count.table <- read.delim(tableWithCounts, check.names = FALSE, quote = "'")
# Because this script is written for RNA-seq:
count.table$gene_id <- with(count.table, paste0(`#chr`, ":", start, "-", end))

samples.plan.df <- data.frame(sample = colnames(count.table)[4:(ncol(count.table) - 1)])
samples.plan.df$time <- sapply(strsplit(samples.plan.df$sample, "_"), head, n = 1)
samples.plan.df$time <- factor(samples.plan.df$time, levels = unique(samples.plan.df$time))
samples.plan.df$rep <- sapply(strsplit(samples.plan.df$sample, "_"), tail, n = 1)
rownames(samples.plan.df) <- samples.plan.df$sample

big.table.fn.long <- "summary_long.txt"

if (!file.exists(file.path(pathForDESeq2, "summary.txt"))) {
  big.annot <- unique(count.table[, c("gene_id", "#chr", "start", "end")])
  rownames(count.table) <- count.table[, "gene_id"]
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
          levels = intersect(
            levels(new.samples.plan[, factorToStudy]),
            unique(new.samples.plan[, factorToStudy])
          )
        )
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
        tail(big.annot)
        all.res.fmt <- subset(
          all.res,
          select = intersect(
            c("gene_id", "gene_short_name", "baseMean", "log2FoldChange", "padj"),
            colnames(all.res)
          )
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

#### General stats ####
signif <- subset(big.annot2, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > log2FC.threshold)
signif$sign <- "up"
signif$sign[signif$log2FoldChange < 0] <- "down"
signif$value <- factor(signif$value, levels = levels(samples.plan.df$time))
signif$comparison.name <- paste0(signif$value, "VS", signif$ref.value)
signif$comparison.name <- factor(signif$comparison.name, levels = unique(signif$comparison.name[order(signif$value)]))

g <- ggplot(signif, aes(x = comparison.name)) +
  geom_bar() +
  geom_text(
    stat = "count", aes(label = after_stat(count)),
    vjust = -1
  ) +
  ggtitle("Number of DE regions") +
  xlab("Comparison") +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(file.path(pathForDESeq2, "DESeq2_nb.pdf"), width = 5, height = 6)

g <- ggplot(signif, aes(x = comparison.name, fill = sign)) +
  geom_bar(position = "dodge") +
  geom_text(
    stat = "count", aes(label = after_stat(count)),
    vjust = -1, position = position_dodge(width = 0.9)
  ) +
  ggtitle("Number of DE regions colored\nby sign of l2fc") +
  xlab("Comparison") +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(file.path(pathForDESeq2, "DESeq2_nb_sign.pdf"), width = 6, height = 8)

#### Output regions to plot ####
df <- subset(signif, value == "96h" & ref.value == "72h")
df <- merge(df, count.table)
# For increased we sort by decreasing signal in 96h
sample <- "96h_WT"
df[, "rep1"] <- df[, paste0(sample, "_rep1")] / sum(df[, paste0(sample, "_rep1")])
df[, "rep2"] <- df[, paste0(sample, "_rep2")] / sum(df[, paste0(sample, "_rep2")])
norm.signal <- apply(
  subset(df, select = c("start", "end", paste0("rep", 1:2))),
  1,
  function(v) {
    width <- v["end"] - v["start"]
    return((v["rep1"] + v["rep2"]) / 2 / width)
  }
)
write.table(subset(df[order(-norm.signal), ], subset = sign == "up", select = c("#chr", "start", "end")),
  file.path(pathForDESeq2, "increased_96h.bed"),
  sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE,
)

# For decreasing we sort by decreasing signal in 72h
sample <- "72h_WT"
df[, "rep1"] <- df[, paste0(sample, "_rep1")] / sum(df[, paste0(sample, "_rep1")])
df[, "rep2"] <- df[, paste0(sample, "_rep2")] / sum(df[, paste0(sample, "_rep2")])
norm.signal <- apply(
  subset(df, select = c("start", "end", paste0("rep", 1:2))),
  1,
  function(v) {
    width <- v["end"] - v["start"]
    return((v["rep1"] + v["rep2"]) / 2 / width)
  }
)
write.table(subset(df[order(-norm.signal), ], subset = sign == "down", select = c("#chr", "start", "end")),
  file.path(pathForDESeq2, "decreased_96h.bed"),
  sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE,
)
