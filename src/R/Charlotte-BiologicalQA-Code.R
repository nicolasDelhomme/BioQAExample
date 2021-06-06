# temporary fix (RSQlite - RStudio issue)
options(connectionObserver = NULL)

library(here)
datadir <- here("data")

meta <- read.delim(here("data/airway/airway_meta.txt"), 
                   header = TRUE, as.is = TRUE)
rownames(meta) <- meta$names

meta$dex <- factor(meta$dex)
meta$cell <- factor(meta$cell)

meta

# Do NOT run
# library(Rsubread)
# fc <- featureCounts(files = files, 
#                     annot.ext = gtf, 
#                     isGTFAnnotationFile = TRUE,
#                     GTF.featureType = "exon", 
#                     GTF.attrType = "gene_id", 
#                     useMetaFeatures = TRUE, 
#                     strandSpecific = 0, 
#                     isPairedEnd = TRUE, 
#                     nthreads = 6)

fc <- readRDS(here("data/airway/featureCounts/star_featureCounts.rds"))
names(fc)
counts_featurecounts <- fc$counts
head(counts_featurecounts)
dim(counts_featurecounts)
fc$stat

suppressPackageStartupMessages({
  library(tximeta)
  library(DESeq2)
  library(org.Hs.eg.db)
  library(SummarizedExperiment)
})

## List all quant.sf output files from Salmon
salmonfiles <- here("data/airway/salmon/", meta$names, "/quant.sf")
names(salmonfiles) <- meta$names
stopifnot(all(file.exists(salmonfiles)))

## Add a column "files" to the metadata table. This table must contain at least
## two columns: "names" and "files"
coldata <- cbind(meta, files = salmonfiles, stringsAsFactors = FALSE)

## Import quantifications on the transcript level
st <- tximeta::tximeta(coldata)

## Summarize quantifications on the gene level
sg <- tximeta::summarizeToGene(st)

## Add gene symbols
sg <- tximeta::addIds(sg, "SYMBOL", gene = TRUE)
sg

counts_salmon <- round(assay(sg, "counts"))

spl <- "SRR1039508"
gns <- rownames(counts_salmon)
quants <- data.frame(featureCounts = counts_featurecounts[gns, spl],
                     salmon = counts_salmon[gns, spl])
pairs(quants)

colData(sg)

colData(sg)$dex

colData(sg)$cell 

colData(sg)$dex <- relevel(colData(sg)$dex, ref = "untrt")
colData(sg)$dex

ds_se <- DESeqDataSet(sg, design = ~ cell + dex)

stopifnot(all(colnames(counts_salmon) == rownames(meta)))
meta$dex <- relevel(meta$dex, ref = "untrt")
ds_matrix <- DESeqDataSetFromMatrix(countData = counts_salmon, 
                                    colData = meta,
                                    design = ~ cell + dex)

suppressPackageStartupMessages({
  library(edgeR)
})
genetable <- data.frame(gene.id = rownames(counts_salmon),
                        stringsAsFactors = FALSE)
stopifnot(all(rownames(meta) == colnames(counts_salmon)))
dge <- DGEList(counts = counts_salmon, 
               samples = meta, 
               genes = genetable)
names(dge)

avetxlengths <- assay(sg, "length")
stopifnot(all(rownames(avetxlengths) == rownames(counts_salmon)))
stopifnot(all(colnames(avetxlengths) == colnames(counts_salmon)))
avetxlengths <- avetxlengths/exp(rowMeans(log(avetxlengths)))
offsets <- log(calcNormFactors(counts_salmon/avetxlengths)) + 
  log(colSums(counts_salmon/avetxlengths))
dge <- scaleOffset(dge, t(t(log(avetxlengths)) + offsets))
names(dge)

dge <- edgeR::calcNormFactors(dge)
dge$samples

vsd <- DESeq2::vst(ds_se)

class(vsd)

head(colData(vsd), 3)

DESeq2::plotPCA(vsd, intgroup = "cell")
DESeq2::plotPCA(vsd, intgroup = "dex")

plotMDS(dge, top = 500, labels = NULL, col = as.numeric(dge$samples$dex), 
        pch = as.numeric(dge$samples$cell), cex = 2, gene.selection = "common")
