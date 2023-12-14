#' ---
#' title: "Zygotic embryogenesis Biological QA"
#' subtitle: "FMG and ZE at stage B4 and B8"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    fig_width: 9
#'    fig_height: 6
#'    toc: true
#'    number_sections: true
#'    toc_depth: 3
#'    toc_float:
#'      collapsed: TRUE
#'      smooth_scroll: TRUE
#'    code_folding: hide
#'    theme: "flatly"
#'    highlight: pygments
#'    includes:
#'      before_body: header.html
#'      after_body: footer.html
#'    css: style.css
#' ---
#' 
#' <hr />
#' &nbsp;
#' 
#' # Setup
#' This section and the next are relevant for reproducibility purposes. For results, please skip to section 3 (Quality Control)
#' 
#' * Libraries
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(emoji)
  library(gplots)
  library(gtools)
  library(here)
  library(hyperSpec)
  library(limma)
  library(magrittr)
  library(parallel)
  library(patchwork)
  library(PCAtools)
  library(pheatmap)
  library(plotly)
  library(pvclust)
  library(RColorBrewer)
  library(tidyverse)
  library(tximport)
  library(vsn)
})

#' * Helper functions
source(here("UPSCb-common/src/R/featureSelection.R"))

#' * Graphics
hpal <- colorRampPalette(c("blue","white","red"))(100)
pal <- brewer.pal(9,"Blues")

#' # Data
#' * Sample information
samples <- read_tsv(here("doc/samples.tsv"),
                    col_types=cols(.default=col_factor()))

#' Read the expression at the gene level
load(here("data/tximport.rda"))
counts <- txi$counts
colnames(counts) <- samples$SampleID

#' 
#' <hr />
#' &nbsp;
#' 
#' # Quality Control
#' * "Not expressed" genes
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

#' ## Sequencing depth
#' * Let us take a look at the sequencing depth, colouring by Stage and splitting by Tissue
dat <- tibble(x=colnames(counts),y=colSums(counts)) %>% 
  bind_cols(samples)

ggplot(dat,aes(x,y,fill=Stage)) + 
  geom_col() + 
  scale_y_continuous(name="reads") +
  facet_grid(~ factor(Tissue), scales = "free") +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=4),
        axis.title.x=element_blank())

#' `r emoji("point_right")` **We observe +/- 20% difference in the raw sequencing depth**
#' 
#' ## per-gene mean expression
#' 
#' _i.e._ the mean raw count of every gene across samples is calculated
#' and displayed on a log10 scale.
#' 

ggplot(data.frame(value=log10(rowMeans(counts))),aes(x=value)) + 
  geom_density(na.rm = TRUE) +
  ggtitle("gene mean raw counts distribution") +
  scale_x_continuous(name="mean raw counts (log10)") + 
  theme_bw()

#' `r emoji("point_right")` **The cumulative gene coverage is as biased towards low level expression.**
#' 
#' ## Per-sample expression

dat <- as.data.frame(log10(counts)) %>% 
  utils::stack() %>% 
  mutate(Tissue=samples$Tissue[match(ind,samples$SampleID)]) %>% 
  mutate(Stage=samples$Stage[match(ind,samples$SampleID)])

ggplot(dat,aes(x=values,group=ind,col=Stage)) + 
  geom_density(na.rm = TRUE) + 
  ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)") + 
  theme_bw() +
  facet_grid(~ factor(Tissue), scales = "free")

ggplot(dat,aes(x=values,group=ind,col=Tissue)) + 
  geom_density(na.rm = TRUE) + 
  ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)") + 
  theme_bw()

#' `r emoji("point_right")` **All samples have the same sequencing depth characteristics and there is a slight deviation when we look at the tissue type. The FMG has more medium expressed genes and less highly expressed ones tham the ZE.**
#' 
#' * Export raw expression data
dir.create(here("data/analysis/salmon"),showWarnings=FALSE,recursive=TRUE)
write.csv(counts,file=here("data/analysis/salmon/raw-unormalised-gene-expression_data.csv"))
#' 
#' <hr />
#' &nbsp;
#' 
#' # Data normalisation 
#' ## Preparation
#' 
#' For visualization, the data is submitted to a variance stabilization
#' transformation using _DESeq2_. The dispersion is estimated independently
#' of the sample tissue and replicate. 
#'  
#'  ```{r CHANGEME7,eval=FALSE,echo=FALSE}
#'  # In the following, we provide the expected expression model, based on the study design.
#'  # It is technically irrelevant here, as we are only doing the quality assessment of the data, 
#'  # but it does not harm setting it correctly for the differential expression analyses that may follow.
#'  ```
dds <- DESeqDataSetFromTximport(
  txi=txi,
  colData = samples,
  design = ~ Stage * Tissue)

colnames(dds) <- samples$SampleID

save(dds,file=here("data/analysis/salmon/dds.rda"))

#' ## size factors 
#' (_i.e._ the sequencing library size effect)
#' 
dds <- estimateSizeFactors(dds) %>% 
  suppressMessages()

boxplot(normalizationFactors(dds),
        main="Sequencing libraries size factor",
        las=2,log="y")
abline(h=1, col = "Red", lty = 3)

#' and without outliers:
boxplot(normalizationFactors(dds),
        main="Sequencing libraries size factor",
        las=2,log="y", outline=FALSE)
abline(h=1, col = "Red", lty = 3)

#' `r emoji("point_right")` **There is some slight differences in the libraries' size factors. They are all within +/- 60% of the average library size. Some library also show a larger spread around the median, which indicates that genes tend to be quite variable in their expression.**
#' 
#' Assess whether there might be a difference in library size linked to a
#' given metadata
boxplot(split(t(normalizationFactors(dds)),dds$Stage),las=2,
        main="Sequencing libraries size factor by Tissue",
        outline=FALSE,notch=TRUE)

boxplot(split(t(normalizationFactors(dds)),dds$Tissue),las=2,
        main="Sequencing libraries size factor by Tissue",
        outline=FALSE,notch=TRUE)

#' `r emoji("point_right")` **The scaling factor distribution is dependent on both variables. This is potentially a caveat for the DE comparison between tissues.**

plot(colMeans(normalizationFactors(dds)),
     log10(colSums(counts(dds))),ylab="log10 raw depth",
     xlab="average scaling factor",
     col=rainbow(n=nlevels(dds$Tissue))[as.integer(dds$Tissue)],
     pch=c(17,19)[as.integer(dds$Stage)])
legend("bottomright",fill=rainbow(n=nlevels(dds$Tissue)),
       legend=levels(dds$Tissue),cex=0.6)

#' `r emoji("point_right")` **The scaling factor appear linearly proportional to the sequencing depth, but is clearly partitioned by tissue.**
#' 
#' ## Variance Stabilising Transformation
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' ## Validation
#' 
#' let's look at standard deviations before and after VST normalization. 
#' This plot is to see whether there is a dependency of SD on the mean. 
#' 
#' Before:  
meanSdPlot(log1p(counts(dds))[rowSums(counts(dds))>0,])

#' After:
meanSdPlot(vst[rowSums(vst)>0,])

#' After VST normalization, the red line is almost horizontal which indicates
#' no dependency of variance on mean (homoscedastic).
#' 
#' `r emoji("point_right")` **We can conclude that the variance stabilization provides only minor gain. This is something to keep in mind in the rest of the analysis. It is most likely due to the difference between tissues and the amount of lowly expressed genes.**
#' 
#' <hr />
#' &nbsp;
#' 
#' ## QC on the normalised data
#' ### PCA
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' Using PCAtools
p <- pca(vst,colData(dds))

#' ### Scree plot
#' 
#' We define the number of variable of the model:
vars <- all.vars(design(dds))
nvar <- length(vars)
nlevel<-reduce(sapply(vars,function(v){nlevels(eval(parse(text=paste0("dds$",v))))}),`*`)

#' We devise the optimal number of components using two methods
horn <- suppressWarnings(parallelPCA(vst))
elbow <- findElbowPoint(p$variance)

#' We plot the percentage explained by different components and try to empirically assess whether
#' the observed number of components would be in agreement with our model's assumptions.
#' 
#' * the red line represents number of variables in the model  
#' * the orange line represents number of variable combinations.
#' * the black dotted, annotate lines represent the optimal number of components 
#' reported by the horn and elbow methods.
#' 
ggplot(tibble(x=1:length(percent),y=cumsum(percent),p=percent),aes(x=x,y=y)) +
  geom_line() + geom_col(aes(x,p)) + scale_y_continuous("variance explained (%)",limits=c(0,100)) +
  scale_x_continuous("Principal component",breaks=1:length(percent),minor_breaks=NULL) + 
  geom_vline(xintercept=nvar,colour="red",linetype="dashed",linewidth=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nvar],colour="red",linetype="dashed",linewidth=0.5) +
  geom_vline(xintercept=nlevel,colour="orange",linetype="dashed",linewidth=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nlevel],colour="orange",linetype="dashed",linewidth=0.5) +
  geom_vline(xintercept=c(horn$n,elbow),colour="black",linetype="dotted",linewidth=0.5) +
  geom_label(aes(x = horn$n + 1, y = cumsum(percent)[horn$n],label = 'Horn', vjust = 1)) +
  geom_label(aes(x = elbow + 1, y = cumsum(percent)[elbow],label = 'Elbow', vjust = 1))

#' `r emoji("point_right")` **The first component explains 29% of the data variance. Both metrics, Horn and Elbow suggest that two or three components (the dark doted lines) are those that are informative. Indeed the slope of the curve is fairly linear past PC3/PC4 and that would indicate that the remaining PCs only capture sample specific noise. While this is only empirical, the scree plot support having only few variables of importance in the dataset.**
#'
#'
#' ### PCA plot
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    PC3=pc$x[,3],
                    as.data.frame(colData(dds)))

#' * PC1 vs PC2
p1 <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=Tissue,shape=Stage,text=SampleID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p1) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC2 (",percent[2],"%)",sep=""))) %>% suppressWarnings()

#' The same as a biplot
biplot(p,
       colby = 'Tissue',
       colLegendTitle = 'Tissue',
       encircle = TRUE,
       encircleFill = TRUE,
       legendPosition = 'top', 
       legendLabSize = 16, legendIconSize = 8.0)

#' `r emoji("point_right")` **The replicates cluster together, the first dimension separates the Tissue while the second one separates the Satges**
#' 
#' * PC1 vs PC3
p2 <- ggplot(pc.dat,aes(x=PC1,y=PC3,col=Tissue,shape=Stage,text=SampleID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p2) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC3 (",percent[3],"%)",sep=""))) %>% suppressWarnings()

#' The same as a biplot
p$metadata$Condition <- paste(samples$Tissue,samples$Stage)
biplot(p,x = 'PC1', y = 'PC3',
       colby = 'Condition',
       colLegendTitle = 'Condition',
       encircle = TRUE,
       encircleFill = TRUE,
       legendPosition = 'top', 
       legendLabSize = 16, legendIconSize = 8.0)

#' `r emoji("point_right")` **The third dimension is interesting as it separates probably some interaction between Tissue and Stage**
#' 
#' ```{r subplot, out.width = '100%'}
#' subplot(style(p1, showlegend = FALSE), p2,
#'         titleX = TRUE, titleY = TRUE, nrows = 1, margin = c(0.05, 0.05, 0, 0))
#' ```
#'
#' ### Pairs plot
#' This allows for looking at more dimensions, five by default
#' 
suppressMessages(pairsplot(p,colby='Tissue',shape='Stage'))

#' `r emoji("point_right")` **While PC 1 to 3 can be linked to the study design, it looks more like PC4 and 5 are linked to individual samples.**
#' 
#' ### Loadings
#' Loadings, _i.e._ the individual effect of every gene in a component can be studied. Here the most important ones are visualized throughout the different PCs
plotloadings(p,
             rangeRetain = 0.01,
             labSize = 4.0,
             title = 'Loadings plot',
             subtitle = 'PC1 to PC5',
             caption = 'Top 1% variables',
             drawConnectors = TRUE)

#' `r emoji("point_right")` **These are genes that have a strong influence in the different components. Probably interesting candidates to look at, in a non-DE fashion, although they are likely to be DE, given the first three PCs are linked to our design.**
#' 
#' ### Correlation
#' This is a plot showing the correlation between the PC and the model variables. Note that while this might be relevant 
#' for a linear variable, it is less so for categorical variables. Sorting categorical variables in a linear order according to the PCs above might help.
#' 
suppressWarnings(eigencorplot(p,metavars=c('Tissue','Stage')))

#' `r emoji("point_right")` **Keep in mind that the two variables are categorical, so they have been assigned an integer value based on that. They are not truly continuous variables. But clearly as we could see PC1 is strongly linked to Tissue and PC2 to Stage.**
#' 
#' ### Samples Distance
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(sampleDistMatrix) <- paste(dds$Tissue,dds$Stage)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=pal)

#' `r emoji("point_right")` **The sample distance clusters samples as expected**
#' 
#' ## Sequencing depth
#' The figures show the number of genes expressed per condition at different expression cutoffs. The scale on the lower plot is the same as on the upper.
#' The first plot is a heatmap showing the number of genes above a given cutoff. The second plot shows it as a ratio of the number of genes expressed for (a)
#' given variable(s) divided by the average number of genes expressed at that cutoff across all variable(s). The latter plot is of course biased at higher cutoff 
#' as the number of genes becomes smaller and smaller.
#' The point of these two plots is to assert whether the number of genes expressed varies between conditions, as this would break some assumptions for normalisation and differential expression.
conds <- factor(paste(dds$Tissue,dds$Stage))
dev.null <- rangeSamplesSummary(counts=vst,
                                conditions=conds,
                                nrep=3)

#' `r emoji("point_right")` **As expected from the raw data QC, the FMG samples have more genes expressed at a higher value than the ZE. This is most likely indicative of a bias in the number of genes expressed with the FMG expressing less genes overall. The library size correction would then "wrongly" inflate the expression of the genes for the FMG. This is because we are breaking the assumption from DESeq2 that there are the same number of genes expressed across samples.**
#' 
#' Plotting the number of genes that are expressed (at any level)
do.call(rbind,split(t(nrow(vst) - colSums(vst==0)),samples$Tissue)) %>% as.data.frame() %>% 
  rownames_to_column("Tissue") %>% pivot_longer(starts_with("V")) %>% 
  ggplot(aes(x=Tissue, y=value, fill=Tissue)) + geom_dotplot(binaxis = "y", stackdir = "center")

#' `r emoji("point_right")` **As suspected, the FMG has about 10000 less genes expressed than the ZE. This is probably due to the low expressed genes (but not necessarily). Here it would be good to do some stronger filtering on the transcripts we are using for salmon.**
#' 
#' ## Heatmap
#' 
#' Here we want to visualise all the informative genes as a heatmap. We first filter the genes to remove those below the selected noise/signal cutoff. 
#' The method employed here is naive, and relies on observing a sharp decrease in the number of genes within the first few low level of expression. 
#' Using an independent filtering method, such as implemented in DESeq2 would be more accurate, but for the purpose of QA validation, a naive approach is sufficient.
#' Note that a sweet spot for computation is around 20 thousand genes, as building the hierarchical clustering for the heatmap scales non-linearly.
#'
sels <- rangeFeatureSelect(counts=vst,
                           conditions=conds,
                           nrep=3) %>% 
  suppressWarnings()

#' `r emoji("point_right")` **Here a cutoff of 5 is applied (reason is 20000 genes can still be classified by a hierarchical clustering in a decent amount of time)**
vst.cutoff <- 5

nn <- nrow(vst[sels[[vst.cutoff+1]],])
tn <- nrow(vst)
pn <- round(nn * 100/ tn, digits=1)

#' `r emoji("warning")` **`r pn`%** (`r nn`) of total `r tn` genes are plotted below:
#'
#'  
mat <- t(scale(t(vst[sels[[vst.cutoff+1]],])))
hm <- pheatmap(mat,
               color = hpal,
               border_color = NA,
               clustering_distance_rows = "correlation",
               clustering_method = "ward.D2",
               show_rownames = FALSE,
               labels_col = conds,
               angle_col = 90,
               legend = FALSE)

#' `r emoji("point_right")` **The samples cluster as expected, note though that now the stage has a stronger importance than the tissue in the grouping!**
#'
#' ## Clustering of samples
#'
#' Below we assess the previous dendrogram's reproducibility and plot the clusters with au and bp where:
#' 
#' * __au (Approximately Unbiased): computed by multiscale bootstrap resampling__ `r emoji("point_left")` a better approximation
#' * __bp (Bootstrap Probability): computed by normal bootstrap resampling__
#' 
#' `r emoji("warning")`Clusters with AU larger than 95% are highlighted by rectangles, which are strongly supported by data
#' 
hm.pvclust <- pvclust(data = t(scale(t(vst[sels[[vst.cutoff+1]],]))),
                       method.hclust = "ward.D2", 
                       nboot = 30, parallel = TRUE)

plot(hm.pvclust, labels = conds)
pvrect(hm.pvclust)

#' `r emoji("point_right")` **The classification here is similar and place the stage together**
#' 
#' <hr />
#' &nbsp;
#' 
#' # Summary
#' `r emoji("star")` **The data is of good quality**
#' 
#' `r emoji("star")` **The replicates group together**
#' 
#' `r emoji("star")` **The PCA first three components are relevant biologically**
#' 
#' `r emoji("star")` **There is good evidence that there will be DE genes as a result of the selected design**
#' 
#' `r emoji("star")` **The FMG has an average 10000 less genes expressed than the ZE, affecting the library size correction and hence the possible DE results.**
#' 
#' `r emoji("star")` **The aforementioned limitation could be addressed by trying to 1. selecting the transcripts more aggressively so as to remove lowly expressed ones - or putative transposable elements or 2. modelling the expression difference and correcting the library size factor correspondingly**
#' 
#' `r emoji("star")` **A third alternative is to run the DE but imposes more stringent cutoffs on the log2 fold changes or to use the independent filtering to help remove uninformative transcripts from the analysis, prior to rerunning the QA and DE.**
#' 
#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#'
#' # Session Info
#' <details><summary>Session Info</summary>
#' ```{r session info}
#' sessionInfo()
#' ```
#' </details>
#'   
#' &nbsp;
#' 
