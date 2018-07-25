source("https://bioconductor.org/biocLite.R")
biocLite(c("Biobase","oligoClasses","knitr","BiocStyle","oligo","geneplotter","ggplot2","dplyr","LSD","gplots","RColorBrewer","ArrayExpress","arrayQualityMetrics","stringr","matrixStats","topGO","genefilter","pd.hugene.1.0.st.v1","hugene10sttranscriptcluster.db","pheatmap","mvtnorm","DAAG","multcomp","limma","ReactomePA","clusterProfiler","openxlsx","devtools","biomaRt","EnrichmentBrowser"))

library(Biobase)
library(oligoClasses)
library(knitr)
library(BiocStyle)
library(oligo)
library(geneplotter)
library(ggplot2)
library(dplyr)
library(LSD)
library(gplots)
library(RColorBrewer)
library(ArrayExpress)
library(arrayQualityMetrics)
library(stringr)
library(matrixStats)
library(topGO)
library(genefilter)
library(pd.hugene.1.0.st.v1)
library(hugene10sttranscriptcluster.db)
library(pheatmap)
library(mvtnorm)
library(DAAG)
library(multcomp)
library(limma)
library(ReactomePA)
library(clusterProfiler)
library(openxlsx)
library(devtools)
library(biomaRt)
library(EnrichmentBrowser)

#Create repository for raw data
if(!dir.exists("raw_data_dir")){
  dir.create("raw_data_dir")
}

#Get datasets (auto)
anno_AE <- getAE("E-MTAB-2967", path="raw_data_dir", type="raw")

#Read SDRF file
#Manual

SDRFtxt <- "E-MTAB-2967.sdrf.txt"
path = "/home/zg267/Documents/Microarray analysis/Affymetrix microarray workflow paper/raw_data_dir/"
SDRF <- read.delim(file=file.path(path,SDRFtxt))
#Auto
SDRF <- read.delim(url("http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2967/E-MTAB-2967.sdrf.txt"))

#Annotate data
rownames(SDRF)<- SDRF$Array.Data.File
SDRF <- AnnotatedDataFrame(SDRF)

#Import CEL files as raw_data
raw_data <- read.celfiles(filenames = file.path("raw_data_dir",SDRF$Array.Data.File),verbose =FALSE, phenoData =SDRF)
validObject(raw_data)
head(pData(raw_data))
head(exprs(raw_data))
stopifnot(validObject(raw_data))
#Retrieve phenotype
pData(raw_data) <- pData(raw_data)[, c("Source.Name",
                                       "Characteristics.individual.",
                                       "Factor.Value.phenotype.",
                                       "Factor.Value.disease."
                                       )]

#Quality control of raw data
#Get log of raw data
exp_raw <- log2(exprs(raw_data))
PCA_raw <- prcomp(t(exp_raw), scale = FALSE)

dataGG <- data.frame(PC1 = PCA_raw$x[,1], 
                     PC2 = PCA_raw$x[,2],
                     Disease = pData(raw_data)$Factor.Value.disease.,
                     Phenotype = pData(raw_data)$Factor.Value.phenotype.,
                     Individual = pData(raw_data)$Characteristics.individual.)

qplot(PC1, PC2, data = dataGG, color = Disease,
      main = "PCA plot of the raw data (log-transformed)", size = I(2),
      asp = 1.0, geom = "text",
      label = Individual)
+ scale_color_brewer(palette = "set2")

boxplot(raw_data, target = "core",
        main = "Boxplots of log2_intensities for the raw data")

#Produce array quality metrics???
ig <- c("Factor.Value.disease.","Factor.Value.phenotype.")
arrayQualityMetrics(expressionset = raw_data,
                    outdir = "Report_for_Palmieri_raw",
                    force = TRUE, do.logtransform = TRUE,
                    intgroup = ig)

#Background adjustment
#Across-array normalization
#Summarization
#One-go preprocessing in oligo (all three step)(RMA)
palmieri_eset <- oligo::rma(raw_data, target="core")


#Generate a heatmap of the sample-to-sample distances
exp_palmieri <- exprs(palmieri_eset)
PCA <- prcomp(t(exp_palmieri), scale = FALSE)

dataGG2 <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                      Disease = pData(palmieri_eset)$Factor.Value.disease.,
                      Phenotype = pData(palmieri_eset)$Factor.Value.phenotype.)

qplot(PC1, PC2, data = dataGG2, color = Disease, shape = Phenotype,
      main = "PCA plot of the calibrated data", size = I(2), asp = 1.0) +scale_color_brewer(palette = "Set2")

dists <- as.matrix(dist(t(exp_palmieri), method = "manhattan"))
colnames(dists) <- NULL
diag(dists) <- NA
rownames(dists) <- pData(palmieri_eset)$Factor.Value.phenotype.
hmcol <- colorRampPalette(rev(brewer.pal(9, "PuOr")))(255)
pheatmap(dists, col = rev(hmcol),clustering_distance_rows = "manhattan",
                       clustering_distance_cols = "manhattan")

#Filter out lowly-expressed genes
#Get number of samples
no_of_samples <- table(paste0(pData(palmieri_eset)$Factor.Value.disease., "_",
                              pData(palmieri_eset)$Factor.Value.phenotype.))
no_of_samples

#Median intensities 
#???
palmieri_medians <- rowMedians(exprs(palmieri_eset))

hist_res <- hist(palmieri_medians, 100, col = "#e7efd8", freq = FALSE,
                 main = "Histogram of the median intensities",
                 xlab = "Median intensities")

emp_mu <- hist_res$breaks[which.max(hist_res$density)]
emp_sd <- mad(palmieri_medians)/2
prop_cental <- 0.50

lines(sort(palmieri_medians), prop_cental * dnorm(sort(palmieri_medians), mean = emp_mu, sd = emp_sd),col = "grey10", lwd = 4)
      
      
cut_val <- 0.05 / prop_cental
thresh_median <- qnorm(0.05 / prop_cental, emp_mu, emp_sd)

samples_cutoff <- min(no_of_samples)

idx_thresh_median <- apply(exprs(palmieri_eset), 1, function(x){
                                 sum(x > thresh_median) >= samples_cutoff})
table(idx_thresh_median)
palmieri_filtered <- subset(palmieri_eset, idx_thresh_median)

#Annotate transcript clusters
#????
anno_palmieri <- AnnotationDbi::select(hugene10sttranscriptcluster.db,
                                       keys = (featureNames(palmieri_filtered)),
                                       columns = c("SYMBOL", "GENENAME"),
                                       keytype = "PROBEID")

#Remove multiple mappings & build custom annotations
probe_stats <- anno_palmieri %>%
  group_by(PROBEID) %>%
  summarise(no_of_matches = n_distinct(SYMBOL)) %>%
  filter(no_of_matches > 1)
probe_stats

#Exclude transcript-clusters
ids_to_exclude <- ((featureNames(palmieri_filtered) %in% probe_stats$PROBEID) |
                     featureNames(palmieri_filtered) %in% subset(anno_palmieri ,
                                                                 is.na(SYMBOL))$PROBEID)
table(ids_to_exclude)
palmieri_final <- subset(palmieri_filtered, !ids_to_exclude)
validObject(palmieri_final)
fData(palmieri_final)$PROBEID <- rownames(fData(palmieri_final))
fData(palmieri_final) <- left_join(fData(palmieri_final), anno_palmieri)

#Restore rownames after left join
rownames(fData(palmieri_final)) <- fData(palmieri_final)$PROBEID
validObject(palmieri_final)


#Fit a linear model
individual <- as.character(pData(palmieri_final)$Characteristics.individual.)
#Group by tissue types
tissue <- str_replace_all(pData(palmieri_final)$Factor.Value.phenotype., " ", "_")
tissue <- ifelse(tissue == "non-inflamed_colonic_mucosa", "nI", "I")
#Group by disease
disease <- str_replace_all(pData(palmieri_final)$Factor.Value.disease., " ", "_")
disease <- ifelse(disease == "Crohn's_disease", "CD", "UC")
#Design matrices
i <- individual[disease == "CD"]
design_palmieri_CD <- model.matrix( ~ 0 + tissue[disease == "CD"] + i)
colnames(design_palmieri_CD)[1:2] <- c("I", "nI")

j<- individual[disease == "UC"]
design_palmieri_UC <- model.matrix( ~ 0 + tissue[disease == "UC"] + j)
colnames(design_palmieri_UC)[1:2] <- c("I", "nI")

head(design_palmieri_CD)
dim(design_palmieri_CD)
min(svd(design_palmieri_CD)$d)

head(design_palmieri_UC)
dim(design_palmieri_UC)
min(svd(design_palmieri_UC)$d)

#Contrasts & fit the linear model
contrast_matrix_CD <- makeContrasts(I-nI, levels = design_palmieri_CD)

palmieri_fit_CD <- eBayes(contrasts.fit(lmFit(palmieri_final[,disease == "CD"],
                                              design = design_palmieri_CD),
                                        contrast_matrix_CD))

contrast_matrix_UC <- makeContrasts(I-nI, levels = design_palmieri_UC)

palmieri_fit_UC <- eBayes(contrasts.fit(lmFit(palmieri_final[,disease == "UC"],
                                              design = design_palmieri_UC),
                                        contrast_matrix_UC))

#Extract results
table_CD <- topTable(palmieri_fit_CD, number = Inf)
head(table_CD)
table(table_CD$adj.P.Val < 0.05)
table(table_CD$P.Value < 0.001)

table_UC <- topTable(palmieri_fit_UC, number = Inf)
head(table_UC)
table(table_UC$adj.P.Val < 0.05)
table(table_UC$P.Value < 0.001)

#Plot p-values
hist(table_CD$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "inflammed vs non-imflamed - Crohn's disease", xlab = "p-values")

hist(table_UC$P.Value, col = brewer.pal(3, name = "Set2")[3],
     main = "inflammed vs non-imflamed - Ulcerative colitis", xlab = "p-values")

#GO enrichment analysis
#Match the background sets of genes
DE_genes_CD <- subset(table_CD, adj.P.Val < 0.1)$PROBEID
back_genes_idx <- genefinder(palmieri_final, as.character(DE_genes_CD),
                             method = "manhattan", scale = "none")
back_genes_idx <- sapply(back_genes_idx, function(x) x$indices)
back_genes <- featureNames(palmieri_final)[back_genes_idx]
back_genes <- setdiff(back_genes, DE_genes_CD)

intersect(back_genes, DE_genes_CD)
length(back_genes)
multidensity(list(
  all= table_CD[,"AveExpr"],
  fore= table_CD[DE_genes_CD, "AveExpr"],
  back= table_CD[rownames(table_CD) %in% back_genes, "AveExpr"]),
  col = c("#e46981","#ae7ee2","#a7ad4a"),
  xlab = "mean expression",
  main = "DE genes for CD - background - matching"
  )

#Run topGO
gene_IDs <- rownames(table_CD)
in_universe <- gene_IDs %in% c(DE_genes_CD, back_genes)
inSelection <- gene_IDs %in% DE_genes_CD
all_genes <- factor(as.integer(inSelection[in_universe]))
names(all_genes) <- gene_IDs[in_universe]

#Initialize topGO
ont <- "BP"

top_GO_data <- new("topGOdata", ontology = ont, allGenes = all_genes,
                   nodeSize = 10, annot=annFUN.db, affyLib = "hugene10sttranscriptcluster.db")

#Test topGOdata
result_top_GO_elim <- runTest(top_GO_data, algorithm = "elim", statistic = "Fisher")
result_top_GO_classic <- runTest(top_GO_data, algorithm = "classic", statistic = "Fisher")

#Look at top 100 GO categories
res_top_GO <- GenTable(top_GO_data, Fisher.elim = result_top_GO_elim,
                       Fisher.classis = result_top_GO_classic,
                       orderBy = "Fisher.elim", topNodes = 100)
genes_top_GO <- printGenes(top_GO_data, whichTerms = res_top_GO$GO.ID,
                           chip = "hugene10sttranscriptcluster.db", geneCutOff = 1000)
res_top_GO$sig_genes <- sapply(genes_top_GO, function(x) {
  str_c(paste0(x[x$'raw p-value' == 2, "Symbol.id"],";"), collapse = "")
})

head(res_top_GO[,1:8], 20)

#Visualize GO analysis results
par(cex = 0.0315)
showSigOfNodes(top_GO_data, score(result_top_GO_elim), firstSigNodes = 3,
               useInfo = 'def')
showSigOfNodes(top_GO_data, score(result_top_GO_classic), firstSigNodes = 3,
               useInfo = 'def')

printGraph(top_GO_data,result_top_GO_elim,firstSigNodes = 3)
printGraph(top_GO_data,result_top_GO_classic,firstSigNodes = 3)

#Pathway enrichment analysis (reactome)
#Map PROBEIDS to entrez 
entrez_ids <- mapIds(hugene10sttranscriptcluster.db,
                     keys = rownames(table_CD),
                     keytype = "PROBEID",
                     column = "ENTREZID")

reactome_enrich <- enrichPathway(gene = entrez_ids[DE_genes_CD],
                                 universe = entrez_ids[c(DE_genes_CD,
                                                         back_genes)],
                                 organism = "human",
                                 pvalueCutoff = 0.01,
                                 qvalueCutoff = 0.8,
                                 readable = TRUE)
reactome_enrich@result$Description <- paste0(str_sub(reactome_enrich@result$Description,
                                                     1, 20),
                                             "...")
head(summary(reactome_enrich), n = 30)[1:6]

#Visualize the results
barplot(reactome_enrich)
enrichMap(reactome_enrich, n = 10, vertex.label.font =2)
emapplot(reactome_enrich, showCategory = 40, color = "p.adjust")

#Save the workspace
save.image(file = "paper_test.RData")
