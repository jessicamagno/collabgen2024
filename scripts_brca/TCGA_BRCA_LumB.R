##########################################################

# ----- TNBC BRCA LumB

# Required packages
library(TCGAbiolinks)
library(SummarizedExperiment)
library(AnnotationHub)
library(DESeq2)
library(ComplexHeatmap)
library(cowplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(org.Hs.eg.db)

# Download GDC BRCA clinical data for TNBC samples
clinical_gdc <- read.delim("ClinicalDataGDC/brca_tcga_pan_can_atlas_2018_clinical_data.tsv",
                           header = T, sep = "\t", dec = ".")

# Subset TNBC samples
LumB_samples <- subset(clinical_gdc, Subtype == "BRCA_LumB")

LumB_barcodes <- LumB_samples$Patient.ID

# Download gene expression data from TCGAbiolinks
query = GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor",
  barcode = LumB_barcodes)

GDCdownload(query)
tcga_brca = GDCprepare(query)
dim(tcga_brca) #60660 197
colnames(colData(tcga_brca))

gexp_LumB <- assay(tcga_brca)

ens <- rownames(gexp_LumB)
ens <- sub('\\.[0-9]*$', '', ens)

rownames(gexp_LumB) <- ens

# Save preprocessed data
save(gexp_LumB, file = "tcgaBRCA_LumB_140624.RData")
save(clinical_gdc, file = "tcgaBRCA_LumB_clinical_140624.RData")

##################################################################

# Load preprocessed data
load("tcgaBRCA_LumB_140624.RData")

# Filtering low counts by edgeR
library(edgeR)
x <- gexp_LumB
keep_x <- filterByExpr(x, min.count = 10)
x <- x[keep_x,] #18133 197
gexp_LumB <- x

# Save data after filtering with edgeR
save(gexp_LumB, file = "tcgaBRCA_LumB_preprocessed.RData")

# Load clinical data
load("tcgaBRCA_LumB_clinical_140624.RData")

samples_BRCA_LumB <- clinical_gdc$Sample.ID 
counts_sample <- colnames(gexp_LumB)
counts_sample <- substr(counts_sample, 1,15)
colnames(gexp_LumB) <- counts_sample
intersect(samples_BRCA_LumB, counts_sample)
length(intersect(samples_BRCA_LumB, counts_sample)) #197

#-- Filtering the samples we have in gene expression matrix
clinical_BRCA_LumB <- clinical_gdc[match(counts_sample,
                                        clinical_gdc$Sample.ID),]
#-- Setting rownames in clinical data
rownames(clinical_BRCA_LumB)
rownames(clinical_BRCA_LumB) <- clinical_BRCA_LumB$Sample.ID

identical(rownames(clinical_BRCA_LumB), colnames(gexp_LumB))

#-- Save counts and clinical data
save(gexp_LumB, file = "BRCA_LumB_final_counts.RData")
save(clinical_BRCA_LumB, file = "BRCA_LumB_final_clin.RData")

#-------------------------------------------------------------

#-- Loading data 
load("BRCA_LumB_final_counts.RData")
load("BRCA_LumB_final_clin.RData")

#-- Checking if counts and clinical data are in the same order
identical(colnames(gexp_LumB), rownames(clinical_BRCA_LumB)) #TRUE

# Normalize data with DESEq2
## Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = gexp_LumB, 
                              colData = clinical_BRCA_LumB,
                              design = ~ 1)
View(counts(dds))
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts_BRCA_LumB <- counts(dds, normalized=TRUE)


################################################################

# Gene expression analysis of PD1 and PDL1
# Selecting ens for these two genes

# for PD1
which(rownames(normalized_counts_BRCA_LumB) == "ENSG00000188389") #12896
gene_exp_PD1 <- normalized_counts_BRCA_LumB[12896,]
gene_exp_PD1 <- as.numeric(gene_exp_PD1)
summary(gene_exp_PD1)

# Adding PD1 expression values as a clinical variable
clinical_BRCA_LumB$PD1 <- gene_exp_PD1


g_PD1 <- ggplot(clinical_BRCA_LumB, aes(x = PD1))
g_PD1 <- g_PD1 + geom_histogram(fill = "lightgreen", color = "white", bins = 20) +
  labs(title = "PD1 expression in BRCA LumB",
       x = "PD1 Expression",
       y = "Frequency") +
  theme_classic() # distribuicao nao normal
g_PD1

# The same for PDL1
# for PDL1
which(rownames(normalized_counts_BRCA_LumB) == "ENSG00000120217") #4371
gene_exp_PDL1 <- normalized_counts_BRCA_LumB[4371,]
gene_exp_PDL1 <- as.numeric(gene_exp_PDL1)
summary(gene_exp_PDL1)

# Adding PDL1 expression values as a clinical variable
clinical_BRCA_LumB$PDL1 <- gene_exp_PDL1

g_PDL1 <- ggplot(clinical_BRCA_LumB, aes(x = PDL1))
g_PDL1 <- g_PDL1 + geom_histogram(fill = "lightpink", color = "white", bins = 20) +
  labs(title = "PDL1 expression in BRCA LumB",
       x = "PDL1 Expression",
       y = "Frequency") +
  theme_classic() # distribuicao nao normal
g_PDL1

# Visualizing both genes
g_PD1_PDL1 <- ggplot(clinical_BRCA_LumB,
                aes(x = PD1, y = PDL1, color = Overall.Survival.Status)) +
  geom_point() + scale_color_manual(values = c("steelblue", "salmon")) +
  theme_classic() +
  labs(title = " PD1 Expression vs PDL1 Expression in BRCA LumB",
       x = "PD1",
       y = "PDL1",
       color = "Overall Survival Status")
g_PD1_PDL1

g_2 <- ggplot(clinical_BRCA_LumB,
                  aes(x = PD1, y = PDL1, fill = Overall.Survival.Status))
g_2 <- g_2 + geom_boxplot() + 
  scale_fill_manual(values = c("salmon", "steelblue")) + 
  theme_minimal() +
  labs("PD1 expression vs PDL1 expression in BRCA LumB",
       x = "PD1 Expression",
       y = "PDL1 Expression",
       fill = "Overall Survival")
g_2


group_PD1_PDL1 <- plot_grid(g_PD1, g_PDL1, g_PD1_PDL1, labels = "AUTO", 
                             nrow = 2,
                             rel_widths = c(2,2,4))

# Exemplo de exportacao da figura em pdf
ggsave2(filename = "PD1_PDL1_BRCA_LumB.pdf", group_PD1_PDL1,
        width = 12, height = 7, units = "in")


###################################################################
########## CIBERSORT ANALYSIS ##################

identical(colnames(normalized_counts_BRCA_LumB), rownames(clinical_BRCA_LumB))

library(IOBR)

eset_brca <- count2tpm(countMat = normalized_counts_BRCA_LumB, idType = "ENSEMBL")
eset_brca_log2 <- log2(eset_brca + 1)

save(eset_brca_log2, file = "eset_BRCA_log2_norm.RData")
load("eset_BRCA_log2.RData")

cibersort_brca_lumB <- deconvo_tme(eset = eset_brca_log2, method = "cibersort", 
                                   arrays = FALSE, perm = 200,
                                   absolute.mode = F)
save(cibersort_brca_lumB, file = "BRCA_LumB_cs_norm.RData")
load("BRCA_LumB_cibersort2406.RData")
View(cibersort_brca_lumB)

# data manipulation for cell bar plot
colnames(cibersort_brca_lumB)[1] <- "Samples"

cols <- hcl.colors(21, palette = "Spectral")
res_BRCA_LumB_norm <- my_cell_bar_plot(input = cibersort_brca_lumB,
                                      title = "CIBERSORT Immune cell fractions",
                                      features = colnames(cibersort_brca_lumB)[-c(1,6,24:26)],
                                      coord_filp = F,
                                      legend.position = "bottom",
                                      cols = cols)



# Getting macrophages values and separating in new data
macrophages_lumB <- cibersort_brca_lumB[,c(1,15:17)]
macrophages_lumB <- as.data.frame(macrophages_lumB)
rownames(macrophages_lumB) <- macrophages_lumB$Samples
# Ordering samples by M2 values
macrophages_lumB_2 <- macrophages_lumB[order(macrophages_lumB$Macrophages_M2_CIBERSORT),]

res_BRCA_LumB_mac <- my_cell_bar_plot(input = macrophages_lumB_2,
                                      title = "CIBERSORT Immune cell fractions",
                                      features = colnames(macrophages_lumB_2)[-1],
                                      coord_filp = F,
                                      legend.position = "bottom",
                                      cols = c("lightblue", "steelblue",
                                               "darkblue"))


# Getting macrophages, linfocytes and neutrophils
macr_linf_neut <- cibersort_brca_lumB[,c(1,5,7,8,15:17,23)]
cols <- hcl.colors(7, palette = "Temps")

res_BRCA_LumB_mln <- my_cell_bar_plot(input = macr_linf_neut,
                                      title = "CIBERSORT Immune cell fractions",
                                      features = colnames(macr_linf_neut)[-1],
                                      coord_filp = F,
                                      legend.position = "bottom",
                                      cols = cols)

##################################################################
