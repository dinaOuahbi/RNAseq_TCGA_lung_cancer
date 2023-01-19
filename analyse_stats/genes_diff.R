



base = read.csv2("/work/shared/ptbc/CNN_Pancreas_V2/projet_RNAseq_Julie/data/base_TCGA_log2.csv", h=T)
base = base[-c(which(is.na(base$Rep))),]

library(biomaRt); library(DESeq2); library(gridExtra)
library("tximport"); library("readr")
library(vsn); library(pheatmap)
library(RColorBrewer)
library(ggplot2)

base_cli = base

base = read.csv2("/work/shared/ptbc/CNN_Pancreas_V2/projet_RNAseq_Julie/data/LUAD_merged.csv", h=T, sep=" ")
base1 = base
base = base[,-c(1:3)] #enlever le status et l'os
base = apply(base, 2, as.numeric) # 2 pour les colonnes
rownames(base) = base1$bcr #indexer sur le nom des patients

#selectionne moi les patients dont la reponse est valable (avec toutes leurs colonnes)
base = base[c(base_cli$Patient[which(!is.na(base_cli$Rep))]),]

# verifier si ya pas de difference entre les deux tables en terme de patients
#base = log2(base+1)
which(rownames(base) != base_cli$Patient)

# on round les RSEM pour avoir des entiers
base = round(base)



# ---- Analyse descriptive des genes ----

# ne prendre que les nom des patients et leurs reponse
samples = data.frame(id_patient = base_cli$Patient, condition = as.factor(base_cli$Rep))


# creer un objet de type DESeqDataSetFromMatrix (sans la conditions, pour reperage des outlayers sur une figure PCA)
ddsTxi <- DESeqDataSetFromMatrix(t(base),
                                 colData = samples,
                                 design =~ 1)

# on ne garde que les genes qui ont plus de 5 comptes
ddsTC <- ddsTxi[ rowSums(counts(ddsTxi)) > 5, ]       

dds <- DESeq(ddsTC)
res <- results(dds)
res

# normalisation
ntd <- normTransform(dds)

#The assay function is used to extract the matrix of normalized values
meanSdPlot(assay(ntd))

# ACP
plotPCA(ntd)
# pas d'outlier

# principal comp analysis of normalized matrix
PCA_LPI <- prcomp(t(assay(ntd)))

# select only PC1 and PC2
ploting_PCA_LPI<-data.frame(PCA_LPI$x[,1:2])

# get the % of variance of each PC
percentVar<-round(100*summary(PCA_LPI)$importance[2,1:2],0)

# plot PCA
ggplot(data=ploting_PCA_LPI)+
  geom_point(aes(x=PC1,y=PC2),size=2)+
  xlab(paste0("PC1: ",percentVar[1],"% variance"))+
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  coord_fixed(ratio=1.0)+
  #ggtitle("Gene Count + HERVs (telescope)")+
  theme(text = element_text(size=15),aspect.ratio=1)

#On a un point compl?tement extr?me:
out=which(PCA_LPI$x[,2]>125)
out

# ---- Analyse différentielle ----
ddsTxi <- DESeqDataSetFromMatrix(t(base),
                                 colData = samples,
                                 design =~ condition)

# on ne garde que les genes qui ont plus de 5 comptes
ddsTC <- ddsTxi[ rowSums(counts(ddsTxi)) > 5, ]       

# on verifie les conditions :
ddsTC$condition
ddsTC$condition = relevel(ddsTC$condition, ref="no_rep") # si selon l'ordre alphanumerique


#Le fold-change est le rapport du niveau moyen d'expression d'un gène dans une condition par rapport à une autre
# si pour un gene, le foldchange >0, ca veut dire que ce gene est sous exprimé chez la reference  
dds <- DESeq(ddsTC)
res <- results(dds)
res

# normalisation
#DESeq2 performs an internal normalization where geometric mean is calculated for each gene across all samples. 
#The counts for a gene in each sample is then divided by this mean.
ntd <- normTransform(dds)

# d'autres types de normalisation existe : rld / vsd 
# le choix de la normalisation se fait visuelement sur le graphique meanSdPlot (on prend la ou c'est homogene)

#Plot row standard deviations versus row means
meanSdPlot(assay(ntd))

# ACP : differencier entre les deux classes de patient
plotPCA(ntd)
# pas d'outlier

# Stat
# ordoner la table en fonction des p value (dans l'ordre croissant)
resOrdered <- res[order(res$pvalue),]


# ce qui nous importe le plus dans les res, c'est la p value et le fold change
    # la p value pour reperer les genes sign
    # le fold change c'est pour voir dans quel sensle gene est different 
summary(res)
summary(res$pvalue)
summary(res$log2FoldChange)
hist(res$log2FoldChange) # si le foldchange est proche de 0 

# selectionner les genes diff avec un taux d'erreur de 5% et un logfold change >2
res[which(res$padj<0.05 & abs(res$log2FoldChange)>2),]


# --- Enregistrement ---
write.csv2(res, "/work/shared/ptbc/CNN_Pancreas_V2/projet_RNAseq_Julie/Resultats/Res_allgenes_RepvsNorep_Deseq2_RSEM.csv", row.names=T)

# ---- Volcano plot ----
# Creation de l'objet pour ggplot : df avec le nom du gene / le log2 / les p values (+adjust)
res_volcano <- as.data.frame(res[,c('log2FoldChange','pvalue','padj')])
res_volcano$Gene = rownames(res_volcano) #indexer sur le nom
res_volcano$color <- ifelse(res_volcano$padj<0.05,'adjusted p-value<0.05','NS')
res_volcano$color <- ifelse(abs(res_volcano$log2FoldChange)>1,'abs(LFC)>1',res_volcano$color)
res_volcano$color <- ifelse(res_volcano$padj<0.05 & abs(res_volcano$log2FoldChange)>1,'ajusted p-value <0.05 & abs(LFC)>1',res_volcano$color)
res_volcano$color = as.factor(res_volcano$color)
table(res_volcano$color)
res_volcano = res_volcano[c(which(!is.na(res_volcano$padj))),]

# Graph
library(ggplot2); library(ggrepel)
ggplot(res_volcano, aes(x = log2FoldChange, y = -log10(padj), colour=color)) + geom_point()+
  scale_colour_manual(values=c("orange", "green", "red", "grey"))+
  theme_classic() +
  geom_text_repel(
    data = subset(res_volcano, (padj < 0.05 & abs(log2FoldChange)>1)),
    aes(label = Gene),
    size = 4,
    colour="black",
    box.padding = unit(0.4, "lines"),
    point.padding = unit(0.2, "lines")
  )+
  theme_bw(base_size = 12) + theme(legend.title = element_blank(),legend.position = "bottom", legend.text = element_text(size=14)) +
  xlab('log2FoldChange')+ ylab('-log10 adjusted p-value')+
  ggtitle("Rep vs No rep")


# --- Heatmap des genes d'interet ---
data = assay(ntd)

# selectionner des genes au choix
#data = data[c(which(rownames(data) %in% c("IL1B", "IL1RN", "IL1R1", "IL1R2", "AIM2", "CD8A", "CXCL10","TXNIP"))),]

# selectionner des genes significative et tres differents (LFC > 3)
names_genes = res_volcano$Gene[which(res_volcano$padj<0.05 & abs(res_volcano$log2FoldChange)>3)]
data = data[c(which(rownames(data) %in% names_genes)),]


df <- data.frame(colData(dds)[,c("condition")]) # select gene name and class
rownames(df) = colData(dds)[,1] ; names(df) = "Condition" #rename
ann_colors = list("Condition" = c(rep = "#FF9933", no_rep = "#666666"))  #choose color between rep and non rep

data = data[, c(which(colnames(data) %in% rownames(df)[which(df$Condition=="rep")]),
                which(colnames(data) %in% rownames(df)[which(df$Condition=="no_rep")])
                )]

pheatmap(data, cluster_rows=T, show_rownames=T,border_color="black",treeheight_col=0,
         color = colorRampPalette(rev(brewer.pal(n=10, name = "RdBu")))(100),
         cluster_cols=F, show_colnames = F, annotation_col=df, annotation_colors = ann_colors)


