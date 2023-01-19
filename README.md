# RNAseq_TCGA_lung_cancer
Analyzing RNA-seq data with DESeq2

Une tâche fondamentale dans l'analyse des données de comptage provenant de RNA-seq est la détection des gènes exprimés de manière différentielle. Les données de comptage sont présentées sous la forme d'un tableau qui rapporte, pour chaque échantillon, le nombre de fragments de séquence qui ont été attribués à chaque gène.


Ici, nous avons des données RNAseq du cancer de poumons (source TCGA)
nous allons comparer l'expression differentielle entre les repondeurs et les non repondeur pour un traitement donnée 

### Données d'entrée : 
- tables des comptes
- table des condition (patients & reponse au trt)

### Analyse descriptive des genes
- si la table comporte des continues, il faut faire une ROUND 
- on selectionne les genes avec plus de 5 comptes
- on fait une normalisation de type dds (http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts)
- on fait une ACP et on affiche les deux composantes principales (reperage des outliers)

![Image of aciduino on protoboard](https://github.com/dinaOuahbi/RNAseq_TCGA_lung_cancer/blob/main/PCA_ntd.png)

### Analyse différentielle
- on relance une nouvelle analyse en integrant cette fois notre condition
- selectionner les genes diff avec un taux d'erreur de 5% et un logfold change >2
Quelques notions importantes : 
    Le fold-change est le rapport du niveau moyen d'expression d'un gène dans une condition par rapport à une autre
    si pour un gene, le foldchange >0, ca veut dire que ce gene est sous exprimé chez la reference
    ntd normalisation = le compte pour chaque gene dans chaque echantillon / moyenne geometrique de chaque gene a travers tous les echantillons
    
  ![Image of aciduino on protoboard](https://github.com/dinaOuahbi/RNAseq_TCGA_lung_cancer/blob/main/heatmap.png)
  ![Image of aciduino on protoboard](https://github.com/dinaOuahbi/RNAseq_TCGA_lung_cancer/blob/main/log2foldchange.png)
    
