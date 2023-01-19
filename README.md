# RNAseq_TCGA_lung_cancer
Analyzing RNA-seq data with DESeq2


A basic task in the analysis of count data from RNA-seq is the detection of differentially expressed genes. The count data are presented as a table which reports, for each sample, the number of sequence fragments that have been assigned to each gene.


ici, nous avons des données RNAseq du cancer de poumons (source TCGA)
nous allons comparer l'expression differentielle entre les repondeurs et les non repondeur pour un traitement donnée 

### données d'entrée : 
- tables des comptes
- table des condition (patients & reponse au trt)

### Analyse descriptive des genes
- si la table comporte des continues, il faut faire une ROUND 
- on selectionne les genes avec plus de 5 comptes
- on fait une normalisation de type dds (http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts)
- on fait une ACP et on affiche les deux composantes principales (reperage des outliers)

![Image of aciduino on protoboard]()

### analyse différentielle
- on relance une nouvelle analyse en integrant cette fois notre condition
- selectionner les genes diff avec un taux d'erreur de 5% et un logfold change >2
Quelques notions importantes : 
    Le fold-change est le rapport du niveau moyen d'expression d'un gène dans une condition par rapport à une autre
    si pour un gene, le foldchange >0, ca veut dire que ce gene est sous exprimé chez la reference
    ntd normalisation = le compte pour chaque gene dans chaque echantillon / moyenne geometrique de chaque gene a travers tous les echantillons
    
  ![Image of aciduino on protoboard]()
  ![Image of aciduino on protoboard]()
    
