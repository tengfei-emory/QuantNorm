# QuantNorm
Quantile Normalization for Batch Effect Adjustment

Our method aims to effectively cluster the cells by modifying the correlation matrix of cell RNA-seq data by quantile normalization. This method has been published in Bioinformatics (https://doi.org/10.1093/bioinformatics/bty117).

### QuantNorm: Vectorization approach
![Vectorization](https://github.com/tengfei-emory/Image/blob/master/f4.png)

### QuantNorm: Row/column iterative approach
![Row/Column](https://github.com/tengfei-emory/Image/blob/master/f6.png)

# Installation Guide
```{r}
library(devtools)
install_github('tengfei-emory/QuantNorm')
library(QuantNorm)
```


# Example 1 - Human and Mouse Brain RNA-seq Data

The dataset used in this example is downloaded from http://web.stanford.edu/group/barres_lab/brainseq2/brainseq2.html. The data consists of brain RNA-seq measurements for human cells and mouse cells. There are 41 human cell samples and 21 mouse cell samples, the types of which is the column name of the data. For each cell sample, there are 15041 gene counts. The dataset can be loaded by data(brain).

```{r}
library(rgl) #for 3D PCA display

data("brain")

#Numbering the cells by cell types
celltype <- c(rep(1,8),rep(7,6),rep(1,12),rep(2,1),rep(3,5),rep(4,3),rep(5,2),rep(6,4),
              rep(1,2),rep(1,4),rep(2,2),rep(3,6),rep(4,2),rep(5,2),rep(6,3))

#Assigning the batch number that the 62 subjects belonging to.
batches <- c(rep(1,41),rep(2,21))

#Plot the 3D PCA for the uncorrected batch effect data
plot3d(princomp(1-cor(brain,method='spearman'))$scores[,1:3], col=celltype, size=10)

#QuantNorm correction
ccc <- QuantNorm(brain,batches,tol=1e-4)
plot3d(princomp(ccc)$scores[,1:3], col=celltype, size=10)

#ComBat
library(sva) #(may need to install package sva from bioconductor)
cleandat <- ComBat(log(brain+1),batches)
ccc.1<-1-cor(cleandat,method="spearman")
plot3d(princomp(ccc)$scores[,1:3], col=celltype, size=10)
```

# Example 2 - ENCODE Data for Human and Mouse Tissues
The dataset used in this example is reproduced according to Gilad et al (2015). The data contains the log-normalized counts for 13 different tissues for human and their counterparts in mouse. Our algorithm obtains clusters by tissues. The dataset can be loaded by data(ENCODE).

```{r}
library(pheatmap) #drawing heatmap
data("ENCODE")

#Before correction, the subjects are clustered by species
pheatmap(cor(ENCODE))

#Assigning the batches based on species
batches <- c(rep(1,13),rep(2,13))

#QuantNorm correction
corrected.distance.matrix <- QuantNorm(ENCODE,batches,method='row/column', cor_method='pearson', logdat=F,standardize = T)
pheatmap(1-corrected.distance.matrix)

#ComBat
library(sva) #(may need to install package sva from bioconductor)
cleandat <- ComBat(ENCODE,batches)
pheatmap(cor(cleandat),clustering_method="average")
```

# References

Fei, Teng, et al. "Mitigating the adverse impact of batch effects in sample pattern detection", Bioinformatics (2018, accepted).

Gilad, Yoav, and Orna Mizrahi-Man. "A reanalysis of mouse ENCODE comparative gene expression data." F1000Research 4 (2015).

Johnson, W. Evan, Cheng Li, and Ariel Rabinovic. "Adjusting batch effects in microarray expression data using empirical Bayes methods." Biostatistics 8.1 (2007): 118-127.

Zhang, Ye, et al. "Purification and characterization of progenitor and mature human astrocytes reveals transcriptional and functional differences with mouse." Neuron 89.1 (2016): 37-53.


