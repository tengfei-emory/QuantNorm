# QuantNorm
Quantile Normalization for Batch Effect Removal 

# Installation Guide
```{r}
library(devtools)
install_github('tengfei-emory/QuantNorm')
```

# Dataset

The dataset used in the package is downloaded from http://web.stanford.edu/group/barres_lab/brainseq2/brainseq2.html. The data consists of brain RNA-seq measurements for human cells and mouse cells. There are 41 human cell samples and 21 mouse cell samples, the types of which is the column name of the data. For each cell sample, there are 15041 gene counts.

The dataset can be loaded by data(humanmouse).

# Example
```{r}
library(rgl) #for 3D PCA display

data("humanmouse")

#Numbering the cells by cell types
celltype <- c(rep(1,8),rep(7,6),rep(1,12),rep(2,1),rep(3,5),rep(4,3),rep(5,2),rep(6,4),
              rep(1,2),rep(1,4),rep(2,2),rep(3,6),rep(4,2),rep(5,2),rep(6,3))

#Assigning the batch number that the 62 subjects belonging to.
batches <- c(rep(1,41),rep(2,21))

#Plot the 3D PCA for the uncorrected batch effect data
plot3d(princomp(1-cor(humanmouse,method='spearman'))$scores[,1:3], col=celltype, size=10)

#QuantNorm correction (one iteration)
ccc <- QuantNorm(humanmouse,batches,iter=1)
plot3d(princomp(ccc)$scores[,1:3], col=celltype, size=10)

#QuantNorm correction (ten iterations)
ccc <- QuantNorm(humanmouse,batches,iter=10)
plot3d(princomp(ccc)$scores[,1:3], col=celltype, size=10)
```


# References
Zhang, Ye, et al. "Purification and characterization of progenitor and mature human astrocytes reveals transcriptional and functional differences with mouse." Neuron 89.1 (2016): 37-53.
