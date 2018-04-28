# QuantNorm
Mitigating the adverse impact of batch effects in sample pattern detection

QuantNorm modifies the distance matrix obtained from data with batch effects, so as to improve the performance of sample pattern detection, such as clustering, dimension reduction, and construction of networks between subjects. The method has been published in Bioinformatics (Fei et al, 2018, https://doi.org/10.1093/bioinformatics/bty117).


### QuantNorm: Vectorization approach
![Vectorization](https://github.com/tengfei-emory/Image/blob/master/f4.png)



### QuantNorm: Row/column iterative approach
![Row/Column](https://github.com/tengfei-emory/Image/blob/master/f5_4.png)



# Installation Guide
```{r}
library(devtools)
install_github('tengfei-emory/QuantNorm')
library(QuantNorm)
```
Currently QuantNorm supports R version >= 3.3.0.

# Preprocessing

Before putting the data in QuantNorm, please feel free to conduct preprocessing, such as cell filtering, log-transformation, or standardization (zero mean and one standard deviation). Currently, QuantNorm only features simple built-in preprocessing task, such as log- transformation and standardization. 

The choice of preprocessing will probably affect the result of running QuantNorm. Since the behavior is pretty much dataset-specific, it is good to try different preprocessing methods.


# Example 1 - Human and Mouse Brain RNA-seq Data

The dataset (Zhang et al, 2016) is downloaded from http://web.stanford.edu/group/barres_lab/brainseq2/brainseq2.html. The data consists of brain RNA-seq measurements for human cells and mouse cells. There are 41 human cell samples and 21 mouse cell samples, the types of which is the column name of the data. For each cell sample, there are 15041 gene counts. The dataset can be loaded by data(brain).

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
corrected.distance.matrix <- QuantNorm(ENCODE,batches,method='row/column', cor_method='pearson', logdat=F,standardize = T,tol=1e-4)
pheatmap(1-corrected.distance.matrix)
```
QuantNorm (left) vs ComBat (right):

![Heatmaps](https://github.com/tengfei-emory/Image/blob/master/f7.png)

# Incorporating with the SC3 method (Under Construction)

(This section is still under construction as I plan to include a real data analysis exmaple of single-cell RNA-seq soon.)

As mentioned in our paper, our method can improve the performance of a current powerful clustering method, Single-Cell Consensus Clustering ([SC3](http://www.bioconductor.org/packages/release/bioc/html/SC3.html), Kiselev VY et al, 2017). The following code shows how we can plug in the corrected distance matrix to the SC3 algorithm in R. For more detailed tutorial about SC3, please refer to [this page](http://www.bioconductor.org/packages/release/bioc/vignettes/SC3/inst/doc/SC3.html) by Vladimir Kiselev.

The following example is based on the single-cell RNA-seq data of mouse neuron (Usoskin et al, 2015). The SingleCellExperiment object 'usoskin.rds' is obtained from [Hemberg-lab's RNA-seq dataset website](https://hemberg-lab.github.io/scRNA.seq.datasets/mouse/brain/).

```{r}
library(SC3)
library(scater)
library(SingleCellExperiment)
library(QuantNorm)
library(Biobase)
library(mclust)

# Import the SCE object data
scenet<-readRDS(url('https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/usoskin.rds'))
counts(scenet) <- normcounts(scenet)

# This is the true cell type
cell.type <- scenet@colData$cell_type1

# This is the expression matrix
dat.sc3 <- exprs(scenet)

# This is batch information
batch <- scenet@colData$Library
```

For the expressoin matrix dat.sc3, we can obtain the 2 corrected distance matrix from QuantNorm, one is based on spearman correlation and one is based on pearson correlation: (In order to save time, I set the maximum iteration as 5. This can be adjusted)

```{r}
pearson.cor <- QuantNorm(dat.sc3[rowSums(dat.sc3 != 0) > 622/4,],batch=as.numeric(as.factor(batch)),cor_method='pearson',logdat=F,max=5)

spearman.cor <- QuantNorm(dat.sc3[rowSums(dat.sc3 != 0) > 622/4,],batch=as.numeric(as.factor(batch)),cor_method='spearman',logdat=F,max=5)
```
Then we could use SC3 in the following way (SC3 codes are borrowed from [SC3 bioconductor manual](http://www.bioconductor.org/packages/release/bioc/vignettes/SC3/inst/doc/SC3.html#singlecellexperiment-qc-and-scater)) 

```{r}
# Run the SC3 algorithm step by step
sce <- sc3_prepare(scenet)
sce <- sc3_estimate_k(sce)
sce <- sc3_calc_dists(sce)

# Replacing the original pearson and spearman distance matrices by the corrected ones.
sce@metadata$sc3$distances$spearman <- spearman.cor
sce@metadata$sc3$distances$pearson <- pearson.cor

# Continue finishing the standard steps
sce <- sc3_calc_transfs(sce)
sce <- sc3_kmeans(sce, ks = 4)
sce <- sc3_calc_consens(sce)

# Check clustering results and ARI
library(mclust)
svm_labels <- colData(sce)$sc3_4_clusters
adjustedRandIndex(cell.type[,1], svm_labels)
```

# References

Fei, Teng, et al. "Mitigating the adverse impact of batch effects in sample pattern detection", Bioinformatics, epub ahead of printing (2018).

Gilad, Yoav, and Orna Mizrahi-Man. "A reanalysis of mouse ENCODE comparative gene expression data." F1000Research 4 (2015).

Hemberg lab's RNA-seq dataset page https://hemberg-lab.github.io/scRNA.seq.datasets/.

Johnson, W. Evan, Cheng Li, and Ariel Rabinovic. "Adjusting batch effects in microarray expression data using empirical Bayes methods." Biostatistics 8.1 (2007): 118-127.

Kiselev, Vladimir Yu, et al. "SC3: consensus clustering of single-cell RNA-seq data." Nature methods 14.5 (2017): 483.

Usoskin, D. et al. Unbiased classification of sensory neuron types by large-scale single-cell RNA sequencing. Nat. Neurosci. 18, 145–153 (2015)

Zhang, Ye, et al. "Purification and characterization of progenitor and mature human astrocytes reveals transcriptional and functional differences with mouse." Neuron 89.1 (2016): 37-53.


