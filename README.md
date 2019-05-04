# QuantNorm
Mitigating the adverse impact of batch effects in sample pattern detection

[![CRAN Status Badge](http://www.r-pkg.org/badges/version/QuantNorm)](http://cran.r-project.org/web/packages/QuantNorm)
[![Downloads badge](https://cranlogs.r-pkg.org/badges/QuantNorm)](https://cranlogs.r-pkg.org/badges/QuantNorm)
[![](https://img.shields.io/badge/doi-10.1093/bioinformatics/bty117-blue.svg)](https://doi.org/10.1093/bioinformatics/bty117)

QuantNorm modifies the distance matrix obtained from data with batch effects, so as to improve the performance of sample pattern detection, such as clustering, dimension reduction, and construction of networks between subjects. The method has been published in Bioinformatics (Fei et al, 2018, https://doi.org/10.1093/bioinformatics/bty117).


### QuantNorm: Vectorization approach
![Vectorization](https://github.com/tengfei-emory/Image/blob/master/f4.png)



### QuantNorm: Row/column iterative approach
![Row/Column](https://github.com/tengfei-emory/Image/blob/master/f5_4.png)

# Study design requirement and running speed

Our method can be applied on data sets with roughly balanced study design, which means the cell populations should be relatively evenly distributed among batches. For data sets with sample size larger than 500, we recommend setting a maximum iteration number (for example max = 10) to save time. We expect up to an hour for the algorithm to complete even if a maximum iteration number is set.

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

# Example - ENCODE Data for Human and Mouse Tissues
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

# Incorporating with the SC3 method

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
Then we could use SC3 in the following way (SC3 codes are borrowed from [SC3 bioconductor manual](http://www.bioconductor.org/packages/release/bioc/vignettes/SC3/inst/doc/SC3.html#singlecellexperiment-qc-and-scater)). Since SC3 is not a deterministic method, the resulting ARI index can vary, as shown in the supplementary fig. S5 in [Fei et al, 2018](https://doi.org/10.1093/bioinformatics/bty117) paper.

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
adjustedRandIndex(cell.type, svm_labels)

# Visualization
sc3_interactive(sce)
```

# References

Fei, Teng, et al. "Mitigating the adverse impact of batch effects in sample pattern detection", Bioinformatics, epub ahead of printing (2018).

Gilad, Yoav, and Orna Mizrahi-Man. "A reanalysis of mouse ENCODE comparative gene expression data." F1000Research 4 (2015).

Hemberg lab's RNA-seq dataset page https://hemberg-lab.github.io/scRNA.seq.datasets/.

Johnson, W. Evan, Cheng Li, and Ariel Rabinovic. "Adjusting batch effects in microarray expression data using empirical Bayes methods." Biostatistics 8.1 (2007): 118-127.

Kiselev, Vladimir Yu, et al. "SC3: consensus clustering of single-cell RNA-seq data." Nature methods 14.5 (2017): 483.

Usoskin, D. et al. Unbiased classification of sensory neuron types by large-scale single-cell RNA sequencing. Nat. Neurosci. 18, 145–153 (2015)


