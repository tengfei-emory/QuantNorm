#' Conducting batch effect correction by quantile normalization
#'
#'
#'
#' @param dat The original p*n batch effect data with n subjects and p RNA-seq measurements.
#' @param batch The vector of length n indicating which batch the subjects belong to.
#' @param method Method for the quantile normalization. There are three options: "refB", "ref1" and "block".
#' @param iter Times of iteration for "refB" and "ref1" methods when no standardization is applied.
#' @param logdat Whether conducting log transformation to data or not.
#' @param standardize Whether conducting standardization [(dat - mean)/var] to data or not.
#' @return The corrected 1-correlation matrix between subjects.
#' @export
#' @examples
#'
#' library(rgl) #for 3D PCA display
#'
#' data("humanmouse")
#'
#' #Numbering the cells by cell types
#' celltype <- c(rep(1,8),rep(7,6),rep(1,12),rep(2,1),rep(3,5),rep(4,3),rep(5,2),rep(6,4),
#'               rep(1,2),rep(1,4),rep(2,2),rep(3,6),rep(4,2),rep(5,2),rep(6,3))
#'
#' #Assigning the batch number that the 62 subjects belonging to.
#' batches <- c(rep(1,41),rep(2,21))
#'
#' #Plot the 3D PCA for the uncorrected batch effect data
#' plot3d(princomp(1-cor(humanmouse,method='spearman'))$scores[,1:3], col=celltype, size=10)
#'
#' #QuantNorm correction
#' ccc <- QuantNorm(humanmouse,batches,iter=10)
#' plot3d(princomp(ccc)$scores[,1:3], col=celltype, size=10)


QuantNorm <-
function(dat,batch,method="refB",iter=1,logdat=TRUE,standardize=FALSE){

  if(standardize==TRUE){

    ccc <- standardization(dat,batch)

    if(method=="ref1"){
      ccc <- qnorm(ccc,batch)
    }else if(method=="refB"){
      ccc <- qnorm1(ccc,batch)
    }else if(method=="block"){
      ccc <- qnorm2(ccc,batch)
    }

  }else{

    if(logdat == FALSE){
      ccc <- 1-cor(dat,method="spearman")
    }else{
      ccc <- 1-cor(log(dat+1),method="spearman")
    }

    if (method=='block'){

      ccc <- qnorm2(ccc,batch)

    }else{

      if (iter < 1){
        cat('iter must be larger than 1')
      }else{
        if (method == "ref1"){
          for (i in 1:floor(iter)){
            ccc <- qnorm(ccc,batch)
          }
        }else if (method == "refB"){
          for (i in 1:floor(iter)){
            ccc <- qnorm1(ccc,batch)
          }
        }else if (method == "block"){
          ccc <- qnorm2(ccc,batch)
        }else{
          cat("method must be 'ref1', 'refB' or 'block'")
        }
      }
    }
  }
  return(ccc)
}
