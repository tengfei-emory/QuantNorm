#' Adjust the distance matrix by quantile normalization
#'
#'
#'
#' @param dat The original p*n batch effect data with n subjects and p RNA-seq measurements.
#' @param batch The vector of length n indicating which batch the subjects belong to.
#' @param method Method for the quantile normalization. There are two options: "row/column" and "vectorize".
#' @param cor_method Method to calculate the correlation matrix, can be 'spearman'(default), 'pearson' or 'kendall'.
#' @param tol The tolerance for the iterative method "row/column", which is the Euclidean distance of the vectorized two dissimilarity matrices before and after each iteration.
#' @param max Maximum number of the iteration if the tolerance is not reached.
#' @param logdat Whether conducting log transformation to data or not.
#' @param standardize Whether conducting standardization [(dat - mean)/sqrt(var)] to data or not.
#' @return The corrected 1-correlation matrix between subjects.
#' @export
#' @examples
#'
#' library(rgl) #for 3D PCA display
#'
#' data("brain")
#'
#' #Numbering the cells by cell types
#' celltype <- c(rep(1,8),rep(7,6),rep(1,12),rep(2,1),rep(3,5),rep(4,3),rep(5,2),rep(6,4),
#'               rep(1,2),rep(1,4),rep(2,2),rep(3,6),rep(4,2),rep(5,2),rep(6,3))
#'
#' #Assigning the batch number that the 62 subjects belonging to.
#' batches <- c(rep(1,41),rep(2,21))
#'
#' #Plot the 3D PCA for the uncorrected batch effect data
#' plot3d(princomp(1-cor(brain,method='spearman'))$scores[,1:3], col=celltype, size=10)
#'
#' #QuantNorm correction
#' ccc <- QuantNorm(brain,batches,tol=1e-4)
#' plot3d(princomp(ccc)$scores[,1:3], col=celltype, size=10)


QuantNorm <- function (dat, batch, method = "row/column", cor_method = 'spearman', tol = 1e-2, max = 50, logdat = TRUE,
                       standardize = FALSE)
{
  dist = 10
  iter = 0

  if (standardize == TRUE) {
    ccc <- standardization(dat, batch, method=cor_method)
  }
  else if (logdat == FALSE) {
    ccc <- 1 - cor(dat, method = cor_method)
  }
  else {
    ccc <- 1 - cor(log(dat + 1), method = cor_method)
  }

  if (method == "vectorize") {
    ccc <- qnorm2(ccc, batch)
  }

  else if (method == "row/column") {
    while (dist > tol && iter < max){
      ccc.0 <- ccc
      ccc <- qnorm1(ccc, batch)
      dist = sqrt(sum((as.vector(ccc)-as.vector(ccc.0))^2))
      iter = iter+1
    }

    if (dist <= tol){
      cat(paste("Algorithm converged after", iter, "iterations. \n"))
    }else{
      cat(paste("The algorithm did not converge after", max,"iteration. \n The difference between the last two iteration is", dist,".\n"))
    }

  }else{
    cat("method must be 'row/column' or 'vectorize'. \n")
  }
  return(ccc)
}
