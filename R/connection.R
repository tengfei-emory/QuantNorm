#' Construct connection matrix for network analysis
#'
#' @param mat n*n dissimilarity (1-correlation) matrix (e.g. obtained by QuantNorm).
#' @param label n-dimension vector for the labels of the n subjects. Replicates share the same label.
#' @param threshold A number between 0 to 1. Two groups will be regarded as connected if average 1-correlation < threshold.
#' @param closest True or False. Whether connect the closest group or not if the closest group cannot satisfy the threshold condition.
#' @export
#' @examples
#'
#' library(network); library(GGally) ;library(ggplot2) #drawing network graph
#'
#' data("ENCODE")
#'
#' #Assigning the batches based on species
#' batches <- c(rep(1,13),rep(2,13))
#'
#' #QuantNorm correction
#' corrected.distance.matrix <- QuantNorm(ENCODE,batches,method='row/column', cor_method='pearson',
#'                                        logdat=FALSE,standardize = TRUE, tol=1e-4)
#'
#' #Constructing connection matrix
#'
#' mat <- connection.matrix(corrected.distance.matrix,colnames(corrected.distance.matrix))
#'
#' #Creating network object and plot
#' ENCODE.net=network(mat, directed=FALSE)
#' ENCODE.net.col <- c("Human" = "#ff69b4", "Mouse" = "#0099ff")
#' ENCODE.net %v% "Species" <- c(rep('Human',13),rep('Mouse',13))
#' p0 <- ggnet2(ENCODE.net,label=TRUE,color = 'Species', palette = "Set2",
#'              size = 3, vjust = -0.6,mode = "kamadakawai",label.size = 3,
#'              color.legend = 'Species')+theme(legend.position = 'bottom')
#' plot(p0)


connection.matrix <- function(mat,label,threshold=0.15,closest=TRUE){

  type.no <- as.numeric(factor(label))
  unique.no <- unique(type.no)
  cell.name <- unique(label)
  m <- length(cell.name)
  connection.matrix <- matrix(0,m,m)

  for (i in unique.no){
    for (j in unique.no){
      if (i == j){
        connection.matrix[unique.no == i,unique.no ==j] <- 0
      }else{
        connection.matrix[unique.no == i,unique.no ==j] <- mean(mat[type.no==i,type.no==j])
      }
    }
  }

  connection.matrix.01 <- connection.matrix
  rownames(connection.matrix.01) <- colnames(connection.matrix.01) <- cell.name
  connection.matrix.01[connection.matrix <= threshold] = 1
  connection.matrix.01[connection.matrix > threshold] = 0
  #rownames(connection.matrix) <- colnames(connection.matrix) <- cell.name

  if (closest == T){
    diag(connection.matrix)<-1
    for (i in 1:m){
      connection.matrix.01[max.col(-connection.matrix)[i],i] = 1
      connection.matrix.01[i,max.col(t(-connection.matrix))[i]] = 1
    }
    diag(connection.matrix)<-0
  }

  connection.matrix.01
}
