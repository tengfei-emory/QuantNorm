qnorm2 <-
function(ccc,batches){
  
  #ccc.0<-ccc<-1-cor(newbatch,method="spearman")
  
  batch_num <- max(batches)
  max_dim <- 0
  
  for (i in 1:batch_num){
    dim <- dim(ccc[batches == i, batches == i])[1]
    if (dim > max_dim){
      max_dim = dim
      big = i
    }
  }
  
  qt.x <- as.vector(ccc[batches == big, batches == big])
  
  for(ba1 in unique(batches))
  {
    for(ba2 in unique(batches))
    {
      qt.y<-ccc[batches == ba1, batches == ba2]
      dy <- dim(qt.y)
      qt.y <- as.vector(qt.y)
      qt.y <- zp.quantile(qt.x,qt.y)
      qt.y <- as.matrix(qt.y,nrow = dy[1],ncol=dy[2])
      ccc[batches == ba1, batches == ba2] <- qt.y
    }
  }
  
  return(ccc)
}


qnorm2s <-
  function(ccc,batches){
    
    #ccc.0<-ccc<-1-cor(newbatch,method="spearman")
    
    batch_num <- max(batches)
    max_dim <- 0
    
    for (i in 1:batch_num){
      dim <- dim(ccc[batches == i, batches == i])[1]
      if (dim > max_dim){
        max_dim = dim
        big = i
      }
    }
    
    
    qt.x <- ccc[batches == big, batches == big]
    qt.xrow <- qt.x
    qt.xcol <- qt.x
    
    for(i in 1:nrow(qt.x)){
      qt.xrow[i,] <- zp.quantile(as.vector(qt.x),qt.x[i,])
    }
    for(i in 1:ncol(qt.x)){
      qt.xcol[,i] <- zp.quantile(as.vector(qt.x),qt.x[,i])
    }
    
    qt.xrow <- as.vector(qt.xrow)
    qt.xcol <- as.vector(qt.xcol)
    
    #qt.x <- as.vector(ccc[batches == big, batches == big])
    
    for(ba1 in unique(batches))
    {
      for(ba2 in unique(batches))
      {
        qt.y<-ccc[batches == ba1, batches == ba2]
        qt.y1<-qt.y
        qt.y2<-qt.y
        dy <- dim(qt.y)
        qt.y1 <- as.vector(qt.y1)
        qt.y1 <- zp.quantile(qt.x,qt.y1)
        qt.y1 <- as.matrix(qt.y1,nrow = dy[1],ncol=dy[2])
        
        qt.y2 <- as.vector(qt.y2)
        qt.y2 <- zp.quantile(qt.xcol,qt.y2)
        qt.y2 <- as.matrix(qt.y2,nrow = dy[1],ncol=dy[2])
        
        ccc[batches == ba1, batches == ba2] <- (qt.y1+qt.y2)/2
      }
    }
    
    return(ccc)
  }
