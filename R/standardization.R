standardization <-
function(dat,batch,method){
  if(is.matrix(dat)==FALSE){
    dat <- as.matrix(dat)
  }
  cat('It may take a long time to standardize the data.')
  for (i in 1:dim(dat)[1]){
    for (j in unique(batch)){
      mean <- mean(dat[i,batch==j])
      var <- var(dat[i,batch==j])
      dat[i,batch==j] <- (dat[i,batch==j] - mean)/sqrt(var)
    }
  }
  dat <- na.omit(dat)
  ccc <- 1-cor(dat,method=method)
  return(ccc)
}
