#Use to optimize alpha hyperparameter in the lasso regression

optimize_alpha_lasso<-function(x,y,family,grid,typemeasure){
  set.seed(100)
  
  mean.error=vector(mode="numeric",length=length(grid))
  for(i in 1:length(grid)){ 
    foldid = sample(rep(seq(10), length = 535)) 
    #for(j in 1:100){
    cv.lasso <- cv.glmnet(as.matrix(x), as.matrix(y), family=family, foldid=foldid, alpha=grid[i], parallel=TRUE, standardize=TRUE, type.measure=typemeasure)
    #cv.lasso <- cv.glmnet(as.matrix(x), as.matrix(y), family=family, alpha=grid[i], parallel=TRUE, standardize=TRUE, type.measure=typemeasure)
    #hundrep[j]=mean(cv.lasso$cvm)
    #} #this inner loop is not needed when foldid is set explicitly
    #mean.error[i]=mean(hundrep)
    mean.error[i]=mean(cv.lasso$cvm)
  }
  
  mindex=which(mean.error==min(mean.error))
  
  best.alpha=grid[mindex]
  return(best.alpha)
}





