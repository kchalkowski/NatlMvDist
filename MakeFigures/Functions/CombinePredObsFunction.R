CombinePredObs<-function(CVstats.random,CVstats.region){
  pred.region.table=data.frame(matrix(nrow=0,ncol=4))
  colnames(pred.region.table)=c("region","CVmethod","test","preds")
  
  for(i in 1:length(CVstats.random[[2]])){
    #combine the random set results 
    predset=CVstats.random[[2]][[i]]
    tempdf=data.frame(matrix(nrow=length(predset$test.set),ncol=4))
    binded=cbind("test"=predset$test.set,"preds"=predset$pred.set)
    pivotdf=as.data.frame(pivot_longer(as.data.frame(binded),cols=c("test","preds")))
    randf=data.frame(matrix(nrow=nrow(binded),ncol=4))
    randf[,1]=rep(names(CVstats.random[[2]])[i],nrow(binded))
    randf[,2]=rep("random",nrow(binded))
    randf[,3]=binded[,1]
    randf[,4]=binded[,2]
    colnames(randf)=c("region","CVmethod","test","preds")
    
    if(missing(CVstats.region)){
    #combine the region set results
    pred.region.table_i=randf
    pred.region.table=rbind(pred.region.table,pred.region.table_i)  
    } else{
      predset=CVstats.region[[2]][[i]]
      tempdf=data.frame(matrix(nrow=length(predset$test.set),ncol=4))
      binded=cbind("test"=predset$test.set,"preds"=predset$pred.set)
      pivotdf=as.data.frame(pivot_longer(as.data.frame(binded),cols=c("test","preds")))
      regdf=data.frame(matrix(nrow=nrow(binded),ncol=4))
      regdf[,1]=rep(names(CVstats.region[[2]])[i],nrow(binded))
      regdf[,2]=rep("region",nrow(binded))
      regdf[,3]=binded[,1]
      regdf[,4]=binded[,2]
      colnames(regdf)=c("region","CVmethod","test","preds")
      
      pred.region.table_i=rbind(randf,regdf)
      pred.region.table=rbind(pred.region.table,pred.region.table_i)
    }
  }
  return(pred.region.table)
}
