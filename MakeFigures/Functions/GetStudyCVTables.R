#The purpose of this script is to define a set of functions to get CV stats for each individual out-of-sample prediction by study
#CVstats_sl.random=GetCVStats_Table(pigsums_sl,X_vec_list$Xsl,"sl_mean",sl.opt.params.kfold,"poisson","random",studydf)



MakeSitePredSets=function(dataset,X.vec,response,opt.params,family){
  ksplits=k_split_group(dataset,"region")
  out.region.list=vector(mode="list",length=length(ksplits))
  region.names=vector(mode="character",length=length(ksplits))
  for(i in 1:length(ksplits)){
    test_i=ksplits[[i]] #held out sample in outer loop
    train_i=anti_join(dataset,test_i) #training sample in outer loop plus test sample in inner loop
    
    gbm_i=gbm.fixed(data=train_i, gbm.x=X.vec, gbm.y=which(colnames(dataset)==response),
                    learning.rate=opt.params$learning.rate, 
                    tree.complexity=opt.params$tree.complexity, 
                    n.trees=opt.params$n.trees,
                    bag.fraction=opt.params$bag.fraction,
                    family=family)
    
    pred_i <- predict.gbm(gbm_i, test_i, opt.params$n.trees, "response")
    
    out.list=list(test_i[,which(colnames(dataset)==response)],pred_i)
    names(out.list)=c("test.set","pred.set")
    out.region.list[[i]]=out.list
    region.names[i]=as.character(ksplits[[i]]$region[1])
  }
  names(out.region.list)=region.names
  return(out.region.list)
}


Get_OutofSite_CVstats<-function(predset,type){
  CVstats=data.frame(matrix(nrow=length(predset),ncol=4))
  colnames(CVstats)=c("region","RMSE","R2","type")
  CVstats[,4]=rep(type,nrow(CVstats))
  
  for(i in 1:length(predset)){
    preds_i=predset[[i]]$pred.set
    test_i=predset[[i]]$test.set
    
    CVstats[i,1]=names(predset)[i]
    CVstats[i,2]=sqrt(mean((test_i - preds_i)^2))
    CVstats[i,3]=cor(test_i,preds_i)^2
    
  }
  return(CVstats)
}

#CVstats_sigma.sl.random=GetCVStats_Table(pigsums_sigmasl,X_vec_list$sigmasl,"sl_disp",sigma.sl.opt.params.kfold,"gaussian","random",studydf)
#pigsums2,X_vec.start,"sl_",sl.opt.params.kfold,"poisson","random",studydf
dataset=pigsums_sigmasl
X.vec=X_vec_list$sigmasl
response="sl_disp"
opt.params=sigma.sl.opt.params.kfold
family="gaussian"
CVtype="random"
studydf=studydf
GetCVStats_Table<-function(dataset,X.vec,response,opt.params,family,CVtype,studydf){
region.pred.set=MakeSitePredSets(dataset,X.vec,response,opt.params,family)
CV_table=Get_OutofSite_CVstats(region.pred.set,CVtype)

return(CV_table)
}



