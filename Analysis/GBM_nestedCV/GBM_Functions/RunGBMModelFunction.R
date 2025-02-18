#Steps for outer wrapper function
#1-get optimal parameters and mean RMSE from k*n nested cross validation 
#2-run the model with the selected optimal parameters
#3-output:
#a-the model object
#b-the mean RMSE of the model
#c-the optimal parameter set

#inputs:
#pigsums,response,X_vec,ko_t,ki_t

#split_type="random"
#split_type=groupname
#dataset=pigsums_sl
#response="sl_mean"
#X_vec=X_vec.start
#family="poisson"
#ntreemax=8000

Run.GBM.Model<-function(dataset,response,X_vec,ko_t,ki_t,family,split_type,ntreemax){
  
  outer.kfold.list=outer_kfold_loop(dataset,response,X_vec,ko_t,ki_t,family,split_type,ntreemax)
  #output.list=list(meanRMSE,best.model.params,meanR2, final.pred.saver, final.test.saver)
  print("finished outer.kfold.list object")
  opt.params.kfold=outer.kfold.list[[2]]
  print("finished assigning outer.kfold.list object")
  gbm.final <- gbm.fixed(data=dataset, gbm.x=X_vec, gbm.y=which(colnames(dataset)==response),
                         learning.rate=opt.params.kfold$learning.rate, 
                         tree.complexity=opt.params.kfold$tree.complexity, 
                         n.trees=opt.params.kfold$n.trees,
                         bag.fraction = opt.params.kfold$bag.fraction,
                         family=family)  
  
  print("finished final gbm model")
  
  final.output.list=list(gbm.final, outer.kfold.list[[1]], outer.kfold.list[[2]], outer.kfold.list[[3]], outer.kfold.list[[4]], outer.kfold.list[[5]])
  
  print("finished final output list")
  
  return(final.output.list)
}
