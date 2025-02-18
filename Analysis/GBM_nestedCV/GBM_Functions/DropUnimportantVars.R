#Remove variables below certain contribution threshold
#outputs vector of column ids remaining in X_vec.start
#use that X_vec to re-run models without unimportant variables
drop.unimportant<-function(gbm.object,X_vec.start,min_importance){
  #gbm.object$contributions
  #which(gbm.object$contributions<min_importance)
  #str(gbm.object$contributions)
  mincontrib.index=which(gbm.object$contributions$rel.inf<0.1)
  mincontrib.vars=gbm.object$contributions$var[mincontrib.index]
  remove.var.index=which(colnames(pigsums)%in%mincontrib.vars)
  X_vec.2=X_vec.start[which(!(X_vec.start%in%remove.var.index))]
  
  return(X_vec.2)
}