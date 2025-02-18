#The purpose of this function is to determine the RMSE value with a given test and training set
hyperparam_model_fit <- function(n.trees, learning.rate, tree.complexity, bag.fraction, innerpigtrain, innerpigtest, response, X_vec, family) {
  set.seed(123)
  m <- gbm.fixed(data=innerpigtrain, 
                 gbm.x=X_vec, 
                 gbm.y=which(colnames(innerpigtrain)==response),
                 learning.rate=learning.rate,
                 tree.complexity=tree.complexity, 
                 n.trees=n.trees,
                 bag.fraction = bag.fraction,
                 family=family)
  
  pred <- predict.gbm(m, innerpigtest, n.trees=n.trees, "response")
  RMSE=sqrt(mean((innerpigtest[,which(colnames(innerpigtest)==response)] - pred)^2))
  
  # compute RMSE
  return(RMSE)
}