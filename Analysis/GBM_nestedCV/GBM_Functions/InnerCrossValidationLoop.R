#Purpose of this function:
#Run the inner cross validation loop

inner_kfold_loop<-function(outerpigtrain,ki_t,response,X_vec,family,ntreemax){

############################################
##Start Inner K-fold function here
#inputs:
#outerpigtrain is the subsetted pigsums, all but the current outer training set
#ki_t is the total number of inner k-folds
#response is the response variable name (i.e. "sl_" or "displacement")
#X_vec is vector of column numbers to use as independent variables
  
#split outer training set into k cross validation subsets
inner_k_fold=k_split(outerpigtrain,ki_t)

#initialize output data frame
opt.params.kfold=data.frame(matrix(nrow=ki_t,ncol=6))
colnames(opt.params.kfold)=c("n.trees","learning.rate","tree.complexity","bag.fraction","rmse.inner","rmse.outer")

##inner for loop start here
for(ki in 1:ki_t){
  
  print(paste("starting inner replicate", ki))
  print("splitting train and test sets")
  #sample=nrow(pigsums)*0.90 #sample 3/4 of the data
  
  #pigtrain <- as.data.frame(pigsums %>% slice_sample(n=sample)%>% ungroup()) #subset for training
  #pigtest <- anti_join(pigsums,pigtrain) #out of sample set
  
  innerpigtest=inner_k_fold[[ki]]
  innerpigtrain=anti_join(outerpigtrain,innerpigtest)
  
  #3-Fix tree hyperparameters and tune learning rate and assess speed vs. performance.
  # create grid search
  
  hyper_grid <- expand.grid(
    learning_rate = c(0.3, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0003, 0.0001),
    trees = NA,
    RMSE = NA,
    time = NA
  )
  
  # execute grid search
  for(i in seq_len(nrow(hyper_grid))) {
    set.seed(123) #for reproducibility
    train.time <- system.time({
      gbm <- gbm.fixed(data=innerpigtrain, gbm.x=X_vec, gbm.y=which(colnames(innerpigtrain)==response),
                       learning.rate=hyper_grid$learning_rate[i], 
                       tree.complexity=5, 
                       n.trees=ntreemax,
                       bag.fraction = 0.5,
                       family=family)
      
      tree.list <- seq(100, ntreemax, by=100)
      
      pred <- predict.gbm(gbm, innerpigtest, n.trees=tree.list, "response")
      
      pig.pred.deviance <- rep(0,length(tree.list))
      for (j in 1:length(tree.list)) {
        #pig.pred.deviance[j] <- calc.deviance(pigtest$sl_,pred[,j],family="poisson", calc.mean=TRUE)
        pig.pred.deviance[j]=sqrt(mean((innerpigtest[,which(colnames(innerpigtest)==response)] - pred[,j])^2))
      }
    })
    #which(pig.pred.deviance==min(pig.pred.deviance))
    hyper_grid$trees[i]=tree.list[which(pig.pred.deviance==min(pig.pred.deviance))]
    hyper_grid$RMSE[i]=min(pig.pred.deviance)
    hyper_grid$time[i]=train.time[["elapsed"]]
  }
  
  #best params:
  opt.params=hyper_grid[which(hyper_grid$RMSE==min(hyper_grid$RMSE)),]
  print("got optimal learning rate and ntrees")
  
  #4-Tune tree-specific parameters for decided learning rate.
  #create search grid
  hyper_grid <- expand.grid(
    n.trees = opt.params$trees,
    learning.rate = opt.params$learning_rate,
    tree.complexity = c(1, 3, 5, 7, 9, 10, 12),
    bag.fraction = c(0.2, 0.5, 0.6, 0.8, 0.9)
  )
  
  #fit the model with all parameter combos in hyper_grid
  hyper_grid$rmse <- purrr::pmap_dbl(
    hyper_grid,
    ~ hyperparam_model_fit(
      n.trees = ..1,
      learning.rate = ..2,
      tree.complexity = ..3,
      bag.fraction = ..4,
      innerpigtrain, 
      innerpigtest, 
      response, 
      X_vec,
      family
    )
  )
  
  opt.params=hyper_grid[which(hyper_grid$rmse==min(hyper_grid$rmse))[1],]
  print("got optimal tree complexity and bag fraction")
  
  
  #5-Once tree-specific parameters have been found, lower the learning rate
  #to assess for any improvements in accuracy.
  hyper_grid <- expand.grid(
    learning_rate = c(0.3, 0.1, 0.05, 0.01, 0.005,0.001,0.0005,0.0003,0.0001),
    trees = NA,
    RMSE = NA,
    time = NA
  )
  
  for(i in seq_len(nrow(hyper_grid))) {
    set.seed(123) #for reproducibility
    train.time <- system.time({
      gbm <- gbm.fixed(data=innerpigtrain, gbm.x=X_vec, gbm.y=which(colnames(innerpigtrain)==response),
                       learning.rate=hyper_grid$learning_rate[i], 
                       tree.complexity=opt.params$tree.complexity, 
                       n.trees=opt.params$n.trees,
                       bag.fraction = opt.params$bag.fraction,
                       family=family)
      
      tree.list <- seq(100, ntreemax, by=100)
      
      pred <- predict.gbm(gbm, innerpigtest, n.trees=tree.list, "response")
      
      pig.pred.deviance <- rep(0,50)
      for (j in 1:50) {
        #pig.pred.deviance[j] <- calc.deviance(pigtest$sl_,pred[,j],family="poisson", calc.mean=TRUE)
        pig.pred.deviance[j]=sqrt(mean((innerpigtest[,which(colnames(innerpigtest)==response)] - pred[,j])^2))
      }
    })
    #which(pig.pred.deviance==min(pig.pred.deviance))
    hyper_grid$trees[i]=tree.list[which(pig.pred.deviance==min(pig.pred.deviance))]
    hyper_grid$RMSE[i]=min(pig.pred.deviance)
    hyper_grid$time[i]=train.time[["elapsed"]]
  }
  
  opt.params$learning.rate=hyper_grid[which(hyper_grid$RMSE==min(hyper_grid$RMSE)),]$learning_rate
  opt.params$rmse.inner=hyper_grid[which(hyper_grid$RMSE==min(hyper_grid$RMSE)),]$RMSE
  opt.params$rmse.outer=NA
  opt.params$rmse=NULL
  
  
  print("got final optimal parameter set")
  print(opt.params)
  
  opt.params.kfold[ki,]=opt.params
  print("added optimal params to kfold parameter dataframe")
  
} #inner k fold loop closing bracket

return(opt.params.kfold)

}