##Run outer loop from k*n-fold nested cross validation

#for testing
#ko_t=2 #inner k-fold cross validations
#ki_t=3 #inner k-fold cross validations
#response="sl_"
#pigsums
#X_vec=c(3:12,18:53)

#outer_kfold_loop(pigsums,response,X_vec,ko_t,ki_t)

outer_kfold_loop<- function(pigsums,response,X_vec,ko_t,ki_t,family,split_type,ntreemax){

####Start outer cross validation function here
#inputs: 
#ko_t
#ki_t
#pigsums
#response
#X_vec

#output: (as list)
#meanRMSE
#best.model.params

#OuterCrossValidatio
#initialize output data frame
outer.opt.params.kfold=data.frame(matrix(nrow=ko_t,ncol=7))
outer.r2.kfold=data.frame(matrix(nrow=ko_t,ncol=1))
colnames(outer.opt.params.kfold)=c("ko", "n.trees", "learning.rate", "tree.complexity", "bag.fraction", "rmse.inner", "rmse.outer")

#get outer k fold split data set
if(split_type=="random"){
outer_k_fold=k_split(pigsums,ko_t)
} else{
outer_k_fold=k_split_group(pigsums,split_type)
}

#gather set of predictions/test set for each top model of each outer replicate
final.pred.saver=vector(mode="list",length=ko_t)
final.test.saver=vector(mode="list",length=ko_t)

for(ko in 1:ko_t){
  print(paste("STARTING OUTER REPLICATE",ko))
  print("#############################")
outerpigtest=outer_k_fold[[ko]] #held out sample in outer loop
outerpigtrain=anti_join(pigsums,outerpigtest) #training sample in outer loop plus test sample in inner loop

inner.r2.kfold=data.frame(matrix(nrow=ki_t,ncol=1))
opt.params.kfold=inner_kfold_loop(outerpigtrain,ki_t,response,X_vec,family,ntreemax)
#after inner loop, we get our little optimal parameters table
#each row of the parameters table gives the optimal parameters for that train/test combo
#the row corresponds to which number in the list was the test set, the train set is all others combined

#now, use each set of optimal hyperparameters to test the train.outer set against the test.outer set
#output needed:
#1-chosen optimal hyperparameters
#2-averaged RMSE of ki*ko nested cross=validation for each outer ko-fold

#steps: 
#1-go through each row of optimal hyperparameters
#2-run a gbm with those hyperparameters using the outer training set
#3-get predictions
#4-get RMSE for outer train vs outer test for each optimal hyperparameter set

#gather predictions from each inner rep
pred.saver=vector(mode="list",length=nrow(opt.params.kfold))
test.saver=vector(mode="list",length=nrow(opt.params.kfold))
for(q in 1:nrow(opt.params.kfold)){
gbm.outer <- gbm.fixed(data=outerpigtrain, gbm.x=X_vec, gbm.y=which(colnames(outerpigtrain)==response),
                 learning.rate=opt.params.kfold$learning.rate[q], 
                 tree.complexity=opt.params.kfold$tree.complexity[q], 
                 n.trees=opt.params.kfold$n.trees[q],
                 bag.fraction = opt.params.kfold$bag.fraction[q],
                 family=family)  
print("predicting outer pred set")
pred.outer <- predict.gbm(gbm.outer, outerpigtest, opt.params.kfold$n.trees[q], "response")
test.saver[[q]]=outerpigtest
pred.saver[[q]]=pred.outer

#Do RMSE calcs and put in data frame
RMSE=sqrt(mean((outerpigtest[,response] - pred.outer)^2))
opt.params.kfold$rmse.outer[q]=RMSE

#Do R2 calcs and put in data frame
#rsq <- function (x, y) cor(x, y) ^ 2
R2=cor(outerpigtest[,response],pred.outer)^2
inner.r2.kfold[q,1]=R2

} #inner opt params loop closing bracket, testing on outer testing set

#get min RMSEindex from those runs against outer validation set
minRMSEindex=which(opt.params.kfold$rmse.outer==min(opt.params.kfold$rmse.outer))

#populate ko'th row of outer.opt.params.kfold set
outer.opt.params.kfold$ko[ko]=ko
outer.opt.params.kfold$n.trees[ko]=opt.params.kfold[minRMSEindex,]$n.trees
outer.opt.params.kfold$learning.rate[ko]=opt.params.kfold[minRMSEindex,]$learning.rate
outer.opt.params.kfold$tree.complexity[ko]=opt.params.kfold[minRMSEindex,]$tree.complexity
outer.opt.params.kfold$bag.fraction[ko]=opt.params.kfold[minRMSEindex,]$bag.fraction
outer.opt.params.kfold$rmse.inner[ko]=opt.params.kfold[minRMSEindex,]$rmse.inner
outer.opt.params.kfold$rmse.outer[ko]=opt.params.kfold[minRMSEindex,]$rmse.outer

outer.r2.kfold[ko,1]=inner.r2.kfold[minRMSEindex,]
print("got outer r2, min RMSE")
#use index to get preds/train values corresponding to the lowest RMSE
final.pred.saver[[ko]]=pred.saver[[minRMSEindex]]
final.test.saver[[ko]]=pred.saver[[minRMSEindex]]
print("assigned final pred/test savers")
} #outer k fold loop closing bracket


#get average RMSE from all the outer ko validation tests
meanRMSE=mean(outer.opt.params.kfold$rmse.outer)
meanR2=mean(outer.r2.kfold[,1])

#get best model hyperparameters
best.model.index=which(outer.opt.params.kfold$rmse.outer==min(outer.opt.params.kfold$rmse.outer))

#subset opt parameter matrix
best.model.params=outer.opt.params.kfold[best.model.index,]

output.list=list(meanRMSE ,best.model.params, meanR2, final.pred.saver, final.test.saver)

return(output.list)

} #function closing bracket


