#run lasso function
runLassofunc<-function(x,y,family,grid,typemeasure,nrep){
#optimizing alpha
  print("optimizing alpha")
best.alpha=optimize_alpha_lasso(x,y,"poisson",grid,"mae") #function here

df_coef_mat=matrix(nrow=(ncol(x)+1),ncol=nrep)

#starting lasso replicates with best alpha
  print("starting lasso reps with best alpha")
for(i in 1:nrep){
  print(paste0("rep ", i))
cv.lasso <- cv.glmnet(as.matrix(x), as.matrix(y), family=family, alpha=best.alpha, parallel=TRUE, standardize=TRUE, type.measure=typemeasure)
#from each rep...
df_coef_mat[,i] <- round(as.matrix(coef(cv.lasso, s=cv.lasso$lambda.min)), 2)
#matrix ncol nrep
}
rownames(df_coef_mat)=rownames(round(as.matrix(coef(cv.lasso, s=cv.lasso$lambda.min)), 2))  
varmeans=rowSums(df_coef_mat)/nrep
varmeans=varmeans[varmeans!=0]
return(varmeans)
}


