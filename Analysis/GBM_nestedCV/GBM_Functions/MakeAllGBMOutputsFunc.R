#vars for input: repname, split_type, response
#response naming conventions:
  #colname from pigsums: (response)
    #sl_mean
    #sigma_sl **check
    #disp
    #sigma_disp **check
  #file paths: (response_str)
    #sl
    #sigma.sl
    #disp
    #sigma.disp
  #pigsums:
    #pigsums_sl
    #pigsums_sigmasl **
    #pigsums_disp
    #pigsums_sigmadisp **
#distribution
    #poisson -sl/disp
    #gaussian -sigmasl/disp

MakeAllGBMOutputs<-function(repname,split_type,pigsums,response,response_str,distribution){

  filestr=paste("4_Outputs",repname,split_type,sep="/")
  if(!dir.exists(file.path("4_Outputs/",repname, fsep = .Platform$file.sep))){dir.create(file.path("4_Outputs/",repname, fsep = .Platform$file.sep))}
  if(!dir.exists(filestr)){dir.create(filestr)}
  
#for testing
ko_t=10 #outer k-fold cross validations
ki_t=3 #inner k-fold cross validations

#X vec is a vector of indices with parameters included in the model
X_vec.start=c(
  which(colnames(pigsums)=="sex"),
  which(colnames(pigsums)=="season"),
  which(colnames(pigsums)=="period"),
  which(colnames(pigsums)=="mean_tc"):
  which(colnames(pigsums)=="var_lc_24")
        )

#Run models
res=Run.GBM.Model(pigsums,response,X_vec.start,ko_t,ki_t,distribution,split_type,ntreemax=8000)

#Export CV stats from first run
#set file.str and make directory for saving
write.csv(res[[2]],file.path(filestr,paste0(response_str,"_meanRMSE.csv")))
write.csv(res[[3]],file.path(filestr,paste0(response_str,"_bestmodelparams.csv")))
write.csv(res[[4]],file.path(filestr,paste0(response_str,"_meanR2.csv")))
saveRDS(res[[5]],file.path(filestr,paste0(response_str,"_preds.rds")))
saveRDS(res[[6]],file.path(filestr,paste0(response_str,"_testobs.rds")))


#Get x_vecs without unimportant variables
#Remove variables that had less than 0.1 influence
X_vec.2=drop.unimportant(res[[1]],X_vec.start,0.1)

#Save new vectors to file
write.csv(X_vec.2,paste0(filestr,"/X_vec_01Drop/",paste0("X",response_str,".csv")))

#Run models without dropped variables
res=Run.GBM.Model(pigsums,response,X_vec.2,ko_t,ki_t,distribution,split_type,ntreemax=8000)

#Export CV stats
write.csv(res[[2]],file.path(filestr,"GBM_01Drop",paste0(response_str,"_meanRMSE.csv")))
write.csv(res[[3]],file.path(filestr,"GBM_01Drop",paste0(response_str,"_bestmodelparams.csv")))
write.csv(res[[4]],file.path(filestr,"GBM_01Drop",paste0(response_str,"_meanR2.csv")))
saveRDS(res[[5]],file.path(filestr,"GBM_01Drop",paste0(response_str,"_preds.rds")))
saveRDS(res[[6]],file.path(filestr,"GBM_01Drop",paste0(response_str,"_testobs.rds")))

#Run lasso function for each model
#get list of remaining important variables to use
nrep=100
typemeasure="mae"

#some tidying of x, cant have characters
x=pigsums[,X_vec.start]
x$sexf=NA
x[x$sex=="Male",]$sexf=1
x[x$sex=="Female",]$sexf=0
x=x[,-2]

#make alpha search grid
grid=seq(0.01,1,by=0.01)

#Run lasso on each model
keepVars=runLassofunc(x,pigsums[,which(colnames(pigsums)==response)],distribution,grid,typemeasure,nrep)
X_vec.lasso=which(colnames(pigsums)%in%names(keepVars))

#Write out the post lasso X vecs
write.csv(X_vec.lasso,file.path(filestr,"X_vec_postlasso",paste0("X",response,"_lasso.csv")))

res=Run.GBM.Model(pigsums,response,X_vec.lasso,ko_t,ki_t,distribution,split_type,ntreemax=8000)
  
#Export CV stats
write.csv(res[[2]],file.path(filestr,"GBM_Lasso",paste0(response_str,"_meanRMSE.csv")))
write.csv(res[[3]],file.path(filestr,"GBM_Lasso",paste0(response_str,"_bestmodelparams.csv")))
write.csv(res[[4]],file.path(filestr,"GBM_Lasso",paste0(response_str,"_meanR2.csv")))
saveRDS(res[[5]],file.path(filestr,"GBM_Lasso",paste0(response_str,"_preds.rds")))
saveRDS(res[[6]],file.path(filestr,"GBM_Lasso",paste0(response_str,"_testobs.rds")))

pigsums$gbm.null=rep(1,nrow(pigsums))
X_null=ncol(pigsums)
  
res=Run.GBM.Model(pigsums,response,X_null,ko_t,ki_t,distribution,split_type,ntreemax=8000)
  
#Export CV stats
write.csv(res[[2]],paste0(filestr,"GBM_Null",paste0(response_str,"_meanRMSE.csv")))
write.csv(res[[3]],paste0(filestr,"GBM_Null",paste0(response_str,"_bestmodelparams.csv")))
write.csv(res[[4]],paste0(filestr,"GBM_Null",paste0(response_str,"_meanR2.csv")))

}
