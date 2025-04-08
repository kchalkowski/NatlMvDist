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
#response=responses[r]
#response_str=response_strings[r]
#distribution=distributions[r]
#pigsums=pigsums_list[[r]]
#repname="21FEB24_Runs"
#split_type="Region"
MakeAllGBMOutputs<-function(path,split_type,pigsums,response,response_str,distribution){

  #make addl folders
  if(!dir.exists(file.path(path,"X_vec_01Drop"))){dir.create(file.path(path,"X_vec_01Drop"))}
  if(!dir.exists(file.path(path,"GBM_01Drop"))){dir.create(file.path(path,"GBM_01Drop"))}
  if(!dir.exists(file.path(path,"GBM_Null"))){dir.create(file.path(path,"GBM_Null"))}
  
#for testing
ko_t=10 #outer k-fold cross validations
ki_t=2 #inner k-fold cross validations

#X vec is a vector of indices with parameters included in the model
X_vec.start=c(
  which(colnames(pigsums)=="sex"),
  which(colnames(pigsums)=="season"),
  which(colnames(pigsums)=="period"),
  which(colnames(pigsums)=="mean_tc"):
  which(colnames(pigsums)=="mean_rd3"),
  which(colnames(pigsums)=="mean_drt"):
  which(colnames(pigsums)=="var_lc_24")
        )

#Run models
res=Run.GBM.Model(pigsums,response,X_vec.start,ko_t,ki_t,distribution,split_type,ntreemax=8000)
#res=Run.GBM.Model(pigsums[1:1000,],response,X_vec.start[1:5],2,2,distribution,split_type,ntreemax=8000)

#Export CV stats from first run
#set file.str and make directory for saving
write.csv(res[[2]],file.path(path,paste0(response_str,"_meanRMSE.csv")))
write.csv(res[[3]],file.path(path,paste0(response_str,"_bestmodelparams.csv")))
write.csv(res[[4]],file.path(path,paste0(response_str,"_meanR2.csv")))
saveRDS(res[[5]],file.path(path,paste0(response_str,"_preds.rds")))
saveRDS(res[[6]],file.path(path,paste0(response_str,"_testobs.rds")))


#Get x_vecs without unimportant variables
#Remove variables that had less than 0.1 influence
X_vec.2=drop.unimportant(res[[1]],X_vec.start,0.1)

#Save new vectors to file
write.csv(X_vec.2,file.path(path,"X_vec_01Drop",paste0("X",response_str,".csv")))

#Run models without dropped variables
res=Run.GBM.Model(pigsums,response,X_vec.2,ko_t,ki_t,distribution,split_type,ntreemax=8000)

#Export CV stats
write.csv(res[[2]],file.path(path,"GBM_01Drop",paste0(response_str,"_meanRMSE.csv")))
write.csv(res[[3]],file.path(path,"GBM_01Drop",paste0(response_str,"_bestmodelparams.csv")))
write.csv(res[[4]],file.path(path,"GBM_01Drop",paste0(response_str,"_meanR2.csv")))
saveRDS(res[[5]],file.path(path,"GBM_01Drop",paste0(response_str,"_preds.rds")))
saveRDS(res[[6]],file.path(path,"GBM_01Drop",paste0(response_str,"_testobs.rds")))

#Run lasso function for each model
#get list of remaining important variables to use
nrep=100
typemeasure="mae"

pigsums$gbm.null=rep(1,nrow(pigsums))
X_null=ncol(pigsums)
  
res=Run.GBM.Model(pigsums,response,X_null,ko_t,ki_t,distribution,split_type,ntreemax=8000)
  
#Export CV stats
write.csv(res[[2]],paste0(path,"GBM_Null",paste0(response_str,"_meanRMSE.csv")))
write.csv(res[[3]],paste0(path,"GBM_Null",paste0(response_str,"_bestmodelparams.csv")))
write.csv(res[[4]],paste0(path,"GBM_Null",paste0(response_str,"_meanR2.csv")))

}
