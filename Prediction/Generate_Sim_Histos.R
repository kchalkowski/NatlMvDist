
# Purpose ----------------------------

# Process ----------------------------

# Setup ----------------------------

#vars
filestr="04APR25_Runs"

#load libraries
library(simstudy)
library(plyr)
library(dplyr)
library(ggplot2)

#set directories
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/2_Projects/StatPigMvmt/Pipeline_R2"
outdir<-file.path(home,"4_Outputs")
objdir<-file.path(home,"2_Data","Objects")
indir<-file.path(home,"2_Data","Input","Env_Cov")

#read geolocation data
pigs<-readRDS(file.path(objdir,"allsteps_daily.rds"))

#read watersheds
wash<-readRDS(file.path(indir,"wash3.rds"))

#read in uncertainty matrices
unc.fs=list.files(file.path(outdir,filestr,"FigTab","UncPredMats"),full.names=TRUE)

#keep quarterly ones
unc.fs=unc.fs[grep("q",unc.fs)]

#source functions, reldiffhistoutil
source(file.path(home,"1_Scripts","Prediction","Functions","reldiffhistos2.R"))

# Bootstrap watershed data ---------------------

#read in uncertainty matrices
#pmat.list=list.files(paste0(getwd(),"/Outputs/PredMats_22APR23/"))
season=c("q1","q2","q3","q4")
response.disp=c("sigmasl","sigmadisp")
response.mean=c("sl","disp")

#make empty mat of nrow=num watersheds, ncol=103 (season,response,objecti and 100 reps)
bootwash=matrix(nrow=0,ncol=105)
colnames(bootwash)<-c("season","response","objecti","u","disp",1:100)

#loop through season
for(s in 1:length(season)){
  print(paste0("season: ",season[s]))
  
  #loop through reasons
  for(r in 1:length(response.mean)){
    
    print(paste0("response: ", response.mean[r]))
    minibootwash=as.data.frame(matrix(nrow=nrow(wash),ncol=105))
    colnames(minibootwash)<-c("season","response","objecti","u","disp",1:100)
    minibootwash[,1]=season[s]
    minibootwash[,2]=response.mean[r]
    minibootwash[,3]=1:nrow(wash)

    print("reading pmat files")
    
    #pull needed season
    str=unc.fs[grep(paste0(s,".rds"),unc.fs)]
    #read pmats
    pmat=readRDS(str[grep(paste0("/",response.mean[r]),str)])
    pmat.disp=readRDS(str[grep(response.disp[r],str)])
    
    #for each watershed (row)
    for(i in 1:nrow(pmat)){
      if(i%%100==0){print(paste0("wash: ", i))}
      
      #make distribution from each mean and disp/var
      minibootwash[i,4]=mean(pmat[i,])
      minibootwash[i,5]=max(mean(pmat.disp[i,]),0) #switch in case disp is <0
      gampams=simstudy::gammaGetShapeRate(minibootwash[i,4],minibootwash[i,5])
      
      #sample values from mean distribution
      minibootwash[i,6:105]=rgamma(100,shape=gampams$shape,rate=gampams$rate)
      
    }
    
    
    bootwash=rbind(bootwash,minibootwash)
    
  }
  
  
  
}

#use for downstream tests
#saveRDS(bootwash,file.path(outdir,filestr,"UncPredMats","bootwash.rds"))
bootwash=readRDS(file.path(outdir,filestr,"UncPredMats","bootwash.rds"))


slq1=bootwash[bootwash$season=="q1"&bootwash$response=="sl",]
#dispq1=bootwash[bootwash$season=="q1"&bootwash$response=="disp",]
#slq2=bootwash[bootwash$season=="q2"&bootwash$response=="sl",]
#dispq2=bootwash[bootwash$season=="q2"&bootwash$response=="disp",]
#slq3=bootwash[bootwash$season=="q3"&bootwash$response=="sl",]
#dispq3=bootwash[bootwash$season=="q3"&bootwash$response=="disp",]
#slq4=bootwash[bootwash$season=="q4"&bootwash$response=="sl",]
#dispq4=bootwash[bootwash$season=="q4"&bootwash$response=="disp",]
df2.q1=data.frame(type="boots",sl=unlist(c(slq1[,6:8])))

#df2.q1=data.frame(type="boots",sl=unlist(c(slq1[,6:8])),disp=unlist(c(dispq1[,6:8])))
#df2.q2=data.frame(type="boots",sl=unlist(c(slq2[,6:8])),disp=unlist(c(dispq2[,6:8])))
#df2.q3=data.frame(type="boots",sl=unlist(c(slq3[,6:8])),disp=unlist(c(dispq3[,6:8])))
#df2.q4=data.frame(type="boots",sl=unlist(c(slq4[,6:8])),disp=unlist(c(dispq4[,6:8])))
#write.csv(df2.q1,paste0(outdir,"df2.q1.csv"))
#write.csv(df2.q2,paste0(outdir,"df2.q2.csv"))
#write.csv(df2.q3,paste0(outdir,"df2.q3.csv"))
#write.csv(df2.q4,paste0(outdir,"df2.q4.csv"))
#write.csv(pigs.q1,paste0(outdir,"pigs.q1.csv"))
#write.csv(pigs.q2,paste0(outdir,"pigs.q2.csv"))
#write.csv(pigs.q3,paste0(outdir,"pigs.q3.csv"))
#write.csv(pigs.q4,paste0(outdir,"pigs.q4.csv"))

#df2.qX are matrices, split up by response/q, with only predicted values

# Make histograms ---------------------
#Make histogram for each season/response [q1:4, sl/disp]
#simplify this

#Make function to make histo for each response/quarter
DoAllHisto<-function(pigs,bootwash,response,ybreaks=NULL){
  coln=colnames(pigs)[grep(response,colnames(pigs))]
  pigs2=pigs[!is.na(pigs[,which(colnames(pigs)==coln)]),]
  
  for(q in 1:4){
    pigsq=pigs2[pigs2$season==as.factor(paste0("q",q)),which(colnames(pigs)==coln)]
    df=bootwash[bootwash$season==paste0("q",q)&bootwash$response==response,]
    df=data.frame(type="boots",response=unlist(c(df[,6:8])))
    
    pq=RelDiffHisto_util(top=df$response,bottom=pigsq,20,response,paste0("q",q),ybreaks)
    histname=paste0("hist_q",q,"_",response,".png")
    ggsave(file.path(outdir,filestr,"FigTab","podiff_histos",histname),plot=pq,width=4,height=8)
    
  }
  
}

#Run histos for all q, response
DoAllHisto(pigs,bootwash,"sl",ybreaks=c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2))
DoAllHisto(pigs,bootwash,"disp",ybreaks=c(-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12))

# Make QQ plots and Run KS tests ---------------------

## QQ Plot functions ---------

#function to scale datasets
fun_range <- function(x,tot) {                              
  (x - min(tot)) / (max(tot) - min(tot))
}

#function to make combined qqplot 
combined_qq<-function(meanobs,meanpred,stepobs,steppred,quantiles = seq(0, 1, 0.01),response_season){
  tots.u=c(meanobs,meanpred)
  tots.s=c(stepobs,steppred)
  
  meanobs.scal=fun_range(meanobs,tots.u)
  meanpred.scal=fun_range(meanpred,tots.u)
  stepobs.scal=fun_range(stepobs,tots.s)
  steppred.scal=fun_range(steppred,tots.s)
  
  ggplot() + 
    geom_abline(aes(slope = 1, intercept = 0), linetype = 2, linewidth=4, color="black")+
    
    geom_point(mapping = aes(x = quantile(stepobs.scal, quantiles), 
                             y = quantile(steppred.scal, quantiles)),color="#db4200",size=6) +
    geom_point(mapping = aes(x = quantile(meanobs.scal, quantiles), 
                             y = quantile(meanpred.scal, quantiles)),color="#5318db",size=6) +
    xlim(0,1)+ylim(0,1)+coord_equal()+theme_linedraw() +
    xlab(response_season) +
    ylab("preds") +
    theme(axis.text.x = element_text(size=50,angle=90,vjust=0.4))+
    theme(axis.text.y = element_text(size=50))+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=5))
  
}

## Combined KS test/qq plot loop -------------

#response="disp"
Do.KS.Tests<-function(bootwash,pigs,response){
  coln=colnames(pigs)[grep(response,colnames(pigs))]
  pigs2=pigs[!is.na(pigs[,which(colnames(pigs)==coln)]),]
  
  sim.step.vec=vector(mode="list",length=4)
  pig.step.vec=vector(mode="list",length=4)
  sim.mu.vec=vector(mode="list",length=4)
  pig.mu.vec=vector(mode="list",length=4)
  ks.res=data.frame("response"=response,"q"=1:4,"mu_D"=NA,"step_D"=NA)
  for(q in 1:4){
    pigsq=pigs2[pigs2$season==as.factor(paste0("q",q)),]
    df=bootwash[bootwash$season==paste0("q",q)&bootwash$response==response,]
    
    df2=data.frame(type="boots",response=unlist(c(df[,6:8])))
    
    #Steps!
    sim.step.vec[[q]]=df2$response
    pig.step.vec[[q]]=pigsq[,which(colnames(pigs)==coln)]
    
    #Means!
    sim.mu.vec[[q]]=df$u
    pig.mu=pigsq %>% dplyr::group_by(animalid) %>% dplyr::summarise(mu=mean(eval(parse(text=coln)))) %>% dplyr::select(mu) %>% as.data.frame()
    pig.mu.vec[[q]]=pig.mu$mu
    
    model.mu=ks.test(pig.mu.vec[[q]],sim.mu.vec[[q]]) #0.28455
    model.step=ks.test(pig.step.vec[[q]],sim.step.vec[[q]]) #0.0796
    
    ks.res[q,3]=model.mu$statistic
    ks.res[q,4]=model.step$statistic
    
      response_season=paste0(response,"_",q)
      qq=combined_qq(pig.mu.vec[[q]],sim.mu.vec[[q]],pig.step.vec[[q]],sim.step.vec[[q]],quantiles = seq(0, 1, 0.01),response_season)
      ggsave(file.path(outdir,"04APR25_Runs","FigTab","qqplots",paste0(response,"_",q,"_qq.png")),qq,width=7,height=6)
    
  }
  
return(list("pig.mu.vec"=pig.mu.vec,
            "sim.mu.vec"=sim.mu.vec,
            "pig.step.vec"=pig.step.vec,
            "sim.step.vec"=sim.step.vec,
            "ks.res"=ks.res))
  
}

sl_ks_list=Do.KS.Tests(bootwash,pigs,"sl")
displ_ks_list=Do.KS.Tests(bootwash,pigs,"disp")

ks_res_out=rbind(sl_ks_list$ks.res,displ_ks_list$ks.res)

write.csv(ks_res_out,file.path(outdir,"04APR25_Runs","FigTab","kstabs","ks_table.csv"))

