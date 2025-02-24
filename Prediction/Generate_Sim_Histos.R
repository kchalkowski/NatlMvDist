
# Purpose ----------------------------

# Process ----------------------------

# Setup ----------------------------

#vars
filestr="TestRuns"

#load libraries
library(simstudy)
library(plyr)
library(dplyr)
library(ggplot2)

#set directories
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/StatPigMvmt/Pipeline_R2"
outdir<-file.path(home,"4_Outputs")
objdir<-file.path(home,"2_Data","Objects")

#read geolocation data
pigs<-readRDS(file.path(objdir,"allsteps_daily.rds"))

#read watersheds
wash<-readRDS(file.path(objdir,"wash_envcov_final.rds"))

#read in uncertainty matrices
unc.fs=list.files(file.path(outdir,filestr,"UncPredMats"),full.names=TRUE)
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
      print(paste0("wash: ", i))
      
      #make distribution from each mean and disp/var
      minibootwash[i,4]=mean(pmat[i,])
      minibootwash[i,5]=mean(pmat.disp[i,])
      gampams=simstudy::gammaGetShapeRate(minibootwash[i,4],minibootwash[i,5])
      
      #sample values from mean distribution
      minibootwash[i,6:105]=rgamma(100,shape=gampams$shape,rate=gampams$rate)
      
    }
    
    bootwash=rbind(bootwash,minibootwash)
    
  }
}


saveRDS(bootwash,paste0(file.path(outdir,filestr,"UncPredMats"),"bootwash.rds"))

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
DoAllHisto<-function(pigs,bootwash,response){
  coln=colnames(pigs)[grep(response,colnames(pigs))]
  pigs2=pigs[!is.na(pigs[,which(colnames(pigs)==coln)]),]
  
  for(q in 1:4){
    pigsq=pigs2[pigs2$season==as.factor(paste0("q",q)),which(colnames(pigs)==coln)]
    df=bootwash[bootwash$season==paste0("q",q)&bootwash$response==response,]
    df=data.frame(type="boots",response=unlist(c(df[,6:8])))
    
    pq=RelDiffHisto_util(top=df$response,bottom=pigsq,20,response,paste0("q",q))
    histname=paste0("hist_q",q,"_",response,".png")
    ggsave(paste0(outdir,filestr,histname),plot=pq,width=4,height=8)
    
  }
  
}

#Run histos for all q, response
DoAllHisto(pigs,bootwash,"sl")
DoAllHisto(pigs,bootwash,"disp")

# Do KS tests ---------------------

#vectors to go through ks stats
#steps
df2.vec.sl=list(df2.q1$sl, df2.q2$sl, df2.q3$sl, df2.q4$sl)
df2.vec.disp=list(df2.q1$disp, df2.q2$disp, df2.q3$disp, df2.q4$disp)
pig.vec.sl=list(pigs.q1$sl_, pigs.q2$sl_, pigs.q3$sl_, pigs.q4$sl_)
pig.vec.disp=list(pigs.q1[which(pigs.q1$displacement<1e6),]$displacement, 
                  pigs.q2[which(pigs.q2$displacement<1e6),]$displacement, 
                  pigs.q3[which(pigs.q3$displacement<1e6),]$displacement, 
                  pigs.q4[which(pigs.q4$displacement<1e6),]$displacement)

#means
pigs.u.sl.vec=list(pigs.q1.u.sl,pigs.q2.u.sl,pigs.q3.u.sl,pigs.q4.u.sl) #obs_u_sl/disp
pigs.u.di.vec=list(pigs.q1.u.di,pigs.q2.u.di,pigs.q3.u.di,pigs.q4.u.di) #obs_u_sl/disp
slq.vec=list(slq1,slq2,slq3,slq4) #sim_u_sl
dispq.vec=list(dispq1,dispq2,dispq3,dispq4) #sim_u_disp

ks.res=data.frame(matrix(nrow=4,ncol=6))
colnames(ks.res)=c("response","level","q1","q2","q3","q4")
ks.res[,1]=c("sl","sl","disp","disp")
ks.res[,2]=c("u","s","u","s")

for(q in 1:4){
model.u.sl=ks.test(pigs.u.sl.vec[[q]]$usl,slq.vec[[q]]$u) #0.28455
model.s.sl=ks.test(pig.vec.sl[[q]],df2.vec.sl[[q]]) #0.0796

model.u.di=ks.test(pigs.u.di.vec[[q]]$udisp,dispq.vec[[q]]$u) #0.28455
model.s.di=ks.test(pig.vec.disp[[q]],df2.vec.disp[[q]]) #0.0796

ks.res[1,(2+q)]=model.u.sl$statistic
ks.res[2,(2+q)]=model.s.sl$statistic
ks.res[3,(2+q)]=model.u.di$statistic
ks.res[4,(2+q)]=model.s.di$statistic

}

write.csv(ks.res,paste0(outdir,"ksres.csv"))

# Make qq plots ---------------------

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

#make combined qq for each season/response
#combined_qq(meanobs,meanpred,stepobs,steppred,quantiles = seq(0, 1, 0.01))
sl_1=combined_qq(pigs.u.sl.vec[[1]]$usl,slq.vec[[1]]$u,pig.vec.sl[[1]],df2.vec.sl[[1]],quantiles = seq(0, 1, 0.01),"sl_1")
sl_2=combined_qq(pigs.u.sl.vec[[2]]$usl,slq.vec[[2]]$u,pig.vec.sl[[2]],df2.vec.sl[[2]],quantiles = seq(0, 1, 0.01),"sl_2")
sl_3=combined_qq(pigs.u.sl.vec[[3]]$usl,slq.vec[[3]]$u,pig.vec.sl[[3]],df2.vec.sl[[3]],quantiles = seq(0, 1, 0.01),"sl_3")
sl_4=combined_qq(pigs.u.sl.vec[[4]]$usl,slq.vec[[4]]$u,pig.vec.sl[[4]],df2.vec.sl[[4]],quantiles = seq(0, 1, 0.01),"sl_4")

di_1=combined_qq(pigs.u.di.vec[[1]]$udisp,dispq.vec[[1]]$u,pig.vec.disp[[1]],df2.vec.disp[[1]],quantiles = seq(0, 1, 0.01),"disp_1")
di_2=combined_qq(pigs.u.di.vec[[2]]$udisp,dispq.vec[[2]]$u,pig.vec.disp[[2]],df2.vec.disp[[2]],quantiles = seq(0, 1, 0.01),"disp_2")
di_3=combined_qq(pigs.u.di.vec[[3]]$udisp,dispq.vec[[3]]$u,pig.vec.disp[[3]],df2.vec.disp[[3]],quantiles = seq(0, 1, 0.01),"disp_3")
di_4=combined_qq(pigs.u.di.vec[[4]]$udisp,dispq.vec[[4]]$u,pig.vec.disp[[4]],df2.vec.disp[[4]],quantiles = seq(0, 1, 0.01),"disp_4")

ggsave(paste0(outdir,"sl_1_qq.png"),sl_1,width=7,height=6)
ggsave(paste0(outdir,"sl_2_qq.png"),sl_2,width=7,height=6)
ggsave(paste0(outdir,"sl_3_qq.png"),sl_3,width=7,height=6)
ggsave(paste0(outdir,"sl_4_qq.png"),sl_4,width=7,height=6)

ggsave(paste0(outdir,"di_1_qq.png"),di_1,width=7,height=6)
ggsave(paste0(outdir,"di_2_qq.png"),di_2,width=7,height=6)
ggsave(paste0(outdir,"di_3_qq.png"),di_3,width=7,height=6)
ggsave(paste0(outdir,"di_4_qq.png"),di_4,width=7,height=6)



