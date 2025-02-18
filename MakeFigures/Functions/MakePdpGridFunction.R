#testing
#numplots=12
#rel.inf=rel.inf.sigma.sl
#gbm.object=gbm.sigma.sl
#data=pigsums2
#opt.params=sigma.sl.opt.params.kfold
#ylims=c(-0.3,0.3)

#####Making above loop into a function for repeated use
MakePdpGrid<-function(numplots,rel.inf,gbm.object,data,opt.params,ylims){
  source(paste0(home,"/Scripts_Polished/Making_Plots_Tables/PlotFunctions/multiplotfunction.R"))
  myplots <- vector('list', length=numplots) #initiate empty plotlist
  vars=rownames(rel.inf)[1:numplots] #set vector of top numplots most important vars
  varneat=rel.inf[1:12,1]
  for (i in seq_along(vars)) {
    message(i) 
    myplots[[i]] <- local({
      pd=partial(gbm.object, pred.var = vars[i], train = data, n.trees = opt.params$n.trees,recursive=FALSE)
      pd$centred=pd[,2]-mean(pd[,2]) 
      i <- i
      if(is.numeric(pd[,1])){
        p1 <- ggplot(data=pd,aes(pd[,1],pd[,3]),inherit.aes=FALSE)+
          geom_line(colour="black")+
          theme_minimal()+
          ylim(ylims)+
          xlab(varneat[i])+
          ylab("centered yhat")+
          geom_rug(data=data,aes(x=data[,which(colnames(data)==vars[i])]),inherit.aes=FALSE)
        print(p1)
      } else{
        p1 <- ggplot(data=pd,aes(pd[,1],pd[,3]),inherit.aes=FALSE)+
          geom_boxplot(colour="black")+
          theme_minimal()+
          ylim(ylims)+
          xlab(varneat[i])+
          ylab("centered yhat")+
          geom_rug(data=data,aes(x=data[,which(colnames(data)==vars[i])]),inherit.aes=FALSE)
        print(p1)
      }
    })
  }
  
  multiplot=multiplot(plotlist = myplots, cols = 4)
  return(multiplot)
}
