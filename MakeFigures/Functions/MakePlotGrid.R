MakePlotGrid<-function(predobs.df2,studydf){
  plotlist=vector(mode="list",length=length(studydf[,1]))
  for(i in 1:length(studydf[,1])){
    plotlist[[i]]=predobs.df2[predobs.df2$region==studydf[i,1],]%>%
      ggplot(aes(x=value, fill=name)) +
      geom_density(alpha=0.3)+ 
      scale_x_log10()+
      labs(x= "Step Length")+
      facet_wrap(vars(CVmethod))+
      theme_minimal()
  }
  
  pg=plot_grid(plotlist[[1]],
            plotlist[[2]],
            plotlist[[3]],
            plotlist[[4]],
            plotlist[[5]],
            plotlist[[6]],
            plotlist[[7]],
            plotlist[[8]],
            plotlist[[9]],
            plotlist[[10]],
            plotlist[[11]],
            plotlist[[12]],
            plotlist[[13]],
            plotlist[[14]],
            plotlist[[15]],
            plotlist[[16]],
            plotlist[[17]],
            plotlist[[18]],
            plotlist[[19]],
            plotlist[[20]],
            plotlist[[21]],
            plotlist[[22]],
            plotlist[[23]],
            plotlist[[24]],
            plotlist[[25]],
            labels = as.character(region.counts[,1]),
            label_size = 10,
            ncol=6)
  
  return(pg)
  
}