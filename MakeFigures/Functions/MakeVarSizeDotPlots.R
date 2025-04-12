#The purpose of this script is to make dotplots for RMSE and R2 values by study
#with two different geom point aesthetics, and including two different legends for each

#first, this function for extracting a legend
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}


#CVstats=CVstats_sl.random
#studycounts.df=studycounts
#response="RMSE"
#lab="sl"
make.varsize.dotplots<-function(CVstats, studycounts.df, response, lab){
if(length(unique(CVstats$type))>1){
colors=c("#039e8c", "#832b9e") 
} else{
  if(unique(CVstats$type)=="region"){
    colors="#832b9e"
  } else{
    if(unique(CVstats$type)=="random"){
      colors="#039e8c"
    }
  }
}
  
if(lab=="sl"){xlab=paste0("Daily step length mean prediction ",response)}
if(lab=="sigmasl"){xlab=paste0("Daily step length dispersion prediction ",response)}
if(lab=="disp"){xlab=paste0("Daily displacement mean prediction ",response)}
if(lab=="sigmadisp"){xlab=paste0("Daily displacement dispersion prediction ",response)}

  #input vars needed: CVstats df, studycounts
  CVstats=left_join(CVstats, studycounts.df, by="region")
  
  ###Make dotplot of prediction RMSE and R2
  p1=ggplot(CVstats, aes(CVstats[,which(colnames(CVstats)==response)], region,color=type)) + 
    geom_point(aes(color=type, size=numpigs))+
    scale_size_continuous(range=c(0.1,0.2))
  p2=ggplot(CVstats, aes(CVstats[,which(colnames(CVstats)==response)], region,color=type)) + 
    geom_point(aes(size=numpigdays), shape=1)+
  scale_size_continuous(range=c(0.1,0.2))
  legend1 <- g_legend(p1)
  legend2 <- g_legend(p2)
  legend.width <- sum(legend2$width)  
  
  p3=ggplot(CVstats, aes(CVstats[,which(colnames(CVstats)==response)],region,color=type)) + 
    geom_point(aes(color=type, size=numpigs)) + 
    scale_size_continuous(range = c(3,20))+
    geom_point(aes(size=numpigdays), shape=1)+
    facet_grid(State ~ ., scales = "free", space = "free") +
    theme(strip.text.y = element_text(angle = 0))+
    theme_minimal()+theme(panel.border = element_rect(color = "gray 50", fill = NA))+
    scale_color_manual(values=colors, name="CV method")+
    labs(y="Study", x=xlab)+theme(text = element_text(size = 35))+
    theme(panel.spacing = unit(1, "lines"))+
    theme(strip.text.y.right = element_text(angle = 0))+
    labs(x=response)
  
  gplot<-grid.arrange(p3 +theme(legend.position = "none"), legend1, legend2,
                        ncol = 2, nrow = 2,
                        layout_matrix = rbind(c(1,2 ),  
                                              c(1,3 )), 
                        widths = unit.c(unit(1, "npc") - legend.width, legend.width))
return(gplot)
}
