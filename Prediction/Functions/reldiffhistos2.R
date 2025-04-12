#Makes top-bottom histogram figures for predictive map figures


#pq=RelDiffHisto_util(top=df$response,bottom=pigsq,20,response,paste0("q",q),ybreaks)
top=df$response
bottom=pigsq
bins=20
response="sl"
quarter=paste0("q",q)
RelDiffHisto_util<-function(top,bottom,bins,response,quarter,ybreaks=NULL){

  #make dataframe format
  topdf=data.frame(name="top",resp=top)
  bottomdf=data.frame(name="bottom",resp=bottom)
  
  #cut values into bins
  topcut<-topdf %>%
    mutate(bin = cut(top, pretty(top, bins)))
  bottomcut<- bottomdf %>%
    mutate(bin = cut(bottom, pretty(bottom, bins)))
  
  #combine top/bottom to single df
  all=rbind(topdf,bottomdf)
  
  #get frequencies of each val, and convert bottom freqs to negative nums
  #negative=plots on bottom/upside-down
  washvec.df2 <-
    all %>%
    mutate(bin = cut(resp, pretty(resp, bins))) %$%
    table(bin, name) %>%
    as.data.frame() %>%
    mutate(plotVal = ifelse(name == "bottom"
                            , -1*Freq
                            , Freq))
  
  #get total values of each frequency set
  bottomtotal<-sum(washvec.df2[washvec.df2$name=="bottom",]$plotVal)
  toptotal<-sum(washvec.df2[washvec.df2$name=="top",]$plotVal)
  
  #convert frequencies to percentages
  washvec.df2$plotVal2<-0
  washvec.df2[washvec.df2$name=="bottom",]$plotVal2=(washvec.df2[washvec.df2$name=="bottom",]$plotVal/bottomtotal)*100
  washvec.df2[washvec.df2$name=="top",]$plotVal2=(washvec.df2[washvec.df2$name=="top",]$plotVal/toptotal)*100

  #trim unneeded interim columns
  washvec.df3<-washvec.df2[,c(1,2,5)]
  
  #pivot wider for plotting
  washvec.df4=as.data.frame(pivot_wider(washvec.df3,names_from=name,values_from=plotVal2))
  
  #bottom needs be negative again
  washvec.df4$bottom=-1*washvec.df4$bottom
  
  #get differences
  washvec.df4$podiff=washvec.df4$top+washvec.df4$bottom
  
  #get key for plotting
  washvec.df4$podiff2=if_else(washvec.df4$podiff>0,"pos","neg")
  
  #make scale if ybreaks NULL
  if(is.null(ybreaks)){
  ybreaks=c(0,max(washvec.df4$podiff))
  }
  
  p<-ggplot(washvec.df4) + 
    geom_col(position="identity",alpha=1,mapping = aes(x=bin, y = podiff,fill=podiff2))+
    scale_fill_manual(values=c("#C0C0C0","#F96545"))+
    theme(legend.position = "none")+
    xlab(paste(response,quarter,sep="_"))+
    scale_y_continuous(breaks = ybreaks, limits=c(min(ybreaks),max(ybreaks)))+
    theme(
      panel.background = element_rect(fill='transparent'),
      panel.border = element_blank(),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.background = element_rect(fill='transparent'),
      legend.box.background = element_rect(fill='transparent'))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme(axis.text.y = element_text(size=20),axis.title.y=element_text(size=20))
  return(p)
}
