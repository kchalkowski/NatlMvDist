k_split <- function(dat,k, seed = 1){
  set.seed(seed)
  grp <- rep(1:k, length.out = nrow(dat))
  klist=dat %>%
    mutate(grp = sample(grp, nrow(dat), replace = F)) |>
    group_split(grp)|>
    map(\(d) dplyr::select(d, -grp))
  klist=lapply(klist,as.data.frame)  
  
  return(klist)
}

k_split_group<- function(dat,grouping_var){
group_var_vec=unique(dat[,which(colnames(dat)==grouping_var)])
grouping_list=vector(mode="list",length=length(group_var_vec))
for(i in 1:length(group_var_vec)){
grouping_list[[i]]=dat[dat[,which(colnames(dat)==grouping_var)]==group_var_vec[i],]
}
return(grouping_list)
}

