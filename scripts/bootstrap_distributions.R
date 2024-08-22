#####################
### Adaptive loci ###
#####################

library(ggpubr)

path=getwd()

allele_freqs<-read.csv(file=file.path(path, "Re_analysis_2023", "PopGen", "AF_preds", "AllFreq.csv"))
Sites<-allele_freqs$Group.1
allele_freqs<-allele_freqs[,c(2:7672)]
allele_freqs<-allele_freqs[, colMeans(allele_freqs, na.rm = TRUE) >= 0.05 &
                             colMeans(allele_freqs, na.rm = TRUE) <= 0.95, drop = FALSE]
allele_freqs$Group.1<-Sites
invasive<-c("WP","HM","CP","TO","RB")
inv_source<-c("WP","HM","CP","TO","RB","WH","CT")
native<-c("DM","GB","WH","CT","OY","DB","BF","HP","FB","SK")

invasive_freqs<-allele_freqs[allele_freqs$Group.1 %in% inv_source,]
ATL_freqs<-allele_freqs[allele_freqs$Group.1 %in% native,]

Env<-read.csv(file=file.path(path, "Re_analysis_2023", "PopGen", "environment_data", "In_situ", "environmental_summary_stats_2023.08.09.csv"))
Env$X<-NULL
# Max
##################################################################################################

### Maximum SST: 
#max_cor<-read.csv(file=file.path(path, "Re_analysis_2023", "PopGen", "AF_preds", "max_corr.csv"))
#max_cor_50<-max_cor[max_cor$correlation >= 0.5,] # N = 54
max_cor_loci<-gsub(":",".",paste0("X",max_comb))
max_cor_AFs<-allele_freqs[colnames(allele_freqs)  %in% max_cor_loci]
max_cor_AFs<-as.data.frame(cbind(allele_freqs$Group.1, max_cor_AFs))
colnames(max_cor_AFs)[1]<-"Pop"
Max<-Env[,c(1,3)]

### T-tests between invasive and source locations: ###
max_list<-list()
for(loc in 1:5){
  site=invasive[loc]
  tframe<-as.data.frame(t(max_cor_AFs[max_cor_AFs$Pop == "CT" |
                        max_cor_AFs$Pop == site,]))
  colnames(tframe)<-c(tframe[1,1], tframe[1,2])
  tframe<-tframe[2:nrow(tframe),]
  test1<-t.test(as.numeric(tframe[,1]), as.numeric(tframe[,2]), paired = T)
  
  tframe<-as.data.frame(t(max_cor_AFs[max_cor_AFs$Pop == "WH" |
                                        max_cor_AFs$Pop == site,]))
  colnames(tframe)<-c(tframe[1,1], tframe[1,2])
  tframe<-tframe[2:nrow(tframe),]
  test2<-t.test(as.numeric(tframe[,1]), as.numeric(tframe[,2]), paired = T)
  
  row1<-as.data.frame(cbind(paste0(site),paste0("CT"), test1$parameter, test1$statistic, test1$p.value))
  row2<-as.data.frame(cbind(paste0(site),paste0("WH"), test2$parameter, test2$statistic, test2$p.value))
  colnames(row1)<-colnames(row2)<-c("Invasive","Source","df","T-stat","P-value")
  comb<-as.data.frame(rbind(row1,row2))
  comb$Set<-"Max_SST"
  
  max_list[[loc]]<-comb
}
max_ttest<-do.call(rbind, max_list)

max_list<-list()
for(locus in 1:length(max_cor_loci)){
  df1<-as.data.frame(max_cor_AFs[max_cor_AFs$Pop %in% native, c(1,locus+1)])
  df1<-merge(df1, Max, by.x = "Pop", by.y = "Site")
  regr<-lm(df1[,2] ~ df1[,3])
  df2<-max_cor_AFs[max_cor_AFs$Pop %in% invasive, c(1,locus+1)]
  df2<-merge(df2, Max, by.x = "Pop", by.y = "Site")
  df2$Residual<-abs(df2[,2] - (regr$coefficients[2] * df2$annual_max + regr$coefficients[1]))
  df2$CT_diff<-df2[,2] - df1[2,2]
  df2$WH_diff<-df2[,2] - df1[9,2]
  colnames(df2)[2]<-"AF"
  max_list[[locus]]<-df2
}
max_loci<-do.call(rbind, max_list)
Residual_agg<-cbind(aggregate(max_loci$Residual ~ max_loci$Pop, FUN = "median"), "Residual")
colnames(Residual_agg)<-c("Pop","Deviation","Measure")
CT_agg<-cbind(aggregate(max_loci$CT_diff ~ max_loci$Pop, FUN = "median"), "CT_diff")
colnames(CT_agg)<-c("Pop","Deviation","Measure")
WH_agg<-cbind(aggregate(max_loci$WH_diff ~ max_loci$Pop, FUN = "median"), "WH_diff")
colnames(WH_agg)<-c("Pop","Deviation","Measure")
max_median<-as.data.frame(rbind(Residual_agg, CT_agg, WH_agg))


medians<-list()
rand_max_ttest<-list()
for(m in 1:1000){
  random_num<-sample(ncol(allele_freqs)-1, length(max_cor_loci), replace = F)
  rand_AFs<-allele_freqs[,c(random_num)]
  rand_AFs<-as.data.frame(cbind(Sites, rand_AFs))
  colnames(rand_AFs)[1]<-"Pop"
  rand_max_ttest2<-list()
  for(b in 1:5){
    site=invasive[b]
    rframe<-as.data.frame(t(rand_AFs[rand_AFs$Pop == "CT" |
                                          rand_AFs$Pop == site,]))
    colnames(rframe)<-c(rframe[1,1], rframe[1,2])
    rframe<-rframe[2:nrow(rframe),]
    test1<-t.test(as.numeric(rframe[,1]), as.numeric(rframe[,2]), paired = T)
    
    rframe<-as.data.frame(t(rand_AFs[rand_AFs$Pop == "WH" |
                                          rand_AFs$Pop == site,]))
    colnames(rframe)<-c(rframe[1,1], rframe[1,2])
    rframe<-rframe[2:nrow(rframe),]
    test2<-t.test(as.numeric(rframe[,1]), as.numeric(rframe[,2]), paired = T)
    
    row1<-as.data.frame(cbind(paste0(site),paste0("CT"), test1$parameter, test1$statistic, test1$p.value))
    row2<-as.data.frame(cbind(paste0(site),paste0("WH"), test2$parameter, test2$statistic, test2$p.value))
    colnames(row1)<-colnames(row2)<-c("Invasive","Source","df","T-stat","P-value")
    comb<-as.data.frame(rbind(row1,row2))
    comb$Set<-"Max_SST"
    
    rand_max_ttest2[[b]]<-comb
  }
  locus_max_ttest<-do.call(rbind, rand_max_ttest2)
  rand_max_ttest[[m]]<-locus_max_ttest
  
  locus_list<-list()
  for(locus in 1:length(max_cor_loci)){
    df1<-as.data.frame(rand_AFs[rand_AFs$Pop %in% native, c(1,locus+1)])
    df1<-merge(df1, Max, by.x = "Pop", by.y = "Site")
    regr<-lm(df1[,2] ~ df1[,3])
    df2<-rand_AFs[rand_AFs$Pop %in% invasive, c(1,locus+1)]
    df2<-merge(df2, Max, by.x = "Pop", by.y = "Site")
    df2$Residual<-abs(df2[,2] - (regr$coefficients[2] * df2$annual_max + regr$coefficients[1]))
    df2$CT_diff<-df2[,2] - df1[2,2]
    df2$WH_diff<-df2[,2] - df1[9,2]
    colnames(df2)[2]<-"AF"
    locus_list[[locus]]<-df2
  }
  comb_loci<-do.call(rbind, locus_list)
  Residual_agg<-cbind(aggregate(comb_loci$Residual ~ comb_loci$Pop, FUN = "median"), "Residual")
  colnames(Residual_agg)<-c("Pop","Deviation","Measure")
  CT_agg<-cbind(aggregate(comb_loci$CT_diff ~ comb_loci$Pop, FUN = "median"), "CT_diff")
  colnames(CT_agg)<-c("Pop","Deviation","Measure")
  WH_agg<-cbind(aggregate(comb_loci$WH_diff ~ comb_loci$Pop, FUN = "median"), "WH_diff")
  colnames(WH_agg)<-c("Pop","Deviation","Measure")
  median<-as.data.frame(rbind(Residual_agg, CT_agg, WH_agg))
  medians[[m]]<-median
  }
max_distributions<-do.call(rbind, medians)
max_ttest_distr<-do.call(rbind, rand_max_ttest)
max_ttest_distr$`P-value`<-as.numeric(max_ttest_distr$`P-value`)

max_ttest_comps1<-list()
for(q in 1:5){
  site=invasive[q]
  popnCT<-max_ttest_distr[max_ttest_distr$Invasive == site &
                           max_ttest_distr$Source == "CT",]
  CT_ttest<-as.numeric(max_ttest[max_ttest$Invasive == site &
                        max_ttest$Source == "CT",4])
  CT_more<-nrow(popnCT[abs(as.numeric(popnCT$`T-stat`)) >= abs(CT_ttest),])
  
  popnWH<-max_ttest_distr[max_ttest_distr$Invasive == site &
                            max_ttest_distr$Source == "WH",]
  WH_ttest<-as.numeric(max_ttest[max_ttest$Invasive == site &
                                   max_ttest$Source == "WH",4])
  WH_more<-nrow(popnWH[abs(as.numeric(popnWH$`T-stat`)) >= abs(WH_ttest),])
  dat<-as.data.frame(rbind(cbind(site, "CT", CT_ttest, CT_more/1000),
                           cbind(site, "WH", WH_ttest, WH_more/1000)))
  colnames(dat)<-c("Invasive", "Source", "T-test P-value", "Prop. Rand. Higher P-value")
  max_ttest_comps1[[q]]<-dat
}
max_ttest_comps<-do.call(rbind, max_ttest_comps1)
max_ttest_comps$Set<-"Max_SST"
max_ttest_distr$`T-stat`<-as.numeric(max_ttest_distr$`T-stat`)

ggdensity(data = max_ttest_distr, x = "`T-stat`", y = "..density..",
          facet.by = "Invasive", fill = "Source")

max_plots<-list()
for(q in 1:5){
  pop<-invasive[q]
  distr<-max_distributions[max_distributions$Pop == pop,]
  outliers<-max_median[max_median$Pop == pop,]
  max_plot<-ggdensity(distr, x = "Deviation", y = "..density..", fill = "Measure",
                      title = paste0(pop)) +
    geom_vline(xintercept =  outliers$Deviation[1], col = "green", linewidth = 1.2) + 
    geom_vline(xintercept =  outliers$Deviation[2], col = "red", linewidth = 1.2, lty = 2) +
    geom_vline(xintercept =  outliers$Deviation[3], col = "blue", linewidth = 1.2, lty = 2)
    
  ggsave(filename = file.path(path, "Re_analysis_2023", "PopGen", "AF_preds",
                              "Distribution_plots", "max", paste0(pop,".png")),
         width = 8, height = 6)
  max_plots[[q]]<-max_plot
}
maxs<-ggarrange(plotlist = max_plots, common.legend = T)
ggsave(maxs, filename = file.path(path, "Re_analysis_2023", "PopGen", "AF_preds",
                                  "Distribution_plots", "max", "All_pops_max.png"),
       width = 16, height = 10)


# Min
##################################################################################################

### minimum SST: 
#min_cor<-read.csv(file=file.path(path, "Re_analysis_2023", "PopGen", "AF_preds", "min_corr.csv"))
#min_cor_50<-min_cor[min_cor$correlation >= 0.5,] # N = 54
min_cor_loci<-gsub(":",".",paste0("X",min_comb))
min_cor_AFs<-allele_freqs[colnames(allele_freqs)  %in% min_cor_loci]
min_cor_AFs<-as.data.frame(cbind(allele_freqs$Group.1, min_cor_AFs))
colnames(min_cor_AFs)[1]<-"Pop"
min<-Env[,c(1,3)]

### T-tests between invasive and source locations: ###
min_list<-list()
for(loc in 1:5){
  site=invasive[loc]
  tframe<-as.data.frame(t(min_cor_AFs[min_cor_AFs$Pop == "CT" |
                                        min_cor_AFs$Pop == site,]))
  colnames(tframe)<-c(tframe[1,1], tframe[1,2])
  tframe<-tframe[2:nrow(tframe),]
  test1<-t.test(as.numeric(tframe[,1]), as.numeric(tframe[,2]), paired = T)
  
  tframe<-as.data.frame(t(min_cor_AFs[min_cor_AFs$Pop == "WH" |
                                        min_cor_AFs$Pop == site,]))
  colnames(tframe)<-c(tframe[1,1], tframe[1,2])
  tframe<-tframe[2:nrow(tframe),]
  test2<-t.test(as.numeric(tframe[,1]), as.numeric(tframe[,2]), paired = T)
  
  row1<-as.data.frame(cbind(paste0(site),paste0("CT"), test1$parameter, test1$statistic, test1$p.value))
  row2<-as.data.frame(cbind(paste0(site),paste0("WH"), test2$parameter, test2$statistic, test2$p.value))
  colnames(row1)<-colnames(row2)<-c("Invasive","Source","df","T-stat","P-value")
  comb<-as.data.frame(rbind(row1,row2))
  comb$Set<-"min_SST"
  
  min_list[[loc]]<-comb
}
min_ttest<-do.call(rbind, min_list)

min_list<-list()
for(locus in 1:length(min_cor_loci)){
  df1<-as.data.frame(min_cor_AFs[min_cor_AFs$Pop %in% native, c(1,locus+1)])
  df1<-merge(df1, min, by.x = "Pop", by.y = "Site")
  regr<-lm(df1[,2] ~ df1[,3])
  df2<-min_cor_AFs[min_cor_AFs$Pop %in% invasive, c(1,locus+1)]
  df2<-merge(df2, min, by.x = "Pop", by.y = "Site")
  df2$Residual<-abs(df2[,2] - (regr$coefficients[2] * df2$annual_min + regr$coefficients[1]))
  df2$CT_diff<-df2[,2] - df1[2,2]
  df2$WH_diff<-df2[,2] - df1[9,2]
  colnames(df2)[2]<-"AF"
  min_list[[locus]]<-df2
}
min_loci<-do.call(rbind, min_list)
Residual_agg<-cbind(aggregate(min_loci$Residual ~ min_loci$Pop, FUN = "median"), "Residual")
colnames(Residual_agg)<-c("Pop","Deviation","Measure")
CT_agg<-cbind(aggregate(min_loci$CT_diff ~ min_loci$Pop, FUN = "median"), "CT_diff")
colnames(CT_agg)<-c("Pop","Deviation","Measure")
WH_agg<-cbind(aggregate(min_loci$WH_diff ~ min_loci$Pop, FUN = "median"), "WH_diff")
colnames(WH_agg)<-c("Pop","Deviation","Measure")
min_median<-as.data.frame(rbind(Residual_agg, CT_agg, WH_agg))


medians<-list()
rand_min_ttest<-list()
for(m in 1:1000){
  random_num<-sample(ncol(allele_freqs)-1, length(min_cor_loci), replace = F)
  rand_AFs<-allele_freqs[,c(random_num)]
  rand_AFs<-as.data.frame(cbind(Sites, rand_AFs))
  colnames(rand_AFs)[1]<-"Pop"
  rand_min_ttest2<-list()
  for(b in 1:5){
    site=invasive[b]
    rframe<-as.data.frame(t(rand_AFs[rand_AFs$Pop == "CT" |
                                       rand_AFs$Pop == site,]))
    colnames(rframe)<-c(rframe[1,1], rframe[1,2])
    rframe<-rframe[2:nrow(rframe),]
    test1<-t.test(as.numeric(rframe[,1]), as.numeric(rframe[,2]), paired = T)
    
    rframe<-as.data.frame(t(rand_AFs[rand_AFs$Pop == "WH" |
                                       rand_AFs$Pop == site,]))
    colnames(rframe)<-c(rframe[1,1], rframe[1,2])
    rframe<-rframe[2:nrow(rframe),]
    test2<-t.test(as.numeric(rframe[,1]), as.numeric(rframe[,2]), paired = T)
    
    row1<-as.data.frame(cbind(paste0(site),paste0("CT"), test1$parameter, test1$statistic, test1$p.value))
    row2<-as.data.frame(cbind(paste0(site),paste0("WH"), test2$parameter, test2$statistic, test2$p.value))
    colnames(row1)<-colnames(row2)<-c("Invasive","Source","df","T-stat","P-value")
    comb<-as.data.frame(rbind(row1,row2))
    comb$Set<-"min_SST"
    
    rand_min_ttest2[[b]]<-comb
  }
  locus_min_ttest<-do.call(rbind, rand_min_ttest2)
  rand_min_ttest[[m]]<-locus_min_ttest
  
  locus_list<-list()
  for(locus in 1:length(min_cor_loci)){
    df1<-as.data.frame(rand_AFs[rand_AFs$Pop %in% native, c(1,locus+1)])
    df1<-merge(df1, min, by.x = "Pop", by.y = "Site")
    regr<-lm(df1[,2] ~ df1[,3])
    df2<-rand_AFs[rand_AFs$Pop %in% invasive, c(1,locus+1)]
    df2<-merge(df2, min, by.x = "Pop", by.y = "Site")
    df2$Residual<-abs(df2[,2] - (regr$coefficients[2] * df2$annual_min + regr$coefficients[1]))
    df2$CT_diff<-df2[,2] - df1[2,2]
    df2$WH_diff<-df2[,2] - df1[9,2]
    colnames(df2)[2]<-"AF"
    locus_list[[locus]]<-df2
  }
  comb_loci<-do.call(rbind, locus_list)
  Residual_agg<-cbind(aggregate(comb_loci$Residual ~ comb_loci$Pop, FUN = "median"), "Residual")
  colnames(Residual_agg)<-c("Pop","Deviation","Measure")
  CT_agg<-cbind(aggregate(comb_loci$CT_diff ~ comb_loci$Pop, FUN = "median"), "CT_diff")
  colnames(CT_agg)<-c("Pop","Deviation","Measure")
  WH_agg<-cbind(aggregate(comb_loci$WH_diff ~ comb_loci$Pop, FUN = "median"), "WH_diff")
  colnames(WH_agg)<-c("Pop","Deviation","Measure")
  median<-as.data.frame(rbind(Residual_agg, CT_agg, WH_agg))
  medians[[m]]<-median
}
min_distributions<-do.call(rbind, medians)
min_ttest_distr<-do.call(rbind, rand_min_ttest)
min_ttest_distr$`P-value`<-as.numeric(min_ttest_distr$`P-value`)

min_ttest_comps1<-list()
for(q in 1:5){
  site=invasive[q]
  popnCT<-min_ttest_distr[min_ttest_distr$Invasive == site &
                            min_ttest_distr$Source == "CT",]
  CT_ttest<-as.numeric(min_ttest[min_ttest$Invasive == site &
                                   min_ttest$Source == "CT",4])
  CT_more<-nrow(popnCT[abs(as.numeric(popnCT$`T-stat`)) >= abs(CT_ttest),])
  
  popnWH<-min_ttest_distr[min_ttest_distr$Invasive == site &
                            min_ttest_distr$Source == "WH",]
  WH_ttest<-as.numeric(min_ttest[min_ttest$Invasive == site &
                                   min_ttest$Source == "WH",4])
  WH_more<-nrow(popnWH[abs(as.numeric(popnWH$`T-stat`)) >= abs(WH_ttest),])
  dat<-as.data.frame(rbind(cbind(site, "CT", CT_ttest, CT_more/1000),
                           cbind(site, "WH", WH_ttest, WH_more/1000)))
  colnames(dat)<-c("Invasive", "Source", "T-test P-value", "Prop. Rand. Higher P-value")
  min_ttest_comps1[[q]]<-dat
}
min_ttest_comps<-do.call(rbind, min_ttest_comps1)
min_ttest_comps$Set<-"min_SST"
min_ttest_distr$`T-stat`<-as.numeric(min_ttest_distr$`T-stat`)

ggdensity(data = min_ttest_distr, x = "`T-stat`", y = "..density..",
          facet.by = "Invasive", fill = "Source")

min_plots<-list()
for(q in 1:5){
  pop<-invasive[q]
  distr<-min_distributions[min_distributions$Pop == pop,]
  outliers<-min_median[min_median$Pop == pop,]
  min_plot<-ggdensity(distr, x = "Deviation", y = "..density..", fill = "Measure",
                      title = paste0(pop)) +
    geom_vline(xintercept =  outliers$Deviation[1], col = "green", linewidth = 1.2) + 
    geom_vline(xintercept =  outliers$Deviation[2], col = "red", linewidth = 1.2, lty = 2) +
    geom_vline(xintercept =  outliers$Deviation[3], col = "blue", linewidth = 1.2, lty = 2)
  
  ggsave(filename = file.path(path, "Re_analysis_2023", "PopGen", "AF_preds",
                              "Distribution_plots", "min", paste0(pop,".png")),
         width = 8, height = 6)
  min_plots[[q]]<-min_plot
}
mins<-ggarrange(plotlist = min_plots, common.legend = T)
ggsave(mins, filename = file.path(path, "Re_analysis_2023", "PopGen", "AF_preds",
                                  "Distribution_plots", "min", "All_pops_min.png"),
       width = 16, height = 10)

# Mean
########################################################################################

### meanimum SST: 
#mean_cor<-read.csv(file=file.path(path, "Re_analysis_2023", "PopGen", "AF_preds", "mean_corr.csv"))
#mean_cor_50<-mean_cor[mean_cor$correlation >= 0.5,] # N = 54
mean_cor_loci<-gsub(":",".",paste0("X",mean_comb))
mean_cor_AFs<-allele_freqs[colnames(allele_freqs)  %in% mean_cor_loci]
mean_cor_AFs<-as.data.frame(cbind(allele_freqs$Group.1, mean_cor_AFs))
colnames(mean_cor_AFs)[1]<-"Pop"
mean<-Env[,c(1,4)]

### T-tests between invasive and source locations: ###
mean_list<-list()
for(loc in 1:5){
  site=invasive[loc]
  tframe<-as.data.frame(t(mean_cor_AFs[mean_cor_AFs$Pop == "CT" |
                                        mean_cor_AFs$Pop == site,]))
  colnames(tframe)<-c(tframe[1,1], tframe[1,2])
  tframe<-tframe[2:nrow(tframe),]
  test1<-t.test(as.numeric(tframe[,1]), as.numeric(tframe[,2]), paired = T)
  
  tframe<-as.data.frame(t(mean_cor_AFs[mean_cor_AFs$Pop == "WH" |
                                        mean_cor_AFs$Pop == site,]))
  colnames(tframe)<-c(tframe[1,1], tframe[1,2])
  tframe<-tframe[2:nrow(tframe),]
  test2<-t.test(as.numeric(tframe[,1]), as.numeric(tframe[,2]), paired = T)
  
  row1<-as.data.frame(cbind(paste0(site),paste0("CT"), test1$parameter, test1$statistic, test1$p.value))
  row2<-as.data.frame(cbind(paste0(site),paste0("WH"), test2$parameter, test2$statistic, test2$p.value))
  colnames(row1)<-colnames(row2)<-c("Invasive","Source","df","T-stat","P-value")
  comb<-as.data.frame(rbind(row1,row2))
  comb$Set<-"mean_SST"
  
  mean_list[[loc]]<-comb
}
mean_ttest<-do.call(rbind, mean_list)

mean_list<-list()
for(locus in 1:length(mean_cor_loci)){
  df1<-as.data.frame(mean_cor_AFs[mean_cor_AFs$Pop %in% native, c(1,locus+1)])
  df1<-merge(df1, mean, by.x = "Pop", by.y = "Site")
  regr<-lm(df1[,2] ~ df1[,3])
  df2<-mean_cor_AFs[mean_cor_AFs$Pop %in% invasive, c(1,locus+1)]
  df2<-merge(df2, mean, by.x = "Pop", by.y = "Site")
  df2$Residual<-abs(df2[,2] - (regr$coefficients[2] * df2$annual_mean + regr$coefficients[1]))
  df2$CT_diff<-df2[,2] - df1[2,2]
  df2$WH_diff<-df2[,2] - df1[9,2]
  colnames(df2)[2]<-"AF"
  mean_list[[locus]]<-df2
}
mean_loci<-do.call(rbind, mean_list)
Residual_agg<-cbind(aggregate(mean_loci$Residual ~ mean_loci$Pop, FUN = "median"), "Residual")
colnames(Residual_agg)<-c("Pop","Deviation","Measure")
CT_agg<-cbind(aggregate(mean_loci$CT_diff ~ mean_loci$Pop, FUN = "median"), "CT_diff")
colnames(CT_agg)<-c("Pop","Deviation","Measure")
WH_agg<-cbind(aggregate(mean_loci$WH_diff ~ mean_loci$Pop, FUN = "median"), "WH_diff")
colnames(WH_agg)<-c("Pop","Deviation","Measure")
mean_median<-as.data.frame(rbind(Residual_agg, CT_agg, WH_agg))


medians<-list()
rand_mean_ttest<-list()
for(m in 1:1000){
  random_num<-sample(ncol(allele_freqs)-1, length(mean_cor_loci), replace = F)
  rand_AFs<-allele_freqs[,c(random_num)]
  rand_AFs<-as.data.frame(cbind(Sites, rand_AFs))
  colnames(rand_AFs)[1]<-"Pop"
  rand_mean_ttest2<-list()
  for(b in 1:5){
    site=invasive[b]
    rframe<-as.data.frame(t(rand_AFs[rand_AFs$Pop == "CT" |
                                       rand_AFs$Pop == site,]))
    colnames(rframe)<-c(rframe[1,1], rframe[1,2])
    rframe<-rframe[2:nrow(rframe),]
    test1<-t.test(as.numeric(rframe[,1]), as.numeric(rframe[,2]), paired = T)
    
    rframe<-as.data.frame(t(rand_AFs[rand_AFs$Pop == "WH" |
                                       rand_AFs$Pop == site,]))
    colnames(rframe)<-c(rframe[1,1], rframe[1,2])
    rframe<-rframe[2:nrow(rframe),]
    test2<-t.test(as.numeric(rframe[,1]), as.numeric(rframe[,2]), paired = T)
    
    row1<-as.data.frame(cbind(paste0(site),paste0("CT"), test1$parameter, test1$statistic, test1$p.value))
    row2<-as.data.frame(cbind(paste0(site),paste0("WH"), test2$parameter, test2$statistic, test2$p.value))
    colnames(row1)<-colnames(row2)<-c("Invasive","Source","df","T-stat","P-value")
    comb<-as.data.frame(rbind(row1,row2))
    comb$Set<-"mean_SST"
    
    rand_mean_ttest2[[b]]<-comb
  }
  locus_mean_ttest<-do.call(rbind, rand_mean_ttest2)
  rand_mean_ttest[[m]]<-locus_mean_ttest
  
  locus_list<-list()
  for(locus in 1:length(mean_cor_loci)){
    df1<-as.data.frame(rand_AFs[rand_AFs$Pop %in% native, c(1,locus+1)])
    df1<-merge(df1, mean, by.x = "Pop", by.y = "Site")
    regr<-lm(df1[,2] ~ df1[,3])
    df2<-rand_AFs[rand_AFs$Pop %in% invasive, c(1,locus+1)]
    df2<-merge(df2, mean, by.x = "Pop", by.y = "Site")
    df2$Residual<-abs(df2[,2] - (regr$coefficients[2] * df2$annual_mean + regr$coefficients[1]))
    df2$CT_diff<-df2[,2] - df1[2,2]
    df2$WH_diff<-df2[,2] - df1[9,2]
    colnames(df2)[2]<-"AF"
    locus_list[[locus]]<-df2
  }
  comb_loci<-do.call(rbind, locus_list)
  Residual_agg<-cbind(aggregate(comb_loci$Residual ~ comb_loci$Pop, FUN = "median"), "Residual")
  colnames(Residual_agg)<-c("Pop","Deviation","Measure")
  CT_agg<-cbind(aggregate(comb_loci$CT_diff ~ comb_loci$Pop, FUN = "median"), "CT_diff")
  colnames(CT_agg)<-c("Pop","Deviation","Measure")
  WH_agg<-cbind(aggregate(comb_loci$WH_diff ~ comb_loci$Pop, FUN = "median"), "WH_diff")
  colnames(WH_agg)<-c("Pop","Deviation","Measure")
  median<-as.data.frame(rbind(Residual_agg, CT_agg, WH_agg))
  medians[[m]]<-median
}
mean_distributions<-do.call(rbind, medians)
mean_ttest_distr<-do.call(rbind, rand_mean_ttest)
mean_ttest_distr$`P-value`<-as.numeric(mean_ttest_distr$`P-value`)

mean_ttest_comps1<-list()
for(q in 1:5){
  site=invasive[q]
  popnCT<-mean_ttest_distr[mean_ttest_distr$Invasive == site &
                            mean_ttest_distr$Source == "CT",]
  CT_ttest<-as.numeric(mean_ttest[mean_ttest$Invasive == site &
                                   mean_ttest$Source == "CT",4])
  CT_more<-nrow(popnCT[abs(as.numeric(popnCT$`T-stat`)) >= abs(CT_ttest),])
  
  popnWH<-mean_ttest_distr[mean_ttest_distr$Invasive == site &
                            mean_ttest_distr$Source == "WH",]
  WH_ttest<-as.numeric(mean_ttest[mean_ttest$Invasive == site &
                                   mean_ttest$Source == "WH",4])
  WH_more<-nrow(popnWH[abs(as.numeric(popnWH$`T-stat`)) >= abs(WH_ttest),])
  dat<-as.data.frame(rbind(cbind(site, "CT", CT_ttest, CT_more/1000),
                           cbind(site, "WH", WH_ttest, WH_more/1000)))
  colnames(dat)<-c("Invasive", "Source", "T-test P-value", "Prop. Rand. Higher P-value")
  mean_ttest_comps1[[q]]<-dat
}
mean_ttest_comps<-do.call(rbind, mean_ttest_comps1)
mean_ttest_comps$Set<-"mean_SST"
mean_ttest_distr$`T-stat`<-as.numeric(mean_ttest_distr$`T-stat`)

ggdensity(data = mean_ttest_distr, x = "`T-stat`", y = "..density..",
          facet.by = "Invasive", fill = "Source")

mean_plots<-list()
for(q in 1:5){
  pop<-invasive[q]
  distr<-mean_distributions[mean_distributions$Pop == pop,]
  outliers<-mean_median[mean_median$Pop == pop,]
  mean_plot<-ggdensity(distr, x = "Deviation", y = "..density..", fill = "Measure",
                      title = paste0(pop)) +
    geom_vline(xintercept =  outliers$Deviation[1], col = "green", linewidth = 1.2) + 
    geom_vline(xintercept =  outliers$Deviation[2], col = "red", linewidth = 1.2, lty = 2) +
    geom_vline(xintercept =  outliers$Deviation[3], col = "blue", linewidth = 1.2, lty = 2)
  
  ggsave(filename = file.path(path, "Re_analysis_2023", "PopGen", "AF_preds",
                              "Distribution_plots", "mean", paste0(pop,".png")),
         width = 8, height = 6)
  mean_plots[[q]]<-mean_plot
}
means<-ggarrange(plotlist = mean_plots, common.legend = T)
ggsave(means, filename = file.path(path, "Re_analysis_2023", "PopGen", "AF_preds",
                                  "Distribution_plots", "mean", "All_pops_mean.png"),
       width = 16, height = 10)

#Range
##############################################################################

### rangeimum SST: 
#range_cor<-read.csv(file=file.path(path, "Re_analysis_2023", "PopGen", "AF_preds", "range_corr.csv"))
#range_cor_50<-range_cor[range_cor$correlation >= 0.5,] # N = 54
range_cor_loci<-gsub(":",".",paste0("X",range_comb))
range_cor_AFs<-allele_freqs[colnames(allele_freqs)  %in% range_cor_loci]
range_cor_AFs<-as.data.frame(cbind(allele_freqs$Group.1, range_cor_AFs))
colnames(range_cor_AFs)[1]<-"Pop"
range<-Env[,c(1,5)]

### T-tests between invasive and source locations: ###
range_list<-list()
for(loc in 1:5){
  site=invasive[loc]
  tframe<-as.data.frame(t(range_cor_AFs[range_cor_AFs$Pop == "CT" |
                                        range_cor_AFs$Pop == site,]))
  colnames(tframe)<-c(tframe[1,1], tframe[1,2])
  tframe<-tframe[2:nrow(tframe),]
  test1<-t.test(as.numeric(tframe[,1]), as.numeric(tframe[,2]), paired = T)
  
  tframe<-as.data.frame(t(range_cor_AFs[range_cor_AFs$Pop == "WH" |
                                        range_cor_AFs$Pop == site,]))
  colnames(tframe)<-c(tframe[1,1], tframe[1,2])
  tframe<-tframe[2:nrow(tframe),]
  test2<-t.test(as.numeric(tframe[,1]), as.numeric(tframe[,2]), paired = T)
  
  row1<-as.data.frame(cbind(paste0(site),paste0("CT"), test1$parameter, test1$statistic, test1$p.value))
  row2<-as.data.frame(cbind(paste0(site),paste0("WH"), test2$parameter, test2$statistic, test2$p.value))
  colnames(row1)<-colnames(row2)<-c("Invasive","Source","df","T-stat","P-value")
  comb<-as.data.frame(rbind(row1,row2))
  comb$Set<-"range_SST"
  
  range_list[[loc]]<-comb
}
range_ttest<-do.call(rbind, range_list)

range_list<-list()
for(locus in 1:length(range_cor_loci)){
  df1<-as.data.frame(range_cor_AFs[range_cor_AFs$Pop %in% native, c(1,locus+1)])
  df1<-merge(df1, range, by.x = "Pop", by.y = "Site")
  regr<-lm(df1[,2] ~ df1[,3])
  df2<-range_cor_AFs[range_cor_AFs$Pop %in% invasive, c(1,locus+1)]
  df2<-merge(df2, range, by.x = "Pop", by.y = "Site")
  df2$Residual<-abs(df2[,2] - (regr$coefficients[2] * df2$annual_range + regr$coefficients[1]))
  df2$CT_diff<-df2[,2] - df1[2,2]
  df2$WH_diff<-df2[,2] - df1[9,2]
  colnames(df2)[2]<-"AF"
  range_list[[locus]]<-df2
}
range_loci<-do.call(rbind, range_list)
Residual_agg<-cbind(aggregate(range_loci$Residual ~ range_loci$Pop, FUN = "median"), "Residual")
colnames(Residual_agg)<-c("Pop","Deviation","Measure")
CT_agg<-cbind(aggregate(range_loci$CT_diff ~ range_loci$Pop, FUN = "median"), "CT_diff")
colnames(CT_agg)<-c("Pop","Deviation","Measure")
WH_agg<-cbind(aggregate(range_loci$WH_diff ~ range_loci$Pop, FUN = "median"), "WH_diff")
colnames(WH_agg)<-c("Pop","Deviation","Measure")
range_median<-as.data.frame(rbind(Residual_agg, CT_agg, WH_agg))


medians<-list()
rand_range_ttest<-list()
for(m in 1:1000){
  random_num<-sample(ncol(allele_freqs)-1, length(range_cor_loci), replace = F)
  rand_AFs<-allele_freqs[,c(random_num)]
  rand_AFs<-as.data.frame(cbind(Sites, rand_AFs))
  colnames(rand_AFs)[1]<-"Pop"
  rand_range_ttest2<-list()
  for(b in 1:5){
    site=invasive[b]
    rframe<-as.data.frame(t(rand_AFs[rand_AFs$Pop == "CT" |
                                       rand_AFs$Pop == site,]))
    colnames(rframe)<-c(rframe[1,1], rframe[1,2])
    rframe<-rframe[2:nrow(rframe),]
    test1<-t.test(as.numeric(rframe[,1]), as.numeric(rframe[,2]), paired = T)
    
    rframe<-as.data.frame(t(rand_AFs[rand_AFs$Pop == "WH" |
                                       rand_AFs$Pop == site,]))
    colnames(rframe)<-c(rframe[1,1], rframe[1,2])
    rframe<-rframe[2:nrow(rframe),]
    test2<-t.test(as.numeric(rframe[,1]), as.numeric(rframe[,2]), paired = T)
    
    row1<-as.data.frame(cbind(paste0(site),paste0("CT"), test1$parameter, test1$statistic, test1$p.value))
    row2<-as.data.frame(cbind(paste0(site),paste0("WH"), test2$parameter, test2$statistic, test2$p.value))
    colnames(row1)<-colnames(row2)<-c("Invasive","Source","df","T-stat","P-value")
    comb<-as.data.frame(rbind(row1,row2))
    comb$Set<-"range_SST"
    
    rand_range_ttest2[[b]]<-comb
  }
  locus_range_ttest<-do.call(rbind, rand_range_ttest2)
  rand_range_ttest[[m]]<-locus_range_ttest
  
  locus_list<-list()
  for(locus in 1:length(range_cor_loci)){
    df1<-as.data.frame(rand_AFs[rand_AFs$Pop %in% native, c(1,locus+1)])
    df1<-merge(df1, range, by.x = "Pop", by.y = "Site")
    regr<-lm(df1[,2] ~ df1[,3])
    df2<-rand_AFs[rand_AFs$Pop %in% invasive, c(1,locus+1)]
    df2<-merge(df2, range, by.x = "Pop", by.y = "Site")
    df2$Residual<-abs(df2[,2] - (regr$coefficients[2] * df2$annual_range + regr$coefficients[1]))
    df2$CT_diff<-df2[,2] - df1[2,2]
    df2$WH_diff<-df2[,2] - df1[9,2]
    colnames(df2)[2]<-"AF"
    locus_list[[locus]]<-df2
  }
  comb_loci<-do.call(rbind, locus_list)
  Residual_agg<-cbind(aggregate(comb_loci$Residual ~ comb_loci$Pop, FUN = "median"), "Residual")
  colnames(Residual_agg)<-c("Pop","Deviation","Measure")
  CT_agg<-cbind(aggregate(comb_loci$CT_diff ~ comb_loci$Pop, FUN = "median"), "CT_diff")
  colnames(CT_agg)<-c("Pop","Deviation","Measure")
  WH_agg<-cbind(aggregate(comb_loci$WH_diff ~ comb_loci$Pop, FUN = "median"), "WH_diff")
  colnames(WH_agg)<-c("Pop","Deviation","Measure")
  median<-as.data.frame(rbind(Residual_agg, CT_agg, WH_agg))
  medians[[m]]<-median
}
range_distributions<-do.call(rbind, medians)
range_ttest_distr<-do.call(rbind, rand_range_ttest)
range_ttest_distr$`P-value`<-as.numeric(range_ttest_distr$`P-value`)

range_ttest_comps1<-list()
for(q in 1:5){
  site=invasive[q]
  popnCT<-range_ttest_distr[range_ttest_distr$Invasive == site &
                            range_ttest_distr$Source == "CT",]
  CT_ttest<-as.numeric(range_ttest[range_ttest$Invasive == site &
                                   range_ttest$Source == "CT",4])
  CT_more<-nrow(popnCT[abs(as.numeric(popnCT$`T-stat`)) >= abs(CT_ttest),])
  
  popnWH<-range_ttest_distr[range_ttest_distr$Invasive == site &
                            range_ttest_distr$Source == "WH",]
  WH_ttest<-as.numeric(range_ttest[range_ttest$Invasive == site &
                                   range_ttest$Source == "WH",4])
  WH_more<-nrow(popnWH[abs(as.numeric(popnWH$`T-stat`)) >= abs(WH_ttest),])
  dat<-as.data.frame(rbind(cbind(site, "CT", CT_ttest, CT_more/1000),
                           cbind(site, "WH", WH_ttest, WH_more/1000)))
  colnames(dat)<-c("Invasive", "Source", "T-test P-value", "Prop. Rand. Higher P-value")
  range_ttest_comps1[[q]]<-dat
}
range_ttest_comps<-do.call(rbind, range_ttest_comps1)
range_ttest_comps$Set<-"range_SST"
range_ttest_distr$`T-stat`<-as.numeric(range_ttest_distr$`T-stat`)

ggdensity(data = range_ttest_distr, x = "`T-stat`", y = "..density..",
          facet.by = "Invasive", fill = "Source")

range_plots<-list()
for(q in 1:5){
  pop<-invasive[q]
  distr<-range_distributions[range_distributions$Pop == pop,]
  outliers<-range_median[range_median$Pop == pop,]
  range_plot<-ggdensity(distr, x = "Deviation", y = "..density..", fill = "Measure",
                      title = paste0(pop)) +
    geom_vline(xintercept =  outliers$Deviation[1], col = "green", linewidth = 1.2) + 
    geom_vline(xintercept =  outliers$Deviation[2], col = "red", linewidth = 1.2, lty = 2) +
    geom_vline(xintercept =  outliers$Deviation[3], col = "blue", linewidth = 1.2, lty = 2)
  
  ggsave(filename = file.path(path, "Re_analysis_2023", "PopGen", "AF_preds",
                              "Distribution_plots", "range", paste0(pop,".png")),
         width = 8, height = 6)
  range_plots[[q]]<-range_plot
}
ranges<-ggarrange(plotlist = range_plots, common.legend = T)
ggsave(ranges, filename = file.path(path, "Re_analysis_2023", "PopGen", "AF_preds",
                                  "Distribution_plots", "range", "All_pops_range.png"),
       width = 16, height = 10)


########
comb_ttest<-as.data.frame(rbind(max_ttest, mean_ttest, range_ttest))
comb_ttest$`T-stat`<-as.numeric(comb_ttest$`T-stat`)
comb_ttest$`P-value`<-as.numeric(comb_ttest$`P-value`)

ttests<-as.data.frame(rbind(max_ttest_distr, mean_ttest_distr,
                            range_ttest_distr))
ttests$`T-stat`<-as.numeric(ttests$`T-stat`)
ttests$`P-value`<-as.numeric(ttests$`P-value`)

pviolins<-ggviolin(ttests, x = "Invasive", y = "`T-stat`", facet.by = "Set", fill = "Source") +
  geom_point(data = comb_ttest, shape = factor(comb_ttest$Source), size = 3) +
  geom_hline(yintercept = 0.0, lty = 2)
ggsave("Re_analysis_2023/Final_plots/2024.08.20_T-stat_deviations_from_source.svg",pviolins,
       height = 5, width = 12)

ggscatter(comb_ttest, x = "Invasive", y = "P-value", facet.by = "Set", shape = "Source") +
  geom_hline(yintercept = 0.05, lty = 2)

########
resid_dist<-list()
for(q in 1:5){
  site=invasive[q]
  max_dist<-max_distributions[max_distributions$Pop == site &
                                max_distributions$Measure == "Residual",]
  max_dist$Set<-"Max_SST"
  max_regr<-max_median[max_median$Pop == site &
                         max_median$Measure == "Residual",2]
  max_rand<-max_dist[max_dist$Deviation <= max_regr,]
  max_num<-nrow(max_dist[max_dist$Deviation <= max_regr,])
  
  #min_dist<-min_distributions[min_distributions$Pop == site &
  #                              min_distributions$Measure == "Residual",]
  #min_dist$Set<-"Min_SST"
  #min_regr<-min_median[min_median$Pop == site &
  #                       min_median$Measure == "Residual",2]
  #min_num<-nrow(min_dist[min_dist$Deviation <= min_regr,])
  
  mean_dist<-mean_distributions[mean_distributions$Pop == site &
                                mean_distributions$Measure == "Residual",]
  mean_dist$Set<-"Mean_SST"
  mean_regr<-mean_median[mean_median$Pop == site &
                         mean_median$Measure == "Residual",2]
  mean_num<-nrow(mean_dist[mean_dist$Deviation <= mean_regr,])
  
  range_dist<-range_distributions[range_distributions$Pop == site &
                                range_distributions$Measure == "Residual",]
  range_dist$Set<-"Range_SST"
  range_regr<-range_median[range_median$Pop == site &
                         range_median$Measure == "Residual",2]
  range_num<-nrow(range_dist[range_dist$Deviation <= range_regr,])
  
  df<-rbind(max_dist, mean_dist, range_dist)
  resid_dist[[q]]<-df
  
}
resid_dist<-do.call(rbind, resid_dist)

max_median$Set<-"Max_SST"
#min_median$Set<-"Min_SST"
mean_median$Set<-"Mean_SST"
range_median$Set<-"Range_SST"

resid<-rbind(max_median, mean_median, range_median)
resid<-resid[resid$Measure == "Residual",]

resid_dens<-ggdensity(resid_dist, x = "Deviation", y = "..density..", facet.by = "Set",
          fill = "Pop") + geom_vline(data = resid, xintercept = resid$Deviation)
ggsave("Re_analysis_2023/Final_plots/2024.08.20_residual_density_plots.svg", resid_dens,
       height = 6, width = 10)

resid_box<-ggboxplot(data = resid_dist, x = "Pop", y = "Deviation", fill = "Pop", facet.by = "Set") +
  geom_point(data = resid, shape = 22, size = 2)
ggsave("Re_analysis_2023/Final_plots/2024.08.20_residual_box_plots.svg", resid_box,
       height = 5, width = 12)
