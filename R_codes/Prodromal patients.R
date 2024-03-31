library(stringr) ; library(ggplot2) ; library(edgeR) ; library(matrixStats) ; library(MatchIt) ; library(reshape) ; library (caret) ; library(statmod)
library(factoextra) ; library(ggfortify) ; library(precrec) ; library(ggrepel) ; library(pROC) ; library(yardstick) ; library(MLeval) 
dir1<-function(x) paste0('Y:/nimrod.madrer/ADPD/PPMI/tRFs/',x) ; dir1('')

## Read files:
# Final diagnosis of part of the Prodromal patients:
dig<-read.csv(dir1('Prodromal_status.csv'),fileEncoding="latin1")

# Read and process tRFs metadata:
ppmiMeta<-read.csv(dir1('metadat_tRF_updated_030222.csv'))
ppmiMeta$len<-unlist(lapply(ppmiMeta$MINTbase.Unique.ID, function(x) strsplit(as.character(x),'-')[[1]][2]))
ppmiMeta$trna<-unlist(lapply(ppmiMeta$Sequence.locations.in.tRNA.space..comma.deliminated., function(x) substr(strsplit(as.character(x),'_')[[1]][2],1,3)))
colnames(ppmiMeta)<-gsub('Sequence.locations.in.tRNA.space..comma.deliminated.','details',colnames(ppmiMeta))

# Define Nuclear and Mitochondrial tRFs:
ppmiMeta$nuclear<-F
for(i in 1:nrow(ppmiMeta)){
  ppmiMeta$nuclear[i]<-sum(str_count(as.character(ppmiMeta$details[i]),as.character('MT')))}
ppmiMeta$nuclear<-factor(ppmiMeta$nuclear,labels=c('Nuclear','Mitochondrial'))

# Find PD-tRFs:
ppmiMeta$motifs<-F
for(i in 1:nrow(ppmiMeta)){
  ppmiMeta$motifs[i]<-sum(str_count(as.character(ppmiMeta$tRF.sequence[i]),as.character('[AG]GTTC[GA]A')))}
ppmiMeta$motifs<-factor(ppmiMeta$motifs,labels=c('No RGTTCRA','RGTTCRA motif'))

# Find qPCR-matching PD-tRFs:
mts<-c('TAACTTAGCATTAACCTTTTAA') ; mtfs<-c('GGTCCCTGGTTCAA') ; ppmiMeta$mtfs3<-0
for(i in 1:nrow(ppmiMeta)){
  for(a in mtfs){if(ppmiMeta$mtfs3[i]==0){ppmiMeta$mtfs3[i]<-sum(str_count(as.character(ppmiMeta$tRF.sequence[i]),as.character(a)))}}
  for(b in mts){if(ppmiMeta$mtfs3[i]==0){ppmiMeta$mtfs3[i]<-sum(str_count(as.character(ppmiMeta$tRF.sequence[i]),as.character(b)))*2}}
}  ; ppmiMeta$mtfs3<-factor(ppmiMeta$mtfs3,labels=c('No motifs','RGTTCRA','MT'))  ; table(ppmiMeta$mtfs3)

# Read tRF counts table:
ctsPPMI<-read.csv(dir1('tRNA_Exclusive_Combined_data.csv')) ; cts1<-ctsPPMI ; rownames(cts1)<-cts1$X ; cts1<-cts1[,-1]
colnames(cts1)<-unlist(lapply(colnames(cts1),function(x) paste0(strsplit(strsplit(x,'IR1.')[[1]][2],'\\.')[[1]][1:2],collapse='_')))

# Read and process Metadata:
col_ppmi<-read.csv(dir1('PPMI_all_cleaned_orgenized_data.csv')) ; col_ppmi<-unique(col_ppmi[,c(1:27,196,199)]) ; xxy<-c("3054","3209","3767","3780","4079_BL")
col_ppmi<-subset(col_ppmi,col_ppmi$f_shortRNA==1 & col_ppmi$f_RSNEXC_term==1 & col_ppmi$f_QC_longRNA=='HighQC')
col_ppmi<-subset(col_ppmi,!col_ppmi$Subject %in% xxy & is.na(col_ppmi$Type_of_surgery_for_Parkinson_disease_term))
col_ppmi<-subset(col_ppmi,col_ppmi$Ethnicity!='Unknown' & col_ppmi$Genetic_background!='Unknown')
srt_cols<-c('s_mus_inflammation','s_psych_anxiety','s_urinary_bladder','s_diabetes')
for(c in srt_cols){
  for(x in(1:nrow(col_ppmi))){
    n<-col_ppmi[col_ppmi$timePoint=='BL' & col_ppmi$Subject==col_ppmi$Subject[x],c]
    if(length(n)>0){col_ppmi[x,c]<-col_ppmi[col_ppmi$timePoint=='BL' & col_ppmi$Subject==col_ppmi$Subject[x],c]}
    if(length(n)==0){col_ppmi[x,c]<-NA}}
  print(paste(length(unique(col_ppmi[is.na(col_ppmi[,c]),'Subject']))==
                length(unique(col_ppmi$Subject))-length(unique(subset(col_ppmi$Subject,col_ppmi$timePoint=='BL'))),c))
  col_ppmi[,c]<-gsub(1,strsplit(as.character(c),'_')[[1]][length(strsplit(as.character(c),'_')[[1]])],col_ppmi[,c])
}
col_ppmi$Other_diseases<-unlist(lapply(1:nrow(col_ppmi),function(x) 
  gsub('_0','',gsub('0_','',paste(unique(as.character(col_ppmi[x,srt_cols])),collapse='_')))))
col_ppmi<-col_ppmi[,!colnames(col_ppmi) %in% c(srt_cols,"f_shortRNA","f_RSNEXC_term","f_QC_longRNA")]
col_ppmi$cond<-factor(col_ppmi$Group,
                      levels=c('Genetic Cohort - Unaffected','Genetic Registry - Unaffected','Healthy Control','Genetic Cohort - PD','Genetic Registry - PD',"Parkinson's disease",'Prodromal','SWEDD'),
                      labels=c(rep('Ctrl',3),rep('PD',3),'Prodromal','SWEDD')) ; 

# Create Prodromal patients table:
f_ppmi<-subset(col_ppmi,col_ppmi$Group %in% c('Prodromal','Healthy Control'))
f_ppmi<-subset(f_ppmi,f_ppmi$PD_meds %in% c('non','UnKnown') & f_ppmi$Age_timepoint>=58 & f_ppmi$Other_diseases %in% c('0','NA') & 
                 f_ppmi$timePoint %in% c('BL','V02','V04')) ; f_ppmi$cond<-droplevels(f_ppmi$cond)

biomarker<-data.frame(sample=colnames(cts1),mt=colSums(cts1[rownames(cts1) %in% subset(ppmiMeta$MINTbase.Unique.ID,ppmiMeta$mtfs3=='MT'),]),
                      mtf=colSums(cts1[rownames(cts1) %in% subset(ppmiMeta$MINTbase.Unique.ID,ppmiMeta$mtfs3=='RGTTCRA'),]))
biomarker$score<-biomarker$mtf/(biomarker$mt+1) ; biomarker<-merge(biomarker,col_ppmi,by='sample')
# write.csv(biomarker,dir1('PPMI_Mtf.MT_scores.csv'))

ggplot(biomarker,aes(cond,(score),col=cond))+theme_classic()+geom_boxplot()+facet_wrap(~score<=15,nrow=2,scales='free')


#
## Matching #####
toKeep<-c()
for(s in unique(f_ppmi$Subject)){
  tmp<-subset(f_ppmi,f_ppmi$Subject==s) ; toKeep<-append(toKeep,as.character(tmp$sample[order(tmp$timePoint,decreasing=F)][1]))}
tmp1<-subset(f_ppmi,f_ppmi$sample %in% toKeep)

# Optimal:
out_opt<-matchit(cond ~ Sex+Age_timepoint+s_Study,data=tmp1,method='optimal',link='probit',distance='glm')
plot(out_opt, type = "density", interactive = FALSE,which.xs = ~Age_timepoint + Sex + s_Study)
m.data<-match.data(out_opt) ; summary(out_opt,un=T) ; plot(summary(out_opt))
table(m.data$subclass,m.data$cond) ; plot(out_opt, type = "jitter", interactive = FALSE)

# Full:
out_opt2<-matchit(cond ~ Sex+Age_timepoint+s_Study,data=tmp1,method='full',link='probit',distance='glm')
plot(out_opt2, type = "density", interactive = FALSE,which.xs = ~Age_timepoint + Sex + s_Study)
m.data2<-match.data(out_opt2) ; summary(out_opt2,un=T) ; plot(summary(out_opt2))
table(m.data2$subclass,m.data2$cond) ; plot(out_opt2, type = "jitter", interactive = FALSE)
# write.csv(m.data2,dir1('ProdromalData_Supp.csv'))

# Create counts tables for the matched data:
trfs0<-ctsPPMI ; rownames(trfs0)<-trfs0$X ; trfs0<-trfs0[,-1]
colnames(trfs0)<-unlist(lapply(colnames(trfs0),
                               function(x) paste0(strsplit(strsplit(x,'IR1.')[[1]][2],'\\.')[[1]][1:2],collapse='_')))
trfs2<-trfs0[,colnames(trfs0) %in% m.data2$sample] ; trfs<-trfs0[,colnames(trfs0) %in% m.data$sample]


## Calculating PD/MT-tRF and clinical ratios for the matched data:
# Fully matched:
f_motif<-data.frame(sample=colnames(trfs2),sum=colSums(trfs2),
                    mtf=colSums(trfs2[rownames(trfs2) %in% subset(ppmiMeta$MINTbase.Unique.ID,ppmiMeta$motifs=='RGTTCRA motif'),]),
                    mt=colSums(trfs2[rownames(trfs2) %in% subset(ppmiMeta$MINTbase.Unique.ID,ppmiMeta$nuclear=='Mitochondrial'),]))
f_motif$ratio<-100*f_motif$mtf/f_motif$sum ;  f_motif$mt_pct<-100*f_motif$mt/f_motif$sum ; f_motif$mtf_mt<-f_motif$mtf/(f_motif$mt+1)
f_motif<-merge(f_motif,m.data2,by='sample')

# Normalizing motifs ratio according to matched samples:
f_motif$mtf_mt2<-unlist(lapply(1:nrow(f_motif),function(n)
  (f_motif$mtf_mt[n]/mean(subset(f_motif$mtf_mt,f_motif$subclass==f_motif$subclass[n])))))#/sd(subset(f_motif$mtf_mt,f_motif$subclass==f_motif$subclass[n]))
f_motif$Genetic_background<-factor(f_motif$Genetic_background,levels=c('LRRK2-/SNCA-/GBA-','GBA+'))
f_motif$nrm_ratio<-unlist(lapply(1:nrow(f_motif),function(n)
  (f_motif$ratio[n]/mean(subset(f_motif$ratio,f_motif$subclass==f_motif$subclass[n])))))#/sd(subset(f_motif$ratio,f_motif$subclass==f_motif$subclass[n]))
tmp<-data.frame(mt=colSums(trfs2[rownames(trfs2) %in% subset(ppmiMeta$MINTbase.Unique.ID,ppmiMeta$nuclear=='Mitochondrial'),]),
                sample=colnames(trfs2),all=colSums(trfs2)) ; tmp$mt_pct<-100*tmp$mt/tmp$all ; f_motif<-merge(f_motif,tmp[,c("sample","mt_pct")])
f_motif$mt_pct<-unlist(lapply(1:nrow(f_motif),function(n)
  (f_motif$mt_pct[n]/mean(subset(f_motif$mt_pct,f_motif$subclass==f_motif$subclass[n])))))#/sd(subset(f_motif$mt_pct,f_motif$subclass==f_motif$subclass[n]))

# Optimally matched:
motif<-data.frame(sample=colnames(trfs),sum=colSums(trfs),
                  mtf=colSums(trfs[rownames(trfs) %in% subset(ppmiMeta$MINTbase.Unique.ID,ppmiMeta$motifs=='RGTTCRA motif'),]),
                  mt=colSums(trfs[rownames(trfs) %in% subset(ppmiMeta$MINTbase.Unique.ID,ppmiMeta$nuclear=='Mitochondrial'),]))
motif$ratio<-100*motif$mtf/motif$sum ; motif$mt_pct<-100*motif$mt/motif$sum ; motif$mtf_mt<-motif$mtf/(motif$mt+1) ; motif<-merge(motif,m.data,by='sample')
motif$Genetic_background<-factor(motif$Genetic_background,levels=c('LRRK2-/SNCA-/GBA-','GBA+'))

# Normalizing motifs ratio according to matched samples:
motif$nrm_ratio<-unlist(lapply(1:nrow(motif),function(n)
  (motif$ratio[n]/mean(subset(motif$ratio,motif$subclass==motif$subclass[n])))))#/sd(subset(motif$ratio,motif$subclass==motif$subclass[n]))
motif$mtf_mt2<-unlist(lapply(1:nrow(motif),function(n)
  (motif$mtf_mt[n]/mean(subset(motif$mtf_mt,motif$subclass==motif$subclass[n])))))#/sd(subset(motif$mtf_mt,motif$subclass==motif$subclass[n]))
motif$mt_pct<-unlist(lapply(1:nrow(motif),function(n)
  (motif$mt_pct[n]/mean(subset(motif$mt_pct,motif$subclass==motif$subclass[n])))))#/sd(subset(motif$mt_pct,motif$subclass==motif$subclass[n]))

# Remove patients without matches:
toRmv<-c() ; for (s in unique(motif$subclass)){tmp<-subset(motif,motif$subclass==s) 
if(sum(tmp$cond=='Ctrl')==0 | sum(tmp$cond=='Prodromal')==0){toRmv<-append(toRmv,s)}} ; motif<-subset(motif,!motif$subclass %in% toRmv)


###
## ROC curves optimal matching #####

# Nrm-PD/MT-tRFs
roc_all0<-motif[,c("cond","mtf_mt2","Ethnicity")]
ctrl0 <- trainControl(method="cv", summaryFunction=twoClassSummary, number=5, classProbs=T,savePredictions = T)
o_fit_all0 <- train(cond ~ .,data=roc_all0,method="gbm",trControl=ctrl0,na.action=na.omit)

# Clinical score
roc_all1<-motif[,c("cond","Ethnicity","Hoehn.and.Yahr.staging","UPDRS.score.III")]
ctrl1 <- trainControl(method="cv", summaryFunction=twoClassSummary, number=5, classProbs=T,savePredictions = T)
o_fit_all1 <- train(cond ~ .,data=roc_all1,method="gbm",trControl=ctrl1,na.action=na.omit)

# Mixed labels
fk_mtf<-motif ; fk_mtf$cond<-'Ctrl' ; fk_mtf$cond[sample(1:nrow(fk_mtf),60,replace=F)]<-'Prodromal'
smps1<-c(sample(subset(fk_mtf$sample,fk_mtf$cond=='Prodromal'),60,replace=F),sample(subset(fk_mtf$sample,fk_mtf$cond=='Ctrl'),60,replace=F))
tmp1<-data.frame(sample=smps1) ; tmp1<-merge(tmp1,fk_mtf,by='sample')

roc_all2<-tmp1[,c("cond","mtf_mt2","Ethnicity","Hoehn.and.Yahr.staging","UPDRS.score.III")]
ctrl2 <- trainControl(method="cv", summaryFunction=twoClassSummary, number=5, classProbs=T,savePredictions = T)
o_fit_all2 <- train(cond ~ .,data=roc_all2,method="gbm",trControl=ctrl2,na.action=na.omit)

# ROC
set.seed(123) ; res<-evalm(list1=list(o_fit_all2,o_fit_all1,o_fit_all0),gnames=c('Mixed lables','Only clinic','Only tRFs'))
# ggsave(dir1('ROC_Prodromal_alltRFs.svg'))

#
# ROC curves full matching

a_all_true<-c() ; a_all_fake<-c() ; a_clnc_true<-c() ; a_tRF_true<-c()
for(i in 1:10000){
  smps<-c(sample(subset(f_motif$sample,f_motif$cond=='Prodromal'),60,replace=F),sample(subset(f_motif$sample,f_motif$cond=='Ctrl'),60,replace=F))
  tmp<-data.frame(sample=smps) ; tmp<-merge(tmp,f_motif,by='sample')
  
  # Nrm-PD/MT-tRFs + clinical scores:
  roc_all<-tmp[,c("cond","mtf_mt2","Ethnicity","Genetic_background","Hoehn.and.Yahr.staging","UPDRS.score.III")]
  ctrl <- trainControl(method="cv", summaryFunction=twoClassSummary, number=5, classProbs=T,savePredictions = T)
  fit_all <- train(cond ~ .,data=roc_all,method="gbm",trControl=ctrl,na.action=na.omit)
  res <- evalm(fit_all) ; a_all_true<-append(a_all_true,as.numeric(res$stdres$`Group 1`['AUC-ROC','Score']))
  
  # Nrm-PD/MT-tRFs
  roc_all0<-tmp[,c("cond","mtf_mt2","Ethnicity","Genetic_background")]
  ctrl0 <- trainControl(method="cv", summaryFunction=twoClassSummary, number=5, classProbs=T,savePredictions = T)
  fit_all0 <- train(cond ~ .,data=roc_all0,method="gbm",trControl=ctrl0,na.action=na.omit)
  res0 <- evalm(fit_all0) ; a_tRF_true<-append(a_tRF_true,as.numeric(res0$stdres$`Group 1`['AUC-ROC','Score']))
  
  # Clinical scores
  roc_all1<-tmp[,c("cond","Ethnicity","Genetic_background","Hoehn.and.Yahr.staging","UPDRS.score.III")]
  ctrl1 <- trainControl(method="cv", summaryFunction=twoClassSummary, number=5, classProbs=T,savePredictions = T)
  fit_all1 <- train(cond ~ .,data=roc_all1,method="gbm",trControl=ctrl1,na.action=na.omit)
  res1 <- evalm(fit_all1) ; a_clnc_true<-append(a_clnc_true,as.numeric(res1$stdres$`Group 1`['AUC-ROC','Score']))
  
  # Mixed labels
  fk_mtf<-f_motif ; fk_mtf$cond<-'Ctrl' ; fk_mtf$cond[sample(1:nrow(fk_mtf),60,replace=F)]<-'Prodromal'
  smps1<-c(sample(subset(fk_mtf$sample,fk_mtf$cond=='Prodromal'),60,replace=F),sample(subset(fk_mtf$sample,fk_mtf$cond=='Ctrl'),60,replace=F))
  tmp1<-data.frame(sample=smps1) ; tmp1<-merge(tmp1,fk_mtf,by='sample')
  
  roc_all2<-tmp1[,c("cond","mtf_mt2","Ethnicity","Genetic_background","Hoehn.and.Yahr.staging","UPDRS.score.III")]
  ctrl2 <- trainControl(method="cv", summaryFunction=twoClassSummary, number=5, classProbs=T,savePredictions = T)
  fit_all2 <- train(cond ~ .,data=roc_all2,method="gbm",trControl=ctrl2,na.action=na.omit)
  res2 <- evalm(fit_all2) ; a_all_fake<-append(a_all_fake,as.numeric(res2$stdres$`Group 1`['AUC-ROC','Score']))
  
  print(paste0('                                                                                ',i))
  print(paste0('                                                                                ',i))
  print(paste0('                                                                                ',i))
}

pv<-data.frame('all_true'=a_all_true,'all_fake'=a_all_fake,'clnc_true'=a_clnc_true,'tRF_true'=a_tRF_true) 
pv$id<-rownames(pv) ; pv<-melt(pv,id.vars='id') ; tmp<-subset(pv,pv$variable!='tRF_true')
pv$variable<-factor(pv$variable,levels=c('all_true','tRF_true','clnc_true','all_fake'),
                        labels=c('tRFs + UPDRS III + H&Y','Only tRFs','UPDRS III + H&Y','Mixed labels'))
ggplot(pv,aes(variable,value,col=variable))+theme_classic()+geom_boxplot() ; TukeyHSD(aov(pv$value~pv$variable))
# ggsave(dir1('Boxplot_prodromal_alltRFs.svg'))

#

## Features selection #####

# Finding the Motif Seq
ppmiMeta$mtf_loc_s<-NA ; ppmiMeta$mtf_loc_e<-NA
for(i in 1:nrow(ppmiMeta)){ppmiMeta[i,c('mtf_loc_s','mtf_loc_e')]<-as.numeric(str_locate(ppmiMeta$tRF.sequence[i],'[AG]GTTC[AG]A'))}
ppmiMeta[is.na(ppmiMeta$mtf_loc_s),'mtf_loc_s']<-0 ; ppmiMeta[is.na(ppmiMeta$mtf_loc_e),'mtf_loc_e']<-0
ppmiMeta$ptrn<-unlist(lapply(1:nrow(ppmiMeta),function(x) substr(ppmiMeta$tRF.sequence[x],ppmiMeta$mtf_loc_s[x]-7,ppmiMeta$mtf_loc_e[x]+0)))
ppmiMeta[ppmiMeta$mtf_loc_e==0,'ptrn']<-'No motif' ; ppmiMeta[ppmiMeta$mtf_loc_s<7,'ptrn']<-'No motif'
View(table(ppmiMeta$ptrn))

chp_tRFs0<-trfs0[rownames(trfs0) %in% subset(ppmiMeta$MINTbase.Unique.ID,ppmiMeta$ptrn!='No motif'),] ; chp_trfs<-chp_tRFs0[,colnames(chp_tRFs0) %in% m.data$sample]
tmp1<-as.data.frame(table(ppmiMeta$ptrn)) ; tmp1<-subset(tmp1,tmp1$Freq<1000 & tmp1$Freq>100) ; pvTab0<-as.data.frame(t(data.frame(row.names=c('mtf','dif','p'))))
for(f in as.character(tmp1$Var1)){
  chip<-data.frame(sample=colnames(chp_trfs),mtf=colSums(chp_trfs[rownames(chp_trfs) %in% subset(ppmiMeta$MINTbase.Unique.ID,ppmiMeta$ptrn==f),]))
  chip<-merge(chip,m.data,by='sample') ; tt<-wilcox.test(chip$mtf~chip$cond)
  md<-as.numeric(median(chip[chip$cond=='Prodromal','mtf'])-median(chip[chip$cond=='Ctrl','mtf']))
  pvTab0<-rbind(pvTab0,data.frame('mtf'=f,'dif'=md,'p'=tt$p.value))}
paste(subset(pvTab0$mtf,pvTab0$dif>0 & pvTab0$p<0.1),collapse="','")

# Mtf
hmap_mtf<-trfs2[rownames(trfs2) %in% subset(ppmiMeta$MINTbase.Unique.ID,ppmiMeta$motifs=='RGTTCRA motif'),] ; hmap_mtf<-hmap_mtf[rowSums(hmap_mtf)>1,]
hmap_mtf<-hmap_mtf[unlist(lapply(1:nrow(hmap_mtf), function(x) sum(hmap_mtf[x,]!=0)>15)),]
for(s in unique(m.data2$subclass)){
  sps<-subset(m.data2$sample,m.data2$subclass==s)
  hmap_mtf[,colnames(hmap_mtf) %in% sps]<-hmap_mtf[,colnames(hmap_mtf) %in% sps]/rowMeans(hmap_mtf[,colnames(hmap_mtf) %in% sps])
} ; hmap_mtf[is.na(hmap_mtf)]<-0 ; hmap_mtf$trf<-rownames(hmap_mtf) 

hmap_mtf$sdCt<-rowSds(as.matrix(hmap_mtf[,colnames(hmap_mtf) %in% subset(m.data2$sample,m.data2$cond=='Ctrl')]))
hmap_mtf$sdPr<-rowSds(as.matrix(hmap_mtf[,colnames(hmap_mtf) %in% subset(m.data2$sample,m.data2$cond=='Prodromal')]))
hmap_mtf$dif<-rowMeans(hmap_mtf[,colnames(hmap_mtf) %in% subset(m.data2$sample,m.data2$cond=='Prodromal')])/
  rowMeans(hmap_mtf[,colnames(hmap_mtf) %in% subset(m.data2$sample,m.data2$cond=='Ctrl')])

hmap_mtf<-melt(hmap_mtf,id.vars=c('trf','dif','sdCt','sdPr'),variable_name='sample')
hmap_mtf<-merge(hmap_mtf,f_motif,by='sample') ; hmap_mtf<-merge(hmap_mtf,ppmiMeta,by.x='trf',by.y='MINTbase.Unique.ID')
hmap_mtf$var<-(hmap_mtf$dif)/(hmap_mtf$sdCt+hmap_mtf$sdPr)
tmp<-unique(hmap_mtf[,c("trf","var")]) ; hmap_mtf$trf<-factor(hmap_mtf$trf,levels=tmp$trf[order(tmp$var)])
hmap_mtf$sample<-factor(hmap_mtf$sample,levels=colnames(trfs2)[order(colSums(trfs2))])

quantile(hmap_mtf$var,probs=0.99)

tmp<-subset(hmap_mtf,hmap_mtf$var>1.025)
ggplot(tmp,aes(sample,trf,fill=(value)))+theme_classic()+facet_wrap(~cond,scales='free')+
  geom_tile()+scale_fill_gradient2(low='red',high='blue')+theme(axis.text=element_blank())

ggplot(tmp,aes(sample,value))+theme_classic()+facet_wrap(~cond,scales='free_x')+
  stat_summary(fun='sum',geom='point')+theme(axis.text=element_blank())

ggplot(tmp,aes(cond,value))+theme_classic()+stat_summary(fun='mean',geom='point')

mtf_Sep<-unique(hmap_mtf[hmap_mtf$var>1.025,colnames(hmap_mtf) %in% c('trf',colnames(ppmiMeta))])


# MT
hmap_mt<-trfs2[rownames(trfs2) %in% subset(ppmiMeta$MINTbase.Unique.ID,ppmiMeta$nuclear=='Mitochondrial'),] ; hmap_mt<-hmap_mt[rowSums(hmap_mt)>1,]
hmap_mt<-hmap_mt[unlist(lapply(1:nrow(hmap_mt), function(x) sum(hmap_mt[x,]!=0)>15)),]
for(s in unique(m.data2$subclass)){
  sps<-subset(m.data2$sample,m.data2$subclass==s)
  hmap_mt[,colnames(hmap_mt) %in% sps]<-hmap_mt[,colnames(hmap_mt) %in% sps]/rowMeans(hmap_mt[,colnames(hmap_mt) %in% sps])
} ; hmap_mt[is.na(hmap_mt)]<-0 ; hmap_mt$trf<-rownames(hmap_mt) 

hmap_mt$sdCt<-rowSds(as.matrix(hmap_mt[,colnames(hmap_mt) %in% subset(m.data2$sample,m.data2$cond=='Ctrl')]))
hmap_mt$sdPr<-rowSds(as.matrix(hmap_mt[,colnames(hmap_mt) %in% subset(m.data2$sample,m.data2$cond=='Prodromal')]))
hmap_mt$dif<-rowMeans(hmap_mt[,colnames(hmap_mt) %in% subset(m.data2$sample,m.data2$cond=='Prodromal')])/
  rowMeans(hmap_mt[,colnames(hmap_mt) %in% subset(m.data2$sample,m.data2$cond=='Ctrl')])
hmap_mt$dif2<-rowMeans(hmap_mt[,colnames(hmap_mt) %in% subset(m.data2$sample,m.data2$cond=='Ctrl')])/
  rowMeans(hmap_mt[,colnames(hmap_mt) %in% subset(m.data2$sample,m.data2$cond=='Prodromal')])

hmap_mt<-melt(hmap_mt,id.vars=c('trf','dif','dif2','sdCt','sdPr'),variable_name='sample')
hmap_mt<-merge(hmap_mt,f_motif,by='sample') ; hmap_mt<-merge(hmap_mt,ppmiMeta,by.x='trf',by.y='MINTbase.Unique.ID')
hmap_mt$sample<-factor(hmap_mt$sample,levels=colnames(trfs2)[order(colSums(trfs2))])
hmap_mt$var<-(hmap_mt$dif)/(hmap_mt$sdCt+hmap_mt$sdPr) ; hmap_mt$var2<-(hmap_mt$dif2)/(hmap_mt$sdCt+hmap_mt$sdPr)
tmp<-unique(hmap_mt[,c("trf","var","var2")]) ; hmap_mt$trf<-factor(hmap_mt$trf,levels=tmp$trf[order(tmp$var2)])

quantile(hmap_mt$var2,probs=0.985)

tmp<-subset(hmap_mt,hmap_mt$var2>2 & hmap_mt$chromosome=='MT')
ggplot(tmp,aes(sample,trf,fill=(value)))+theme_classic()+facet_wrap(~cond,scales='free')+
  geom_tile()+scale_fill_gradient2(low='red',high='blue')+theme(axis.text=element_blank())

ggplot(tmp,aes(cond,value))+theme_classic()+stat_summary(fun='mean',geom='point')

mt_Sep<-unique(hmap_mt[hmap_mt$var2>2 & hmap_mt$chromosome=='MT',colnames(hmap_mt) %in% c('trf',colnames(ppmiMeta))])

#
## qPCR-matching tRFs based on the feature selection process #####
#creating counts data for fully and optimally matched:
chp_tRFs0<-trfs0[rownames(trfs0) %in% subset(ppmiMeta$MINTbase.Unique.ID,ppmiMeta$mtfs3!='No motifs'),]
chp_trfs2<-chp_tRFs0[,colnames(chp_tRFs0) %in% m.data2$sample] ; chp_trfs<-chp_tRFs0[,colnames(chp_tRFs0) %in% m.data$sample]

# Fully matched:
f_chip<-data.frame(sample=colnames(chp_trfs2),sum=colSums(chp_trfs2),
                   mtf=colSums(chp_trfs2[rownames(chp_trfs2) %in% subset(ppmiMeta$MINTbase.Unique.ID,ppmiMeta$mtfs3=='RGTTCRA'),]))
f_chip$ratio<-100*f_chip$mtf/f_chip$sum ; f_chip<-merge(f_chip,m.data2,by='sample')
f_chip$nrm_ratio<-unlist(lapply(1:nrow(f_chip),function(n)
  (f_chip$ratio[n]/mean(subset(f_chip$ratio,f_chip$subclass==f_chip$subclass[n])))))#/sd(subset(f_chip$ratio,f_chip$subclass==f_chip$subclass[n]))
tmp<-data.frame(mt=colSums(chp_trfs2[rownames(chp_trfs2) %in% subset(ppmiMeta$MINTbase.Unique.ID,ppmiMeta$mtfs3=='MT'),]),
                sample=colnames(chp_trfs2),all=colSums(chp_trfs2)) ; tmp$mt_pct<-100*tmp$mt/tmp$all ; f_chip<-merge(f_chip,tmp[,c("sample","mt","mt_pct")])
f_chip$mt_pct<-unlist(lapply(1:nrow(f_chip),function(n)
  (f_chip$mt_pct[n]/mean(subset(f_chip$mt_pct,f_chip$subclass==f_chip$subclass[n])))))#/sd(subset(f_chip$mt_pct,f_chip$subclass==f_chip$subclass[n]))
f_chip$mtf_mt<-f_chip$mtf/(f_chip$mt+1) ; f_chip$mtf_mt3<-log2(f_chip$mtf_mt)
f_chip$mtf_mt2<-unlist(lapply(1:nrow(f_chip),function(n)
  (f_chip$mtf_mt[n]/mean(subset(f_chip$mtf_mt,f_chip$subclass==f_chip$subclass[n])))))#/sd(subset(f_chip$mtf_mt,f_chip$subclass==f_chip$subclass[n]))
f_chip$race<-f_chip$Ethnicity ; for(g in c("Am_Indian_or_Alaska_Native","Ashkenazi_Jewish")){f_chip$race<-gsub(g,'Other',f_chip$race)}
f_chip$race<-factor(f_chip$race,levels=c('White','Hispanic_or_Latino','Asian','Other','Black_or_African_American'),labels=c('White','Hispanic','Other','Other','African American'))

# Optimally matched:
chip<-data.frame(sample=colnames(chp_trfs),sum=colSums(chp_trfs),
                 mtf=colSums(chp_trfs[rownames(chp_trfs) %in% subset(ppmiMeta$MINTbase.Unique.ID,ppmiMeta$mtfs3=='RGTTCRA'),]))
chip$ratio<-100*chip$mtf/chip$sum ; chip<-merge(chip,m.data,by='sample')
chip$nrm_ratio<-unlist(lapply(1:nrow(chip),function(n)
  (chip$ratio[n]/mean(subset(chip$ratio,chip$subclass==chip$subclass[n])))))#/sd(subset(chip$ratio,chip$subclass==chip$subclass[n]))
tmp<-data.frame(mt=colSums(chp_trfs[rownames(chp_trfs) %in% subset(ppmiMeta$MINTbase.Unique.ID,ppmiMeta$mtfs3=='MT'),]),
                sample=colnames(chp_trfs),all=colSums(chp_trfs)) ; tmp$mt_pct<-100*tmp$mt/tmp$all ; chip<-merge(chip,tmp[,c("sample","mt","mt_pct")])
chip$mt_pct<-unlist(lapply(1:nrow(chip),function(n)
  (chip$mt_pct[n]/mean(subset(chip$mt_pct,chip$subclass==chip$subclass[n])))))#/sd(subset(chip$mt_pct,chip$subclass==chip$subclass[n]))
toRmv<-c() ; for (s in unique(chip$subclass)){tmp<-subset(chip,chip$subclass==s) 
if(sum(tmp$cond=='Ctrl')==0 | sum(tmp$cond=='Prodromal')==0){toRmv<-append(toRmv,s)}} ; chip<-subset(chip,!chip$subclass %in% toRmv)
chip$mtf_mt<-chip$mtf/(chip$mt+1) ; chip$mtf_mt3<-log2(chip$mtf_mt)
chip$mtf_mt2<-unlist(lapply(1:nrow(chip),function(n)
  (chip$mtf_mt[n]/mean(subset(chip$mtf_mt,chip$subclass==chip$subclass[n])))))#/sd(subset(chip$mtf_mt,chip$subclass==chip$subclass[n]))
chip$mtf2<-unlist(lapply(1:nrow(chip),function(n)
  (chip$mtf[n]/mean(subset(chip$mtf,chip$subclass==chip$subclass[n])))))#/sd(subset(chip$mtf,chip$subclass==chip$subclass[n]))
chip$mt2<-unlist(lapply(1:nrow(chip),function(n)
  (chip$mt[n]/mean(subset(chip$mt,chip$subclass==chip$subclass[n])))))#/sd(subset(chip$mt,chip$subclass==chip$subclass[n]))
chip$race<-chip$Ethnicity ; for(g in c("Am_Indian_or_Alaska_Native","Ashkenazi_Jewish")){chip$race<-gsub(g,'Other',chip$race)}
chip$race<-factor(chip$race,levels=c('White','Hispanic_or_Latino','Other','Black_or_African_American'),labels=c('White','Hispanic','Other','African American'))
chip<-merge(chip,dig,by.x='Subject',by.y='PATNO')

# Visualise data:
ggplot(chip,aes(Genetic_background,(mtf_mt2),group=interaction(cond,Genetic_background),col=cond,fill=cond))+theme_classic()+
  geom_boxplot(outlier.shape=NA,position = position_dodge2(preserve = "single"))+
  geom_point(position=position_jitterdodge(jitter.width=0.3,jitter.height=0))+facet_grid(~race,space='free',scales='free_x')+
  stat_summary(fun='mean',geom='point',shape=18,size=3,position=position_dodge(width=0.75),col='black')+#ylim(0,100)+
  scale_fill_manual(values=c('grey75','lightblue'))+scale_color_manual(values=c('grey3','darkblue'))+
  ylab('RGTTCRA/MT tRFs')+theme(axis.title.x=element_blank(),legend.title=element_blank())#+geom_hline(yintercept = 5)
# ggsave(dir1('Prodromal_mt.mtf.motifs_boxplots.svg'),width=9.5,height=5)
TukeyHSD(aov(chip$mtf_mt2~chip$cond*chip$Genetic_background))$`chip$cond`

ggplot(chip,aes(Genetic_background,(UPDRS.score.III),group=interaction(cond,Genetic_background),col=cond,fill=cond))+theme_classic()+
  geom_boxplot(outlier.shape=NA,position = position_dodge2(preserve = "single"))+
  geom_point(position=position_jitterdodge(jitter.width=0.3,jitter.height=0.2))+facet_grid(~race,space='free',scales='free_x')+
  stat_summary(fun='mean',geom='point',shape=18,size=3,position=position_dodge(width=0.75),col='black')+#ylim(0,500)+
  scale_fill_manual(values=c('grey75','lightblue'))+scale_color_manual(values=c('grey3','darkblue'))+
  ylab('UPDRSIII')+theme(axis.title.x=element_blank(),legend.title=element_blank())
# ggsave(dir1('Prodromal_UPDRS_boxplots.svg'),width=9.5,height=5)
TukeyHSD(aov(chip$UPDRS.score.III~chip$cond*chip$Genetic_background))

tmp<-subset(chip,chip$DIAG1=='PD')
ggplot(tmp,aes(UPDRS.score.III,(mtf_mt2)))+theme_classic()+geom_point()+ylim(0,2)+geom_abline(slope=2/15)+
  xlab('UPDRS score')+ylab('tRFs Score')
  geom_hline(yintercept=1)+geom_vline(xintercept=2.1)
ggsave(dir1('Dagnosed_tRFS_vs_UPDRS.svg'),width=3,height=3)

tmp<-subset(chip,chip$DIAG1=='PD' | (chip$cond=='Ctrl' & !is.na(chip$UPDRS.score.III)))
ggplot(tmp,aes(cond,mtf_mt2,col=cond))+theme_classic()+geom_boxplot(outlier.shape=NA)+geom_point(position=position_jitter(width=0.3,height=0))+
  stat_summary(fun='mean',geom='point',shape=18,col='black',size=3)+theme(legend.position='none',axis.title.x=element_blank())+ylab('tRFs score')
# ggsave(dir1('Diagnosed_tRFs2.svg'),width=1.85,height=3.1)
t.test(tmp$mtf_mt2~tmp$cond)

ggplot(tmp,aes(cond,UPDRS.score.III,col=cond))+theme_classic()+geom_boxplot(outlier.shape=NA)+geom_point(position=position_jitter(width=0.3,height=0))+
  stat_summary(fun='mean',geom='point',shape=18,col='black',size=3)+theme(legend.position='none',axis.title.x=element_blank())+ylab('UPDRS III')
# ggsave(dir1('Diagnosed_UPDRS.svg'),width=1.85,height=3.1)
t.test(tmp$UPDRS.score.III~tmp$cond)

p.adjust(c(5.715e-06,0.006626),'fdr')
  
ggplot(f_chip,aes(Genetic_background,log10(mtf_mt),group=interaction(cond,Genetic_background),col=cond,fill=cond))+theme_classic()+
  geom_boxplot(outlier.shape=NA,position = position_dodge2(preserve = "single"))+
  geom_point(position=position_jitterdodge(jitter.width=0.3,jitter.height=0))+facet_grid(~race,space='free',scales='free_x')+
  stat_summary(fun='mean',geom='point',shape=18,size=3,position=position_dodge(width=0.75),col='black')+#ylim(0,100)+
  scale_fill_manual(values=c('grey75','lightblue'))+scale_color_manual(values=c('grey3','darkblue'))+
  ylab('Log10(RGTTCRA/MT tRFs)')+theme(axis.title.x=element_blank(),legend.title=element_blank())#+geom_hline(yintercept = 5)
# ggsave(dir1('Prodromal_mt.mtf.motifs_boxplots_unNormalized_FullMatched.svg'),width=9.5,height=5)
TukeyHSD(aov(log10(f_chip$mtf_mt)~f_chip$cond*f_chip$Genetic_background))$`f_chip$cond`
tmp<-chisq.test(table(f_chip$cond,f_chip$mtf_mt>=6)) ; tmp$p.value ; tmp$observed ; tmp$expected 

###

## ROC curves optimal matching #####

# Nrm-PD/MT-tRFs
roc_chp0<-chip[,c("cond","mtf_mt2","Ethnicity")]
ctrl0 <- trainControl(method="cv", summaryFunction=twoClassSummary, number=5, classProbs=T,savePredictions = T)
c_fit_all0 <- train(cond ~ .,data=roc_chp0,method="gbm",trControl=ctrl0,na.action=na.omit)

# Clinical scores
roc_chp1<-chip[,c("cond","Ethnicity","Hoehn.and.Yahr.staging","UPDRS.score.III")]
ctrl1 <- trainControl(method="cv", summaryFunction=twoClassSummary, number=5, classProbs=T,savePredictions = T)
c_fit_all1 <- train(cond ~ .,data=roc_chp1,method="gbm",trControl=ctrl1,na.action=na.omit)

# Mixed lables
fk_mtf<-chip ; fk_mtf$cond<-'Ctrl' ; fk_mtf$cond[sample(1:nrow(fk_mtf),60,replace=F)]<-'Prodromal'
smps1<-c(sample(subset(fk_mtf$sample,fk_mtf$cond=='Prodromal'),60,replace=F),sample(subset(fk_mtf$sample,fk_mtf$cond=='Ctrl'),60,replace=F))
tmp1<-data.frame(sample=smps1) ; tmp1<-merge(tmp1,fk_mtf,by='sample')

roc_chp2<-tmp1[,c("cond","mtf_mt2","Ethnicity","Hoehn.and.Yahr.staging","UPDRS.score.III")]
ctrl2 <- trainControl(method="cv", summaryFunction=twoClassSummary, number=5, classProbs=T,savePredictions = T)
c_fit_all2 <- train(cond ~ .,data=roc_chp2,method="gbm",trControl=ctrl2,na.action=na.omit)

# ROC
set.seed(123) ; res<-evalm(list1=list(c_fit_all2,c_fit_all1,c_fit_all0),gnames=c('Mixed lables','Only clinic','Only tRFs'))
# ggsave(dir1('ROC_Prodromal_chip.svg'))


## Training on optimal k=1 and testing on full #####

#tRFs
set.seed(123)
l0<-c("cond","mtf_mt2","Ethnicity") ; roc_chp0<-chip[,l0]
ctrl0 <- trainControl(method="none", summaryFunction=twoClassSummary,number=1,classProbs=T,savePredictions = T)
c_fit_all0_k1 <- train(cond ~ .,data=roc_chp0,method="gbm",trControl=ctrl0,na.action=na.omit)

tmp1<-subset(f_chip,f_chip$Ethnicity!='Asian')
chp_mtfs<-data.frame(predict(c_fit_all0_k1, newdata = tmp1[,l0], type = "prob"),
                     pred=predict(c_fit_all0_k1, newdata = tmp1[,l0]),true=tmp1$cond,PATNO=tmp1$Subject)
chp_mtfs<-merge(chp_mtfs,dig,by='PATNO')
ggplot(chp_mtfs,aes(true,Prodromal,col=true))+theme_classic()+ylim(0,1)+geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitter(width=0.1,height=0.01))+geom_hline(yintercept=0.5,col='black',linetype='dotted')
# ggsave(dir1('Prodromal_TestOnFullMatched_tRFs.svg'),width=3.6,height=7)

ggplot(chp_clnc,aes(DIAG1=='PD',Prodromal,col=true))+theme_classic()+ylim(0,1)+geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=0.1,jitter.height=0.01))+geom_hline(yintercept=0.5,col='black',linetype='dotted')
table(tolower(chp_mtfs$pred),chp_mtfs$true)*100/nrow(chp_mtfs)
#           Ctrl  Prodromal
# ctrl       89      14
# prodromal  21      46

#Clinical scores
l1<-c("cond","Ethnicity","Hoehn.and.Yahr.staging","UPDRS.score.III") ; roc_chp1<-chip[,l1]
ctrl1 <- trainControl(method="none", summaryFunction=twoClassSummary,number=1,classProbs=T,savePredictions = T)
c_fit_all1_k1 <- train(cond ~ .,data=roc_chp1,method="gbm",trControl=ctrl1,na.action=na.omit)

tmp1<-subset(f_chip,!is.na(f_chip$UPDRS.score.III) & f_chip$Ethnicity!='Asian')
chp_clnc<-data.frame(predict(c_fit_all1_k1, newdata = tmp1[,l1], type = "prob"),
                     pred=predict(c_fit_all1_k1, newdata = tmp1[,l1]),true=tmp1$cond,PATNO=tmp1$Subject) 
chp_clnc<-merge(chp_clnc,dig,by='PATNO')

ggplot(chp_clnc,aes(true,Prodromal,col=true))+theme_classic()+ylim(0,1)+geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitter(width=0.1,height=0.01))+geom_hline(yintercept=0.5,col='black',linetype='dotted')
# ggsave(dir1('Prodromal_TestOnFullMatched_Clinics.svg'),width=3.6,height=7)

ggplot(chp_clnc,aes(DIAG1=='PD',Prodromal,col=true))+theme_classic()+ylim(0,1)+geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width=0.1,jitter.height=0.01))+geom_hline(yintercept=0.5,col='black',linetype='dotted')

table(tolower(chp_clnc$pred),chp_clnc$true)*100/nrow(chp_clnc)
#           Ctrl    Prodromal
# ctrl        75        15
# prodromal   35        45

# Compare predicted score to true diagnosis:
tmp<-subset(chp_mtfs,chp_mtfs$DIAG1=='PD') ; colnames(tmp)<-gsub('Prodromal','tRFs',colnames(tmp))
diagnosed<-subset(chp_clnc,chp_clnc$DIAG1=='PD') ; colnames(diagnosed)<-gsub('Prodromal','Clinical',colnames(diagnosed))
diagnosed<-merge(diagnosed[,c("PATNO","Clinical")],tmp,by='PATNO')

ggplot(diagnosed,aes(Clinical,tRFs,col=DIAG1VIS))+theme_classic()+geom_point()
tmp<-diagnosed[,c("PATNO","tRFs","Clinical")]; diag2<-melt(tmp,id.vars='PATNO',variable_name='method')
diag2$method<-factor(diag2$method,levels=c('Clinical','tRFs'))

ggplot(diag2,aes(method,value,col=method))+theme_classic()+geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitter(width=0.3,height=0))+stat_summary(fun='mean',geom='point',shape=18,col='black',size=3)+
  theme(legend.position='none',axis.title.x=element_blank())+ylab('Prodromal score')+ylim(0,1)
# ggsave(dir1('Diagnosed.svg'),width=2,height=6.5)
table(diag2$value>0.55,diag2$method)
  
#Null
fk_mtf<-chip ; fk_mtf$cond<-'Ctrl' ; fk_mtf$cond[sample(1:nrow(fk_mtf),60,replace=F)]<-'Prodromal'
smps1<-c(sample(subset(fk_mtf$sample,fk_mtf$cond=='Prodromal'),60,replace=F),sample(subset(fk_mtf$sample,fk_mtf$cond=='Ctrl'),60,replace=F))
tmp0<-data.frame(sample=smps1) ; tmp0<-merge(tmp0,fk_mtf,by='sample')

l2<-c("cond","mtf_mt2","Ethnicity","Hoehn.and.Yahr.staging","UPDRS.score.III") ; roc_chp2<-tmp0[,l2]
ctrl2 <- trainControl(method="none", summaryFunction=twoClassSummary,number=1,classProbs=T,savePredictions = T)
c_fit_all2_k1 <- train(cond ~ .,data=roc_chp2,method="gbm",trControl=ctrl2,na.action=na.omit)

tmp1<-subset(f_chip,!is.na(f_chip$UPDRS.score.III) & f_chip$Ethnicity!='Asian')
chp_ctrl<-data.frame(predict(c_fit_all2_k1, newdata = tmp1[,l2], type = "prob"),
                     pred=predict(c_fit_all2_k1, newdata = tmp1[,l2]),true=tmp1$cond,PATNO=tmp1$Subject) 
chp_ctrl<-merge(chp_ctrl,dig,by='PATNO')

ggplot(chp_ctrl,aes(true,Prodromal,col=true))+theme_classic()+ylim(0,1)+geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitter(width=0.1,height=0.01))+geom_hline(yintercept=0.5,col='black',linetype='dotted')
# ggsave(dir1('Prodromal_TestOnFullMatched_Null.svg'),width=3.6,height=7)
table(tolower(chp_ctrl$pred),chp_ctrl$true)*100/nrow(chp_ctrl)
#             Ctrl Prodromal
# ctrl        67        37
# prodromal   43        23

tmp<-chisq.test(table(tolower(chp_mtfs$pred),chp_mtfs$true))
chisq.test(table(tolower(chp_clnc$pred),chp_clnc$true))
p.adjust(c(2.51e-14,2.693e-07),'fdr')

(table(tolower(chp_mtfs$pred),chp_mtfs$true))/(table(tolower(chp_clnc$pred),chp_clnc$true))

chisq.test(as.numeric(table(tolower(chp_mtfs$pred),chp_mtfs$true)),
           p=as.numeric(table(tolower(chp_clnc$pred),chp_clnc$true)/nrow(chp_clnc)))

#
# ROC curves full matching in Prodromals #####

c_all_true<-c() ; c_all_fake<-c() ; c_clnc_true<-c() ; c_tRF_true<-c()
for(i in 1:10000){
  smps<-c(sample(subset(f_chip$sample,f_chip$cond=='Prodromal'),60,replace=F),sample(subset(f_chip$sample,f_chip$cond=='Ctrl'),60,replace=F))
  tmp<-data.frame(sample=smps) ; tmp<-merge(tmp,f_chip,by='sample')
  
  # Nrm-PD/MT-tRFs + Clinical scores
  roc_chp<-tmp[,c("cond","mtf_mt2","Ethnicity","Hoehn.and.Yahr.staging","UPDRS.score.III")]
  ctrl <- trainControl(method="cv", summaryFunction=twoClassSummary, number=5, classProbs=T,savePredictions = T)
  fit_all <- train(cond ~ .,data=roc_chp,method="gbm",trControl=ctrl,na.action=na.omit)
  res <- evalm(fit_all) ; c_all_true<-append(c_all_true,as.numeric(res$stdres$`Group 1`['AUC-ROC','Score']))
  
  # Nrm-PD/MT-tRFs
  roc_chp0<-tmp[,c("cond","mtf_mt2","Ethnicity")]
  ctrl0 <- trainControl(method="cv", summaryFunction=twoClassSummary, number=5, classProbs=T,savePredictions = T)
  fit_all0 <- train(cond ~ .,data=roc_chp0,method="gbm",trControl=ctrl0,na.action=na.omit)
  res0 <- evalm(fit_all0) ; c_tRF_true<-append(c_tRF_true,as.numeric(res0$stdres$`Group 1`['AUC-ROC','Score']))
  
  # Clinical scores
  roc_chp1<-tmp[,c("cond","Ethnicity","Hoehn.and.Yahr.staging","UPDRS.score.III")]
  ctrl1 <- trainControl(method="cv", summaryFunction=twoClassSummary, number=5, classProbs=T,savePredictions = T)
  fit_all1 <- train(cond ~ .,data=roc_chp1,method="gbm",trControl=ctrl1,na.action=na.omit)
  res1 <- evalm(fit_all1) ; c_clnc_true<-append(c_clnc_true,as.numeric(res1$stdres$`Group 1`['AUC-ROC','Score']))
  
  # Mixed lables
  fk_mtf<-f_chip ; fk_mtf$cond<-'Ctrl' ; fk_mtf$cond[sample(1:nrow(fk_mtf),60,replace=F)]<-'Prodromal'
  smps1<-c(sample(subset(fk_mtf$sample,fk_mtf$cond=='Prodromal'),60,replace=F),sample(subset(fk_mtf$sample,fk_mtf$cond=='Ctrl'),40,replace=F))
  tmp1<-data.frame(sample=smps1) ; tmp1<-merge(tmp1,fk_mtf,by='sample')
  
  roc_chp2<-tmp1[,c("cond","mtf_mt2","Ethnicity","Hoehn.and.Yahr.staging","UPDRS.score.III")]
  ctrl2 <- trainControl(method="cv", summaryFunction=twoClassSummary, number=5, classProbs=T,savePredictions = T)
  fit_all2 <- train(cond ~ .,data=roc_chp2,method="gbm",trControl=ctrl2,na.action=na.omit)
  res2 <- evalm(fit_all2) ; c_all_fake<-append(c_all_fake,as.numeric(res2$stdres$`Group 1`['AUC-ROC','Score']))
  
  print(paste0('                                                                                ',i))
  print(paste0('                                                                                ',i))
  print(paste0('                                                                                ',i))
}

chp_pv<-data.frame('all_true'=c_all_true,'all_fake'=c_all_fake,'clnc_true'=c_clnc_true,'tRF_true'=c_tRF_true) 
chp_pv$id<-rownames(chp_pv) ; chp_pv<-melt(chp_pv,id.vars='id')
chp_pv$variable<-factor(chp_pv$variable,levels=c('all_true','tRF_true','clnc_true','all_fake'),
                        labels=c('tRFs + UPDRS III + H&Y','Only tRFs','UPDRS III + H&Y','Mixed labels'))
ggplot(chp_pv,aes(variable,value,col=variable))+theme_classic()+geom_boxplot()+ylab('ROC-AUC')+
  theme(legend.title=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank())
TukeyHSD(aov(chp_pv$value~chp_pv$variable)) ; #ggsave(dir1('Prodromal_ROC_Chip_Boxplot.svg'))

# Training and Testing on Prodromals #####
set.seed(123)

tmp<-chip[,c("cond","mtf_mt2","Ethnicity","Hoehn.and.Yahr.staging","UPDRS.score.III")]
tmp$Ethnicity<-gsub('Am_Indian_or_Alaska_Native','Other',tmp$Ethnicity) ; tmp$Ethnicity<-gsub('Ashkenazi_Jewish','Other',tmp$Ethnicity)
ids<-createDataPartition(tmp$cond, p = 0.65, list = FALSE)
chp_trn<-tmp[ids,] ; chp_tst<-tmp[-ids,]

roc_chp0<-chp_trn[,c("cond","mtf_mt2","Ethnicity")]
ctrl0 <- trainControl(method="cv", summaryFunction=twoClassSummary, number=5, classProbs=T,savePredictions = T)
c_fit_all0 <- train(cond ~ .,data=roc_chp0,method="gbm",trControl=ctrl0,na.action=na.omit)

roc_chp1<-chp_trn[,c("cond","Ethnicity","Hoehn.and.Yahr.staging","UPDRS.score.III")]
ctrl1 <- trainControl(method="cv", summaryFunction=twoClassSummary, number=5, classProbs=T,savePredictions = T)
c_fit_all1 <- train(cond ~ .,data=roc_chp1,method="gbm",trControl=ctrl1,na.action=na.omit)

n<-as.numeric(table(chp_trn$cond)['Prodromal']) ; fk_mtf<-chp_trn ; fk_mtf$cond<-'Ctrl' ; fk_mtf$cond[sample(1:nrow(fk_mtf),n,replace=F)]<-'Prodromal'
roc_chp2<-fk_mtf[,c("cond","mtf_mt2","Ethnicity","Hoehn.and.Yahr.staging","UPDRS.score.III")]
ctrl2 <- trainControl(method="cv", summaryFunction=twoClassSummary, number=5, classProbs=T,savePredictions = T)
c_fit_all2 <- train(cond ~ .,data=roc_chp2,method="gbm",trControl=ctrl2,na.action=na.omit)

tmp1<-subset(chp_tst,!is.na(chp_tst$UPDRS.score.III))
chp_ctrl<-data.frame(predict(c_fit_all2, newdata = tmp1[,c("cond","mtf_mt2","Ethnicity","Hoehn.and.Yahr.staging","UPDRS.score.III")], type = "prob"),
                     pred=predict(c_fit_all2, newdata = tmp1[,c("cond","mtf_mt2","Ethnicity","Hoehn.and.Yahr.staging","UPDRS.score.III")]),true=tmp1$cond) 
ggplot(chp_ctrl,aes(true,Prodromal,col=true))+theme_classic()+geom_boxplot(outlier.shape=NA)+geom_point(position=position_jitter(width=0.1,height=0.01))
table(tolower(chp_ctrl$pred),chp_ctrl$true)*100/nrow(chp_ctrl)

tmp1<-subset(chp_tst,!is.na(chp_tst$UPDRS.score.III))
chp_clnc<-data.frame(predict(c_fit_all1, newdata = tmp1[,c("cond","Ethnicity","Hoehn.and.Yahr.staging","UPDRS.score.III")], type = "prob"),
                     pred=predict(c_fit_all1, newdata = tmp1[,c("cond","Ethnicity","Hoehn.and.Yahr.staging","UPDRS.score.III")]),true=tmp1$cond) 
ggplot(chp_clnc,aes(true,Prodromal,col=true))+theme_classic()+geom_boxplot(outlier.shape=NA)+geom_point(position=position_jitter(width=0.1,height=0.01))
table(tolower(chp_clnc$pred),chp_clnc$true)*100/nrow(chp_clnc)

chp_mtfs<-data.frame(predict(c_fit_all0, newdata = chp_tst[,c("cond","mtf_mt2","Ethnicity")], type = "prob"),
                     pred=predict(c_fit_all0, newdata = chp_tst[,c("cond","mtf_mt2","Ethnicity")]),true=chp_tst$cond) 
ggplot(chp_mtfs,aes(true,Prodromal,col=true))+theme_classic()+geom_boxplot(outlier.shape=NA)+geom_point(position=position_jitter(width=0.1,height=0.01))
table(tolower(chp_mtfs$pred),chp_mtfs$true)*100/nrow(chp_mtfs)


# Checking qPCR-matching tRFs in PDs #####

v_ppmi<-subset(col_ppmi,!col_ppmi$Group %in% c('Prodromal','SWEDD')) ; v_ppmi$cond<-droplevels(v_ppmi$cond)
v_ppmi<-subset(v_ppmi,v_ppmi$Other_diseases %in% c('0','NA') & v_ppmi$timePoint %in% c('BL') & v_ppmi$PD_meds %in% c('non','UnKnown') &
                 v_ppmi$Ethnicity %in% c('Ashkenazi_Jewish','Basque','Black_or_African_American','Hispanic_or_Latino','White') &
                 (v_ppmi$PD_duration==0 | is.na(v_ppmi$PD_duration)) & v_ppmi$s_Study=='phase1') 

# write.csv(v_ppmi,dir1('Supp_ROC_PDs.csv'))

toKeep<-c()
for(s in unique(v_ppmi$Subject)){
  tmp<-subset(v_ppmi,v_ppmi$Subject==s) ; toKeep<-append(toKeep,as.character(tmp$sample[order(tmp$timePoint,decreasing=T)][1]))}
tmp1<-subset(v_ppmi,v_ppmi$sample %in% toKeep)

out_opt<-matchit(cond ~ Sex+Age_timepoint,data=tmp1,method='optimal',link='probit',distance='glm')
plot(out_opt, type = "density", interactive = FALSE,which.xs = ~Age_timepoint + Sex)
v_m.data<-match.data(out_opt) ; summary(out_opt,un=T) ; plot(summary(out_opt))
table(v_m.data$subclass,v_m.data$cond) ; plot(out_opt, type = "jitter", interactive = FALSE)

out_opt2<-matchit(cond ~ Sex+Age_timepoint,data=tmp1,method='full',link='probit',distance='glm')
plot(out_opt2, type = "density", interactive = FALSE,which.xs = ~Age_timepoint + Sex)
v_m.data2<-match.data(out_opt2) ; summary(out_opt2,un=T) ; plot(summary(out_opt2))
table(v_m.data2$subclass,v_m.data2$cond) ; plot(out_opt2, type = "jitter", interactive = FALSE)


# Motifs
v_trfs0<-ctsPPMI ; rownames(v_trfs0)<-v_trfs0$X ; v_trfs0<-v_trfs0[,-1]
v_trfs0<-v_trfs0[rownames(v_trfs0) %in% subset(ppmiMeta$MINTbase.Unique.ID,ppmiMeta$mtfs3!='No motifs'),]
colnames(v_trfs0)<-unlist(lapply(colnames(v_trfs0),function(x) paste0(strsplit(strsplit(x,'IR1.')[[1]][2],'\\.')[[1]][1:2],collapse='_')))
v_trfs2<-v_trfs0[,colnames(v_trfs0) %in% v_m.data2$sample] ; v_trfs<-v_trfs0[,colnames(v_trfs0) %in% v_m.data$sample]

v_motif<-data.frame(sample=colnames(v_trfs2),sum=colSums(v_trfs2),
                    mtf=colSums(v_trfs2[rownames(v_trfs2) %in% subset(ppmiMeta$MINTbase.Unique.ID,ppmiMeta$mtfs3=='RGTTCRA'),]),
                    mt=colSums(v_trfs2[rownames(v_trfs2) %in% subset(ppmiMeta$MINTbase.Unique.ID,ppmiMeta$mtfs3=='MT'),]))
v_motif$mtf_mt<-v_motif$mtf/(v_motif$mt+1) ; v_motif<-merge(v_motif,v_m.data2,by='sample')
v_motif$mtf_mt2<-unlist(lapply(1:nrow(v_motif),function(n)
  (v_motif$mtf_mt[n]/mean(subset(v_motif$mtf_mt,v_motif$subclass==v_motif$subclass[n])))))#/sd(subset(v_motif$mtf_mt,v_motif$subclass==v_motif$subclass[n]))

v_motif1<-data.frame(sample=colnames(v_trfs),sum=colSums(v_trfs),
                     mtf=colSums(v_trfs[rownames(v_trfs) %in% subset(ppmiMeta$MINTbase.Unique.ID,ppmiMeta$mtfs3=='RGTTCRA'),]),
                     mt=colSums(v_trfs[rownames(v_trfs) %in% subset(ppmiMeta$MINTbase.Unique.ID,ppmiMeta$mtfs3=='MT'),]))
v_motif1$mtf_mt<-v_motif1$mtf/(v_motif1$mt+1) ; v_motif1<-merge(v_motif1,v_m.data,by='sample')
v_motif1$mtf_mt2<-unlist(lapply(1:nrow(v_motif1),function(n)
  (v_motif1$mtf_mt[n]/mean(subset(v_motif1$mtf_mt,v_motif1$subclass==v_motif1$subclass[n])))))#/sd(subset(v_motif1$mtf_mt,v_motif1$subclass==v_motif1$subclass[n]))
toRmv<-c() ; for (s in unique(v_motif1$subclass)){tmp<-subset(v_motif1,v_motif1$subclass==s) 
if(sum(tmp$cond=='Ctrl')==0 | sum(tmp$cond=='PD')==0){toRmv<-append(toRmv,s)}} ; v_motif1<-subset(v_motif1,!v_motif1$subclass %in% toRmv)

ggplot(v_motif1,aes(Genetic_background,mtf_mt2,col=cond))+theme_classic()+geom_boxplot()+geom_point(position=position_jitterdodge(jitter.height=0))#+facet_wrap(~Ethnicity)
ggplot(v_motif1,aes(Genetic_background,log10(mtf),col=cond))+theme_classic()+geom_boxplot()+geom_point(position=position_jitterdodge(jitter.height=0))#+facet_wrap(~Ethnicity)
ggplot(v_motif1,aes(Genetic_background,mt,col=cond))+theme_classic()+geom_boxplot()+geom_point(position=position_jitterdodge(jitter.height=0))#+facet_wrap(~Ethnicity)


roc_all0<-v_motif1[,c("cond","mtf_mt2","Genetic_background","Ethnicity")]
ctrl0 <- trainControl(method="cv", summaryFunction=twoClassSummary, number=5, classProbs=T,savePredictions = T)
v_fit_all0 <- train(cond ~ .,data=roc_all0,method="gbm",trControl=ctrl0,na.action=na.omit)

roc_all1<-v_motif1[,c("cond","Ethnicity","Genetic_background","Hoehn.and.Yahr.staging","UPDRS.score.III")]
ctrl1 <- trainControl(method="cv", summaryFunction=twoClassSummary, number=5, classProbs=T,savePredictions = T)
v_fit_all1 <- train(cond ~ .,data=roc_all1,method="gbm",trControl=ctrl1,na.action=na.omit)

fk_mtf<-v_motif1 ; fk_mtf$cond<-'Ctrl' ; fk_mtf$cond[sample(1:nrow(fk_mtf),193,replace=F)]<-'PD'
roc_all2<-fk_mtf[,c("cond","mtf_mt2","Genetic_background","Ethnicity","Hoehn.and.Yahr.staging","UPDRS.score.III")]
ctrl2 <- trainControl(method="cv", summaryFunction=twoClassSummary, number=5, classProbs=T,savePredictions = T)
v_fit_all2 <- train(cond ~ .,data=roc_all2,method="gbm",trControl=ctrl2,na.action=na.omit)

roc_all<-v_motif1[,c("cond","mtf_mt2","Genetic_background","Ethnicity","Hoehn.and.Yahr.staging","UPDRS.score.III")]
ctrl <- trainControl(method="cv", summaryFunction=twoClassSummary, number=5, classProbs=T,savePredictions = T)
v_fit_all <- train(cond ~ .,data=roc_all,method="gbm",trControl=ctrl,na.action=na.omit)

set.seed(123) ; res<-evalm(list1=list(v_fit_all,v_fit_all2,v_fit_all1,v_fit_all0),gnames=c('Clinic & tRFs','Mixed lables','Only clinic','Only tRFs'))
# ggsave('C:/Users/Nimrod/Desktop/ROC_BL.svg')


# ROC curves full matching

v_all_true<-c() ; v_all_fake<-c() ; v_clnc_true<-c() ; v_tRF_true<-c()
for(i in 1:10000){
  smps<-c(sample(subset(v_motif$sample,v_motif$cond=='PD'),192,replace=F),sample(subset(v_motif$sample,v_motif$cond=='Ctrl'),192,replace=F))
  tmp<-data.frame(sample=smps) ; tmp<-merge(tmp,v_motif,by='sample')
  
  roc_val<-tmp[,c("cond","nrm_ratio","mt_pct","Ethnicity","Genetic_background","Hoehn.and.Yahr.staging","UPDRS.score.III")]
  ctrl <- trainControl(method="cv", summaryFunction=twoClassSummary, number=5, classProbs=T,savePredictions = T)
  fit_all <- train(cond ~ .,data=roc_val,method="gbm",trControl=ctrl,na.action=na.omit)
  res <- evalm(fit_all) ; v_all_true<-append(v_all_true,as.numeric(res$stdres$`Group 1`['AUC-ROC','Score']))
  
  roc_val0<-tmp[,c("cond","nrm_ratio","mt_pct","Ethnicity","Genetic_background")]
  ctrl0 <- trainControl(method="cv", summaryFunction=twoClassSummary, number=5, classProbs=T,savePredictions = T)
  fit_all0 <- train(cond ~ .,data=roc_val0,method="gbm",trControl=ctrl0,na.action=na.omit)
  res0 <- evalm(fit_all0) ; v_tRF_true<-append(v_tRF_true,as.numeric(res0$stdres$`Group 1`['AUC-ROC','Score']))
  
  roc_val1<-tmp[,c("cond","Ethnicity","Genetic_background","Hoehn.and.Yahr.staging","UPDRS.score.III")]
  ctrl1 <- trainControl(method="cv", summaryFunction=twoClassSummary, number=5, classProbs=T,savePredictions = T)
  fit_all1 <- train(cond ~ .,data=roc_val1,method="gbm",trControl=ctrl1,na.action=na.omit)
  res1 <- evalm(fit_all1) ; v_clnc_true<-append(v_clnc_true,as.numeric(res1$stdres$`Group 1`['AUC-ROC','Score']))
  
  fk_mtf<-v_motif ; fk_mtf$cond<-'Ctrl' ; fk_mtf$cond[sample(1:nrow(fk_mtf),193,replace=F)]<-'PD'
  smps1<-c(sample(subset(fk_mtf$sample,fk_mtf$cond=='PD'),192,replace=F),sample(subset(fk_mtf$sample,fk_mtf$cond=='Ctrl'),192,replace=F))
  tmp1<-data.frame(sample=smps1) ; tmp1<-merge(tmp1,fk_mtf,by='sample')
  
  roc_val2<-tmp1[,c("cond","nrm_ratio","mt_pct","Ethnicity","Genetic_background","Hoehn.and.Yahr.staging","UPDRS.score.III")]
  ctrl2 <- trainControl(method="cv", summaryFunction=twoClassSummary, number=5, classProbs=T,savePredictions = T)
  fit_all2 <- train(cond ~ .,data=roc_val2,method="gbm",trControl=ctrl2,na.action=na.omit)
  res2 <- evalm(fit_all2) ; v_all_fake<-append(v_all_fake,as.numeric(res2$stdres$`Group 1`['AUC-ROC','Score']))
  
  print(paste0('                                                                                ',i))
  print(paste0('                                                                                ',i))
  print(paste0('                                                                                ',i))
}

val_pv<-data.frame('all_true'=v_all_true,'all_fake'=v_all_fake,'clnc_true'=v_clnc_true,'tRF_true'=v_tRF_true) 
val_pv$id<-rownames(val_pv) ; val_pv<-melt(val_pv,id.vars='id') ; tmp<-subset(val_pv,val_pv$variable!='tRF_true')
ggplot(val_pv,aes(variable,value,col=variable))+theme_classic()+geom_boxplot() ; TukeyHSD(aov(val_pv$value~val_pv$variable))



