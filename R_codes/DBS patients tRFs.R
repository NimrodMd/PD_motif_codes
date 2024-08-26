library(ggplot2); library(ggfortify); library(dplyr); library(matrixStats)

ndf<-function(x) as.data.frame(t(data.frame(row.names = c(x)))) ; I<-function(x,y) x %in% y
dirin<-function(x) paste0(getwd(),x)
dirout<-function(x) paste0(getwd(),x)

DBS_data<-read.csv(dirin('PPMI_data.2.csv')) ; 
colnames(DBS_data)<-gsub('Type_of_surgery_for_Parkinson_disease_term','DBS',colnames(DBS_data))
DBS_data<-subset(DBS_data,DBS_data$cond!='Prodromal' & 
                !(DBS_data$Genetic_background=='LRRK2-/SNCA-/GBA-' & !is.na(DBS_data$DBS)))
DBS_data[is.na(DBS_data$DBS),'DBS']<-'a' ; DBS_data$DBS<-factor(as.factor(DBS_data$DBS),labels=c('No DBS','DBS'))
tmp1<-subset(DBS_data,DBS_data$s_Study=='phase2' & DBS_data$Genetic_background %in% c('LRRK2+','GBA+'))
table(tmp1$DBS,tmp1$cond) # ctrl=9, PD = 12, DBS=18 (13+5)

ggplot(tmp1,aes(interaction(DBS,cond),rt_motif,fill=cond,col=DBS))+theme_classic()+
  geom_boxplot(outlier.shape=NA)+stat_summary(fun='mean',geom='point',shape=18,stroke=3)+
  geom_point(position=position_jitterdodge(jitter.width=0.3,jitter.height=0))+ylab('% RGTTCTA-tRFs')+
  scale_color_manual(values=c('grey13','darkred'))+scale_fill_manual(values=c('grey33','#00A5FF','grey','lightblue'))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.title=element_blank()) 
# ggsave(dirout('PPMI_DBS_motif_ratio.svg'),device='svg',width=3.5,height=3.35)
TukeyHSD(aov(tmp1$rt_motif~interaction(tmp1$DBS,tmp1$cond))) 
# PD(noDBS)-ctrl p<0.027 ; PD(noDBS)-PD(DBS) p<0.0953 ; PD(DBS)-ctrl p<0.6

ggplot(tmp1,aes(interaction(DBS,cond),rt_mt,fill=cond,col=DBS))+theme_classic()+
  geom_boxplot(outlier.shape=NA)+stat_summary(fun='mean',geom='point',shape=18,stroke=3)+
  geom_point(position=position_jitterdodge(jitter.width=0.3,jitter.height=0))+ylab('% RGTTCTA-tRFs')+
  scale_color_manual(values=c('grey13','darkred'))+scale_fill_manual(values=c('grey33','#00A5FF','grey','lightblue'))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.title=element_blank()) 
# ggsave(dirout('PPMI_DBS_mt_ratio.svg'),device='svg',width=3.5,height=3.35)
TukeyHSD(aov(tmp1$rt_mt~interaction(tmp1$DBS,tmp1$cond))) 
# PD(noDBS)-ctrl p<0.88 ; PD(noDBS)-PD(DBS) p<0.86 ; PD(DBS)-ctrl p<0.47
