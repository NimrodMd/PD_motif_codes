library(ggplot2) ; library(matrixStats) ; library('PMCMRplus')
ndf<-function(x) as.data.frame(t(data.frame(row.names=x))) ; Ttmp<-function(x) as.data.frame(t(x))
Ftmp<-function(df,col='exp',id='id') {
  assign('tmp',as.data.frame(t(df))); colnames(tmp)<-col; tmp[,id]<-rownames(tmp) ; return(tmp)}
reform<-function(x){rownames(x)<-unlist(lapply(x[,1], function(y) strsplit(as.character(y),'\\.')[[1]][1])) ; x<-x[,2:ncol(x)]}
dir<-function(x) paste0(c(getwd(),x),collapse='/') ; dir('')

#
          # Brains #####
brn_qpcr<-read.csv(dir('Excels/brains_brn_qpcrs.csv'))
brn_qpcr$dif<-brn_qpcr$mtfBr_avg-brn_qpcr$mt_avg

brn_qpcr$exp_mtf<-0 ; brn_qpcr$exp_mt<-0 ; brn_qpcr$exp_chp<-0
for(i in 1:nrow(brn_qpcr)){
  brn_qpcr$exp_chp[i]<-2^(brn_qpcr$dif[i]-mean(subset(brn_qpcr$dif,brn_qpcr$sex==brn_qpcr$sex[i] & brn_qpcr$group=='Cont')))
}
  
ggplot(brn_qpcr,aes(cnd,exp_chp,col=sex,group=cnd))+theme_classic()+#facet_grid(~sex)+
  geom_boxplot(outlier.shape=NA)+geom_point(position=position_jitter(width=0.15,height=0))+
  stat_summary(fun='mean',geom='point',col='black',shape=18,size=3)+ylab('Relative expression')+
  theme(legend.position='none',axis.title.x=element_blank())
# ggsave(dir('SVGs & PNGs/Brain_qpcr.svg'),width=1.5,height=2.5)
t.test(brn_qpcr$exp_chp~brn_qpcr$cnd)

#
          # Bloods #####
bld_qpcr<-read.csv(dir('Excels/David _384 PD motif and MT motif on blood samples 2024-02-13 15-41-18_384.csv'))
updrs<-read.csv(dir('Excels/PD Trauma Control Blood RNA.csv')) ; colnames(updrs)<-gsub('Running.number','Samp',colnames(updrs)) 
updrs$Samp<-gsub('C0','C',updrs$Samp) ; updrs$Samp<-gsub('P0','PD',updrs$Samp) ; updrs$Samp<-gsub('T0','T',updrs$Samp)
updrs$Samp<-gsub('P10','PD10',updrs$Samp) ; updrs$Samp<-gsub('P11','PD11',updrs$Samp) ; bld_qpcr<-merge(bld_qpcr,updrs,by='Samp')

bld_qpcr$ddCq<-0 
for(i in 1:nrow(bld_qpcr)){
  bld_qpcr$ddCq[i]<-bld_qpcr$dif[i]-mean(subset(bld_qpcr$dif,bld_qpcr$Sex==bld_qpcr$Sex[i] & bld_qpcr$Diagnosis=='Control'))
} ; bld_qpcr$exp_chp<-2^bld_qpcr$ddCq ; bld_qpcr$cnd<-factor(bld_qpcr$Diagnosis,levels=c('PD','Control','Trauma'))

ggplot(bld_qpcr,aes(Diagnosis,exp_chp,col=Sex,group=Diagnosis))+theme_classic()+#facet_grid(~Sex)+
  geom_boxplot(outlier.shape=NA)+geom_point(position=position_jitter(width=0.15,height=0),size=3)+
  stat_summary(fun='mean',geom='point',col='black',shape=18,size=3)+ylab('Relative expression')+
  theme(legend.position='none',axis.title.x=element_blank())
# ggsave(dir('SVGs & PNGs/Blood_qpcr.svg'),width=2,height=2.5)
dunnettTest(bld_qpcr$exp_chp,as.factor(bld_qpcr$cnd))
