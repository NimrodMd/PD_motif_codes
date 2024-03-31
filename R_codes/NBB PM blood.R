library(ggplot2) ; library(matrixStats) ; library(stringr) ; library(edgeR) ; library(ggh4x)

ndf<-function(x) as.data.frame(t(data.frame(row.names=x))) ; Ttmp<-function(x) as.data.frame(t(x))
Ftmp<-function(df,col='exp',id='id') {
  assign('tmp',as.data.frame(t(df))); colnames(tmp)<-col; tmp[,id]<-rownames(tmp) ; return(tmp)}
se<-function(x) sd(x)/sqrt(length(x)) ; errN<-function(x) mean(x)-se(x) ; errP<-function(x) mean(x)+se(x)

reform<-function(x){
  rownames(x)<-unlist(lapply(x[,1], function(y) strsplit(as.character(y),'\\.')[[1]][1])) ; x<-x[,2:ncol(x)]}

dir_in<-function(x) paste0(c(getwd(),'Excels',x),collapse='/') ; dir_in('')
dir_out<-function(x) paste0(c(getwd(),'SVGs & PNGs',x),collapse='/') ; dir_out('')

## Read NBB files #####
sn_PD<-read.csv('C:/Users/nmisr/OneDrive/lab/5. DataSets/SN_rawData.csv') ; sn_PD<-reform(sn_PD)
col_PD<-read.csv('C:/Users/nmisr/OneDrive/lab/5. DataSets/BBB_ColData.csv')
tmp<-as.data.frame(t(sn_PD['ENSG00000180176',])) ; colnames(tmp)<-'exp' ; tmp$Sample<-rownames(tmp)
tmp<-merge(tmp,col_PD,by='Sample') ; TH<-subset(tmp,tmp$Condition=='Parkinson' & tmp$exp>5)

coldata<-read.csv(dir_in('NBB PM blood ColData.csv')) ; colnames(coldata)<-gsub('Sample.Name','id',colnames(coldata))
trfs<-read.csv(dir_in('NBB PM blood CTS.csv')) ; rownames(trfs)<-trfs$X ; trfs<-trfs[,-1]
colnames(trfs)<-unlist(lapply(colnames(trfs),function(x) paste0(strsplit(as.character(x),'_')[[1]][1:2],collapse='/')))
nrm_trfs<-as.data.frame(cpm(trfs))

ptients<-ndf(c('Sex','Group','Age','LewyBodies','Duration','PMD'))
for(s in unique(coldata$Sex)){for(g in unique(coldata$full_Diagnosis)){
  tmp<-subset(coldata,coldata$Sex==s & coldata$full_Diagnosis==g)
  ptients<-rbind(ptients,data.frame('Sex'=s,'Group'=g,'Age'=mean(as.numeric(tmp$Age)),'LewyBodies'=mean(as.numeric(tmp$braaklb),na.rm=T),
                                    'Duration'=mean(as.numeric(tmp$PD_duration),na.rm=T),'PMD'=mean(as.numeric(tmp$PMD_blood),na.rm=T)))
}} #; write.csv(ptients,dir_out('NBB blood patients.csv'))

durs<-read.csv(dir_in('PPMI and NBB durations.csv'))
ggplot(durs,aes(Dur))+theme_classic()+facet_wrap(~group,scales='free_y')+xlab('PD duration')+
  geom_histogram(aes(y=after_stat(count/sum(count)),fill=group),binwidth=2)+geom_density(col='black')+
  theme(axis.title.y=element_blank(),axis.ticks=element_blank(),axis.text.y=element_blank(),legend.position='none')
# ggsave(dir_out('PPMI and NBB durations.svg'),width=4.5,height=1.7)
min(subset(durs$Dur,durs$group=='NBB Blood')) ; max(subset(durs$Dur,durs$group=='NBB Blood'))
min(subset(durs$Dur,durs$group=='NBB SN')) ; max(subset(durs$Dur,durs$group=='NBB SN'))
min(subset(durs$Dur,durs$group=='PPMI Blood')) ; max(subset(durs$Dur,durs$group=='PPMI Blood'))

# Read and process tRFs metadata:
metadata<-read.csv(dir_in('NBB PM blood tRF meta.csv')) ; metadata<-metadata[,-1]
colnames(metadata)<-gsub('MINTbase.Unique.ID','trf',colnames(metadata)) ; metadata$motifs<-F
for(i in 1:nrow(metadata)){
  metadata$motifs[i]<-sum(str_count(as.character(metadata$tRF.sequence[i]),as.character('[AG]GTTC[GA]A')))}
metadata$motifs<-factor(metadata$motifs,labels=c('No RGTTCRA','RGTTCRA motif'))
metadata$nuclear<-factor(unlist(lapply(metadata$Sequence.locations.in.tRNA.space..comma.deliminated.,function(x) 
  sum(str_count(as.character(as.character(x)),as.character('trnaMT')))!=0)),labels=c('Nuclear','Mitochondrial'))

tmp0<-as.data.frame(t(trfs)) ; tmp0$id<-rownames(tmp0) ; tmp0<-merge(coldata[,c("id","Sex","full_Diagnosis")],tmp0)
onlyMotif<-tmp0[,colnames(tmp0) %in% c("id","Sex","full_Diagnosis",subset(metadata$trf,metadata$motifs=='RGTTCRA motif'))]
noMotif<-tmp0[,colnames(tmp0) %in% c("id","Sex","full_Diagnosis",subset(metadata$trf,metadata$motifs!='RGTTCRA motif'))]
# write.csv(tmp0,dir_in('All_tRFs.csv')) ; write.csv(onlyMotif,dir_in('onlyMotif.csv')) ; write.csv(noMotif,dir_in('noMotif.csv'))

motif<-data.frame(
  motif=colSums(nrm_trfs[rownames(nrm_trfs) %in% subset(metadata$trf,metadata$motifs=='RGTTCRA motif'),]),
  mt=colSums(nrm_trfs[rownames(nrm_trfs) %in% subset(metadata$trf,metadata$nuclear=='Mitochondrial' & metadata$motifs!='RGTTCRA motif'),]),
  mean=colMeans(nrm_trfs[rownames(nrm_trfs) %in% subset(metadata$trf,metadata$motifs=='RGTTCRA motif'),]),
  med=colMedians(as.matrix(nrm_trfs[rownames(nrm_trfs) %in% subset(metadata$trf,metadata$motifs=='RGTTCRA motif'),])),
  all=colSums(nrm_trfs),id=colnames(nrm_trfs))
motif$ratio<-100*motif$motif/motif$all ; motif$mt_pct<-100*motif$mt/motif$all ; motif<-merge(motif,coldata,by='id')
motif$group<-factor(as.numeric(motif$id %in% TH$BrainID)+as.numeric(motif$PD_y_n=='PD'),labels=c('Ctrl','PD','SWEDD'))
motif$rm<-!motif$braaklb %in% c(1,3)
motif$full_Diagnosis<-factor(motif$full_Diagnosis,levels=c('NDC','NDPD','DPD'))

# Create datafrem enabling showing outliers as separeted from the non-outlier data:
rpMotif<-motif ; rpMotif$ratio<-as.numeric(gsub(max(rpMotif$ratio),sort(as.numeric(rpMotif$ratio),decreasing=T)[2]+1,rpMotif$ratio))
rpMotif$mtf2<-as.numeric(gsub(max(rpMotif$motif),sort(as.numeric(rpMotif$motif),decreasing=T)[2]+10,rpMotif$motif))
rpMotif$lb2<-gsub(3,4,rpMotif$braaklb) ; rpMotif$lb2<-gsub(1,0,rpMotif$lb2) ;rpMotif[is.na(rpMotif$lb2),'lb2']<-0 ; 
rpMotif$lb2<-factor(rpMotif$lb2,levels=c('0','4','5','6'),labels=c('0-1','3-5','3-5','6'))

#
## Read PPMI data #####
dirin<-function(x) paste0('C:/Users/nmisr/OneDrive/lab/ADPD/1. PD SGs tRFs/Excels/',x)
dirout<-function(x) paste0('C:/Users/nmisr/OneDrive/lab/ADPD/1. PD SGs tRFs/SVGs & PNGs/',x)

PPMI_data<-read.csv(dirin('PPMI_data.2.csv'))
PPMI_data<-subset(PPMI_data,!(PPMI_data$Genetic_background=='LRRK2-/SNCA-/GBA-' & PPMI_data$surgery=='DBS'))
PPMI_data$surgery<-factor(PPMI_data$surgery,levels=c('No','DBS'),labels=c('No DBS','DBS'))
PPMI_data$bioGroup<-interaction(PPMI_data$surgery,PPMI_data$cond)

## Combine NBB and PPMI data:
allblood<-motif[,c("Sex","Age","PD_y_n","ratio","mt_pct")] ; colnames(allblood)<-c("Sex","Age","group","ratio","mt_pct")
allblood$group<-gsub('NDC','Ctrl',allblood$group) ; allblood$genetic<-'Idiopathic' ; allblood$phase<-'Advanced'
allblood$r2<-allblood$ratio/mean(subset(allblood$ratio,allblood$group=='Ctrl'))
allblood$m2<-allblood$mt_pct/mean(subset(allblood$mt_pct,allblood$group=='Ctrl'))

tmp<-PPMI_data[PPMI_data$surgery=='No DBS',c("Sex","Age_timepoint","cond","rt_motif",'rt_mt2',"Genetic_background")]
colnames(tmp)<-c("Sex","Age","group","ratio","mt_pct","genetic") ; tmp$group<-gsub('Control','Ctrl',tmp$group)
tmp$genetic<-as.character(factor(tmp$genetic!='LRRK2-/SNCA-/GBA-', labels=c('Idiopathic','Genetic')))
tmp$phase<-'Early' ; tmp$r2<-tmp$ratio/mean(subset(tmp$ratio,tmp$group=='Ctrl'))
tmp$m2<-tmp$mt_pct/mean(subset(tmp$mt_pct,tmp$group=='Ctrl'))

allblood<-rbind(allblood,tmp) ; allblood$Sex<-str_to_title(allblood$Sex)
allblood$phase<-factor(allblood$phase,levels=c('Early','Advanced'),labels=c('Alive','PM'))
allblood$genetic<-factor(allblood$genetic,levels=c('Idiopathic','Genetic'))

# Remove Prodromal patients:
f_allblood<-subset(allblood,allblood$group!='Prodromal')
f_allblood$ratio2<-as.numeric(gsub(max(f_allblood$ratio),sort(as.numeric(f_allblood$ratio),decreasing=T)[2]+2,f_allblood$ratio))
f_allblood$class<-interaction(f_allblood$genetic,f_allblood$phase)

# Calculate differences significance between groups:
pvTab<-ndf(c('group','p_mtf','p_mt','p_mtf.mt'))
for(c in unique(f_allblood$class)){
  tmp<-subset(f_allblood,f_allblood$class==c)
  tt1<-wilcox.test(tmp$ratio~tmp$group) ; tt2<-wilcox.test(tmp$mt_pct~tmp$group) 
  tt3<-wilcox.test(tmp$ratio/tmp$mt_pct~tmp$group)
  pvTab<-rbind(pvTab,data.frame('group'=c,p_mtf=tt1$p.value,p_mt=tt2$p.value,p_mtf.mt=tt3$p.value))
} ; pvTab$FDR_mtf<-p.adjust(pvTab$p_mtf,'fdr') ; pvTab$FDR_mt<-p.adjust(pvTab$p_mt,'fdr') ; pvTab$FDR_mtf.mt<-p.adjust(pvTab$p_mtf.mt,'fdr')

# Plot results:
ggplot(f_allblood,aes(genetic,ratio2,col=group,fill=group))+theme_classic()+
  facet_grid2(~phase,scales='free',space='free',independent='y')+
  geom_boxplot(outlier.shape=NA)+geom_point(position=position_jitterdodge(jitter.width=0.3,jitter.height=0))+
  ylab('% PD-tRFs')+theme(axis.title.x=element_blank(),legend.title=element_blank())+
  scale_color_manual(values=c('grey35','darkblue'))+scale_fill_manual(values=c('grey','lightblue'))
# ggsave(dir_out('NBB blood nrm ratio.svg'),width=4.5,height=5)

ggplot(f_allblood,aes(genetic,mt_pct,col=group,fill=group))+theme_classic()+
  facet_grid2(~phase,scales='free',space='free',independent='y')+
  geom_boxplot(outlier.shape=NA)+geom_point(position=position_jitterdodge(jitter.width=0.3,jitter.height=0))+
  ylab('% MT-tRFs')+theme(axis.title.x=element_blank(),legend.title=element_blank())+
  scale_color_manual(values=c('grey35','darkblue'))+scale_fill_manual(values=c('grey','lightblue'))
# ggsave(dir_out('NBB blood nrm MT.svg'),width=4.5,height=5)

ggplot(f_allblood,aes(genetic,ratio/mt_pct,col=group,fill=group))+theme_classic()+
  facet_grid(~phase,scales='free',space='free')+
  geom_boxplot(outlier.shape=NA)+geom_point(position=position_jitterdodge(jitter.width=0.3,jitter.height=0))+
  ylab('RGTTCRA-tRFs / MT-tRFs')+theme(axis.title.x=element_blank(),legend.title=element_blank())+
  scale_color_manual(values=c('grey35','darkblue'))+scale_fill_manual(values=c('grey','lightblue'))
# ggsave(dir_out('NBB blood nrm Motif.MT.svg'),width=4.5,height=5)