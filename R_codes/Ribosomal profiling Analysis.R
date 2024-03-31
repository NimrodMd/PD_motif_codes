library(edgeR); library(ggplot2); library(ggfortify); library('biomaRt'); library(matrixStats); library(reshape); library(stringr) ; library(M3C)

reform<-function(x){
  rownames(x)<-unlist(lapply(x[,1], function(y) strsplit(as.character(y),'\\.')[[1]][1])) ; x<-x[,2:ncol(x)]}
dir<-function(x) paste0(c(getwd(),x),collapse='/') ; dir('')
ndf<-function(x) as.data.frame(t(data.frame(row.names=x)))

## Read files #####
# Read tRF counts table:
trfs<-read.csv(dir('tRNA_Exclusive_Combined_data.csv')) ; rownames(trfs)<-trfs$X ; trfs<-trfs[,-1]
tcpm<-as.data.frame(cpm(trfs)) #create CPM counts table

# Read and process tRFs metadata:
meta_trfs<-read.csv(dir('tRF_meta.csv')); rownames(meta_trfs)<-meta_trfs$X; meta_trfs<-meta_trfs[,-1]; meta_trfs<-meta_trfs[rownames(trfs),]
meta_trfs$motifs<-F
for(i in 1:nrow(meta_trfs)){
  meta_trfs$motifs[i]<-sum(str_count(as.character(meta_trfs$tRF.sequence[i]),as.character('[AG]GTTC[GA]A')))}
meta_trfs$motifs<-factor(meta_trfs$motifs,labels=c('No RGTTCRA','RGTTCRA motif'))
meta_trfs$nuclear<-factor(unlist(lapply(meta_trfs$details,function(x) sum(str_count(as.character(as.character(x)),as.character('trnaMT')))!=0)),
                          labels=c('Nuclear','Mitochondrial'))

# Read Metadata
cold<-read.csv(dir('ColData_new.csv')) ; colnames(cold)<-gsub('Run','id',colnames(cold)) ; cold<-subset(cold,cold$id %in% colnames(trfs))
cold$molecule_subtype<-droplevels(cold$molecule_subtype)
cold$group<-factor(cold$group,levels=c('Resting','Depol 0 hours','Depol 2 hours'))

tmp<-data.frame(trf=rownames(tcpm),small=rowSums(tcpm[,as.character(subset(cold$id,cold$molecule_subtype=='small RNA'))]),
                ribo=rowSums(tcpm[,as.character(subset(cold$id,cold$molecule_subtype=='ribosome-protected mRNA'))]))
tmp<-melt(tmp,id.vars='trf') ; tmp<-merge(tmp,meta_trfs,by.x='trf',by.y=colnames(meta_trfs)[1])
ggplot(tmp,aes(len,log10(value),col=variable))+theme_classic()+geom_point(position=position_jitterdodge(jitter.width=0.3,jitter.height=0))+
  facet_wrap(~nuclear,scales='free',ncol=1)


## Analysis
# PD-tRFs boxplots:
anl<-data.frame(
  id=colnames(trfs),
  mt=colSums(trfs[rownames(trfs) %in% subset(meta_trfs$MINTbase.Unique.ID,meta_trfs$nuclear=='Mitochondrial'),]),
  mtf=colSums(trfs[rownames(trfs) %in% subset(meta_trfs$MINTbase.Unique.ID,meta_trfs$motifs=='RGTTCRA motif'),]),
  all=colSums(trfs))
anl$score<-anl$mtf/(anl$mt+1) ; anl$rtMtf<-100*anl$mtf/anl$all ; anl$rtMT<-100*anl$mt/anl$all ; anl$othr<-100*(anl$all-anl$mt-anl$mtf)/anl$all ; anl<-merge(anl,cold,by='id')
anl$molecule_subtype<-factor(anl$molecule_subtype,levels=c("small RNA","ribosome-protected mRNA"))

ggplot(anl,aes(group,rtMtf,col=group))+theme_classic()+facet_wrap(~molecule_subtype,scales='free')+
  geom_boxplot(outlier.shape=NA)+geom_point(position=position_jitter(height=0))+
  scale_color_manual(values=c('orange','pink3','red'))
# ggsave(dir('Depol_Boxplot_rtMtf.svg'),width=4.5,height=3)
TukeyHSD(aov(anl$rtMtf~anl$group*anl$molecule_subtype))

# Barplots:
p2Tab<-data.frame(
  id=colnames(tcpm),
  motif=colSums(tcpm[subset(meta_trfs$MINTbase.Unique.ID,meta_trfs$motifs=='RGTTCRA motif'),]),
  mt=colSums(tcpm[subset(meta_trfs$MINTbase.Unique.ID,meta_trfs$nuclear=='Mitochondrial'),]),
  noMt=colSums(tcpm[subset(meta_trfs$MINTbase.Unique.ID,meta_trfs$nuclear!='Mitochondrial'),]),
  noMotif=colSums(tcpm[subset(meta_trfs$MINTbase.Unique.ID,meta_trfs$motifs!='RGTTCRA motif'),]))
p2Tab$ratio<-p2Tab$motif/p2Tab$noMotif ; p2Tab$ratio2<-p2Tab$mt/p2Tab$noMt ; p2Tab<-merge(p2Tab,cold,by='id')
p2Tab$molecule_subtype<-factor(p2Tab$molecule_subtype,levels=c('small RNA','ribosome-protected mRNA'),labels=c('Small RNA','Ribosome profiling'))

tmp<-melt(p2Tab[,c("id","motif","noMotif","group","molecule_subtype")],id.vars=c("id","group","molecule_subtype"),variable_name='type')
tmp$type<-factor(tmp$type,levels=c('noMotif','motif'),labels=c('No RGTTCRA','RGTTCRA-tRFs'))
tmp$group<-factor(tmp$group,labels=c('Resting','Dep','2h pDP'))
tmp$mn<-unlist(lapply(1:nrow(tmp),function(x) 
  mean(subset(tmp$value,tmp$type=='RGTTCRA-tRFs' & tmp$molecule_subtype==tmp$molecule_subtype[x] & tmp$group==tmp$group[x]))/1000000))
tmp<-tmp[order(interaction(tmp$group,tmp$molecule_subtype,tmp$type)),] ; tmp$id2<-rep(c(1,2,3,4),12)

ggplot(tmp,aes((id),value,fill=type))+theme_classic()+facet_wrap(molecule_subtype~group,nrow=1,scales='free_x')+
  geom_bar(stat="identity",position='fill')+ylab('% of tRFs')+stat_summary(data=tmp,aes(id2,mn),fun='mean',geom='line',col='black',lwd=1.5)+
  scale_fill_manual(values=c('green3','purple3'))+theme(legend.position='bottom',legend.title=element_blank(),
                                                        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
# ggsave(dir('riboProfiling_bars.svg'),device='svg',width=10,height=5)  

tmp<-melt(p2Tab[,c("id","mt","noMt","group","molecule_subtype")],id.vars=c("id","group","molecule_subtype"),variable_name='type')
tmp$type<-factor(tmp$type,levels=c('noMt','mt'),labels=c('NtRFs','MT-tRFs'))
tmp$group<-factor(tmp$group,labels=c('Resting','Dep','2h pDP'))
tmp$mn<-unlist(lapply(1:nrow(tmp),function(x) 
  mean(subset(tmp$value,tmp$type=='MT-tRFs' & tmp$molecule_subtype==tmp$molecule_subtype[x] & tmp$group==tmp$group[x]))/1000000))
tmp<-tmp[order(interaction(tmp$group,tmp$molecule_subtype,tmp$type)),] ; tmp$id2<-rep(c(1,2,3,4),12)

ggplot(tmp,aes((id),value,fill=type))+theme_classic()+facet_wrap(molecule_subtype~group,nrow=1,scales='free_x')+
  geom_bar(stat="identity",position='fill')+ylab('% of tRFs')+stat_summary(data=tmp,aes(id2,mn),fun='mean',geom='line',col='black',lwd=1.5)+
  scale_fill_manual(values=c('green3','purple3'))+theme(legend.position='bottom',legend.title=element_blank(),
                                                        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
# ggsave(dir('riboProfiling_bars_MT-tRFs.svg'),device='svg',width=10,height=5)  
#p<0.0032 for small RNA 2hpDP vs Dep

ggplot(p2Tab)+theme_classic()+facet_wrap(~molecule_subtype,scales='fixed')+
  geom_point(aes(group,motif),position=position_jitter(width=0.1,height=0),col='red')+
  geom_smooth(mapping=aes(as.numeric(group),motif),method='loess',formula='y~x',col='darkred')+
  geom_point(aes(group,noMotif),position=position_jitter(width=0.1,height=0),col='blue')+
  geom_smooth(mapping=aes(as.numeric(group),noMotif),method='loess',formula='y~x',col='darkblue')+
  geom_point(aes(group,noMotif+motif),position=position_jitter(width=0.1,height=0),col='grey2')+
  geom_smooth(mapping=aes(as.numeric(group),noMotif+motif),method='loess',formula='y~x',col='black')+
  ylab('Normalized counts (CPM)')+theme(axis.title.x=element_blank())
# ggsave(dir('Motif_DePol_CytoVsRibo.svg'),width=10,height=5)
TukeyHSD(aov(p2Tab$motif~p2Tab$molecule_subtype*p2Tab$group))


# Only 3'- and i-tRFs:
pB2Tab<-data.frame(
  id=colnames(tcpm),
  motif=colSums(tcpm[subset(meta_trfs$MINTbase.Unique.ID,meta_trfs$motifs=='RGTTCRA motif' & 
                              meta_trfs$nuclear=='Nuclear' & meta_trfs$tRF.type.s. %in% c("3'-tRF","i-tRF")),]),
  noMotif=colSums(tcpm[subset(meta_trfs$MINTbase.Unique.ID,meta_trfs$motifs!='RGTTCRA motif' & 
                              meta_trfs$nuclear=='Nuclear' & meta_trfs$tRF.type.s. %in% c("3'-tRF","i-tRF")),]))
pB2Tab$ratio<-pB2Tab$motif/pB2Tab$noMotif ; pB2Tab<-merge(pB2Tab,cold,by='id')

ggplot(pB2Tab)+theme_classic()+facet_wrap(~molecule_subtype,scales='fixed')+
  geom_point(aes(group,motif),position=position_jitter(width=0.1,height=0),col='red')+
  geom_smooth(mapping=aes(as.numeric(group),motif),method='loess',formula='y~x',col='darkred')+
  geom_point(aes(group,noMotif),position=position_jitter(width=0.1,height=0),col='blue')+
  geom_smooth(mapping=aes(as.numeric(group),noMotif),method='loess',formula='y~x',col='darkblue')+
  ylab('Normalized counts (CPM)')+theme(axis.title.x=element_blank())
# ggsave(dir("Motif_DePol_CytoVsRibo_Supplemental_Only.3.i'tRFs.svg"),width=10,height=5)
TukeyHSD(aov(pB2Tab$noMotif~pB2Tab$molecule_subtype*pB2Tab$group))

