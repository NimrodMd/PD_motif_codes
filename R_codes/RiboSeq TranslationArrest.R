library(edgeR); library(ggplot2); library(ggfortify); library('biomaRt'); library(matrixStats); library(reshape); library(stringr) ; library(M3C)
library(DescTools)

reform<-function(x){rownames(x)<-unlist(lapply(x[,1], function(y) strsplit(as.character(y),'\\.')[[1]][1])) ; x<-x[,2:ncol(x)]}
dir<-function(x) paste0(c(getwd(),x),collapse='/') ; dir('')
ndf<-function(x) as.data.frame(t(data.frame(row.names=x)))

# Read files #####
cold<-read.csv(dir('ColData_new.csv')) ; colnames(cold)<-gsub('Run','id',colnames(cold))
cold$treatment_medium<-gsub(' \\+ 10% dialyzed FBS','',cold$treatment_medium)
cold$treatment_medium<-gsub('DMEM ','',cold$treatment_medium) ; cold$treatment_medium<-gsub('DMEM','Ctrl',cold$treatment_medium)
cold1<-subset(cold,cold$Genotype %in% c('AAVS1-iCAG:hrGFP','wild type') & cold$treatment_medium %in% c('Ctrl','-Arginine','-Leucine'))
cold1$treatment_medium<-factor(cold1$treatment_medium, levels=c('Ctrl','-Arginine','-Leucine'), labels=c('Ctrl','-Arg','-Leu'))


trfs<-read.csv(dir('tRNA_Exclusive_Combined_data.csv')) ; rownames(trfs)<-trfs$X ; trfs<-trfs[,-1]

meta_trfs<-read.csv(dir('tRF_meta.csv')); rownames(meta_trfs)<-meta_trfs$X; meta_trfs<-meta_trfs[,-1]; meta_trfs<-meta_trfs[rownames(trfs),]
meta_trfs$G4s<-as.factor(unlist(lapply(meta_trfs$tRF.sequence, function(x) as.numeric(substr(as.character(x),1,4)=='GGGG'))))
meta_trfs$motifs<-F
for(i in 1:nrow(meta_trfs)){
  meta_trfs$motifs[i]<-sum(str_count(as.character(meta_trfs$tRF.sequence[i]),as.character('[AG]GTTC[GA]A')))}
meta_trfs$motifs<-factor(meta_trfs$motifs,labels=c('No RGTTCRA','RGTTCRA motif'))
meta_trfs$nuclear<-factor(unlist(lapply(meta_trfs$details,function(x) sum(str_count(as.character(as.character(x)),as.character('trnaMT')))!=0)),
                          labels=c('Nuclear','Mitochondrial'))

tcpm<-as.data.frame(cpm(trfs)) ; tcpm<-tcpm[rowMedians(as.matrix(tcpm))>0,]

pTab<-data.frame(id=colnames(trfs),
                 motif=colSums(trfs[subset(meta_trfs$trf,meta_trfs$motifs=='RGTTCRA motif'),]),
                 noMotif=colSums(trfs[subset(meta_trfs$trf,meta_trfs$motifs!='RGTTCRA motif'),]))
p2Tab<-pTab ; p2Tab$ratio<-p2Tab$motif/(p2Tab$noMotif+p2Tab$motif) ; p2Tab<-merge(p2Tab,cold,by='id')
pTab<-melt(pTab,id.vars='id') ; pTab<-merge(pTab,cold,by='id')
p2Tab$sum<-p2Tab$motif+p2Tab$noMotif
p2Tab$Cell_Line<-factor(p2Tab$Cell_Line,levels=c('HEK293T','HCT116','HeLa'))
p2Tab$id<-factor(p2Tab$id,levels=p2Tab$id[order(p2Tab$Cell_Line)])

tmp1<-subset(p2Tab,p2Tab$Genotype %in% c('wild type','AAVS1-iCAG:hrGFP'))
ggplot(tmp1,aes(x=treatment_medium))+theme_classic()+facet_grid(~treatment_length,scales='free_x',space='free_x')+
  theme(axis.title.x=element_blank(),legend.position='bottom',legend.title=element_blank())+
  geom_hline(yintercept=p[2]*100,lwd=1,col='grey',linetype='dashed')+
  geom_point(aes(y=ratio*100,col=Cell_Line,group=id),stat='identity',position=position_dodge(width=.8),lwd=3)+
  # geom_bar(aes(y=sum/5,group=id),stat='identity',position=position_dodge(width=.8),width=0.7,fill='grey')+
  scale_y_continuous(name = "% of RGTTCRA motif tRFs", sec.axis = sec_axis(~.*5, name="# of reads"))+
  scale_color_manual(values=c('#3C5488B2','#00A087B2','#DC0000B2'))
# ggsave(dir('RiboSeq Motif enrichment.svg'),device='svg')


pvTab<-ndf(c('id','fc','p','pChance'))
for(s in subset(cold$id,cold$Genotype %in% c('wild type','AAVS1-iCAG:hrGFP'))){
  tmp<-p2Tab[p2Tab$id==s,c("noMotif","motif")] ; csq<-chisq.test(tmp,p=p) ; chnce<-chisq.test(tmp)
  pvTab<-rbind(pvTab,data.frame('id'=s,'fc'=csq$observed[2]/csq$expected[2],'p'=csq$p.value,'pChance'=chnce$p.value))
} ; pvTab$fdr<-p.adjust(pvTab$p,'fdr'); pvTab$fdrChance<-p.adjust(pvTab$pChance,'fdr') ; pvTab<-merge(pvTab,cold,by='id')
# write.csv(pvTab[,1:10],dir('Enrichment pValues.csv'))

