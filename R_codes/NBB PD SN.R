library(edgeR); library(ggplot2); library(ggfortify); library(matrixStats); library(reshape); library(stringr) ; library(gg3D)

reform<-function(x){
  rownames(x)<-unlist(lapply(x[,1], function(y) strsplit(as.character(y),'\\.')[[1]][1])) ; x<-x[,2:ncol(x)]}

dir<-function(x) paste0(c(getwd(),x),collapse='/') ; dir('')
ndf<-function(x) as.data.frame(t(data.frame(row.names=x)))

# Function to calculate distance of points from a line
dst<-function(i,j,a,b) (abs(i-(j-b)/a)*abs(j-a*i+b))/sqrt((i-(j-b)/a)^2+(j-a*i+b)^2)
clcDist<-function(x,y,cfs=coef(lm(y~x))){unlist(lapply(1:length(x),function(c) dst(x[c],y[c],cfs[2],cfs[1])))}


# Read files #####
# Read Meatadata:
cold<-read.csv(dir('ColData.csv')) 
cold$pmd2<-unlist(lapply(cold$pmd,function(x) as.numeric(strsplit(x,':')[[1]][1])+as.numeric(strsplit(x,':')[[1]][2])/60))

ggplot(cold,aes(PD.length,braaklb,col=to.remove))+theme_classic()+geom_point()
cor.test(cold$PD.length,cold$braaklb,method='pearson')

tmp<-cold[,c("id","age","pmd2")] ; tmp<-melt(tmp,id.vars='id')
ggplot(tmp,aes(value))+theme_classic()+geom_density()+facet_wrap(~variable,scales='free')
#ggsave(dir('duration_pmi.svg'),width=4.5,height=2.5)

# Read tRF counts table:
trfs<-read.csv(dir('tRNA_Exclusive_Combined_data.csv')) ; rownames(trfs)<-trfs$X ; trfs<-trfs[,-1]
colnames(trfs)<-unlist(lapply(colnames(trfs),function(x) paste0('S',paste(strsplit(x,'_')[[1]][2:3],collapse='/'))))
for(c in colnames(trfs)){colnames(trfs)<-gsub(c,subset(cold$id,cold$BrainID==c),colnames(trfs))}
tcpm<-as.data.frame(cpm(trfs)) ; tcpm<-tcpm[rowMedians(as.matrix(tcpm))>0,]

# Read and process tRFs Metadata:
meta_trfs<-read.csv(dir('tRF_meta.csv')); rownames(meta_trfs)<-meta_trfs$X; meta_trfs<-meta_trfs[,-1]; meta_trfs<-meta_trfs[rownames(trfs),]
meta_trfs$len<-unlist(lapply(meta_trfs$trf, function(x) strsplit(as.character(x),'-')[[1]][2]))
meta_trfs$motifs<-F
for(i in 1:nrow(meta_trfs)){
  meta_trfs$motifs[i]<-sum(str_count(as.character(meta_trfs$tRF.sequence[i]),as.character('[AG]GTTC[GA]A')))}
meta_trfs$motifs<-factor(meta_trfs$motifs,labels=c('No RGTTCRA','RGTTCRA motif'))

# Read long dopaminergic-transcripts expression:
long<-read.csv(dir('SN_longRNA.csv')) ; rownames(long)<-long$X ; long<-long[,-1]

#
# Analysis #####

# Check no differenence in Nuclear noMotif tRFs
nuc<-data.frame(
  id=colnames(tcpm),
  Mean=colMeans(tcpm[rownames(tcpm) %in% subset(meta_trfs$trf,meta_trfs$motifs!='RGTTCRA motif' & meta_trfs$nuclear=='Nuclear'),]),
  Sum=colSums(tcpm[rownames(tcpm) %in% subset(meta_trfs$trf,meta_trfs$motifs!='RGTTCRA motif' & meta_trfs$nuclear=='Nuclear'),]),
  Med=colMedians(as.matrix(tcpm[rownames(tcpm) %in% subset(meta_trfs$trf,meta_trfs$motifs!='RGTTCRA motif' & meta_trfs$nuclear=='Nuclear'),])),
  All=colSums(tcpm[rownames(tcpm) %in% subset(meta_trfs$trf,meta_trfs$nuclear=='Nuclear'),]))
nuc$ratio<-100*nuc$Sum/nuc$All ; nuc<-merge(nuc,cold,by='id')

# Motif tRFs
pd<-data.frame(
  id=colnames(trfs),
  Mean=colMeans(trfs[rownames(trfs) %in% subset(meta_trfs$trf,meta_trfs$motifs=='RGTTCRA motif'),]),
  Sum=colSums(trfs[rownames(trfs) %in% subset(meta_trfs$trf,meta_trfs$motifs=='RGTTCRA motif'),]),
  MT=colSums(trfs[rownames(trfs) %in% subset(meta_trfs$trf,meta_trfs$nuclear=='Mitochondrial'),]),
  Med=colMedians(as.matrix(trfs[rownames(trfs) %in% subset(meta_trfs$trf,meta_trfs$motifs=='RGTTCRA motif'),])),
  All=colSums(trfs),Nuc=colSums(trfs[rownames(trfs) %in% subset(meta_trfs$trf,meta_trfs$nuclear=='Nuclear'),]))
pd$ratio<-100*pd$Sum/pd$All ; pd$mt_pct<-100*pd$MT/pd$All ; pd$nucRatio<-100*pd$Sum/pd$Nuc ; pd<-merge(pd,cold,by='id')

# Find SWEDD-like patients based on DA-gene expression
tmp<-as.data.frame(t(long)) ; tmp$id<-rownames(tmp) ; pd<-merge(pd,tmp,by='id')
pd$dis<-clcDist(pd$PD.length,pd$ratio,coef(lm(pd$ratio[pd$to.remove=='no']~pd$PD.length[pd$to.remove=='no'])))
pd$DA<-log10((pd$TH+1)*(pd$DAT+1)) ; tmp<-subset(pd,pd$to.remove=='no') 
pd$grp<-factor(pd$to.remove,levels=c('no','high TH'),labels=c('PD','SWEDD-like'))

pd$cng_all<-100*pd$All/pd[pd$PD.length==min(pd$PD.length),'All']
pd$cng_mtf<-100*pd$Sum/pd[pd$PD.length==min(pd$PD.length),'Sum']
pd$cng_mt<-100*pd$MT/pd[pd$PD.length==min(pd$PD.length),'MT']
pd$highLB<-factor(pd$braaklb>5,labels=c('LB 4/5','LB 6')) ; pd$mtf.mt<-pd$Sum/(1+pd$MT)

## Plot results
ggplot(pd)+theme_classic()+geom_smooth(aes(PD.length,cng_all),method='lm',formula='y~x',col='black',se=F)+geom_point(aes(PD.length,cng_all,shape=grp),col='black')+
  geom_smooth(aes(PD.length,cng_mtf),method='lm',formula='y~x',col='purple',se=F)+geom_point(aes(PD.length,cng_mtf,shape=grp),col='purple')+
  geom_smooth(aes(PD.length,cng_mt),method='lm',formula='y~x',col='green',se=F)+geom_point(aes(PD.length,cng_mt,shape=grp),col='green')+facet_grid(~grp,space='free')
# ggsave(dir('changes_with_time_pct.svg'))

ggplot(pd,aes(PD.length,ratio,col=to.remove))+theme_classic()+geom_smooth(method='lm',formula='y~x')+geom_point()+
  theme(legend.position='none')+ylab('% RGTTCRA-tRFs')+xlab('PD duration (years)')+scale_color_manual(values=c('orange3','purple'))
# ggsave(dir('sSNPD_MotifDuration.svg'),device='svg',width=4,height=4.75)

tmp<-subset(pd,pd$to.remove=='no')
cor.test(tmp$ratio,tmp$PD.length) # r=0.82, p<0.024
cor.test(tmp$ratio,tmp$pmd2) # r=-0.7, p<0.0745

ggplot(pd,aes(PD.length,mt_pct,col=to.remove))+theme_classic()+geom_smooth(method='lm',formula='y~x')+geom_point()+
  scale_color_manual(values=c('orange3','purple'))+theme(legend.position='none')+ylab('% RGTTCRA-tRFs')+xlab('PD duration (years)')
# ggsave(dir('sSNPD_MT_Duration.svg'),device='svg',width=4,height=4.75)
cor.test(tmp$mt_pct,tmp$PD.length) # r=0.85, p<0.015
cor.test(tmp$mt_pct,tmp$pmd2) # r=-0.53, p<0.22

ggplot(pd,aes(age,ratio,col=to.remove))+theme_classic()+geom_smooth(method='lm',formula='y~x')+geom_point()+
  scale_color_manual(values=c('orange3','purple'))+theme(legend.position='none')+ylab('% RGTTCRA-tRFs')+xlab("Patient's age (years)")
# ggsave(dir('sSNPD_MotifAge.svg'),device='svg',width=4,height=4.75)
cor.test(tmp$age,tmp$ratio) # r=-1.5 p<0.744

ggplot(pd,aes(age,mt_pct,col=to.remove))+theme_classic()+geom_smooth(method='lm',formula='y~x')+geom_point()+
  scale_color_manual(values=c('orange3','purple'))+theme(legend.position='none')+ylab('% MT-tRFs')+xlab("Patient's age (years)")
# ggsave(dir('sSNPD_MTAge.svg'),device='svg',width=4,height=4.75)
cor.test(tmp$age,tmp$mt_pct) # r=-0.24 p<0.6

tmp<-pd[,c("to.remove","TH","DAT","ratio","PD.length","dis","id")]; tmp<-melt(tmp,id.vars=c("to.remove","PD.length","ratio","dis","id"))
tmp$to.remove<-gsub('no','PD',tmp$to.remove) ; tmp$to.remove<-gsub('high TH','SWEDD-like',tmp$to.remove)
ggplot(tmp,aes((dis),value,col=to.remove))+theme_classic()+facet_wrap(~variable)+geom_point(lwd=3)+
  scale_color_manual(values=c('purple','orange3'))+ylab('Counts per million')+xlab('Distance from fit')+
  theme(legend.title=element_blank(),legend.position='bottom')
# ggsave(dir('sSNPD_DAT.TH.svg'),device='svg',width=4,height=4.75)
         
#
          # DEseq #####
# Remove expression-outlier tRFs (Mean higher than the 85th quantile):
toRemove<-c() ; for(g in rownames(tcpm)){
  if(quantile(as.numeric(tcpm[g,]),prob=0.85)<mean(as.numeric(tcpm[g,]))){toRemove<-append(toRemove,g)}}

# Run DESeq analysis based on PD duration, Braak LB score and DA gene expression:
for(f in c('PD.length','braaklb','DAT')){
  cts<-trfs ; cold1<-pd ; colnames(cold1)<-gsub(f,'fct',colnames(cold1)) ; cold1<-subset(cold1,cold1$id %in% colnames(cts))
  cts<-cts[,as.character(cold1$id)] ; cts<-cts[,order(cold1$id)] ; cold1<-cold1[order(cold1$id),] ; nrow(cold1)==sum(cold1$id==colnames(cts))
  y <- DGEList(counts=cts,group=cold1$fct) # change the condition to the name of your condition column
  keep <- filterByExpr(y) ; y <- y[keep, , keep.lib.sizes=FALSE] ; y$samples$lib.size <- colSums(y$counts) ; y <- calcNormFactors(y)
  cts1<-as.data.frame(cpm(y,log=F)) ; dsgn<-model.matrix(~fct,data=cold1) ; y<-estimateDisp(y,dsgn,robust = T) ; head(dsgn,8)
  fit <- glmQLFit(y, dsgn, robust = T) ; lrt <- glmLRT(fit,coef = 2) 
  
  sgGens<-as.data.frame(topTags(lrt,adjust.method = 'fdr',n = nrow(cts1))) ; sgGens$transcript<-rownames(sgGens)
  sgGens$isSg<-factor(as.numeric(sgGens$FDR<0.051),labels = c('N.S.','Significant')) ; sgGens<-subset(sgGens,!sgGens$transcript %in% toRemove)
  sgGens<-merge(sgGens,meta_trfs,by.x='transcript',by.y=colnames(meta_trfs)[1])
  assign(f,sgGens)
}
sgCsf<-read.csv('Y:/nimrod.madrer/ADPD/CSF_ADPD/phs000727_ADPD_Short/Analysis/sgGenes.csv')
tmp<-subset(sgCsf,sgCsf$FDR<0.051); tmp<-subset(braaklb,braaklb$transcript %in% tmp$transcript & braaklb$motifs=='RGTTCRA motif')

tmp<-subset(braaklb,!(braaklb$nuclear=='Mitochondrial' & braaklb$motifs=='RGTTCRA motif'))
ggplot(tmp,aes(log2(2.72^logFC),-log10(FDR),col=as.numeric(len),alpha=isSg))+theme_classic()+facet_wrap(motifs~nuclear)+geom_point()+
  geom_vline(xintercept=c(-1,1),linetype='dashed',col='grey')+geom_hline(yintercept=-log10(0.05),linetype='dashed',col='grey')+
  scale_color_gradient2(low='green',mid='lightblue',high='purple',midpoint=35,name='tRF length')+xlab('Log2(FC)')
# ggsave(dir('sSNPD_DE_MotifLB.svg'),device='svg',width=5.5,height=4.75)
