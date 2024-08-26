# Load Data and Libraries #####

# Load necessary libraries
library(edgeR)                   # For differential expression analysis
library(ggplot2)                 # For creating visualizations
library(ggfortify)               # For visualizations
library('biomaRt')               # For accessing biological data
library(matrixStats)             # For statistical functions
library(reshape)                 # For data manipulation
library(stringr)                 # For string manipulation
library(M3C)                     # For clustering analysis
library(DescTools)               # For descriptive statistics
library(MatchIt)                 # For propensity score matching
library(caret)                   # For machine learning
library(factoextra)              # For clustering analysis
library(precrec)                 # For precision-recall curves
library(ggrepel)                 # For text labels in plots
library(pROC)                    # For ROC analysis
library(yardstick)               # For model evaluation
library(MLeval)                  # For model evaluation

reform<-function(x){
  rownames(x)<-unlist(lapply(x[,1], function(y) strsplit(as.character(y),'\\.')[[1]][1])) ; x<-x[,2:ncol(x)]}

dir<-function(x) paste0('Y:/nimrod.madrer/ADPD/CSF_ADPD/phs000727_ADPD_Short/Analysis/',x)
ndf<-function(x) as.data.frame(t(data.frame(row.names=x)))

# CSF data:
ColData<-read.csv(dir('MetaData.csv')) # Read Metadata
ColData<-ColData[-18,]
# patient 02_01 has two occurences:
# - one with ApoE 3;4 and higher BRAAK (row 18)
# - one with ApoE 3;3 and lower BRAAK (row 17)

ColData$group<-factor(ColData$group,levels=c('Ctrl','AD','PD'))

# Data preprocessing:
ColData$BraakScore<-gsub('IV',4,ColData$BraakScore) ; ColData$BraakScore<-gsub('VI',6,ColData$BraakScore)
ColData$BraakScore<-gsub('III',3,ColData$BraakScore) ; ColData$BraakScore<-gsub('II',2,ColData$BraakScore)
ColData$BraakScore<-gsub('I',1,ColData$BraakScore) ; ColData$BraakScore<-gsub('V',5,ColData$BraakScore) 
ColData$PlaqueScore<-as.numeric(as.character(factor(ColData$PlaqueDensity,labels=c(4,3,2,1))))*ColData$PlaqueTotal
ColData$NIA_ReaganScore<-factor(ColData$NIA_ReaganScore,levels=c("criteria not met","not AD","low","intermediate","high"),labels=c("Not AD","Not AD","Low","Intermediate","High"))
ColData$LewyBody_Stage<-factor(ColData$LewyBody_Stage,levels=c("No Lewy bodies","","LB pathology unspecified or not further assessed","Brainstem predominant type","Limbic type","Neocortical type"),
                               labels=c("No Lewy bodies","No Lewy bodies","No Lewy bodies","Brainstem","Limbic","Neocortical"))
ColData$sn_depigmentation<-factor(ColData$sn_depigmentation,levels=c("none","","mild","moderate","severe"),labels=c("None","Unknown","Mild","Moderate","Severe"))
ColData$AD_score<-(as.numeric(ColData$BraakScore)*ColData$TangleTotal*ColData$PlaqueScore)/400

# Identify ADPD patients
ADPD<-as.character(subset(ColData$ID,(ColData$group=='AD' & ColData$sn_depigmentation=='Severe') | 
                            (ColData$group=='PD' & ColData$NIA_ReaganScore %in% c('Intermediate','High') & ColData$AD_score>=2)))
ColData$ADPD<-ColData$ID %in% ADPD

# Remove ADPD and other outlier patients
ToRemove1<-c('SRR1568665','SRR1568516','SRR1568612','SRR1568554','SRR1568520','SRR1568402','SRR1568626','SRR1568616','SRR1568556',
             'SRR1568527','SRR1568441','SRR1568593','SRR1568452','SRR1568726','SRR1568546','SRR1568552') ; p<-ToRemove1
ToRemove1<-append(ToRemove1,as.character(subset(ColData$ID,ColData$group=='Ctrl' & (ColData$AD_score>5 | ColData$sn_depigmentation=='Moderate'))))
ColData$keep_trf<-!ColData$ID %in% ToRemove1 

t_cold<-subset(ColData,ColData$keep_trf) ; cold1<-subset(t_cold,!t_cold$ADPD & t_cold$age>60)

# Generate summary statistics for the subjects:
sbjs<-ndf(c('Group','Sex','N','Avg. age','Avg. SN depigmentation','Avg. Lewy body stage','Avg. Braak score','Avg. Tangle','Avg. plaque','PM interval'))
for(g in unique(cold1$group)){for(s in unique(cold1$sex)){
  tmp<-subset(cold1,cold1$group==g &cold1$sex==s)
  sbjs<-rbind(sbjs,data.frame('Group'=g,'Sex'=s,'N'=nrow(tmp),'Avg. age'=round(mean(tmp$age),1),
                              'Avg. SN depigmentation'=levels(tmp$sn_depigmentation)[round(mean(as.numeric(tmp$sn_depigmentation)),0)],
                              'Avg. Lewy body stage'=levels(tmp$LewyBody_Stage)[round(mean(as.numeric(tmp$LewyBody_Stage)),0)],
                              'Avg. Braak score'=round(mean(as.numeric(tmp$BraakScore)),1),
                              'Avg. Tangle'=round(mean(tmp$TangleTotal),1),'Avg. plaque'=round(mean(tmp$PlaqueTotal),1),
                              'PM interval'=round(mean(tmp$PostMortemInterval),1)))
}}
TukeyHSD(aov(cold1$age~cold1$group*cold1$sex))$`cold1$group:cold1$sex`
TukeyHSD(aov(cold1$PostMortemInterval~cold1$group*cold1$sex))$`cold1$group:cold1$sex`
# write.csv(sbjs,dir('sup CSF Subjects info.csv'))

tmp<-subset(ColData,(ColData$ID %in% ToRemove1 | ColData$age<=60 | ColData$ADPD) & ColData$biofluid=='CSF')
tmp$exclusion<-factor(as.numeric(tmp$ADPD)+2*as.numeric(!tmp$keep_trf)+as.numeric(tmp$ID %in% p),
                      labels=c('Age < 60','ADPD patients','Controls With PD/AD phenotypes','Excluded by PCA'))
tmp<-tmp[,c("ID","group","age","sex","sn_depigmentation","LewyBody_Stage","PlaqueTotal","PlaqueDensity","TangleTotal","BraakScore","exclusion")]
# write.csv(tmp,dir('excluded CSF patients.csv'))

# Read CPM counts table:
cpm_trfs<-read.csv(dir('trf_counts_simulated.csv')); rownames(cpm_trfs)<-cpm_trfs$X; cpm_trfs<-cpm_trfs[,-1] ; cpm_trfs<-as.data.frame(cpm(cpm_trfs))
cpm_trfs<-cpm_trfs[rowMedians(as.matrix(cpm_trfs[,colnames(cpm_trfs) %in% subset(ColData$ID,ColData$biofluid=='CSF')]))>10,]

# Read tRF counts table:
trfs<-read.csv(dir('trf_counts_simulated.csv')); 
rownames(trfs)<-trfs$X; trfs<-trfs[,-1] ; colnames(trfs)<-gsub('_dbGaP.28067','',colnames(trfs))
trfs<-trfs[rownames(cpm_trfs),]

# Remove patients with too low counts:
tmp<-data.frame(ID=colnames(trfs),sum=colSums(trfs)) ; cold1<-merge(cold1,tmp,by='ID') ; cold1<-subset(cold1,cold1$sum>10000)

# Read and process tRFs metadata:
meta_trfs<-read.csv(dir('tRF_meta_simulated.csv')); rownames(meta_trfs)<-meta_trfs$X; meta_trfs<-meta_trfs[,-1]; meta_trfs<-meta_trfs[rownames(trfs),]
meta_trfs$len<-unlist(lapply(meta_trfs$MINTbase.Unique.ID, function(x) strsplit(as.character(x),'-')[[1]][2]))
meta_trfs$trna<-unlist(lapply(meta_trfs$Sequence.locations.in.tRNA.space..comma.deliminated., function(x) substr(strsplit(as.character(x),'_')[[1]][2],1,3)))
colnames(meta_trfs)<-gsub('Sequence.locations.in.tRNA.space..comma.deliminated.','details',colnames(meta_trfs))

# Define Nuclear and Mitochondrial tRFs:
meta_trfs$nuclear<-factor(unlist(lapply(meta_trfs$details,function(x)
  as.numeric(sum(c('trnaMT','trnalookalike2','trnalookalike8') %in% strsplit(as.character(x),'_')[[1]])!=0))),labels=c('Nuclear','Mitochondrial'))

# Find PD-tRFs
meta_trfs$motifs<-F
for(i in 1:nrow(meta_trfs)){
  meta_trfs$motifs[i]<-sum(str_count(as.character(meta_trfs$tRF.sequence[i]),as.character('[AG]GTTC[GA]A')))}
meta_trfs$motifs<-factor(meta_trfs$motifs,labels=c('No RGTTCRA','RGTTCRA motif'))

# Create updated CPM counts table:
cts<-trfs; cld<-subset(cold1,cold1$group!='AD'); cld<-subset(cld,cld$ID %in% colnames(cts))
cts<-cts[,as.character(cld$ID)] ; cts<-cts[,order(cld$ID)] ; cld<-cld[order(cld$ID),]
y <- DGEList(counts=cts,group=cld$group) ; keep<-filterByExpr(y) ; y<-y[keep, , keep.lib.sizes=FALSE]
y$samples$lib.size<-colSums(y$counts) ; y<-calcNormFactors(y) ; cpm_trfs<-as.data.frame(cpm(y,log = F))

table(cold1$group,cold1$sex) ; print('Sample sizes:'); table(subset(cold1,cold1$age>60)$group)

#
# DESeq: #####
cts<-trfs ; cld<-subset(cold1,cold1$ID %in% colnames(cts)) ; cts<-cts[,as.character(cld$ID)] 
cts<-cts[,order(cld$ID)] ; cld<-cld[order(cld$ID),] ; nrow(cld)==sum(cld$ID==colnames(cts))
y <- DGEList(counts=cts,group=cld$group) # change the condition to the name of your condition column
keep <- filterByExpr(y) ; y <- y[keep, , keep.lib.sizes=FALSE] ; y$samples$lib.size <- colSums(y$counts) ; y <- calcNormFactors(y)

cts1<-as.data.frame(cpm(y,log = F)) # creates the normalized counts matrix

# Define design matrix:
dsgn <- model.matrix(~PostMortemInterval+age+sex+group, data = cld)
y <- estimateDisp(y, dsgn, robust = T) ; head(dsgn,3)

# you change the coef to the coeficient that intersts you (the number of the column in the dsgn matrix)
fit <- glmQLFit(y, dsgn, robust = T) ; lrt <- glmLRT(fit,coef = 6) 

# Remove tRFs with expression outliers (Mean higher than the 85th quantile):
toRemove<-c() ; for(g in rownames(cpm_trfs)){
  if(quantile(as.numeric(cpm_trfs[g,]),prob=0.85)<mean(as.numeric(cpm_trfs[g,]))){toRemove<-append(toRemove,g)}}

# Generate significat tRFs table for PD:
t_sgGens<-as.data.frame(topTags(lrt,adjust.method = 'fdr',n = nrow(cts1))) ; t_sgGens$transcript<-rownames(t_sgGens)
t_sgGens$isSg<-factor(as.numeric(t_sgGens$FDR<0.051),labels = c('N.S.','Sgnificant'))
t_sgGens<-merge(t_sgGens,meta_trfs,by.x='transcript',by.y=colnames(meta_trfs)[1])
t_sgGens<-subset(t_sgGens,t_sgGens$logCPM>2.5 & !t_sgGens$transcript %in% toRemove)
t_sgGens$tRF.type.s.<-factor(t_sgGens$tRF.type.s.,levels = c("5'-half","5'-tRF","i-tRF","3'-tRF","3'-half"))

# Generate significat tRFs table for AD:
lrt2 <- glmLRT(fit,coef = 5)
ad_sgGens<-as.data.frame(topTags(lrt2,adjust.method = 'fdr',n = nrow(cts1))) ; ad_sgGens$transcript<-rownames(ad_sgGens)
ad_sgGens$isSg<-factor(as.numeric(ad_sgGens$FDR<0.051),labels = c('N.S.','Sgnificant'))
ad_sgGens<-merge(ad_sgGens,meta_trfs,by.x='transcript',by.y=colnames(meta_trfs)[1])
ad_sgGens$tRF.type.s.<-factor(ad_sgGens$tRF.type.s.,levels = c("5'-half","5'-tRF","i-tRF","3'-tRF","3'-half"))


## Visualise results:
# Density plots:
ggplot(cld,aes(age,col=group))+theme_classic()+geom_density(lwd=1.5)+facet_wrap(~sex)
# ggsave(dir('density_age.svg'))
TukeyHSD(aov(cld$age~cld$group))

ggplot(cld,aes(PostMortemInterval,col=group))+theme_classic()+geom_density(lwd=1.5)+facet_wrap(~sex)
# ggsave(dir('density_PMI.svg'))
TukeyHSD(aov(cld$PostMortemInterval~cld$group))
# AD:
ggplot(ad_sgGens,aes(log2(2.72^(logFC)),-log10(FDR),alpha=FDR<0.05,col=as.numeric(len)))+geom_point(size=3)+theme_classic()+
  scale_color_gradient2(low='green',mid='lightblue',high='purple',midpoint=35,name='tRF length')+
  geom_hline(yintercept=-log10(0.05),col='grey',linetype='dotted',lwd=1.15)+xlab('log2(FC)')+ylab('-log10(p.adjusted)')+
  geom_vline(xintercept=c(-1,1),col='grey',linetype='dotted',lwd=1.15)+scale_alpha_discrete(guide='none')+facet_wrap(~nuclear)
# ggsave(dir('Volcano_AD_tRFs_byLengthAndGenome.svg'),device='svg')

tmp<-subset(ad_sgGens,ad_sgGens$nuclear=='Nuclear')
ggplot(tmp,aes(log2(2.72^logFC),-log10(FDR),col=as.numeric(len)))+geom_point(size=3)+theme_classic()+
  facet_grid(~motifs)+geom_hline(yintercept=-log10(0.05),lwd=1,col='blue')+geom_hline(yintercept=0.6,col='grey')+
  ylab('-Log10(FDR)')+xlab('Log2(FoldChange)')+scale_color_gradient2(low='green',mid='lightblue',high='purple',midpoint=35,name='Length')
# ggsave(dir('motif AD DE.svg'),device = 'svg')

# PD:
ggplot(t_sgGens,aes(log2(2.72^(logFC)),-log10(FDR),alpha=FDR<0.05,col=as.numeric(len)))+geom_point(size=3)+theme_classic()+
  scale_color_gradient2(low='green',mid='lightblue',high='purple',midpoint=35,name='tRF length')+
  geom_hline(yintercept=-log10(0.05),col='grey',linetype='dotted',lwd=1.15)+xlab('log2(FC)')+ylab('-log10(p.adjusted)')+
  geom_vline(xintercept=c(-1,1),col='grey',linetype='dotted',lwd=1.15)+scale_alpha_discrete(guide='none')+facet_wrap(~nuclear)
# ggsave(dir('Volcano_tRFs_byLengthAndGenome.svg'),device='svg')

ggplot(t_sgGens,aes(log2(2.72^logFC),-log10(FDR),alpha=FDR<0.05))+theme_classic()+facet_grid(nuclear~tRF.type.s.)+
  geom_point(col='#AC88FF',size=2)+ylab('-log10(p.adjusted)')+scale_alpha_discrete(guide='none')
# ggsave(dir('Volcano_Nuclear_tRFs_by_Type.svg'),device='svg')

tmp<-subset(t_sgGens,t_sgGens$nuclear=='Nuclear')
ggplot(tmp,aes(log2(2.72^logFC),-log10(FDR),col=as.numeric(len)))+geom_point(size=3)+theme_classic()+
  facet_grid(~motifs)+geom_hline(yintercept=-log10(0.05),lwd=1,col='blue')+geom_hline(yintercept=0.6,col='grey')+
  ylab('-Log10(FDR)')+xlab('Log2(FoldChange)')+scale_color_gradient2(low='green',mid='lightblue',high='purple',midpoint=35,name='Length')
# ggsave(dir('motif DE.svg'),device = 'svg')

# binomial tests:     https://www.statology.org/binomial-test-r/ #
t_sgGens$int <- interaction(t_sgGens$motifs,t_sgGens$nuclear,sep='_')
binom_df<-ndf(c('group','number_decrease','total','p_value','CI_low','CI_high','prob_of_success'))
for (i in unique(t_sgGens$int)){
  tmp <- subset(t_sgGens,t_sgGens$int==i)
  binom_res <- binom.test(x = length(subset(tmp$transcript,tmp$logFC>0)),n = nrow(tmp),p = 0.5)
  binom_df<-rbind(binom_df,data.frame(group=i,number_decrease=binom_res[[1]],total=binom_res[[2]],p_value=binom_res[[3]],
                                      CI_low=binom_res[[4]][1],CI_high=binom_res[[4]][2],prob_of_success=binom_res[[5]]))
} ; binom_df$FDR <- p.adjust(binom_df$p_value, method = "fdr")

tmp<-table(t_sgGens$motifs,t_sgGens$logFC>0) ; colnames(tmp)<-c('Down','Up')
ch1<-chisq.test(tmp) ; ch1$p.value ; ch1$observed/ch1$expected #; ch1$expected ; ch1$observed 

tmp<-subset(t_sgGens,t_sgGens$tRF.type.s. %in% c("i-tRF","3'-tRF") & t_sgGens$nuclear=='Nuclear') ; tmp<-table(tmp$motifs,tmp$logFC>0)
colnames(tmp)<-c('Down','Up') ; ch1<-chisq.test(tmp) ; ch1$p.value ; ch1$observed/ch1$expected #; ch1$observed ; ch1$expected 

tmp<-data.frame(Ctrl=c(mean(colSums(cpm_trfs[rownames(cpm_trfs) %in% subset(meta_trfs$MINTbase.Unique.ID,meta_trfs$motifs=='No RGTTCRA'),
                                     colnames(cpm_trfs) %in% subset(cold1$ID,cold1$group=='Ctrl')])),
                       mean(colSums(cpm_trfs[rownames(cpm_trfs) %in% subset(meta_trfs$MINTbase.Unique.ID,meta_trfs$motifs=='RGTTCRA motif'),
                                            colnames(cpm_trfs) %in% subset(cold1$ID,cold1$group=='Ctrl')]))),
                PD=c(mean(colSums(cpm_trfs[rownames(cpm_trfs) %in% subset(meta_trfs$MINTbase.Unique.ID,meta_trfs$motifs=='No RGTTCRA'),
                                            colnames(cpm_trfs) %in% subset(cold1$ID,cold1$group=='PD')])),
                     mean(colSums(cpm_trfs[rownames(cpm_trfs) %in% subset(meta_trfs$MINTbase.Unique.ID,meta_trfs$motifs=='RGTTCRA motif'),
                                            colnames(cpm_trfs) %in% subset(cold1$ID,cold1$group=='PD')]))))
rownames(tmp)<-c('No RGTTCRA','RGTTCRA motif')
ch1<-chisq.test(tmp) ; ch1$p.value ; ch1$observed/ch1$expected #; ch1$observed ; ch1$expected
