library(ggplot2) ; library(matrixStats) ; library("biomaRt") ; library(reshape)

ndf<-function(x) as.data.frame(t(data.frame(row.names=x))) ; Ttmp<-function(x) as.data.frame(t(x))
Ftmp<-function(df,col='exp',id='id') {
  assign('tmp',as.data.frame(t(df))); colnames(tmp)<-col; tmp[,id]<-rownames(tmp) ; return(tmp)}
reform<-function(x){rownames(x)<-unlist(lapply(x[,1], function(y) strsplit(as.character(y),'\\.')[[1]][1])) ; x<-x[,2:ncol(x)]}
dir<-function(x) paste0(c(getwd(),x),collapse='/') ; dir('')

cts<-read.table(dir('GSE23676.txt.gz'),header=T,skip=68,nrows=1411399)
cold<-read.csv(dir('coldata.csv'))
for(c in colnames(cts)[2:ncol(cts)]){colnames(cts)<-gsub(c,subset(cold$id,cold$sample==c),colnames(cts))}

mart<-useEnsembl(biomart = "ensembl",dataset = "hsapiens_gene_ensembl", mirror = "useast") # View(attributes(mart)$attributes)
G_list<-getBM(filters="hgnc_symbol",attributes=c("affy_huex_1_0_st_v2","ensembl_gene_id","hgnc_symbol"),values='Ang',mart= mart)
G_list$ID_REF<-unlist(lapply(G_list$affy_huex_1_0_st_v2,function(x) strsplit(as.character(x),'_')[[1]][1]))

ang<-ndf(c('id','gene','sum','mean','mx','sd'))
tmp<-cts[cts$ID_REF %in% G_list$ID_REF,] ; tmp<-tmp[,-1]
tmp1<-data.frame(id=colnames(tmp),gene=g,sum=colSums(tmp),mean=colMeans(tmp),
                   mx=as.numeric(tmp[rowSums(tmp)==max(rowSums(tmp)),]),
                   sd=as.numeric(tmp[rowSds(as.matrix(tmp))==max(rowSds(as.matrix(tmp))),]))
ang<-rbind(ang,tmp1) ; ang<-merge(ang,cold,by='id') ; ang$cnd<-factor(ang$cnd,levels=c('Ctrl','pre-DBS','DBS-on','DBS-off'))
ang<-subset(ang,ang$cnd!='DBS-off')

ggplot(ang,aes(cnd,sum,col=group))+theme_classic()+geom_boxplot(outlier.shape=NA)+geom_point()+
  geom_line(aes(group=subject),col='black')+ylab('Ang expression')+#facet_wrap(~gene,scales='free_y')+
  theme(legend.title=element_blank(),axis.title.x=element_blank())
# ggsave(dir('Ang Expression DBS.svg'),width=2.5,height=2.5)
TukeyHSD(aov(ang$sum~ang$cnd))



