library(edgeR); library(ggplot2); library(ggfortify); library('biomaRt'); library(matrixStats); library(reshape); library(stringr) ; library(M3C)

reform<-function(x){rownames(x)<-unlist(lapply(x[,1], function(y) strsplit(as.character(y),'\\.')[[1]][1])) ; x<-x[,2:ncol(x)]}
dir<-function(x) paste0(c(getwd(),x),collapse='/') ; dir('')
ndf<-function(x) as.data.frame(t(data.frame(row.names=x)))

# Read files #####
cold<-read.csv(dir('ColData_new.csv')) ; colnames(cold)<-gsub('Run','id',colnames(cold))
trfs<-read.csv(dir('tRNA_Exclusive_Combined_data.csv')) ; rownames(trfs)<-trfs$X ; trfs<-trfs[,-1]

meta_trfs<-read.csv(dir('tRF_meta.csv')); rownames(meta_trfs)<-meta_trfs$X; meta_trfs<-meta_trfs[,-1]; meta_trfs<-meta_trfs[rownames(trfs),]
meta_trfs$G4s<-as.factor(unlist(lapply(meta_trfs$tRF.sequence, function(x) as.numeric(substr(as.character(x),1,4)=='GGGG'))))
meta_trfs$motifs<-F
for(i in 1:nrow(meta_trfs)){
  meta_trfs$motifs[i]<-sum(str_count(as.character(meta_trfs$tRF.sequence[i]),as.character('[AG]GTTC[GA]A')))}
meta_trfs$motifs<-factor(meta_trfs$motifs,labels=c('No RGTTCRA','RGTTCRA motif'))
meta_trfs$nuclear<-factor(unlist(lapply(meta_trfs$details,function(x) sum(str_count(as.character(as.character(x)),as.character('trnaMT')))!=0)),
                          labels=c('Nuclear','Mitochondrial'))
meta_trfs$loc<-unlist(lapply(meta_trfs$details,function(x) strsplit(as.character(x),'[-,+]')[[1]][1]))

tcpm<-as.data.frame(cpm(trfs))

pd<-data.frame(id=colnames(trfs),sum=colSums(trfs),
               pd=colSums(trfs[rownames(trfs) %in% subset(meta_trfs$MINTbase.Unique.ID,meta_trfs$motifs=='RGTTCRA motif'),]),
               mt=colSums(trfs[rownames(trfs) %in% subset(meta_trfs$MINTbase.Unique.ID,meta_trfs$nuclear=='Mitochondrial'),]),
               ctrl=colSums(trfs[rownames(trfs) %in% subset(meta_trfs$MINTbase.Unique.ID,meta_trfs$tRF.type.s. %in% c("3'-tRF","i-tRF")),]))
pd$r_pd<-100*pd$pd/pd$sum ; pd$r_mt<-100*pd$mt/pd$sum ; pd$ratio<-pd$pd/pd$mt ; pd$ctrl<-100*pd$ctrl/pd$sum ; pd<-merge(pd,cold,by='id')
pd$source_name<-factor(pd$source_name,levels=c("HEK293T EV","HEK293T ANG overexpression",
                                "U2OS WT untreated","U2OS ANG KO untreated","U2OS WT sodium arsenite","U2OS ANG KO sodium arsenite"),
                       labels=c("HEK EV","HEK ANG OE","U2OS WT UT","U2OS ANG KO UT","U2OS WT SA","U2OS ANG KO SA"))

ggplot(pd,aes(source_name,ctrl,fill=source_name))+theme_classic()+facet_wrap(~Cell_Line,scales='free')+
  geom_boxplot(outlier.shape=NA)+geom_point(size=3)+theme(legend.title=element_blank(),axis.title.x=element_blank())+
  scale_fill_manual(values=c('grey33','lightblue','grey33','lightgreen','red','violet'))
# ggsave(dir('PD-tRF Ang.svg'),width=6,height=2.5)

tmp<-pd[pd$Cell_Line=='HEK293T',] ; t.test(tmp$r_pd~tmp$source_name)
tmp<-pd[pd$Cell_Line!='HEK293T',] ; TukeyHSD(aov(tmp$r_pd~tmp$Genotype*tmp$stress))


ggplot(pd,aes(source_name,r_mt,fill=source_name))+theme_classic()+facet_wrap(~Cell_Line,scales='free')+
  geom_boxplot(outlier.shape=NA)+geom_point(size=3)+theme(legend.title=element_blank(),axis.title.x=element_blank())+
  scale_fill_manual(values=c('grey33','lightblue','grey33','lightgreen','red','violet'))
# ggsave(dir('MT-tRF Ang.svg'),width=6,height=2.5)

tmp<-pd[pd$Cell_Line=='HEK293T',] ; t.test(tmp$r_mt~tmp$source_name)
tmp<-pd[pd$Cell_Line!='HEK293T',] ;TukeyHSD(aov(tmp$r_mt~tmp$Genotype*tmp$stress))


#