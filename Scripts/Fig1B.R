#updated by E.S. de Camargo Magalhaes on 20/04/2024 using R v4.3.2 "Eye Holes"

####Mackay et al. 2017 Dataset analysis----
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Rmisc))
suppressPackageStartupMessages(library(survminer))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(rstatix))

setwd('data/Datasets/Fig1B')#custom folder for importing files
pfkfb2_status<-as.data.frame(read.delim('PFKFB2_status.txt',header=T,sep="\t",dec = "."))
pfkfb2_sv<-as.data.frame(read.delim('PFKFB2_sv.txt',header=T,sep="\t",dec = "."))
pfkfb3_status<-as.data.frame(read.delim('PFKFB3_status.txt',header=T,sep="\t",dec = "."))
pfkfb3_sv<-as.data.frame(read.delim('PFKFB3_sv.txt',header=T,sep="\t",dec = "."))
pfkfb4_status<-as.data.frame(read.delim('PFKFB4_status.txt',header=T,sep="\t",dec = "."))
pfkfb4_sv<-as.data.frame(read.delim('PFKFB4_sv.txt',header=T,sep="\t",dec = "."))
hk2_status<-as.data.frame(read.delim('HK2_status.txt',header=T,sep="\t",dec = "."))
hk2_sv<-as.data.frame(read.delim('HK2_sv.txt',header=T,sep="\t",dec = "."))
pfkp_status<-as.data.frame(read.delim('PFKP_status.txt',header=T,sep="\t",dec = "."))
pfkp_sv<-as.data.frame(read.delim('PFKP_sv.txt',header=T,sep="\t",dec = "."))
pgk1_status<-as.data.frame(read.delim('PGK1_status.txt',header=T,sep="\t",dec = "."))
pgk1_sv<-as.data.frame(read.delim('PGK1_sv.txt',header=T,sep="\t",dec = "."))

#Create color object
mycolors1<-c('red3','blue3','green4','darkmagenta','goldenrod4','darkorange','deeppink',
             'gray60','darkcyan','darkturquoise')

#Treating the data for analysis
raw_table<-as.data.frame(join_all(list(pfkfb2_status,pfkfb3_status,pfkfb4_status,hk2_status,pfkp_status,pgk1_status
                                       ,pfkfb2_sv,pfkfb3_sv,pfkfb4_sv,hk2_sv,pfkp_sv,pgk1_sv),by='Sample.Id',type='full'))

colnames(raw_table)<-c('id','histone','PFKFB2','status','PFKFB3','PFKFB4','HK2','PFKP','PGK1','time')
raw_table$time<-round(as.numeric(raw_table$time)/12,digits=2)
raw_table$status<-as.numeric(ifelse(raw_table$status=='1:DECEASED',1,ifelse(raw_table$status=='0:LIVING',0,NA)))
raw_table<-raw_table %>% mutate(across(3:10,function(x) as.numeric(as.character(x))))
raw_table$histone<-factor(ifelse(raw_table$histone=='wt','H3.3-WT',
                                 ifelse(raw_table$histone=='H3.3_G34RV','H3.3-G34RV',
                                        ifelse(raw_table$histone=='H3.3_K27M','H3.3K27M','H3.1-K27M'))))
raw_table$histone<-factor(raw_table$histone,levels=c('H3.3-WT','H3.3-G34RV','H3.3K27M','H3.1-K27M'))
raw_table$id<-gsub("pHGG_META_", "", raw_table$id)
raw_table$id<-as.integer(raw_table$id)
raw_table<-as.data.frame(raw_table[,c('id','histone','status','time','PFKFB2','PFKFB3','PFKFB4','HK2','PFKP','PGK1')])
boxplot(raw_table[,c('PFKFB2','PFKFB3','PFKFB4','HK2','PFKP','PGK1')])
abline(h=0, col='blue')
#mRNA expression data consists in z-scores (considering all samples, the median expression of each gene is zero)

#Survival analysis and gene expression analysis
dim(raw_table)
surv_data<-na.omit(raw_table)
surv_data<-surv_data[surv_data$histone=='H3.3-WT'| surv_data$histone=='H3.3K27M',]
surv_data$histone<-droplevels(surv_data$histone)
dim(surv_data)
rownames(surv_data)<-c(1:158)
#Select only patients with wild type WT or H3.3K27M

summary(surv_data$histone)
dim(surv_data)

#Normality tests
summary(surv_data$histone)
#Numbers are higher than 50, so use KS test
ks.test(surv_data[surv_data$histone=='H3.3-WT',]$PFKFB2,"pnorm")
ks.test(surv_data[surv_data$histone=='H3.3K27M',]$PFKFB2,"pnorm")
ks.test(surv_data[surv_data$histone=='H3.3-WT',]$PFKFB3,"pnorm")
ks.test(surv_data[surv_data$histone=='H3.3K27M',]$PFKFB3,"pnorm")
ks.test(surv_data[surv_data$histone=='H3.3-WT',]$PFKFB4,"pnorm")
ks.test(surv_data[surv_data$histone=='H3.3K27M',]$PFKFB4,"pnorm")
ks.test(surv_data[surv_data$histone=='H3.3-WT',]$HK2,"pnorm")
ks.test(surv_data[surv_data$histone=='H3.3K27M',]$HK2,"pnorm")
ks.test(surv_data[surv_data$histone=='H3.3-WT',]$PGK1,"pnorm")
ks.test(surv_data[surv_data$histone=='H3.3K27M',]$PGK1,"pnorm")
ks.test(surv_data[surv_data$histone=='H3.3-WT',]$PFKP,"pnorm")
ks.test(surv_data[surv_data$histone=='H3.3K27M',]$PFKP,"pnorm")
#Distribution is NOT normal in all cases(p<0.05 in all tests)
#Use non-parametric test for all

surv_melt<-surv_data[,c(2,5:10)]
surv_melt<-as.data.frame(reshape2::melt(surv_melt, id=c('histone')))
colnames(surv_melt)<-c('histone','gene','z_score')
surv_melt<-as.data.frame(surv_melt[,c('gene','z_score','histone')])

#Compare means
test<-surv_melt %>% group_by(gene) %>% pairwise_wilcox_test(z_score~histone)%>% adjust_pvalue(method="fdr")%>%
  add_significance() %>% add_xy_position(x='gene')

theme_set(theme_classic(base_size = 11,base_line_size = 0.05,base_rect_size = 0.05))

mackay_dpl<-ggplot(surv_melt, aes(x=gene, y=z_score))+ggtitle('Glycolysis-related genes')+
  geom_dotplot(aes(fill=histone),binaxis='y',stackdir='center',position=position_dodge(0.8),dotsize=2.5,binwidth=0.06,color='transparent')+
  stat_summary(aes(group=histone),fun.data=mean_se,geom="errorbar",position = position_dodge(width=0.8), color="black", linewidth=0.3,width=0.15)+
  stat_summary(aes(group=histone),fun=mean, geom="crossbar",position = position_dodge(width=0.8),linewidth=0.15,width=0.3, color="black")+
  stat_pvalue_manual(test, label = "p.adj.signif",size=7,face="bold")+
  scale_fill_manual(values=mycolors1[2:1])+labs(y='mRNA expression (Z-scores)')+
  scale_y_continuous(limits = c(-3,6.5),breaks=seq(-2,6,by=2))+
  theme(legend.position="bottom",legend.title=element_blank(),legend.text=element_text(hjust=0.5,size=16,face='bold')
        ,axis.title.y=element_text(size=20,face="bold",margin=margin(t=0,r=15,b=0,l=0)),axis.title.x= element_blank()
        ,axis.text=element_text(size=18,face="bold"),plot.title = element_text(size=20,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=20, color="black",face="bold"),legend.key.size = unit(3.5,"line"))

setwd('/results')#custom destination folder for file generation
ggexport(mackay_dpl,filename='Fig1B.pdf',width=11.69, height=8.27, pointsize=8, res=250) #A4 landscape
