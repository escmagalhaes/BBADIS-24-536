#updated by E.S. de Camargo Magalhaes on 20/04/2024 using R v4.3.2 "Eye Holes"

####Seahorse analysis 1----
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Rmisc))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(rstatix))

###Import excel file generated from wave (overview) and create OCR and ECAR datasets
setwd('data/Datasets/Fig1C_D')#custom folder for importing files
sh_res3<-as.data.frame(read_xlsx('seahorse_assay1.xlsx'))
sh_res3[sh_res3<0]<-NA

#OCR
ocr3<-as.data.frame(sh_res3[c(1:24),c(1:52)])
ocr3<-ocr3[,colSums(is.na(ocr3))<nrow(ocr3)]
#Remove NA rows between samples
colnames(ocr3)<-c('time',paste0('H3.3K27M_rep_',c(1:22)),paste0('H3.3K27M_KO_rep_',c(1:24)))
row.names(ocr3)<-NULL
ocr3<-ocr3[c(-1:-3),]
row.names(ocr3)<-c(1:21)
ocr3<-mutate_all(ocr3, function(x) as.numeric(as.character(x)))
#Remove unused rows, rename rows and convert values to numeric

#ECAR
ecar3<-as.data.frame(sh_res3[c(28:51),c(1:52)])
ecar3<-ecar3[,colSums(is.na(ecar3))<nrow(ecar3)]
colnames(ecar3)<-c('time',paste0('H3.3K27M_rep_',c(1:22)),paste0('H3.3K27M_KO_rep_',c(1:24)))
row.names(ecar3)<-NULL
ecar3<-ecar3[c(-1:-3),]
row.names(ecar3)<-c(1:21)
ecar3<-mutate_all(ecar3, function(x) as.numeric(as.character(x)))

#Create color object
mycolors1<-c('red3','blue3','green4','darkmagenta','goldenrod4','darkorange','deeppink',
             'gray60','darkcyan','darkturquoise')

##Non-normalized data
#OCR
ocr3_tab<-ocr3[,-1]
#remove variable time for calculations and figure generation
non_mit3<-apply(ocr3_tab[c(12:21),],2,function(x){min(x[!is.na(x)])})
non_mit3<-ifelse(non_mit3==Inf,0,non_mit3)
basal_resp3<-ocr3_tab[3,]-non_mit3
max_resp3<-apply(ocr3_tab[c(8:11),],2,function(x){max(x[!is.na(x)])})-non_mit3
proton_leak3<-apply(ocr3_tab[c(4:7),],2,function(x){max(x[!is.na(x)])})-non_mit3
spare_cap_ocr3<-max_resp3-basal_resp3
spare_cap_ocr3_fd<-max_resp3/basal_resp3
spare_cap_ocr3_pct<-((max_resp3-basal_resp3)/max_resp3)*100
atp_prod<-basal_resp3-proton_leak3

ocr3_tab_grp<-rbind(basal_resp3,max_resp3,spare_cap_ocr3,spare_cap_ocr3_pct,spare_cap_ocr3_fd,atp_prod,non_mit3)
ocr3_tab_grp<-as.data.frame(t(ocr3_tab_grp))
ocr3_tab_grp$cell_type3<-factor(c(rep('H3.3K27M',22),rep('H3.3K27M-KO',24)),levels=c('H3.3K27M','H3.3K27M-KO'))
colnames(ocr3_tab_grp)<-c('basal_resp3','max_resp3','spare_cap_ocr3','spare_cap_ocr3_pct','spare_cap_ocr3_fd','atp_prod','non_mit3'
                          ,'cell_type3')

test_basal_resp3<-ocr3_tab_grp %>% wilcox_test(basal_resp3~cell_type3)%>% 
  add_significance() %>% add_xy_position(x='cell_type3') %>% mutate(y.position = max(basal_resp3)+sd(basal_resp3)*2)
test_max_resp3<-ocr3_tab_grp %>% wilcox_test(max_resp3~cell_type3)%>% 
  add_significance() %>% add_xy_position(x='cell_type3') %>% mutate(y.position = max(max_resp3)+sd(max_resp3))
test_spare_cap_ocr3_pct<-ocr3_tab_grp %>% wilcox_test(spare_cap_ocr3_pct~cell_type3)%>% 
  add_significance() %>% add_xy_position(x='cell_type3') %>% mutate(y.position = max(spare_cap_ocr3_pct)+sd(spare_cap_ocr3_pct)*2)

#Linegraph
ocr3_melt<-as.data.frame(reshape2::melt(ocr3, id='time'))
ocr3_melt$time<-round(ocr3_melt$time,digits = 2)
ocr3_melt$time<-factor(ocr3_melt$time)
ocr3_melt$cell_type3<-factor(c(rep(rep('H3.3K27M',22),21),rep(rep('H3.3K27M-KO',24),21)),levels=c('H3.3K27M','H3.3K27M-KO'))
dim(ocr3_melt)
ocr3_melt<-na.omit(ocr3_melt)
dim(ocr3_melt)
ocr3_sum<-summarySE(ocr3_melt, measurevar="value", groupvars=c('time','cell_type3'))
ocr3_sum$time<-as.numeric(levels(ocr3_sum$time))[ocr3_sum$time]

ocr3_tab1<-ggplot(ocr3_tab_grp, aes(x=cell_type3, y=basal_resp3)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=5,binwidth = 1.5,aes(fill=cell_type3))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  stat_pvalue_manual(test_basal_resp3, label = "p.signif",size=6)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='pmol/min')+ggtitle('Basal Respiration')+
  scale_y_continuous(limits = c(0, 120),breaks=seq(0,120,by=20))+
  theme(legend.position="none",axis.title.y=element_text(size=24,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=24,face="bold"),plot.title = element_text(size=30,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=26, color="black",face="bold"))

ocr3_tab2<-ggplot(ocr3_tab_grp, aes(x=cell_type3, y=max_resp3)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=5,binwidth = 1.5,aes(fill=cell_type3))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='pmol/min')+ggtitle('Maximum Respiration')+
  stat_pvalue_manual(test_max_resp3, label = "p.signif",size=6)+
  scale_y_continuous(limits = c(0, 120),breaks=seq(0,120,by=20))+
  theme(legend.position="none",axis.title.y=element_text(size=24,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=24,face="bold"),plot.title = element_text(size=30,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=26, color="black",face="bold"))

ocr3_tab3<-ggplot(ocr3_tab_grp, aes(x=cell_type3, y=spare_cap_ocr3_pct)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=5,binwidth = 1.5,aes(fill=cell_type3))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='Percentage(%)')+ggtitle('Spare Respiratory Capacity')+
  stat_pvalue_manual(test_spare_cap_ocr3_pct, label = "p.signif",size=6)+
  scale_y_continuous(limits = c(0,100),breaks=seq(0,100,by=20))+
  theme(legend.position="none",axis.title.y=element_text(size=24,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=24,face="bold"),plot.title = element_text(size=30,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=26, color="black",face="bold"))

ocr3_tab4<-ggplot(ocr3_sum, aes(x=time, y=value, group=cell_type3,color=cell_type3))+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=1.5,size=1)+geom_line(size=1)+geom_point(size=2)+
  scale_linetype_manual(values=c(rep(c('solid'),2)))+labs(x='Time (min)',y='pmol/min')+
  scale_color_manual(values=mycolors1[1:2])+ggtitle('Oxygen Consumption Rate (OCR)')+
  theme(legend.position="bottom",legend.title=element_blank(),legend.text=element_text(hjust=0.5,size=20,face='bold')
        ,axis.title=element_text(size=24,face="bold"),axis.text = element_text(size=22)
        ,plot.title = element_text(size= 26,hjust = 0.5,face='bold'),legend.key.size = unit(3.4,"line",'point'))

resp3_grp<-ggarrange(ocr3_tab1,NULL,ocr3_tab2,NULL,NULL,NULL,ocr3_tab3,NULL,ocr3_tab4,ncol=3,nrow=3
                     ,widths=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1),heights=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1))+theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))

#ECAR
ecar3_tab<-ecar3[,-1]
basal_ecar3<-ecar3_tab[3,]
max_ecar3<-ecar3_tab[4,]
min_ecar3<-ecar3_tab[21,]
basal_glyco3<-basal_ecar3-min_ecar3
max_glyco3<-max_ecar3-min_ecar3
spare_cap_glyco3<-max_glyco3-basal_glyco3
spare_cap_glyco3_fd<-max_glyco3/basal_glyco3
spare_cap_glyco3_pct<-((max_glyco3-basal_glyco3)/max_glyco3)*100

glyco3_tab_grp<-rbind(basal_glyco3,max_glyco3,spare_cap_glyco3,spare_cap_glyco3_pct,spare_cap_glyco3_fd)
glyco3_tab_grp<-as.data.frame(t(glyco3_tab_grp))
glyco3_tab_grp$cell_type3<-factor(c(rep('H3.3K27M',22),rep('H3.3K27M-KO',24)),levels=c('H3.3K27M','H3.3K27M-KO'))
glyco3_tab_grp<-na.omit(glyco3_tab_grp)
colnames(glyco3_tab_grp)<-c('basal_glyco3','max_glyco3','spare_cap_glyco3','spare_cap_glyco3_pct','spare_cap_glyco3_fd','cell_type3')

#Linegraph
ecar3_melt<-as.data.frame(reshape2::melt(ecar3, id='time'))
ecar3_melt$time<-round(ecar3_melt$time,digits = 2)
ecar3_melt$time<-factor(ecar3_melt$time)
ecar3_melt$cell_type3<-factor(c(rep(rep('H3.3K27M',22),21),rep(rep('H3.3K27M-KO',24),21)),levels=c('H3.3K27M','H3.3K27M-KO'))
dim(ecar3_melt)
ecar3_melt<-na.omit(ecar3_melt)
dim(ecar3_melt)
ecar3_sum<-summarySE(ecar3_melt, measurevar="value", groupvars=c('time','cell_type3'))
ecar3_sum$time<-as.numeric(levels(ecar3_sum$time))[ecar3_sum$time]

test_basal_glyco3<-glyco3_tab_grp %>% wilcox_test(basal_glyco3~cell_type3)%>% 
  add_significance() %>% add_xy_position(x='cell_type3') %>% mutate(y.position = max(glyco3_tab_grp$basal_glyco3)+sd(glyco3_tab_grp$basal_glyco3)*2)
test_max_glyco3<-glyco3_tab_grp %>% wilcox_test(max_glyco3~cell_type3)%>% 
  add_significance() %>% add_xy_position(x='cell_type3') %>% mutate(y.position = max(glyco3_tab_grp$max_glyco3)+sd(glyco3_tab_grp$max_glyco3))
test_spare_cap_glyco3_pct<-glyco3_tab_grp %>% wilcox_test(spare_cap_glyco3_pct~cell_type3)%>% 
  add_significance() %>% add_xy_position(x='cell_type3') %>% mutate(y.position = max(glyco3_tab_grp$spare_cap_glyco3_pct)+sd(glyco3_tab_grp$spare_cap_glyco3_pct)*2)


glyco3_tab1<-ggplot(glyco3_tab_grp, aes(x=cell_type3, y=basal_glyco3)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=2.5,binwidth = 2.5,aes(fill=cell_type3))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='mpH/min')+ggtitle('Basal Glycolysis')+
  stat_pvalue_manual(test_basal_glyco3, label = "p.signif",size=6)+
  scale_y_continuous(limits = c(0,120),breaks=seq(0,100,by=20))+
  theme(legend.position="none",axis.title.y=element_text(size=24,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=24,face="bold"),plot.title = element_text(size=30,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=26, color="black",face="bold"))

glyco3_tab2<-ggplot(glyco3_tab_grp, aes(x=cell_type3, y=max_glyco3)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=2.5,binwidth = 2.5,aes(fill=cell_type3))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='mpH/min')+ggtitle('Maximum Glycolysis')+
  stat_pvalue_manual(test_max_glyco3, label = "p.signif",size=6)+
  scale_y_continuous(limits = c(0,120),breaks=seq(0,100,by=20))+
  theme(legend.position="none",axis.title.y=element_text(size=24,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=24,face="bold"),plot.title = element_text(size=30,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=26, color="black",face="bold"))

glyco3_tab3<-ggplot(glyco3_tab_grp, aes(x=cell_type3, y=spare_cap_glyco3_pct)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=5,binwidth = 1,aes(fill=cell_type3))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='Percentage(%)')+ggtitle('Spare Glycolytic Capacity')+
  stat_pvalue_manual(test_spare_cap_glyco3_pct, label = "p.signif",size=6)+
  scale_y_continuous(limits = c(0,100),breaks=seq(0,100,by=20))+
  theme(legend.position="none",axis.title.y=element_text(size=24,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=24,face="bold"),plot.title = element_text(size=30,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=26, color="black",face="bold"))

glyco3_tab4<-ggplot(ecar3_sum, aes(x=time, y=value, group=cell_type3,color=cell_type3))+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=1.5,size=1)+geom_line(size=1)+geom_point(size=2)+
  scale_linetype_manual(values=c(rep(c('solid'),2)))+labs(x='Time (min)',y='mpH/min')+
  scale_color_manual(values=mycolors1[1:2])+ggtitle('Extracellular Acidification Rate (ECAR)')+
  theme(legend.position="bottom",legend.title=element_blank(),legend.text=element_text(hjust=0.5,size=20,face='bold')
        ,axis.title=element_text(size=24,face="bold"),axis.text = element_text(size=22)
        ,plot.title = element_text(size= 26,hjust = 0.5,face='bold'),legend.key.size = unit(3.4,"line",'point'))

glyco3_grp<-ggarrange(glyco3_tab1,NULL,glyco3_tab2,NULL,NULL,NULL,glyco3_tab3,NULL,glyco3_tab4,ncol=3,nrow=3
                      ,widths=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1),heights=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1))+theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))

###Data Normalized by protein
bca_res_bt3<-as.data.frame(read_xlsx('BCA_seahorse1.xlsx'))
bca_res_bt3<-bca_res_bt3[c(22,23),c(2:49)]
#Create column with time variable for further calculations
time<-as.matrix(c('time',1))
ptn3<-as.data.frame(cbind(time,bca_res_bt3))
ptn3<-as.data.frame(t(ptn3[,colSums(is.na(ptn3))<nrow(ptn3)]))
rownames(ptn3)<-NULL
colnames(ptn3)<-c('names','protein')

ocr3_norm<-as.data.frame(t(ocr3))
ocr3_norm<-cbind(ocr3_norm,ptn3)
ocr3_norm<-ocr3_norm %>% dplyr::select('names','protein',everything())
#check if names match
ocr3_norm<-ocr3_norm[,-1]
ocr3_norm<-mutate_all(ocr3_norm, function(x) as.numeric(as.character(x)))
ocr3_norm<- ocr3_norm %>% mutate_at(vars(-protein),funs(./protein))%>% dplyr::select(c(1:22))
ocr3_norm<-ocr3_norm[,-1]
ocr3_norm<-as.data.frame(t(ocr3_norm))

ecar3_norm<-as.data.frame(t(ecar3))
ecar3_norm<-cbind(ecar3_norm,ptn3)
ecar3_norm<-ecar3_norm %>% dplyr::select('names','protein',everything())
#View(ecar3_norm) #check if names match
ecar3_norm<-ecar3_norm[,-1]
ecar3_norm<-mutate_all(ecar3_norm, function(x) as.numeric(as.character(x)))
ecar3_norm<- ecar3_norm %>% mutate_at(vars(-protein),funs(./protein))%>% dplyr::select(c(1:22))
ecar3_norm<-ecar3_norm[,-1]
ecar3_norm<-as.data.frame(t(ecar3_norm))

#Create subsets, treat data and generate figures
#Normalized OCR
ocr3_norm_tab<-ocr3_norm[,-1]
non_mit3_norm<-apply(ocr3_norm_tab[c(12:21),],2,function(x){min(x[!is.na(x)])})
non_mit3_norm<-ifelse(non_mit3_norm==Inf,0,non_mit3_norm)
basal_resp3_norm<-ocr3_norm_tab[3,]-non_mit3_norm
max_resp3_norm<-apply(ocr3_norm_tab[c(8:11),],2,function(x){max(x[!is.na(x)])})-non_mit3_norm
proton_leak3_norm<-apply(ocr3_norm_tab[c(4:7),],2,function(x){max(x[!is.na(x)])})-non_mit3_norm
spare_cap_ocr3_norm<-max_resp3_norm-basal_resp3_norm
spare_cap_ocr3_norm_fd<-max_resp3_norm/basal_resp3_norm
spare_cap_ocr3_norm_pct<-((max_resp3_norm-basal_resp3_norm)/max_resp3_norm)*100
atp_prod_norm<-basal_resp3_norm-proton_leak3_norm

ocr3_norm_tab_grp<-rbind(basal_resp3_norm,max_resp3_norm,spare_cap_ocr3_norm,spare_cap_ocr3_norm_pct
                         ,spare_cap_ocr3_norm_fd,atp_prod_norm,non_mit3_norm)
ocr3_norm_tab_grp<-as.data.frame(t(ocr3_norm_tab_grp))
ocr3_norm_tab_grp$cell_type3<-factor(c(rep('H3.3K27M',22),rep('H3.3K27M-KO',24)),levels=c('H3.3K27M','H3.3K27M-KO'))
colnames(ocr3_norm_tab_grp)<-c('basal_resp3_norm','max_resp3_norm','spare_cap_ocr3_norm','spare_cap_ocr3_norm_pct'
                               ,'spare_cap_ocr3_norm_fd','atp_prod_norm','non_mit3_norm','cell_type3')

test_basal_resp3_norm<-ocr3_norm_tab_grp %>% wilcox_test(basal_resp3_norm~cell_type3)%>% 
  add_significance() %>% add_xy_position(x='cell_type3') %>% mutate(y.position = max(basal_resp3_norm)+sd(basal_resp3_norm)*2)
test_max_resp3_norm<-ocr3_norm_tab_grp %>% wilcox_test(max_resp3_norm~cell_type3)%>% 
  add_significance() %>% add_xy_position(x='cell_type3') %>% mutate(y.position = max(max_resp3_norm)+sd(max_resp3_norm))
test_spare_cap_ocr3_norm_pct<-ocr3_norm_tab_grp %>% wilcox_test(spare_cap_ocr3_norm_pct~cell_type3)%>% 
  add_significance() %>% add_xy_position(x='cell_type3') %>% mutate(y.position = max(spare_cap_ocr3_norm_pct)+sd(spare_cap_ocr3_norm_pct)*2)

#Linegraph
ocr3_norm_melt<-as.data.frame(reshape2::melt(ocr3_norm, id='time'))
ocr3_norm_melt$time<-round(ocr3_norm_melt$time,digits = 2)
ocr3_norm_melt$time<-factor(ocr3_norm_melt$time)
ocr3_norm_melt$cell_type3<-factor(c(rep(rep('H3.3K27M',22),21),rep(rep('H3.3K27M-KO',24),21)),levels=c('H3.3K27M','H3.3K27M-KO'))
dim(ocr3_norm_melt)
ocr3_norm_melt<-na.omit(ocr3_norm_melt)
dim(ocr3_norm_melt)
ocr3_norm_sum<-summarySE(ocr3_norm_melt, measurevar="value", groupvars=c('time','cell_type3'))
ocr3_norm_sum$time<-as.numeric(levels(ocr3_norm_sum$time))[ocr3_norm_sum$time]

ocr3_norm_tab1<-ggplot(ocr3_norm_tab_grp, aes(x=cell_type3, y=basal_resp3_norm)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=5,binwidth = 0.01,aes(fill=cell_type3))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  stat_pvalue_manual(test_basal_resp3_norm, label = "p.signif",size=6)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='(pmol/min)/(ug/mL)')+ggtitle('Normalized Basal Respiration')+
  scale_y_continuous(limits = c(0,1.4),breaks=seq(0,1.4,by=0.2))+
  theme(legend.position="none",axis.title.y=element_text(size=14,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=12,face="bold"),plot.title = element_text(size= 14,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=12, color="black",face="bold"))

ocr3_norm_tab2<-ggplot(ocr3_norm_tab_grp, aes(x=cell_type3, y=max_resp3_norm)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=5,binwidth = 0.01,aes(fill=cell_type3))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  stat_pvalue_manual(test_max_resp3_norm, label = "p.signif",size=6)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='(pmol/min)/(ug/mL)')+ggtitle('Normalized Maximum Respiration')+
  scale_y_continuous(limits = c(0,2.2),breaks=seq(0,2,by=0.2))+
  theme(legend.position="none",axis.title.y=element_text(size=14,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=12,face="bold"),plot.title = element_text(size= 14,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=12, color="black",face="bold"))

ocr3_norm_tab3<-ggplot(ocr3_norm_tab_grp, aes(x=cell_type3, y=spare_cap_ocr3_norm_pct)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=4.5,binwidth = 1,aes(fill=cell_type3))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  stat_pvalue_manual(test_spare_cap_ocr3_norm_pct, label = "p.signif",size=6)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='Percentage(%)')+ggtitle('Normalized Spare Respiratory Capacity')+
  scale_y_continuous(limits = c(0,120),breaks=seq(0,120,by=20))+
  theme(legend.position="none",axis.title.y=element_text(size=14,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=12,face="bold"),plot.title = element_text(size= 14,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=12, color="black",face="bold"))

ocr3_norm_tab4<-ggplot(ocr3_norm_sum, aes(x=time, y=value,color=cell_type3,group=cell_type3))+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=1.5,size=1)+geom_line(size=1)+geom_point(size=2)+
  scale_linetype_manual(values=c(rep(c('solid'),2)))+labs(x='Time (min)',y='(pmol/min)/(ug/mL)')+
  scale_color_manual(values=c(rep(mycolors1[1:2]),rep(mycolors1[2],2)))+ggtitle('Normalized Oxygen Consumption Rate (OCR)')+
  theme(legend.position="bottom",legend.title=element_blank(),legend.text=element_text(hjust=0.5,size=12,face='bold')
        ,axis.title=element_text(size=16,face="bold"),axis.text = element_text(size=14)
        ,plot.title = element_text(size= 14,hjust = 0.5,face='bold'))

resp3_norm_grp<-ggarrange(ocr3_norm_tab1,NULL,ocr3_norm_tab2,NULL,NULL,NULL,ocr3_norm_tab3,NULL,ocr3_norm_tab4,ncol=3,nrow=3
                          ,widths=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1),heights=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1))+theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))

#Normalized ECAR
ecar3_norm_tab<-ecar3_norm[,-1]
basal_ecar3_norm<-ecar3_norm_tab[3,]
max_ecar3_norm<-ecar3_norm_tab[4,]
min_ecar3_norm<-ecar3_norm_tab[21,]
basal_glyco3_norm<-basal_ecar3_norm-min_ecar3_norm
max_glyco3_norm<-max_ecar3_norm-min_ecar3_norm
spare_cap_glyco3_norm<-max_glyco3_norm-basal_glyco3_norm
spare_cap_glyco3_norm_fd<-max_glyco3_norm/basal_glyco3_norm
spare_cap_glyco3_norm_pct<-((max_glyco3_norm-basal_glyco3_norm)/max_glyco3_norm)*100

glyco3_norm_tab_grp<-rbind(basal_glyco3_norm,max_glyco3_norm,spare_cap_glyco3_norm,spare_cap_glyco3_norm_pct,spare_cap_glyco3_norm_fd)
glyco3_norm_tab_grp<-as.data.frame(t(glyco3_norm_tab_grp))
glyco3_norm_tab_grp$cell_type3<-factor(c(rep('H3.3K27M',22),rep('H3.3K27M-KO',24)),levels=c('H3.3K27M','H3.3K27M-KO'))
glyco3_norm_tab_grp<-na.omit(glyco3_norm_tab_grp)
colnames(glyco3_norm_tab_grp)<-c('basal_glyco3_norm','max_glyco3_norm','spare_cap_glyco3_norm','spare_cap_glyco3_norm_pct','spare_cap_glyco3_norm_fd'
                                 ,'cell_type3')

test_basal_glyco3_norm<-glyco3_norm_tab_grp %>% wilcox_test(basal_glyco3_norm~cell_type3)%>% 
  add_significance() %>% add_xy_position(x='cell_type3') %>% 
  mutate(y.position = max(glyco3_norm_tab_grp$basal_glyco3_norm)+sd(glyco3_norm_tab_grp$basal_glyco3_norm)*2)
test_max_glyco3_norm<-glyco3_norm_tab_grp %>% wilcox_test(max_glyco3_norm~cell_type3)%>% 
  add_significance() %>% add_xy_position(x='cell_type3') %>% 
  mutate(y.position = max(glyco3_norm_tab_grp$max_glyco3_norm)+sd(glyco3_norm_tab_grp$max_glyco3_norm))
test_spare_cap_glyco3_norm_pct<-glyco3_norm_tab_grp %>% wilcox_test(spare_cap_glyco3_norm_pct~cell_type3)%>% 
  add_significance() %>% add_xy_position(x='cell_type3') %>% 
  mutate(y.position = max(glyco3_norm_tab_grp$spare_cap_glyco3_norm_pct)+sd(glyco3_norm_tab_grp$spare_cap_glyco3_norm_pct)*2)

#Linegraph
ecar3_norm_melt<-as.data.frame(reshape2::melt(ecar3_norm, id='time'))
ecar3_norm_melt$time<-round(ecar3_norm_melt$time,digits = 2)
ecar3_norm_melt$time<-factor(ecar3_norm_melt$time)
ecar3_norm_melt$cell_type3<-factor(c(rep(rep('H3.3K27M',22),21),rep(rep('H3.3K27M-KO',24),21)),levels=c('H3.3K27M','H3.3K27M-KO'))
dim(ecar3_norm_melt)
ecar3_norm_melt<-na.omit(ecar3_norm_melt)
dim(ecar3_norm_melt)
ecar3_norm_sum<-summarySE(ecar3_norm_melt, measurevar="value", groupvars=c('time','cell_type3'))
ecar3_norm_sum$time<-as.numeric(levels(ecar3_norm_sum$time))[ecar3_norm_sum$time]

glyco3_norm_tab1<-ggplot(glyco3_norm_tab_grp, aes(x=cell_type3, y=basal_glyco3_norm)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=2.5,binwidth = 0.012,aes(fill=cell_type3))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  stat_pvalue_manual(test_basal_glyco3_norm, label = "p.signif",size=6)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='(mpH/min)/(ug/mL)')+ggtitle('Basal Glycolysis')+
  scale_y_continuous(limits = c(0,max(max_glyco3_norm)+sd(max_glyco3_norm)),breaks=seq(0,1,by=0.2))+
  theme(legend.position="none",axis.title.y=element_text(size=14,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=12,face="bold"),plot.title = element_text(size= 14,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=12, color="black",face="bold"))

glyco3_norm_tab2<-ggplot(glyco3_norm_tab_grp, aes(x=cell_type3, y=max_glyco3_norm)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=2.5,binwidth = 0.012,aes(fill=cell_type3))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  stat_pvalue_manual(test_max_glyco3_norm, label = "p.signif",size=6)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='(mpH/min)/(ug/mL)')+ggtitle('Normalized Maximum Glycolysis')+
  scale_y_continuous(limits = c(0,max(max_glyco3_norm)+sd(max_glyco3_norm)),breaks=seq(0,1,by=0.2))+
  theme(legend.position="none",axis.title.y=element_text(size=14,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=12,face="bold"),plot.title = element_text(size= 14,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=12, color="black",face="bold"))

glyco3_norm_tab3<-ggplot(glyco3_norm_tab_grp, aes(x=cell_type3, y=spare_cap_glyco3_norm_pct)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=4.5,binwidth = 1,aes(fill=cell_type3))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  stat_pvalue_manual(test_spare_cap_glyco3_norm_pct, label = "p.signif",size=6)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='Percentage(%)')+ggtitle('Normalized Spare Glycolytic Capacity')+
  scale_y_continuous(limits = c(0,120),breaks=seq(0,120,by=20))+
  theme(legend.position="none",axis.title.y=element_text(size=14,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=12,face="bold"),plot.title = element_text(size= 14,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=12, color="black",face="bold"))

glyco3_norm_tab4<-ggplot(ecar3_norm_sum, aes(x=time,y=value,group=cell_type3,color=cell_type3))+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=1.5,size=1)+geom_line(size=1)+geom_point(size=2)+
  scale_linetype_manual(values=c(rep(c('solid'),2)))+labs(x='Time (min)',y='(mpH/min)/(ug/mL)')+
  scale_color_manual(values=mycolors1[1:2])+ggtitle('Normalized Extracellular Acidification Rate (ECAR)')+
  theme(legend.position="bottom",legend.title=element_blank(),legend.text=element_text(hjust=0.5,size=12,face='bold')
        ,axis.title=element_text(size=16,face="bold"),axis.text = element_text(size=14)
        ,plot.title = element_text(size= 14,hjust = 0.5,face='bold'))

glyco3_norm_grp<-ggarrange(glyco3_norm_tab1,NULL,glyco3_norm_tab2,NULL,NULL,NULL,glyco3_norm_tab3,NULL,glyco3_norm_tab4,ncol=3,nrow=3
                           ,widths=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1),heights=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1))+theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))

####Seahorse analysis 2----
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Rmisc))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(rstatix))

###Import excel file generated from wave (overview) and create OCR and ECAR datasets
sh_res4<-as.data.frame(read_xlsx('seahorse_assay2.xlsx'))
sh_res4[sh_res4<0]<-NA

#OCR
ocr4<-as.data.frame(sh_res4[c(1:24),c(1:52)])
ocr4<-ocr4[,colSums(is.na(ocr4))<nrow(ocr4)]
#Remove NA rows between samples
colnames(ocr4)<-c('time',paste0('H3.3K27M_rep_',c(1:22)),paste0('H3.3K27M_KO_rep_',c(1:24)))
row.names(ocr4)<-NULL
ocr4<-ocr4[c(-1:-3),]
row.names(ocr4)<-c(1:21)
ocr4<-mutate_all(ocr4, function(x) as.numeric(as.character(x)))
#Remove unused rows, rename rows and convert values to numeric
#ECAR
ecar4<-as.data.frame(sh_res4[c(28:51),c(1:52)])
ecar4<-ecar4[,colSums(is.na(ecar4))<nrow(ecar4)]
colnames(ecar4)<-c('time',paste0('H3.3K27M_rep_',c(1:22)),paste0('H3.3K27M_KO_rep_',c(1:24)))
ecar4<-ecar4[,-25] #Remove bad sample
row.names(ecar4)<-NULL
ecar4<-ecar4[c(-1:-3),]
row.names(ecar4)<-c(1:21)
ecar4<-mutate_all(ecar4, function(x) as.numeric(as.character(x)))

#Create color object
mycolors1<-c('red3','blue3','green4','darkmagenta','goldenrod4','darkorange','deeppink',
             'gray60','darkcyan','darkturquoise')

##Non-normalized data
#OCR
ocr4_tab<-ocr4[,-1]
#remove variable time for calculations and figure generation
non_mit4<-apply(ocr4_tab[c(12:21),],2,function(x){min(x[!is.na(x)])})
non_mit4<-ifelse(non_mit4==Inf,0,non_mit4)
basal_resp4<-ocr4_tab[3,]-non_mit4
max_resp4<-apply(ocr4_tab[c(8:11),],2,function(x){max(x[!is.na(x)])})-non_mit4
proton_leak<-apply(ocr4_tab[c(4:7),],2,function(x){max(x[!is.na(x)])})-non_mit4
spare_cap_ocr4<-max_resp4-basal_resp4
spare_cap_ocr4_fd<-max_resp4/basal_resp4
spare_cap_ocr4_pct<-((max_resp4-basal_resp4)/max_resp4)*100
atp_prod<-basal_resp4-proton_leak

ocr4_tab_grp<-rbind(basal_resp4,max_resp4,spare_cap_ocr4,spare_cap_ocr4_pct,spare_cap_ocr4_fd,atp_prod,non_mit4)
ocr4_tab_grp<-as.data.frame(t(ocr4_tab_grp))
ocr4_tab_grp$cell_type<-factor(c(rep('H3.3K27M',22),rep('H3.3K27M-KO',24)),levels=c('H3.3K27M','H3.3K27M-KO'))
colnames(ocr4_tab_grp)<-c('basal_resp4','max_resp4','spare_cap_ocr4','spare_cap_ocr4_pct','spare_cap_ocr4_fd','atp_prod','non_mit4'
                          ,'cell_type')

test_basal_resp4<-ocr4_tab_grp %>% wilcox_test(basal_resp4~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% mutate(y.position = max(basal_resp4)+sd(basal_resp4)*2)
test_max_resp4<-ocr4_tab_grp %>% wilcox_test(max_resp4~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% mutate(y.position = max(max_resp4)+sd(max_resp4))
test_spare_cap_ocr4_pct<-ocr4_tab_grp %>% wilcox_test(spare_cap_ocr4_pct~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% mutate(y.position = max(spare_cap_ocr4_pct)+sd(spare_cap_ocr4_pct)*2)

#Linegraph
ocr4_melt<-as.data.frame(reshape2::melt(ocr4, id='time'))
ocr4_melt$time<-round(ocr4_melt$time,digits = 2)
ocr4_melt$time<-factor(ocr4_melt$time)
ocr4_melt$cell_type<-factor(c(rep(rep('H3.3K27M',22),21),rep(rep('H3.3K27M-KO',24),21)),levels=c('H3.3K27M','H3.3K27M-KO'))
dim(ocr4_melt)
ocr4_melt<-na.omit(ocr4_melt)
dim(ocr4_melt)
ocr4_sum<-summarySE(ocr4_melt, measurevar="value", groupvars=c('time','cell_type'))
ocr4_sum$time<-as.numeric(levels(ocr4_sum$time))[ocr4_sum$time]

ocr4_tab1<-ggplot(ocr4_tab_grp, aes(x=cell_type, y=basal_resp4)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=5,binwidth = 1.5,aes(fill=cell_type))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  stat_pvalue_manual(test_basal_resp4, label = "p.signif",size=6)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='pmol/min')+ggtitle('Basal Respiration')+
  scale_y_continuous(limits = c(0, 120),breaks=seq(0,120,by=20))+
  theme(legend.position="none",axis.title.y=element_text(size=24,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=24,face="bold"),plot.title = element_text(size=30,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=26, color="black",face="bold"))

ocr4_tab2<-ggplot(ocr4_tab_grp, aes(x=cell_type, y=max_resp4)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=5,binwidth = 1.5,aes(fill=cell_type))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='pmol/min')+ggtitle('Maximum Respiration')+
  stat_pvalue_manual(test_max_resp4, label = "p.signif",size=6)+
  scale_y_continuous(limits = c(0, 120),breaks=seq(0,120,by=20))+
  theme(legend.position="none",axis.title.y=element_text(size=24,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=24,face="bold"),plot.title = element_text(size=30,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=26, color="black",face="bold"))

ocr4_tab3<-ggplot(ocr4_tab_grp, aes(x=cell_type, y=spare_cap_ocr4_pct)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=5,binwidth = 1.5,aes(fill=cell_type))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='Percentage(%)')+ggtitle('Spare Respiratory Capacity')+
  stat_pvalue_manual(test_spare_cap_ocr4_pct, label = "p.signif",size=6)+
  scale_y_continuous(limits = c(0,100),breaks=seq(0,100,by=20))+
  theme(legend.position="none",axis.title.y=element_text(size=24,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=24,face="bold"),plot.title = element_text(size=30,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=26, color="black",face="bold"))

ocr4_tab4<-ggplot(ocr4_sum, aes(x=time, y=value, group=cell_type,color=cell_type))+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=1.5,size=1)+geom_line(size=1)+geom_point(size=2)+
  scale_linetype_manual(values=c(rep(c('solid'),2)))+labs(x='Time (min)',y='pmol/min')+
  scale_color_manual(values=mycolors1[1:2])+ggtitle('Oxygen Consumption Rate (OCR)')+
  theme(legend.position="bottom",legend.title=element_blank(),legend.text=element_text(hjust=0.5,size=20,face='bold')
        ,axis.title=element_text(size=24,face="bold"),axis.text = element_text(size=22)
        ,plot.title = element_text(size= 26,hjust = 0.5,face='bold'),legend.key.size = unit(3.4,"line",'point'))

resp4_grp<-ggarrange(ocr4_tab1,NULL,ocr4_tab2,NULL,NULL,NULL,ocr4_tab3,NULL,ocr4_tab4,ncol=3,nrow=3
                     ,widths=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1),heights=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1))+theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))

#ECAR
ecar4_tab<-ecar4[,-1]
basal_ecar4<-ecar4_tab[3,]
max_ecar4<-ecar4_tab[4,]
min_ecar4<-ecar4_tab[21,]
basal_glyco4<-basal_ecar4-min_ecar4
max_glyco4<-max_ecar4-min_ecar4
spare_cap_glyco4<-max_glyco4-basal_glyco4
spare_cap_glyco4_fd<-max_glyco4/basal_glyco4
spare_cap_glyco4_pct<-((max_glyco4-basal_glyco4)/max_glyco4)*100

glyco4_tab_grp<-rbind(basal_glyco4,max_glyco4,spare_cap_glyco4,spare_cap_glyco4_pct,spare_cap_glyco4_fd)
glyco4_tab_grp<-as.data.frame(t(glyco4_tab_grp))
glyco4_tab_grp$cell_type<-factor(c(rep('H3.3K27M',22),rep('H3.3K27M-KO',23)),levels=c('H3.3K27M','H3.3K27M-KO'))
glyco4_tab_grp<-na.omit(glyco4_tab_grp)
colnames(glyco4_tab_grp)<-c('basal_glyco4','max_glyco4','spare_cap_glyco4','spare_cap_glyco4_pct','spare_cap_glyco4_fd','cell_type')

#Linegraph
ecar4_melt<-as.data.frame(reshape2::melt(ecar4, id='time'))
ecar4_melt$time<-round(ecar4_melt$time,digits = 2)
ecar4_melt$time<-factor(ecar4_melt$time)
ecar4_melt$cell_type<-factor(c(rep(rep('H3.3K27M',22),21),rep(rep('H3.3K27M-KO',23),21)),levels=c('H3.3K27M','H3.3K27M-KO'))
dim(ecar4_melt)
ecar4_melt<-na.omit(ecar4_melt)
dim(ecar4_melt)
ecar4_sum<-summarySE(ecar4_melt, measurevar="value", groupvars=c('time','cell_type'))
ecar4_sum$time<-as.numeric(levels(ecar4_sum$time))[ecar4_sum$time]

test_basal_glyco4<-glyco4_tab_grp %>% wilcox_test(basal_glyco4~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% mutate(y.position = max(glyco4_tab_grp$basal_glyco4)+sd(glyco4_tab_grp$basal_glyco4)*2)
test_max_glyco4<-glyco4_tab_grp %>% wilcox_test(max_glyco4~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% mutate(y.position = max(glyco4_tab_grp$max_glyco4)+sd(glyco4_tab_grp$max_glyco4))
test_spare_cap_glyco4_pct<-glyco4_tab_grp %>% wilcox_test(spare_cap_glyco4_pct~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% mutate(y.position = max(glyco4_tab_grp$spare_cap_glyco4_pct)+sd(glyco4_tab_grp$spare_cap_glyco4_pct)*2)


glyco4_tab1<-ggplot(glyco4_tab_grp, aes(x=cell_type, y=basal_glyco4)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=2.5,binwidth = 2.5,aes(fill=cell_type))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='mpH/min')+ggtitle('Basal Glycolysis')+
  stat_pvalue_manual(test_basal_glyco4, label = "p.signif",size=6)+
  scale_y_continuous(limits = c(0,120),breaks=seq(0,100,by=20))+
  theme(legend.position="none",axis.title.y=element_text(size=24,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=24,face="bold"),plot.title = element_text(size=30,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=26, color="black",face="bold"))

glyco4_tab2<-ggplot(glyco4_tab_grp, aes(x=cell_type, y=max_glyco4)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=2.5,binwidth = 2.5,aes(fill=cell_type))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='mpH/min')+ggtitle('Maximum Glycolysis')+
  stat_pvalue_manual(test_max_glyco4, label = "p.signif",size=6)+
  scale_y_continuous(limits = c(0,120),breaks=seq(0,100,by=20))+
  theme(legend.position="none",axis.title.y=element_text(size=24,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=24,face="bold"),plot.title = element_text(size=30,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=26, color="black",face="bold"))

glyco4_tab3<-ggplot(glyco4_tab_grp, aes(x=cell_type, y=spare_cap_glyco4_pct)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=5,binwidth = 1,aes(fill=cell_type))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='Percentage(%)')+ggtitle('Spare Glycolytic Capacity')+
  stat_pvalue_manual(test_spare_cap_glyco4_pct, label = "p.signif",size=6)+
  scale_y_continuous(limits = c(0,100),breaks=seq(0,100,by=20))+
  theme(legend.position="none",axis.title.y=element_text(size=24,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=24,face="bold"),plot.title = element_text(size=30,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=26, color="black",face="bold"))

glyco4_tab4<-ggplot(ecar4_sum, aes(x=time, y=value, group=cell_type,color=cell_type))+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=1.5,size=1)+geom_line(size=1)+geom_point(size=2)+
  scale_linetype_manual(values=c(rep(c('solid'),2)))+labs(x='Time (min)',y='mpH/min')+
  scale_color_manual(values=mycolors1[1:2])+ggtitle('Extracellular Acidification Rate (ECAR)')+
  theme(legend.position="bottom",legend.title=element_blank(),legend.text=element_text(hjust=0.5,size=20,face='bold')
        ,axis.title=element_text(size=24,face="bold"),axis.text = element_text(size=22)
        ,plot.title = element_text(size= 26,hjust = 0.5,face='bold'),legend.key.size = unit(3.4,"line",'point'))

glyco4_grp<-ggarrange(glyco4_tab1,NULL,glyco4_tab2,NULL,NULL,NULL,glyco4_tab3,NULL,glyco4_tab4,ncol=3,nrow=3
                      ,widths=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1),heights=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1))+theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))

###Data Normalized by protein
bca_res_bt4<-as.data.frame(read_xlsx('BCA_seahorse2.xlsx'))
bca_res_bt4<-bca_res_bt4[c(22,23),c(2:49)]
#Create column with time variable for further calculations
time<-as.matrix(c('time',1))
ptn4<-as.data.frame(cbind(time,bca_res_bt4))
ptn4<-as.data.frame(t(ptn4[,colSums(is.na(ptn4))<nrow(ptn4)]))
rownames(ptn4)<-NULL
colnames(ptn4)<-c('names','protein')
ptn4_ocr4<-ptn4
ptn4_ecar4<-ptn4[-25,] #Remove bad sample

ocr4_norm<-as.data.frame(t(ocr4))
ocr4_norm<-cbind(ocr4_norm,ptn4_ocr4)
ocr4_norm<-ocr4_norm %>% dplyr::select('names','protein',everything())
#check if names match
ocr4_norm<-ocr4_norm[,-1]
ocr4_norm<-mutate_all(ocr4_norm, function(x) as.numeric(as.character(x)))
ocr4_norm<- ocr4_norm %>% mutate_at(vars(-protein),funs(./protein))%>% dplyr::select(c(1:22))
ocr4_norm<-ocr4_norm[,-1]
ocr4_norm<-as.data.frame(t(ocr4_norm))

ecar4_norm<-as.data.frame(t(ecar4))
ecar4_norm<-cbind(ecar4_norm,ptn4_ecar4)
ecar4_norm<-ecar4_norm %>% dplyr::select('names','protein',everything())
#View(ecar4_norm) #check if names match
ecar4_norm<-ecar4_norm[,-1]
ecar4_norm<-mutate_all(ecar4_norm, function(x) as.numeric(as.character(x)))
ecar4_norm<- ecar4_norm %>% mutate_at(vars(-protein),funs(./protein))%>% dplyr::select(c(1:22))
ecar4_norm<-ecar4_norm[,-1]
ecar4_norm<-as.data.frame(t(ecar4_norm))

#Create subsets, treat data and generate figures
#Normalized OCR
ocr4_norm_tab<-ocr4_norm[,-1]
non_mit4_norm<-apply(ocr4_norm_tab[c(12:21),],2,function(x){min(x[!is.na(x)])})
non_mit4_norm<-ifelse(non_mit4_norm==Inf,0,non_mit4_norm)
basal_resp4_norm<-ocr4_norm_tab[3,]-non_mit4_norm
max_resp4_norm<-apply(ocr4_norm_tab[c(8:11),],2,function(x){max(x[!is.na(x)])})-non_mit4_norm
proton_leak_norm<-apply(ocr4_norm_tab[c(4:7),],2,function(x){max(x[!is.na(x)])})-non_mit4_norm
spare_cap_ocr4_norm<-max_resp4_norm-basal_resp4_norm
spare_cap_ocr4_norm_fd<-max_resp4_norm/basal_resp4_norm
spare_cap_ocr4_norm_pct<-((max_resp4_norm-basal_resp4_norm)/max_resp4_norm)*100
atp_prod_norm<-basal_resp4_norm-proton_leak_norm

ocr4_norm_tab_grp<-rbind(basal_resp4_norm,max_resp4_norm,spare_cap_ocr4_norm,spare_cap_ocr4_norm_pct
                         ,spare_cap_ocr4_norm_fd,atp_prod_norm,non_mit4_norm)
ocr4_norm_tab_grp<-as.data.frame(t(ocr4_norm_tab_grp))
ocr4_norm_tab_grp$cell_type<-factor(c(rep('H3.3K27M',22),rep('H3.3K27M-KO',24)),levels=c('H3.3K27M','H3.3K27M-KO'))
colnames(ocr4_norm_tab_grp)<-c('basal_resp4_norm','max_resp4_norm','spare_cap_ocr4_norm','spare_cap_ocr4_norm_pct'
                               ,'spare_cap_ocr4_norm_fd','atp_prod_norm','non_mit4_norm','cell_type')

test_basal_resp4_norm<-ocr4_norm_tab_grp %>% wilcox_test(basal_resp4_norm~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% mutate(y.position = max(basal_resp4_norm)+sd(basal_resp4_norm)*2)
test_max_resp4_norm<-ocr4_norm_tab_grp %>% wilcox_test(max_resp4_norm~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% mutate(y.position = max(max_resp4_norm)+sd(max_resp4_norm))
test_spare_cap_ocr4_norm_pct<-ocr4_norm_tab_grp %>% wilcox_test(spare_cap_ocr4_norm_pct~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% mutate(y.position = max(spare_cap_ocr4_norm_pct)+sd(spare_cap_ocr4_norm_pct)*2)

#Linegraph
ocr4_norm_melt<-as.data.frame(reshape2::melt(ocr4_norm, id='time'))
ocr4_norm_melt$time<-round(ocr4_norm_melt$time,digits = 2)
ocr4_norm_melt$time<-factor(ocr4_norm_melt$time)
ocr4_norm_melt$cell_type<-factor(c(rep(rep('H3.3K27M',22),21),rep(rep('H3.3K27M-KO',24),21)),levels=c('H3.3K27M','H3.3K27M-KO'))
dim(ocr4_norm_melt)
ocr4_norm_melt<-na.omit(ocr4_norm_melt)
dim(ocr4_norm_melt)
ocr4_norm_sum<-summarySE(ocr4_norm_melt, measurevar="value", groupvars=c('time','cell_type'))
ocr4_norm_sum$time<-as.numeric(levels(ocr4_norm_sum$time))[ocr4_norm_sum$time]

ocr4_norm_tab1<-ggplot(ocr4_norm_tab_grp, aes(x=cell_type, y=basal_resp4_norm)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=5,binwidth = 0.01,aes(fill=cell_type))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  stat_pvalue_manual(test_basal_resp4_norm, label = "p.signif",size=6)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='(pmol/min)/(ug/mL)')+ggtitle('Normalized Basal Respiration')+
  scale_y_continuous(limits = c(0,1.4),breaks=seq(0,1.4,by=0.2))+
  theme(legend.position="none",axis.title.y=element_text(size=14,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=12,face="bold"),plot.title = element_text(size= 14,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=12, color="black",face="bold"))

ocr4_norm_tab2<-ggplot(ocr4_norm_tab_grp, aes(x=cell_type, y=max_resp4_norm)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=5,binwidth = 0.01,aes(fill=cell_type))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  stat_pvalue_manual(test_max_resp4_norm, label = "p.signif",size=6)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='(pmol/min)/(ug/mL)')+ggtitle('Normalized Maximum Respiration')+
  scale_y_continuous(limits = c(0,2.2),breaks=seq(0,2,by=0.2))+
  theme(legend.position="none",axis.title.y=element_text(size=14,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=12,face="bold"),plot.title = element_text(size= 14,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=12, color="black",face="bold"))

ocr4_norm_tab3<-ggplot(ocr4_norm_tab_grp, aes(x=cell_type, y=spare_cap_ocr4_norm_pct)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=4.5,binwidth = 1,aes(fill=cell_type))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  stat_pvalue_manual(test_spare_cap_ocr4_norm_pct, label = "p.signif",size=6)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='Percentage(%)')+ggtitle('Normalized Spare Respiratory Capacity')+
  scale_y_continuous(limits = c(0,120),breaks=seq(0,120,by=20))+
  theme(legend.position="none",axis.title.y=element_text(size=14,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=12,face="bold"),plot.title = element_text(size= 14,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=12, color="black",face="bold"))

ocr4_norm_tab4<-ggplot(ocr4_norm_sum, aes(x=time, y=value,color=cell_type,group=cell_type))+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=1.5,size=1)+geom_line(size=1)+geom_point(size=2)+
  scale_linetype_manual(values=c(rep(c('solid'),2)))+labs(x='Time (min)',y='(pmol/min)/(ug/mL)')+
  scale_color_manual(values=c(rep(mycolors1[1:2]),rep(mycolors1[2],2)))+ggtitle('Normalized Oxygen Consumption Rate (OCR)')+
  theme(legend.position="bottom",legend.title=element_blank(),legend.text=element_text(hjust=0.5,size=12,face='bold')
        ,axis.title=element_text(size=16,face="bold"),axis.text = element_text(size=14)
        ,plot.title = element_text(size= 14,hjust = 0.5,face='bold'))

resp4_norm_grp<-ggarrange(ocr4_norm_tab1,NULL,ocr4_norm_tab2,NULL,NULL,NULL,ocr4_norm_tab3,NULL,ocr4_norm_tab4,ncol=3,nrow=3
                          ,widths=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1),heights=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1))+theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))

#Normalized ECAR
ecar4_norm_tab<-ecar4_norm[,-1]
basal_ecar4_norm<-ecar4_norm_tab[3,]
max_ecar4_norm<-ecar4_norm_tab[4,]
min_ecar4_norm<-ecar4_norm_tab[21,]
basal_glyco4_norm<-basal_ecar4_norm-min_ecar4_norm
max_glyco4_norm<-max_ecar4_norm-min_ecar4_norm
spare_cap_glyco4_norm<-max_glyco4_norm-basal_glyco4_norm
spare_cap_glyco4_norm_fd<-max_glyco4_norm/basal_glyco4_norm
spare_cap_glyco4_norm_pct<-((max_glyco4_norm-basal_glyco4_norm)/max_glyco4_norm)*100

glyco4_norm_tab_grp<-rbind(basal_glyco4_norm,max_glyco4_norm,spare_cap_glyco4_norm,spare_cap_glyco4_norm_pct,spare_cap_glyco4_norm_fd)
glyco4_norm_tab_grp<-as.data.frame(t(glyco4_norm_tab_grp))
glyco4_norm_tab_grp$cell_type<-factor(c(rep('H3.3K27M',22),rep('H3.3K27M-KO',23)),levels=c('H3.3K27M','H3.3K27M-KO'))
glyco4_norm_tab_grp<-na.omit(glyco4_norm_tab_grp)
colnames(glyco4_norm_tab_grp)<-c('basal_glyco4_norm','max_glyco4_norm','spare_cap_glyco4_norm','spare_cap_glyco4_norm_pct','spare_cap_glyco4_norm_fd'
                                 ,'cell_type')

test_basal_glyco4_norm<-glyco4_norm_tab_grp %>% wilcox_test(basal_glyco4_norm~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% 
  mutate(y.position = max(glyco4_norm_tab_grp$basal_glyco4_norm)+sd(glyco4_norm_tab_grp$basal_glyco4_norm)*2)
test_max_glyco4_norm<-glyco4_norm_tab_grp %>% wilcox_test(max_glyco4_norm~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% 
  mutate(y.position = max(glyco4_norm_tab_grp$max_glyco4_norm)+sd(glyco4_norm_tab_grp$max_glyco4_norm))
test_spare_cap_glyco4_norm_pct<-glyco4_norm_tab_grp %>% wilcox_test(spare_cap_glyco4_norm_pct~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% 
  mutate(y.position = max(glyco4_norm_tab_grp$spare_cap_glyco4_norm_pct)+sd(glyco4_norm_tab_grp$spare_cap_glyco4_norm_pct)*2)

#Linegraph
ecar4_norm_melt<-as.data.frame(reshape2::melt(ecar4_norm, id='time'))
ecar4_norm_melt$time<-round(ecar4_norm_melt$time,digits = 2)
ecar4_norm_melt$time<-factor(ecar4_norm_melt$time)
ecar4_norm_melt$cell_type<-factor(c(rep(rep('H3.3K27M',22),21),rep(rep('H3.3K27M-KO',23),21)),levels=c('H3.3K27M','H3.3K27M-KO'))
dim(ecar4_norm_melt)
ecar4_norm_melt<-na.omit(ecar4_norm_melt)
dim(ecar4_norm_melt)
ecar4_norm_sum<-summarySE(ecar4_norm_melt, measurevar="value", groupvars=c('time','cell_type'))
ecar4_norm_sum$time<-as.numeric(levels(ecar4_norm_sum$time))[ecar4_norm_sum$time]

glyco4_norm_tab1<-ggplot(glyco4_norm_tab_grp, aes(x=cell_type, y=basal_glyco4_norm)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=2.5,binwidth = 0.012,aes(fill=cell_type))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  stat_pvalue_manual(test_basal_glyco4_norm, label = "p.signif",size=6)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='(mpH/min)/(ug/mL)')+ggtitle('Basal Glycolysis')+
  scale_y_continuous(limits = c(0,max(max_glyco4_norm)+sd(max_glyco4_norm)),breaks=seq(0,1,by=0.2))+
  theme(legend.position="none",axis.title.y=element_text(size=14,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=12,face="bold"),plot.title = element_text(size= 14,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=12, color="black",face="bold"))

glyco4_norm_tab2<-ggplot(glyco4_norm_tab_grp, aes(x=cell_type, y=max_glyco4_norm)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=2.5,binwidth = 0.012,aes(fill=cell_type))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  stat_pvalue_manual(test_max_glyco4_norm, label = "p.signif",size=6)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='(mpH/min)/(ug/mL)')+ggtitle('Normalized Maximum Glycolysis')+
  scale_y_continuous(limits = c(0,max(max_glyco4_norm)+sd(max_glyco4_norm)),breaks=seq(0,1,by=0.2))+
  theme(legend.position="none",axis.title.y=element_text(size=14,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=12,face="bold"),plot.title = element_text(size= 14,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=12, color="black",face="bold"))

glyco4_norm_tab3<-ggplot(glyco4_norm_tab_grp, aes(x=cell_type, y=spare_cap_glyco4_norm_pct)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=4.5,binwidth = 1,aes(fill=cell_type))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  stat_pvalue_manual(test_spare_cap_glyco4_norm_pct, label = "p.signif",size=6)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='Percentage(%)')+ggtitle('Normalized Spare Glycolytic Capacity')+
  scale_y_continuous(limits = c(0,120),breaks=seq(0,120,by=20))+
  theme(legend.position="none",axis.title.y=element_text(size=14,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=12,face="bold"),plot.title = element_text(size= 14,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=12, color="black",face="bold"))

glyco4_norm_tab4<-ggplot(ecar4_norm_sum, aes(x=time,y=value,group=cell_type,color=cell_type))+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=1.5,size=1)+geom_line(size=1)+geom_point(size=2)+
  scale_linetype_manual(values=c(rep(c('solid'),2)))+labs(x='Time (min)',y='(mpH/min)/(ug/mL)')+
  scale_color_manual(values=mycolors1[1:2])+ggtitle('Normalized Extracellular Acidification Rate (ECAR)')+
  theme(legend.position="bottom",legend.title=element_blank(),legend.text=element_text(hjust=0.5,size=12,face='bold')
        ,axis.title=element_text(size=16,face="bold"),axis.text = element_text(size=14)
        ,plot.title = element_text(size= 14,hjust = 0.5,face='bold'))

glyco4_norm_grp<-ggarrange(glyco4_norm_tab1,NULL,glyco4_norm_tab2,NULL,NULL,NULL,glyco4_norm_tab3,NULL,glyco4_norm_tab4,ncol=3,nrow=3
                           ,widths=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1),heights=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1))+theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))

####Group analysis of normalized results----

#Non-normalized data
ocr5<-ocr3
colnames(ocr5)<-c('time',paste0('H3.3K27M_rep_',c(23:44)),paste0('H3.3K27M-KO_rep_',c(25:48)))
ocr5<-cbind(ocr4,ocr5)
ocr5<-ocr5[,-48]

#Create subsets, treat data and generate figures
#OCR
ocr5_tab<-ocr5[,-1]
non_mit5<-apply(ocr5_tab[c(12:21),],2,function(x){min(x[!is.na(x)])})
non_mit5<-ifelse(non_mit5==Inf,0,non_mit5)
basal_resp5<-ocr5_tab[3,]-non_mit5
max_resp5<-apply(ocr5_tab[c(8:11),],2,function(x){max(x[!is.na(x)])})-non_mit5
proton_leak<-apply(ocr5_tab[c(4:7),],2,function(x){max(x[!is.na(x)])})-non_mit5
spare_cap_ocr5<-max_resp5-basal_resp5
spare_cap_ocr5_fd<-max_resp5/basal_resp5
spare_cap_ocr5_pct<-((max_resp5-basal_resp5)/max_resp5)*100
atp_prod<-basal_resp5-proton_leak

ocr5_tab_grp<-rbind(basal_resp5,max_resp5,spare_cap_ocr5,spare_cap_ocr5_pct
                    ,spare_cap_ocr5_fd,atp_prod,non_mit5)
ocr5_tab_grp<-as.data.frame(t(ocr5_tab_grp))
ocr5_tab_grp$cell_type<-factor(c(rep('H3.3K27M',22),rep('H3.3K27M-KO',24),rep('H3.3K27M',22),rep('H3.3K27M-KO',24)),levels=c('H3.3K27M','H3.3K27M-KO'))
colnames(ocr5_tab_grp)<-c('basal_resp5','max_resp5','spare_cap_ocr5','spare_cap_ocr5_pct'
                          ,'spare_cap_ocr5_fd','atp_prod','non_mit5','cell_type')

test_basal_resp5<-ocr5_tab_grp %>% wilcox_test(basal_resp5~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% mutate(y.position = 115)
test_max_resp5<-ocr5_tab_grp %>% wilcox_test(max_resp5~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% mutate(y.position = 115)
test_spare_cap_ocr5_pct<-ocr5_tab_grp %>% wilcox_test(spare_cap_ocr5_pct~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% mutate(y.position = 115)

#Linegraph
ocr5_melt<-as.data.frame(reshape2::melt(ocr5, id='time'))
ocr5_melt$time<-round(ocr5_melt$time,digits = 2)
ocr5_melt$time<-factor(ocr5_melt$time)
ocr5_melt$cell_type<-factor(c(rep(rep('H3.3K27M',22),21),rep(rep('H3.3K27M-KO',24),21),rep(rep('H3.3K27M',22),21)
                              ,rep(rep('H3.3K27M-KO',24),21)),levels=c('H3.3K27M','H3.3K27M-KO'))
dim(ocr5_melt)
ocr5_melt<-na.omit(ocr5_melt)
dim(ocr5_melt)
ocr5_sum<-summarySE(ocr5_melt, measurevar="value", groupvars=c('time','cell_type'))
ocr5_sum$time<-as.numeric(levels(ocr5_sum$time))[ocr5_sum$time]

ocr5_tab1<-ggplot(ocr5_tab_grp, aes(x=cell_type, y=basal_resp5)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=4.5,binwidth = 1,aes(fill=cell_type))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  stat_pvalue_manual(test_basal_resp5, label = "p.signif",size=6)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='pmol/min')+ggtitle('Basal Respiration')+
  scale_y_continuous(limits = c(0,120),breaks=seq(0,110,by=20))+
  theme(legend.position="none",axis.title.y=element_text(size=14,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=12,face="bold"),plot.title = element_text(size= 14,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=12, color="black",face="bold"))

ocr5_tab2<-ggplot(ocr5_tab_grp, aes(x=cell_type, y=max_resp5)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=4.5,binwidth = 1,aes(fill=cell_type))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  stat_pvalue_manual(test_max_resp5, label = "p.signif",size=6)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='pmol/min')+ggtitle('Maximum Respiration')+
  scale_y_continuous(limits = c(0,120),breaks=seq(0,110,by=20))+
  theme(legend.position="none",axis.title.y=element_text(size=14,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=12,face="bold"),plot.title = element_text(size= 14,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=12, color="black",face="bold"))

ocr5_tab3<-ggplot(ocr5_tab_grp, aes(x=cell_type, y=spare_cap_ocr5_pct)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=4.5,binwidth = 1,aes(fill=cell_type))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  stat_pvalue_manual(test_spare_cap_ocr5_pct, label = "p.signif",size=6)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='Percentage(%)')+ggtitle('Spare Respiratory Capacity')+
  scale_y_continuous(limits = c(0,120),breaks=seq(0,100,by=25))+
  theme(legend.position="none",axis.title.y=element_text(size=14,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=12,face="bold"),plot.title = element_text(size= 14,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=12, color="black",face="bold"))

ocr5_tab4<-ggplot(ocr5_sum, aes(x=time, y=value,color=cell_type,group=cell_type))+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=1.5,size=1)+geom_line(size=1)+geom_point(size=2)+
  scale_linetype_manual(values=c(rep(c('solid'),2)))+labs(x='Time (min)',y='pmol/min')+
  scale_y_continuous(limits = c(0,75),breaks=seq(0,75,by=25))+
  scale_color_manual(values=c(rep(mycolors1[1:2]),rep(mycolors1[2],2)))+ggtitle('Oxygen Consumption Rate (OCR)')+
  theme(legend.position="bottom",legend.title=element_blank(),legend.text=element_text(hjust=0.5,size=12,face='bold')
        ,axis.title=element_text(size=16,face="bold"),axis.text = element_text(size=14)
        ,plot.title = element_text(size= 14,hjust = 0.5,face='bold'))

resp5_grp<-ggarrange(ocr5_tab1,NULL,ocr5_tab2,NULL,NULL,NULL,ocr5_tab3,NULL,ocr5_tab4,ncol=3,nrow=3
                     ,widths=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1),heights=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1))+theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))

#ECAR
ecar5<-ecar3
colnames(ecar5)<-c('time',paste0('H3.3K27M_rep_',c(23:44)),paste0('H3.3K27M-KO_rep_',c(25:48)))
ecar5<-cbind(ecar4,ecar5)
ecar5<-ecar5[,-47]

ecar5_tab<-ecar5[,-1]
basal_ecar5<-ecar5_tab[3,]
max_ecar5<-ecar5_tab[4,]
min_ecar5<-ecar5_tab[21,]
basal_glyco5<-basal_ecar5-min_ecar5
max_glyco5<-max_ecar5-min_ecar5
spare_cap_glyco5<-max_glyco5-basal_glyco5
spare_cap_glyco5_fd<-max_glyco5/basal_glyco5
spare_cap_glyco5_pct<-((max_glyco5-basal_glyco5)/max_glyco5)*100

glyco5_tab_grp<-rbind(basal_glyco5,max_glyco5,spare_cap_glyco5,spare_cap_glyco5_pct,spare_cap_glyco5_fd)
glyco5_tab_grp<-as.data.frame(t(glyco5_tab_grp))
glyco5_tab_grp$cell_type<-factor(c(rep('H3.3K27M',22),rep('H3.3K27M-KO',23),rep('H3.3K27M',22),rep('H3.3K27M-KO',24)),levels=c('H3.3K27M','H3.3K27M-KO'))
glyco5_tab_grp<-na.omit(glyco5_tab_grp)
colnames(glyco5_tab_grp)<-c('basal_glyco5','max_glyco5','spare_cap_glyco5','spare_cap_glyco5_pct','spare_cap_glyco5_fd'
                            ,'cell_type')

test_basal_glyco5<-glyco5_tab_grp %>% wilcox_test(basal_glyco5~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% mutate(y.position = 115)
test_max_glyco5<-glyco5_tab_grp %>% wilcox_test(max_glyco5~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% mutate(y.position = 115)
test_spare_cap_glyco5_pct<-glyco5_tab_grp %>% wilcox_test(spare_cap_glyco5_pct~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% mutate(y.position = 115)

#Linegraph
ecar5_melt<-as.data.frame(reshape2::melt(ecar5, id='time'))
ecar5_melt$time<-round(ecar5_melt$time,digits = 2)
ecar5_melt$time<-factor(ecar5_melt$time)
ecar5_melt$cell_type<-factor(c(rep(rep('H3.3K27M',22),21),rep(rep('H3.3K27M-KO',24),21),rep(rep('H3.3K27M',22),21),rep(rep('H3.3K27M-KO',23),21)),levels=c('H3.3K27M','H3.3K27M-KO'))
dim(ecar5_melt)
ecar5_melt<-na.omit(ecar5_melt)
dim(ecar5_melt)
ecar5_sum<-summarySE(ecar5_melt, measurevar="value", groupvars=c('time','cell_type'))
ecar5_sum$time<-as.numeric(levels(ecar5_sum$time))[ecar5_sum$time]

glyco5_tab1<-ggplot(glyco5_tab_grp, aes(x=cell_type, y=basal_glyco5)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=4,binwidth = 1,aes(fill=cell_type))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  stat_pvalue_manual(test_basal_glyco5, label = "p.signif",size=6)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='mpH/min')+ggtitle('Basal Glycolysis')+
  scale_y_continuous(limits = c(0,120),breaks=seq(0,100,by=25))+
  theme(legend.position="none",axis.title.y=element_text(size=14,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=12,face="bold"),plot.title = element_text(size= 14,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=12, color="black",face="bold"))

glyco5_tab2<-ggplot(glyco5_tab_grp, aes(x=cell_type, y=max_glyco5)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=4,binwidth = 1,aes(fill=cell_type))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  stat_pvalue_manual(test_max_glyco5, label = "p.signif",size=6)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='mpH/min')+ggtitle(' Maximum Glycolysis')+
  scale_y_continuous(limits = c(0,120),breaks=seq(0,100,by=25))+
  theme(legend.position="none",axis.title.y=element_text(size=14,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=12,face="bold"),plot.title = element_text(size= 14,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=12, color="black",face="bold"))

glyco5_tab3<-ggplot(glyco5_tab_grp, aes(x=cell_type, y=spare_cap_glyco5_pct)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=4,binwidth = 1,aes(fill=cell_type))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  stat_pvalue_manual(test_spare_cap_glyco5_pct, label = "p.signif",size=6)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='Percentage(%)')+ggtitle('Spare Glycolytic Capacity')+
  scale_y_continuous(limits = c(0,120),breaks=seq(0,100,by=25))+
  theme(legend.position="none",axis.title.y=element_text(size=14,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=12,face="bold"),plot.title = element_text(size= 14,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=12, color="black",face="bold"))

glyco5_tab4<-ggplot(ecar5_sum, aes(x=time,y=value,group=cell_type,color=cell_type))+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=1.5,size=1)+geom_line(size=1)+geom_point(size=2)+
  scale_linetype_manual(values=c(rep(c('solid'),2)))+labs(x='Time (min)',y='mpH/min')+
  scale_color_manual(values=mycolors1[1:2])+ggtitle('Extracellular Acidification Rate (ECAR)')+
  scale_y_continuous(limits = c(0,75),breaks=seq(0,75,by=25))+
  theme(legend.position="bottom",legend.title=element_blank(),legend.text=element_text(hjust=0.5,size=12,face='bold')
        ,axis.title=element_text(size=16,face="bold"),axis.text = element_text(size=14)
        ,plot.title = element_text(size= 14,hjust = 0.5,face='bold'))

glyco5_grp<-ggarrange(glyco5_tab1,NULL,glyco5_tab2,NULL,NULL,NULL,glyco5_tab3,NULL,glyco5_tab4,ncol=3,nrow=3
                      ,widths=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1),heights=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1))+theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))

#Normalized data
ocr5_norm<-ocr3_norm
colnames(ocr5_norm)<-c('time',paste0('H3.3K27M_rep_',c(23:44)),paste0('H3.3K27M-KO_rep_',c(25:48)))
ocr5_norm<-cbind(ocr4_norm,ocr5_norm)
ocr5_norm<-ocr5_norm[,-48]

#Create subsets, treat data and generate figures
#Normalized OCR
ocr5_norm_tab<-ocr5_norm[,-1]
non_mit5_norm<-apply(ocr5_norm_tab[c(12:21),],2,function(x){min(x[!is.na(x)])})
non_mit5_norm<-ifelse(non_mit5_norm==Inf,0,non_mit5_norm)
basal_resp5_norm<-ocr5_norm_tab[3,]-non_mit5_norm
max_resp5_norm<-apply(ocr5_norm_tab[c(8:11),],2,function(x){max(x[!is.na(x)])})-non_mit5_norm
proton_leak_norm<-apply(ocr5_norm_tab[c(4:7),],2,function(x){max(x[!is.na(x)])})-non_mit5_norm
spare_cap_ocr5_norm<-max_resp5_norm-basal_resp5_norm
spare_cap_ocr5_norm_fd<-max_resp5_norm/basal_resp5_norm
spare_cap_ocr5_norm_pct<-((max_resp5_norm-basal_resp5_norm)/max_resp5_norm)*100
atp_prod_norm<-basal_resp5_norm-proton_leak_norm

ocr5_norm_tab_grp<-rbind(basal_resp5_norm,max_resp5_norm,spare_cap_ocr5_norm,spare_cap_ocr5_norm_pct
                         ,spare_cap_ocr5_norm_fd,atp_prod_norm,non_mit5_norm)
ocr5_norm_tab_grp<-as.data.frame(t(ocr5_norm_tab_grp))
ocr5_norm_tab_grp$cell_type<-factor(c(rep('H3.3K27M',22),rep('H3.3K27M-KO',24),rep('H3.3K27M',22),rep('H3.3K27M-KO',24)),levels=c('H3.3K27M','H3.3K27M-KO'))
colnames(ocr5_norm_tab_grp)<-c('basal_resp5_norm','max_resp5_norm','spare_cap_ocr5_norm','spare_cap_ocr5_norm_pct'
                               ,'spare_cap_ocr5_norm_fd','atp_prod_norm','non_mit5_norm','cell_type')

test_basal_resp5_norm<-ocr5_norm_tab_grp %>% wilcox_test(basal_resp5_norm~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% mutate(y.position = 2)
test_max_resp5_norm<-ocr5_norm_tab_grp %>% wilcox_test(max_resp5_norm~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% mutate(y.position = 2)
test_spare_cap_ocr5_norm_pct<-ocr5_norm_tab_grp %>% wilcox_test(spare_cap_ocr5_norm_pct~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% mutate(y.position = 95)

#Linegraph
ocr5_norm_melt<-as.data.frame(reshape2::melt(ocr5_norm, id='time'))
ocr5_norm_melt$time<-round(ocr5_norm_melt$time,digits = 2)
ocr5_norm_melt$time<-factor(ocr5_norm_melt$time)
ocr5_norm_melt$cell_type<-factor(c(rep(rep('H3.3K27M',22),21),rep(rep('H3.3K27M-KO',24),21),rep(rep('H3.3K27M',22),21)
                                   ,rep(rep('H3.3K27M-KO',24),21)),levels=c('H3.3K27M','H3.3K27M-KO'))
dim(ocr5_norm_melt)
ocr5_norm_melt<-na.omit(ocr5_norm_melt)
dim(ocr5_norm_melt)
ocr5_norm_sum<-summarySE(ocr5_norm_melt, measurevar="value", groupvars=c('time','cell_type'))
ocr5_norm_sum$time<-as.numeric(levels(ocr5_norm_sum$time))[ocr5_norm_sum$time]

theme_set(theme_classic(base_size = 11,base_line_size = 0.05,base_rect_size = 0.05))

ocr5_norm_tab1<-ggplot(ocr5_norm_tab_grp, aes(x=cell_type, y=basal_resp5_norm)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=3,binwidth = 0.02,aes(fill=cell_type),color='transparent')+
  stat_summary(fun=mean, geom="crossbar",size=0.15,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.3)+
  stat_pvalue_manual(test_basal_resp5_norm, label = "p.signif",size=8,face="bold")+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='(pmol/min)/(ug/mL)')+ggtitle('Basal Respiration')+
  scale_y_continuous(limits = c(0,2.1),breaks=seq(0,2,by=0.4))+
  theme(legend.position="none",axis.title.y=element_text(size=20,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=20,face="bold"),plot.title = element_text(size= 24,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=24, color="black",face="bold"))

ocr5_norm_tab2<-ggplot(ocr5_norm_tab_grp, aes(x=cell_type, y=max_resp5_norm)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=3,binwidth = 0.02,aes(fill=cell_type),color='transparent')+
  stat_summary(fun=mean, geom="crossbar",size=0.15,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.3)+
  stat_pvalue_manual(test_max_resp5_norm, label = "p.signif",size=8,face="bold")+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='(pmol/min)/(ug/mL)')+ggtitle('Maximum Respiration')+
  scale_y_continuous(limits = c(0,2.1),breaks=seq(0,2,by=0.4))+
  theme(legend.position="none",axis.title.y=element_text(size=20,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=22,face="bold"),plot.title = element_text(size= 24,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=24, color="black",face="bold"))

ocr5_norm_tab3<-ggplot(ocr5_norm_tab_grp, aes(x=cell_type, y=spare_cap_ocr5_norm_pct)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=3,binwidth = 1,aes(fill=cell_type),color='transparent')+
  stat_summary(fun=mean, geom="crossbar",size=0.15,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.3)+
  stat_pvalue_manual(test_spare_cap_ocr5_norm_pct, label = "p.signif",size=8,face="bold")+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='Percentage(%)')+ggtitle('Spare Respiratory Capacity')+
  scale_y_continuous(limits = c(0,100),breaks=seq(0,100,by=20))+
  theme(legend.position="none",axis.title.y=element_text(size=20,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=20,face="bold"),plot.title = element_text(size= 24,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=24, color="black",face="bold"))

ocr5_norm_tab4<-ggplot(ocr5_norm_sum, aes(x=time, y=value,color=cell_type,group=cell_type))+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=1,size=0.15)+geom_line(size=0.2)+geom_point(size=0.15)+
  scale_linetype_manual(values=c(rep(c('solid'),2)))+labs(x='Time (min)',y='(pmol/min)/(ug/mL)')+
  scale_y_continuous(limits = c(0,1),breaks=seq(0,1,by=0.2))+
  scale_color_manual(values=c(rep(mycolors1[1:2]),rep(mycolors1[2],2)))+ggtitle('OCR over time')+
  theme(legend.position="bottom",legend.title=element_blank(),legend.text=element_text(hjust=0.5,size=22,face='bold')
        ,axis.title=element_text(size=18,face="bold"),axis.text = element_text(size=20,face='bold')
        ,plot.title = element_text(size= 24,hjust = 0.5,face='bold'),legend.key.size = unit(1,"cm"))+
  guides(linetype = guide_legend(override.aes = list(size = 0.3)))

resp5_norm_grp<-ggarrange(ocr5_norm_tab1,NULL,ocr5_norm_tab2,NULL,NULL,NULL,ocr5_norm_tab3,NULL,ocr5_norm_tab4,ncol=3,nrow=3
                          ,widths=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1),heights=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1))+theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))

#Normalized ECAR
ecar5_norm<-ecar3_norm
colnames(ecar5_norm)<-c('time',paste0('H3.3K27M_rep_',c(23:44)),paste0('H3.3K27M-KO_rep_',c(25:48)))
ecar5_norm<-cbind(ecar4_norm,ecar5_norm)
ecar5_norm<-ecar5_norm[,-47]

ecar5_norm_tab<-ecar5_norm[,-1]
basal_ecar5_norm<-ecar5_norm_tab[3,]
max_ecar5_norm<-ecar5_norm_tab[4,]
min_ecar5_norm<-ecar5_norm_tab[21,]
basal_glyco5_norm<-basal_ecar5_norm-min_ecar5_norm
max_glyco5_norm<-max_ecar5_norm-min_ecar5_norm
spare_cap_glyco5_norm<-max_glyco5_norm-basal_glyco5_norm
spare_cap_glyco5_norm_fd<-max_glyco5_norm/basal_glyco5_norm
spare_cap_glyco5_norm_pct<-((max_glyco5_norm-basal_glyco5_norm)/max_glyco5_norm)*100

glyco5_norm_tab_grp<-rbind(basal_glyco5_norm,max_glyco5_norm,spare_cap_glyco5_norm,spare_cap_glyco5_norm_pct,spare_cap_glyco5_norm_fd)
glyco5_norm_tab_grp<-as.data.frame(t(glyco5_norm_tab_grp))
glyco5_norm_tab_grp$cell_type<-factor(c(rep('H3.3K27M',22),rep('H3.3K27M-KO',23),rep('H3.3K27M',22),rep('H3.3K27M-KO',24)),levels=c('H3.3K27M','H3.3K27M-KO'))
glyco5_norm_tab_grp<-na.omit(glyco5_norm_tab_grp)
colnames(glyco5_norm_tab_grp)<-c('basal_glyco5_norm','max_glyco5_norm','spare_cap_glyco5_norm','spare_cap_glyco5_norm_pct','spare_cap_glyco5_norm_fd'
                                 ,'cell_type')

test_basal_glyco5_norm<-glyco5_norm_tab_grp %>% wilcox_test(basal_glyco5_norm~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% 
  mutate(y.position = 1.75)
test_max_glyco5_norm<-glyco5_norm_tab_grp %>% wilcox_test(max_glyco5_norm~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% 
  mutate(y.position = 1.75)
test_spare_cap_glyco5_norm_pct<-glyco5_norm_tab_grp %>% wilcox_test(spare_cap_glyco5_norm_pct~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% 
  mutate(y.position = 95)

#Linegraph
ecar5_norm_melt<-as.data.frame(reshape2::melt(ecar5_norm, id='time'))
ecar5_norm_melt$time<-round(ecar5_norm_melt$time,digits = 2)
ecar5_norm_melt$time<-factor(ecar5_norm_melt$time)
ecar5_norm_melt$cell_type<-factor(c(rep(rep('H3.3K27M',22),21),rep(rep('H3.3K27M-KO',24),21),rep(rep('H3.3K27M',22),21),rep(rep('H3.3K27M-KO',23),21)),levels=c('H3.3K27M','H3.3K27M-KO'))
dim(ecar5_norm_melt)
ecar5_norm_melt<-na.omit(ecar5_norm_melt)
dim(ecar5_norm_melt)
ecar5_norm_sum<-summarySE(ecar5_norm_melt, measurevar="value", groupvars=c('time','cell_type'))
ecar5_norm_sum$time<-as.numeric(levels(ecar5_norm_sum$time))[ecar5_norm_sum$time]

theme_set(theme_classic(base_size = 11,base_line_size = 0.05,base_rect_size = 0.05))

glyco5_norm_tab1<-ggplot(glyco5_norm_tab_grp, aes(x=cell_type, y=basal_glyco5_norm)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=2.5,binwidth = 0.02,aes(fill=cell_type),color='transparent')+
  stat_summary(fun=mean, geom="crossbar",size=0.15,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.3)+
  stat_pvalue_manual(test_basal_glyco5_norm, label = "p.signif",size=7,face="bold")+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='(mpH/min)/(ug/mL)')+ggtitle('Basal Glycolysis')+
  scale_y_continuous(limits = c(0,1.8),breaks=seq(0,1.5,by=0.5))+
  theme(legend.position="none",axis.title.y=element_text(size=24,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=20,face="bold"),plot.title = element_text(size= 24,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=24, color="black",face="bold"))

glyco5_norm_tab2<-ggplot(glyco5_norm_tab_grp, aes(x=cell_type, y=max_glyco5_norm)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=2.5,binwidth = 0.02,aes(fill=cell_type),color='transparent')+
  stat_summary(fun=mean, geom="crossbar",size=0.15,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.3)+
  stat_pvalue_manual(test_max_glyco5_norm, label = "p.signif",size=7,face="bold")+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='(mpH/min)/(ug/mL)')+ggtitle('Maximum Glycolysis')+
  scale_y_continuous(limits = c(0,1.8),breaks=seq(0,1.5,by=0.5))+
  theme(legend.position="none",axis.title.y=element_text(size=24,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=20,face="bold"),plot.title = element_text(size= 24,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=24, color="black",face="bold"))

glyco5_norm_tab3<-ggplot(glyco5_norm_tab_grp, aes(x=cell_type, y=spare_cap_glyco5_norm_pct)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=3,binwidth = 1,aes(fill=cell_type),color='transparent')+
  stat_summary(fun=mean, geom="crossbar",size=0.15,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.3)+
  stat_pvalue_manual(test_spare_cap_glyco5_norm_pct, label = "p.signif",size=8,face="bold")+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='Percentage(%)')+ggtitle('Spare Glycolytic Capacity')+
  scale_y_continuous(limits = c(0,100),breaks=seq(0,100,by=20))+
  theme(legend.position="none",axis.title.y=element_text(size=24,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=20,face="bold"),plot.title = element_text(size= 24,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=24, color="black",face="bold"))

glyco5_norm_tab4<-ggplot(ecar5_norm_sum, aes(x=time,y=value,group=cell_type,color=cell_type))+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=1,size=0.15)+geom_line(size=0.2)+geom_point(size=0.15)+
  scale_linetype_manual(values=c(rep(c('solid'),2)))+labs(x='Time (min)',y='(mpH/min)/(ug/mL)')+
  scale_color_manual(values=mycolors1[1:2])+ggtitle('ECAR over time')+
  scale_y_continuous(limits = c(0,1),breaks=seq(0,1,by=0.2))+
  theme(legend.position="bottom",legend.title=element_blank(),legend.text=element_text(hjust=0.5,size=24,face='bold')
        ,axis.title=element_text(size=24,face="bold"),axis.text = element_text(size=20,face='bold')
        ,plot.title = element_text(size= 24,hjust = 0.5,face='bold'),legend.key.size = unit(1,"cm"))+
  guides(linetype = guide_legend(override.aes = list(size = 0.3)))

glyco5_norm_grp<-ggarrange(glyco5_norm_tab1,NULL,glyco5_norm_tab2,NULL,NULL,NULL,glyco5_norm_tab3,NULL,glyco5_norm_tab4,ncol=3,nrow=3
                           ,widths=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1),heights=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1))+theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))


####Export Figures----

setwd('/results')#custom destination folder for file generation
ggexport(glyco5_norm_grp,filename='Fig1C.pdf',width=11.69, height=8.27,res=250) #A4 landscape

setwd('/results')#custom destination folder for file generation
ggexport(resp5_norm_grp,filename='Fig1D.pdf',width=11.69, height=8.27, res=250) #A4 landscape
