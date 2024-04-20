#updated by E.S. de Camargo Magalhaes on 20/04/2024 using R v4.3.2 "Eye Holes"

####Seahorse analysis----
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Rmisc))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(rstatix))

###Import excel file generated from wave (overview) and create OCR and ECAR datasets
setwd('data/Datasets/Fig1H_I')#custom folder for importing files
sh_res<-as.data.frame(read_xlsx('seahorse_assay3.xlsx'))
sh_res[sh_res<0]<-NA

ocr<-as.data.frame(sh_res[c(1:24),c(1:104)])
ocr<-ocr[,colSums(is.na(ocr))<nrow(ocr)]
#Remove NA rows between samples
colnames(ocr)<-c('time',paste0('H3.3K27M_DMSO_rep_',c(1:22)),paste0('H3.3K27M_PFK15_5\u00b5M_rep_',c(1:24))
                 ,paste0('H3.3K27M-KO_PFK15_5\u00b5M_rep_',c(1:24)),paste0('H3.3K27M-KO_DMSO_rep_',c(1:22)))
ocr<-ocr[,-51]
#Remove bad samples
row.names(ocr)<-NULL
ocr<-ocr[c(-1:-3),]
row.names(ocr)<-c(1:21)
ocr<-mutate_all(ocr, function(x) as.numeric(as.character(x)))
#Remove unused rows, rename rows and convert values to numeric

#ECAR
ecar<-as.data.frame(sh_res[c(28:51),c(1:104)])
ecar<-ecar[,colSums(is.na(ecar))<nrow(ecar)]
colnames(ecar)<-c('time',paste0('H3.3K27M_DMSO_rep_',c(1:22)),paste0('H3.3K27M_PFK15_5\u00b5M_rep_',c(1:24))
                  ,paste0('H3.3K27M-KO_PFK15_5\u00b5M_rep_',c(1:24)),paste0('H3.3K27M-KO_DMSO_rep_',c(1:22)))
ecar<-ecar[,-51]
row.names(ecar)<-NULL
ecar<-ecar[c(-1:-3),]
row.names(ecar)<-c(1:21)
ecar<-mutate_all(ecar, function(x) as.numeric(as.character(x)))

#Create color object
mycolors1<-c('red3','blue3','green4','darkmagenta','goldenrod4','darkorange','deeppink',
             'gray60','darkcyan','darkturquoise')

##Non-normalized data
#OCR
ocr_tab<-ocr[,-1]
#remove variable time for calculations and figure generation
non_mit<-apply(ocr_tab[c(12:21),],2,function(x){min(x[!is.na(x)])})
non_mit<-ifelse(non_mit==Inf,0,non_mit)
basal_resp<-ocr_tab[3,]-non_mit
max_resp<-apply(ocr_tab[c(8:11),],2,function(x){max(x[!is.na(x)])})-non_mit
proton_leak<-apply(ocr_tab[c(4:7),],2,function(x){max(x[!is.na(x)])})-non_mit
spare_cap_ocr<-max_resp-basal_resp
spare_cap_ocr_fd<-max_resp/basal_resp
spare_cap_ocr_pct<-((max_resp-basal_resp)/max_resp)*100
atp_prod<-basal_resp-proton_leak
#non_mit = minimal value after rotenone/anti-mycin A injection (non-mitochondrial respiration)
#basal_resp = OCR value before first injection (oligomycin) minus non-mit
#max_resp = maximum value after DNP injection minus non-mit
#proton_leak = minimum value after oligomycin injection minus non-mit
#spare_cap = spare respiratory capacity of cells
#spare_cap_fd = fold increase in respiration (baseline is basal respiration)
#spare_cap_pct = spare capacity relative to total respiratory capacity (max resp) in percentage
#atp_prod = Respiration that effectively led to ATP production (basal resp - proton leak)

ocr_tab_grp<-rbind(basal_resp,max_resp,spare_cap_ocr,spare_cap_ocr_pct,spare_cap_ocr_fd,atp_prod,non_mit)
ocr_tab_grp<-as.data.frame(t(ocr_tab_grp))
ocr_tab_grp$cell_type<-factor(c(rep('H3.3K27M',46),rep('H3.3K27M-KO',45)),levels=c('H3.3K27M','H3.3K27M-KO'))
ocr_tab_grp$treatment<-factor(c(rep('DMSO',22),rep('PFK15\n5\u00b5M',24),rep('PFK15\n5\u00b5M',23),rep('DMSO',22)),levels=c('DMSO','PFK15\n5\u00b5M'))
ocr_tab_grp$group<-factor(c(rep('H3.3K27M\nDMSO',22),rep('H3.3K27M\nPFK15\n5\u00b5M',24),rep('H3.3K27M-KO\nPFK15\n5\u00b5M',23)
                            ,rep('H3.3K27M-KO\nDMSO',22)),levels=c('H3.3K27M\nDMSO','H3.3K27M\nPFK15\n5\u00b5M','H3.3K27M-KO\nDMSO','H3.3K27M-KO\nPFK15\n5\u00b5M'))
#Create objects to separate samples by groups and add to dataframe
colnames(ocr_tab_grp)<-c('basal_resp','max_resp','spare_cap_ocr','spare_cap_ocr_pct','spare_cap_ocr_fd','atp_prod','non_mit'
                         ,'cell_type','treatment','group')

#Normality tests
norm_tests_basal_resp<-ocr_tab_grp %>% group_by(group) %>% shapiro_test(basal_resp)
norm_tests_max_resp<-ocr_tab_grp %>% group_by(group) %>% shapiro_test(max_resp)
norm_tests_spare_cap_ocr_pct<-ocr_tab_grp %>% group_by(group) %>% shapiro_test(spare_cap_ocr_pct)
#No test is significant, data is normal

#Compare Means
test_basal_resp<-ocr_tab_grp %>% group_by(cell_type) %>% t_test(basal_resp~treatment)%>% 
  add_significance() %>% add_xy_position(x='treatment') %>% mutate(y.position = max(ocr_tab_grp$max_resp)+sd(ocr_tab_grp$max_resp))
test_max_resp<-ocr_tab_grp %>% group_by(cell_type) %>% t_test(max_resp~treatment)%>% 
  add_significance() %>% add_xy_position(x='treatment') %>% mutate(y.position = max(ocr_tab_grp$max_resp)+sd(ocr_tab_grp$max_resp))
test_spare_cap_ocr_pct<-ocr_tab_grp %>% group_by(cell_type) %>% t_test(spare_cap_ocr_pct~treatment)%>% 
  add_significance() %>% add_xy_position(x='treatment') %>% mutate(y.position=115)

#Linegraph
ocr_melt<-as.data.frame(reshape2::melt(ocr, id='time'))
ocr_melt$time<-round(ocr_melt$time,digits = 2)
ocr_melt$time<-factor(ocr_melt$time)
ocr_melt$cell_type<-factor(c(rep(rep('H3.3K27M',46),21),rep(rep('H3.3K27M-KO',45),21)),levels=c('H3.3K27M','H3.3K27M-KO'))
ocr_melt$treatment<-factor(c(rep(rep('DMSO',22),21),rep(rep('PFK15\n5\u00b5M',24),21),rep(rep('PFK15\n5\u00b5M',23),21)
                             ,rep(rep('DMSO',22),21)),levels=c('DMSO','PFK15\n5\u00b5M'))
ocr_melt$group<-factor(c(rep(rep('H3.3K27M\nDMSO',22),21),rep(rep('H3.3K27M\nPFK15\n5\u00b5M',24),21),rep(rep('H3.3K27M-KO\nPFK15\n5\u00b5M',23),21)
                         ,rep(rep('H3.3K27M-KO\nDMSO',22),21)),levels=c('H3.3K27M\nDMSO','H3.3K27M\nPFK15\n5\u00b5M','H3.3K27M-KO\nDMSO','H3.3K27M-KO\nPFK15\n5\u00b5M'))
ocr_melt<-na.omit(ocr_melt)
dim(ocr_melt)
ocr_sum<-summarySE(ocr_melt, measurevar="value", groupvars=c('group','time','cell_type','treatment'))
ocr_sum$time<-as.numeric(levels(ocr_sum$time))[ocr_sum$time]

#Create Figures
ocr_tab1<-ggplot(ocr_tab_grp, aes(x=treatment,y=basal_resp)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=3,binwidth = 4.5,aes(fill=cell_type))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  stat_pvalue_manual(test_basal_resp, label = "p.signif",size=7)+facet_grid(.~ cell_type)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='pmol/min')+ggtitle('Basal Respiration')+
  scale_y_continuous(limits = c(0,400),breaks=seq(0,400,by=100))+
  theme(legend.position="none",axis.title.y=element_text(size=24,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=24,face="bold"),plot.title = element_text(size= 30,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=26, color="black",face="bold"))

ocr_tab2<-ggplot(ocr_tab_grp, aes(x=treatment, y=max_resp)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=3,binwidth = 4.5,aes(fill=cell_type))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  stat_pvalue_manual(test_max_resp, label = "p.signif",size=7)+facet_grid(.~ cell_type)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='pmol/min')+ggtitle('Maximum Respiration')+
  scale_y_continuous(limits = c(0,400),breaks=seq(0,400,by=100))+
  theme(legend.position="none",axis.title.y=element_text(size=24,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=24,face="bold"),plot.title = element_text(size= 30,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=26, color="black",face="bold"))

ocr_tab3<-ggplot(ocr_tab_grp, aes(x=treatment, y=spare_cap_ocr_pct)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=5.5,binwidth = 0.8,aes(fill=cell_type))+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  stat_pvalue_manual(test_spare_cap_ocr_pct, label = "p.signif",size=7)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='pmol/min')+ggtitle('Spare Respiratory Capacity')+
  scale_y_continuous(limits = c(0, 120),breaks=seq(0,120,by=25))+facet_grid(.~ cell_type)+
  theme(legend.position="none",axis.title.y=element_text(size=24,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=24,face="bold"),plot.title = element_text(size= 30,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=26, color="black",face="bold"))

ocr_tab4<-ggplot(ocr_sum, aes(x=time, y=value, group=group,color=group,linetype=group))+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=1.5,size=1)+geom_line(size=1)+geom_point(size=2)+
  scale_linetype_manual(values=c(rep(c('solid','dashed'),2)))+labs(x='Time (min)',y='pmol/min')+
  scale_color_manual(values=c(rep(mycolors1[1],2),rep(mycolors1[2],2)))+ggtitle('Oxygen Consumption Rate (OCR)')+
  theme(legend.position="bottom",legend.title=element_blank(),legend.text=element_text(hjust=0.5,size=20,face='bold')
        ,axis.title=element_text(size=24,face="bold"),axis.text = element_text(size=22)
        ,plot.title = element_text(size= 26,hjust = 0.5,face='bold'),legend.key.size = unit(3.4,"line",'point'))

resp_grp<-ggarrange(ocr_tab1,NULL,ocr_tab2,NULL,NULL,NULL,ocr_tab3,NULL,ocr_tab4,ncol=3,nrow=3
                    ,widths=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1),heights=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1))+theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))

#ECAR
ecar_tab<-ecar[,-1]
#remove variable time for calculations and figure generation
basal_ecar<-ecar_tab[3,]
max_ecar<-ecar_tab[4,]
min_ecar<-ecar_tab[21,]
basal_glyco<-basal_ecar-min_ecar
max_glyco<-max_ecar-min_ecar
spare_cap_glyco<-max_glyco-basal_glyco
spare_cap_glyco_fd<-max_glyco/basal_glyco
spare_cap_glyco_pct<-((max_glyco-basal_glyco)/max_glyco)*100

glyco_tab_grp<-rbind(basal_glyco,max_glyco,spare_cap_glyco,spare_cap_glyco_pct,spare_cap_glyco_fd)
glyco_tab_grp<-as.data.frame(t(glyco_tab_grp))
glyco_tab_grp$cell_type<-factor(c(rep('H3.3K27M',46),rep('H3.3K27M-KO',45)),levels=c('H3.3K27M','H3.3K27M-KO'))
glyco_tab_grp$treatment<-factor(c(rep('DMSO',22),rep('PFK15\n5\u00b5M',24),rep('PFK15\n5\u00b5M',23),rep('DMSO',22)),levels=c('DMSO','PFK15\n5\u00b5M'))
glyco_tab_grp$group<-factor(c(rep('H3.3K27M\nDMSO',22),rep('H3.3K27M\nPFK15\n5\u00b5M',24),rep('H3.3K27M-KO\nPFK15\n5\u00b5M',23)
                              ,rep('H3.3K27M-KO\nDMSO',22)),levels=c('H3.3K27M\nDMSO','H3.3K27M\nPFK15\n5\u00b5M','H3.3K27M-KO\nDMSO','H3.3K27M-KO\nPFK15\n5\u00b5M'))
#Create objects to separate samples by groups and add to dataframe
colnames(glyco_tab_grp)<-c('basal_glyco','max_glyco','spare_cap_glyco','spare_cap_glyco_pct','spare_cap_glyco_fd'
                           ,'cell_type','treatment','group')

#Normality tests
norm_tests_basal_glyco<-glyco_tab_grp %>% group_by(group) %>% shapiro_test(basal_glyco)
norm_tests_max_glyco<-glyco_tab_grp %>% group_by(group) %>% shapiro_test(max_glyco)
norm_tests_spare_cap_glyco_pct<-glyco_tab_grp %>% group_by(group) %>% shapiro_test(spare_cap_glyco_pct)
#Only H3.3K27M\nPFK15\n5\u00b5M in max_glyco is not normal 

#Compare Means
test_basal_glyco<-glyco_tab_grp %>% group_by(cell_type) %>% t_test(basal_glyco~treatment)%>% 
  add_significance() %>% add_xy_position(x='treatment') %>% mutate(y.position = max(glyco_tab_grp$max_glyco)+sd(glyco_tab_grp$max_glyco)*1.5)
test_max_glyco<-glyco_tab_grp %>% group_by(cell_type) %>% wilcox_test(max_glyco~treatment)%>% 
  add_significance() %>% add_xy_position(x='treatment') %>% mutate(y.position = max(glyco_tab_grp$max_glyco)+sd(glyco_tab_grp$max_glyco)*1.5)
test_spare_cap_glyco_pct<-glyco_tab_grp %>% group_by(cell_type) %>% t_test(spare_cap_glyco_pct~treatment)%>% 
  add_significance() %>% add_xy_position(x='treatment') %>% mutate(y.position=115)

#Linegraph
ecar_melt<-as.data.frame(reshape2::melt(ecar, id='time'))
ecar_melt$time<-round(ecar_melt$time,digits = 2)
ecar_melt$time<-factor(ecar_melt$time)
ecar_melt$cell_type<-factor(c(rep(rep('H3.3K27M',46),21),rep(rep('H3.3K27M-KO',45),21)),levels=c('H3.3K27M','H3.3K27M-KO'))
ecar_melt$treatment<-factor(c(rep(rep('DMSO',22),21),rep(rep('PFK15\n5\u00b5M',24),21),rep(rep('PFK15\n5\u00b5M',23),21)
                              ,rep(rep('DMSO',22),21)),levels=c('DMSO','PFK15\n5\u00b5M'))
ecar_melt$group<-factor(c(rep(rep('H3.3K27M\nDMSO',22),21),rep(rep('H3.3K27M\nPFK15\n5\u00b5M',24),21),rep(rep('H3.3K27M-KO\nPFK15\n5\u00b5M',23),21)
                          ,rep(rep('H3.3K27M-KO\nDMSO',22),21)),levels=c('H3.3K27M\nDMSO','H3.3K27M\nPFK15\n5\u00b5M','H3.3K27M-KO\nDMSO','H3.3K27M-KO\nPFK15\n5\u00b5M'))
ecar_melt<-na.omit(ecar_melt)
dim(ecar_melt)
ecar_sum<-summarySE(ecar_melt, measurevar="value", groupvars=c('group','time','cell_type','treatment'))
ecar_sum$time<-as.numeric(levels(ecar_sum$time))[ecar_sum$time]

glyco_tab1<-ggplot(glyco_tab_grp, aes(x=treatment, y=basal_glyco, fill=cell_type)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=2.5,binwidth = 2.5)+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='mpH/min')+ggtitle('Basal Glycolysis')+
  stat_pvalue_manual(test_basal_glyco, label = "p.signif",size=6)+facet_grid(.~ cell_type)+
  scale_y_continuous(limits = c(0, 165),breaks=seq(0,150,by=25))+
  theme(legend.position="none",axis.title.y=element_text(size=14,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=12,face="bold"),plot.title = element_text(size= 14,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=12, color="black",face="bold"))

glyco_tab2<-ggplot(glyco_tab_grp, aes(x=treatment, y=max_glyco, fill=cell_type)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=2.5,binwidth = 2.5)+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='mpH/min')+ggtitle('Maximum Glycolysis')+
  stat_pvalue_manual(test_max_glyco, label = "p.signif",size=6)+facet_grid(.~ cell_type)+
  scale_y_continuous(limits = c(0, 165),breaks=seq(0,150,by=25))+
  theme(legend.position="none",axis.title.y=element_text(size=14,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=12,face="bold"),plot.title = element_text(size= 14,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=12, color="black",face="bold"))

glyco_tab3<-ggplot(glyco_tab_grp, aes(x=treatment, y=spare_cap_glyco, fill=cell_type)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=2.5,binwidth = 2.5)+
  stat_summary(fun=mean, geom="crossbar",size=0.4,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.8,)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='mpH/min')+ggtitle('Spare Glycolytic Capacity')+
  stat_pvalue_manual(test_spare_cap_glyco_pct, label = "p.signif",size=6)+facet_grid(.~ cell_type)+
  scale_y_continuous(limits = c(0,120),breaks=seq(0,100,by=25))+
  theme(legend.position="none",axis.title.y=element_text(size=14,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=12,face="bold"),plot.title = element_text(size= 14,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=12, color="black",face="bold"))

glyco_tab4<-ggplot(ecar_sum, aes(x=time, y=value, group=group,color=group,linetype=group))+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=1.5,size=1)+geom_line(size=1)+geom_point(size=2)+
  scale_linetype_manual(values=c(rep(c('solid','dashed'),2)))+labs(x='Time (min)',y='mpH/min')+
  scale_y_continuous(limits = c(0,150),breaks=seq(0,150,by=50))+
  scale_color_manual(values=c(rep(mycolors1[1],2),rep(mycolors1[2],2)))+ggtitle('Extracellular Acidification Rate (ECAR)')+
  theme(legend.position="bottom",legend.title=element_blank(),legend.text=element_text(hjust=0.5,size=12,face='bold')
        ,axis.title=element_text(size=16,face="bold"),axis.text = element_text(size=14)
        ,plot.title = element_text(size= 14,hjust = 0.5,face='bold'))

glyco_grp<-ggarrange(glyco_tab1,NULL,glyco_tab2,NULL,NULL,NULL,glyco_tab3,NULL,glyco_tab4,ncol=3,nrow=3
                     ,widths=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1),heights=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1))+theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))

###Data Normalized by protein
bca_res_btwt<-as.data.frame(read_xlsx('BCA_H3.3K27M.xlsx'))
bca_res_btwt<-bca_res_btwt[c(22,23),c(2:49)]
time<-as.matrix(c('time',1))
bca_res_btwt<-cbind(time,bca_res_btwt)
bca_res_btwt<-bca_res_btwt[,colSums(is.na(bca_res_btwt))<nrow(bca_res_btwt)]
#Create column with time variable for further calculations
bca_res_btko<-as.data.frame(read_xlsx('BCA_H3.3K27M-KO.xlsx'))
bca_res_btko<-bca_res_btko[c(22,23),c(2:49)]
bca_res_btko<-bca_res_btko[,colSums(is.na(bca_res_btko))<nrow(bca_res_btko)]
bca_res_btko<-bca_res_btko[,-4]
ptn<-as.data.frame(t(cbind(bca_res_btwt,bca_res_btko)))
rownames(ptn)<-NULL
colnames(ptn)<-c('names','protein')

ocr_norm<-t(ocr)
ocr_norm<-cbind(ptn,ocr_norm)
#View(ocr_norm) #check if names match
ocr_norm<-ocr_norm[,-1]
ocr_norm<-mutate_all(ocr_norm, function(x) as.numeric(as.character(x)))
ocr_norm<- ocr_norm %>% dplyr::mutate_at(vars(-protein),list(~./protein))%>% dplyr::select(c(1:22))
ocr_norm<-ocr_norm[,-1]
ocr_norm<-as.data.frame(t(ocr_norm))

ecar_norm<-t(ecar)
ecar_norm<-cbind(ptn,ecar_norm)
#View(ecar_norm) #check if names match
ecar_norm<-ecar_norm[,-1]
ecar_norm<-mutate_all(ecar_norm, function(x) as.numeric(as.character(x)))
ecar_norm<- ecar_norm %>% dplyr::mutate_at(vars(-protein),funs(./protein))%>% dplyr::select(c(1:22))
ecar_norm<-ecar_norm[,-1]
ecar_norm<-as.data.frame(t(ecar_norm))

#Create subsets, treat data and generate figures
#Normalized OCR
ocr_norm_tab<-ocr_norm[,-1]
non_mit_norm<-apply(ocr_norm_tab[c(12:21),],2,function(x){min(x[!is.na(x)])})
non_mit_norm<-ifelse(non_mit_norm==Inf,0,non_mit_norm)
basal_resp_norm<-ocr_norm_tab[3,]-non_mit_norm
max_resp_norm<-apply(ocr_norm_tab[c(8:11),],2,function(x){max(x[!is.na(x)])})-non_mit_norm
proton_leak_norm<-apply(ocr_norm_tab[c(4:7),],2,function(x){max(x[!is.na(x)])})-non_mit_norm
spare_cap_ocr_norm<-max_resp_norm-basal_resp_norm
spare_cap_ocr_norm_fd<-max_resp_norm/basal_resp_norm
spare_cap_ocr_norm_pct<-((max_resp_norm-basal_resp_norm)/max_resp_norm)*100
atp_prod_norm<-basal_resp_norm-proton_leak_norm

ocr_norm_tab_grp<-rbind(basal_resp_norm,max_resp_norm,spare_cap_ocr_norm,spare_cap_ocr_norm_pct
                        ,spare_cap_ocr_norm_fd,atp_prod_norm,non_mit_norm)
ocr_norm_tab_grp<-as.data.frame(t(ocr_norm_tab_grp))
ocr_norm_tab_grp$cell_type<-factor(c(rep('H3.3K27M',46),rep('H3.3K27M-KO',45)),levels=c('H3.3K27M','H3.3K27M-KO'))
ocr_norm_tab_grp$treatment<-factor(c(rep('DMSO',22),rep('PFK15 5\u00b5M',24),rep('PFK15 5\u00b5M',23),rep('DMSO',22)),levels=c('DMSO','PFK15 5\u00b5M'))
ocr_norm_tab_grp$group<-factor(c(rep('H3.3K27M\nDMSO',22),rep('H3.3K27M\nPFK15 5\u00b5M',24),rep('H3.3K27M-KO\nPFK15 5\u00b5M',23)
                                 ,rep('H3.3K27M-KO\nDMSO',22)),levels=c('H3.3K27M\nDMSO','H3.3K27M\nPFK15 5\u00b5M','H3.3K27M-KO\nDMSO','H3.3K27M-KO\nPFK15 5\u00b5M'))
colnames(ocr_norm_tab_grp)<-c('basal_resp_norm','max_resp_norm','spare_cap_ocr_norm','spare_cap_ocr_norm_pct'
                              ,'spare_cap_ocr_norm_fd','atp_prod_norm','non_mit_norm','cell_type','treatment','group')

#Normality tests
norm_tests_basal_resp_norm<-ocr_norm_tab_grp %>% group_by(group) %>% shapiro_test(basal_resp_norm)
norm_tests_max_resp_norm<-ocr_norm_tab_grp %>% group_by(group) %>% shapiro_test(max_resp_norm)
norm_tests_spare_cap_ocr_norm_pct<-ocr_norm_tab_grp %>% group_by(group) %>% shapiro_test(spare_cap_ocr_norm_pct)
#H3.3K27M-KO\nPFK15 5\u00b5M is significant for basal_resp and H3.3K27M\nDMSO is signif for max_resp, so use Wilcoxon

#Compare Means
test_basal_resp_norm<-ocr_norm_tab_grp %>% group_by(cell_type) %>% wilcox_test(basal_resp_norm~treatment)%>% 
  add_significance() %>% add_xy_position(x='treatment') %>% mutate(y.position = 1.7)
test_max_resp_norm<-ocr_norm_tab_grp %>% group_by(cell_type) %>% wilcox_test(max_resp_norm~treatment)%>% 
  add_significance() %>% add_xy_position(x='treatment') %>% mutate(y.position = 1.7)
test_spare_cap_ocr_norm_pct<-ocr_norm_tab_grp %>% group_by(treatment) %>% t_test(spare_cap_ocr_norm_pct~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% mutate(y.position=90)

#Linegraph
ocr_norm_melt<-as.data.frame(reshape2::melt(ocr_norm, id='time'))
ocr_norm_melt$time<-round(ocr_norm_melt$time,digits = 2)
ocr_norm_melt$time<-factor(ocr_norm_melt$time)
ocr_norm_melt$cell_type<-factor(c(rep(rep('H3.3K27M',46),21),rep(rep('H3.3K27M-KO',45),21)),levels=c('H3.3K27M','H3.3K27M-KO'))
ocr_norm_melt$treatment<-factor(c(rep(rep('DMSO',22),21),rep(rep('PFK15 5\u00b5M',24),21),rep(rep('PFK15 5\u00b5M',23),21)
                                  ,rep(rep('DMSO',22),21)),levels=c('DMSO','PFK15 5\u00b5M'))
ocr_norm_melt$group<-factor(c(rep(rep('H3.3K27M\nDMSO',22),21),rep(rep('H3.3K27M\nPFK15 5\u00b5M',24),21),rep(rep('H3.3K27M-KO\nPFK15 5\u00b5M',23),21)
                              ,rep(rep('H3.3K27M-KO\nDMSO',22),21)),levels=c('H3.3K27M\nDMSO','H3.3K27M\nPFK15 5\u00b5M','H3.3K27M-KO\nDMSO','H3.3K27M-KO\nPFK15 5\u00b5M'))
ocr_norm_melt<-na.omit(ocr_norm_melt)
dim(ocr_norm_melt)
ocr_norm_sum<-summarySE(ocr_norm_melt, measurevar="value", groupvars=c('group','time','cell_type','treatment'))
ocr_norm_sum$time<-as.numeric(levels(ocr_norm_sum$time))[ocr_norm_sum$time]

theme_set(theme_classic(base_size = 11,base_line_size = 0.05,base_rect_size = 0.05))

ocr_norm_tab1<-ggplot(ocr_norm_tab_grp, aes(x=treatment, y=basal_resp_norm)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=3,binwidth = 0.02,aes(fill=cell_type),color='transparent')+
  stat_summary(fun=mean, geom="crossbar",size=0.15,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.3)+
  stat_pvalue_manual(test_basal_resp_norm, label = "p.signif",size=7,face="bold")+facet_grid(.~ cell_type)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='(pmol/min)/(ug/mL)')+ggtitle('Basal Respiration')+
  scale_y_continuous(limits = c(0,1.8),breaks=seq(0,1.6,by=0.4))+
  theme(legend.position="none",axis.title.y=element_text(size=18,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=16,face="bold"),plot.title = element_text(size= 20,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=20, color="black",face="bold"))

ocr_norm_tab2<-ggplot(ocr_norm_tab_grp, aes(x=treatment, y=max_resp_norm)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=3,binwidth = 0.02,aes(fill=cell_type),color='transparent')+
  stat_summary(fun=mean, geom="crossbar",size=0.15,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.3)+
  stat_pvalue_manual(test_max_resp_norm, label = "p.signif",size=7,face="bold")+
  facet_grid(.~ cell_type)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='(pmol/min)/(ug/mL)')+ggtitle('Maximum Respiration')+
  scale_y_continuous(limits = c(0,1.8),breaks=seq(0,1.6,by=0.4))+
  theme(legend.position="none",axis.title.y=element_text(size=18,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=16,face="bold"),plot.title = element_text(size= 20,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=20, color="black",face="bold"))

ocr_norm_tab3<-ggplot(ocr_norm_tab_grp, aes(x=cell_type, y=spare_cap_ocr_norm_pct)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=3,binwidth = 1,aes(fill=cell_type),color='transparent')+
  stat_summary(fun=mean, geom="crossbar",size=0.15,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.3)+
  stat_pvalue_manual(test_spare_cap_ocr_norm_pct, label = "p.signif",size=7,face="bold")+facet_grid(.~ treatment)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='Percentage (%)')+ggtitle('Spare Respiratory Capacity')+
  scale_y_continuous(limits = c(0,100),expand = c(0,0),breaks=seq(0,100,by=20))+
  theme(legend.position="none",axis.title.y=element_text(size=18,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=16,face="bold"),plot.title = element_text(size= 20,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=20, color="black",face="bold"))

ocr_norm_tab4<-ggplot(ocr_norm_sum, aes(x=time, y=value, group=group,color=group,linetype=group))+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=1,size=0.15)+geom_line(size=0.2)+geom_point(size=0.15)+
  scale_linetype_manual(values=c(rep(c('solid','dotted'),2)))+labs(x='Time (min)',y='(pmol/min)/(ug/mL)')+
  scale_y_continuous(limits = c(0,1.2),breaks=seq(0,1.2,by=0.4))+
  scale_color_manual(values=c(rep(mycolors1[1],2),rep(mycolors1[2],2)))+ggtitle('OCR over time')+
  theme(legend.position="bottom",legend.title=element_blank(),legend.text=element_text(hjust=0.5,size=14,face='bold')
        ,axis.title=element_text(size=16,face="bold"),axis.text = element_text(size=16,face='bold')
        ,plot.title = element_text(size= 20,hjust = 0.5,face='bold'),legend.key.size = unit(2,"line",'point'))+
  guides(linetype = guide_legend(override.aes = list(size = 0.3)))

resp_norm_grp<-ggarrange(ocr_norm_tab1,NULL,ocr_norm_tab2,NULL,NULL,NULL,ocr_norm_tab3,NULL,ocr_norm_tab4,ncol=3,nrow=3
                         ,widths=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1),heights=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1))+theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))

#Normalized ECAR
ecar_norm_tab<-ecar_norm[,-1]
basal_ecar_norm<-ecar_norm_tab[3,]
max_ecar_norm<-ecar_norm_tab[4,]
min_ecar_norm<-ecar_norm_tab[21,]
basal_glyco_norm<-basal_ecar_norm-min_ecar_norm
max_glyco_norm<-max_ecar_norm-min_ecar_norm
spare_cap_glyco_norm<-max_glyco_norm-basal_glyco_norm
spare_cap_glyco_norm_fd<-max_glyco_norm/basal_glyco_norm
spare_cap_glyco_norm_pct<-((max_glyco_norm-basal_glyco_norm)/max_glyco_norm)*100

glyco_norm_tab_grp<-rbind(basal_glyco_norm,max_glyco_norm,spare_cap_glyco_norm,spare_cap_glyco_norm_pct,spare_cap_glyco_norm_fd)
glyco_norm_tab_grp<-as.data.frame(t(glyco_norm_tab_grp))
glyco_norm_tab_grp$cell_type<-factor(c(rep('H3.3K27M',46),rep('H3.3K27M-KO',45)),levels=c('H3.3K27M','H3.3K27M-KO'))
glyco_norm_tab_grp$treatment<-factor(c(rep('DMSO',22),rep('PFK15 5\u00b5M',24),rep('PFK15 5\u00b5M',23),rep('DMSO',22)),levels=c('DMSO','PFK15 5\u00b5M'))
glyco_norm_tab_grp$group<-factor(c(rep('H3.3K27M\nDMSO',22),rep('H3.3K27M\nPFK15 5\u00b5M',24),rep('H3.3K27M-KO\nPFK15 5\u00b5M',23)
                                   ,rep('H3.3K27M-KO\nDMSO',22)),levels=c('H3.3K27M\nDMSO','H3.3K27M\nPFK15 5\u00b5M','H3.3K27M-KO\nDMSO','H3.3K27M-KO\nPFK15 5\u00b5M'))
colnames(glyco_norm_tab_grp)<-c('basal_glyco_norm','max_glyco_norm','spare_cap_glyco_norm','spare_cap_glyco_norm_pct','spare_cap_glyco_norm_fd'
                                ,'cell_type','treatment','group')

#Normality tests
norm_tests_basal_glyco_norm<-glyco_norm_tab_grp %>% group_by(group) %>% shapiro_test(basal_glyco_norm)
norm_tests_max_glyco_norm<-glyco_norm_tab_grp %>% group_by(group) %>% shapiro_test(max_glyco_norm)
norm_tests_spare_cap_glyco_norm_pct<-glyco_norm_tab_grp %>% group_by(group) %>% shapiro_test(spare_cap_glyco_norm_pct)
#Only H3.3K27M\nPFK15 5\u00b5M in max_glyco is not normal 

#Compare Means
test_basal_glyco_norm<-glyco_norm_tab_grp %>% group_by(cell_type) %>% wilcox_test(basal_glyco_norm~treatment)%>% 
  add_significance() %>% add_xy_position(x='treatment') %>% mutate(y.position = 0.9)
test_max_glyco_norm<-glyco_norm_tab_grp %>% group_by(cell_type) %>% wilcox_test(max_glyco_norm~treatment)%>% 
  add_significance() %>% add_xy_position(x='treatment') %>% mutate(y.position = 0.9)
test_spare_cap_glyco_norm_pct<-glyco_norm_tab_grp %>% group_by(treatment) %>% t_test(spare_cap_glyco_norm_pct~cell_type)%>% 
  add_significance() %>% add_xy_position(x='cell_type') %>% mutate(y.position=90)

#Linegraph
ecar_norm_melt<-as.data.frame(reshape2::melt(ecar_norm, id='time'))
ecar_norm_melt$time<-round(ecar_norm_melt$time,digits = 2)
ecar_norm_melt$time<-factor(ecar_norm_melt$time)
ecar_norm_melt$cell_type<-factor(c(rep(rep('H3.3K27M',46),21),rep(rep('H3.3K27M-KO',45),21)),levels=c('H3.3K27M','H3.3K27M-KO'))
ecar_norm_melt$treatment<-factor(c(rep(rep('DMSO',22),21),rep(rep('PFK15 5\u00b5M',24),21),rep(rep('PFK15 5\u00b5M',23),21)
                                   ,rep(rep('DMSO',22),21)),levels=c('DMSO','PFK15 5\u00b5M'))
ecar_norm_melt$group<-factor(c(rep(rep('H3.3K27M\nDMSO',22),21),rep(rep('H3.3K27M\nPFK15 5\u00b5M',24),21),rep(rep('H3.3K27M-KO\nPFK15 5\u00b5M',23),21)
                               ,rep(rep('H3.3K27M-KO\nDMSO',22),21)),levels=c('H3.3K27M\nDMSO','H3.3K27M\nPFK15 5\u00b5M','H3.3K27M-KO\nDMSO','H3.3K27M-KO\nPFK15 5\u00b5M'))
ecar_norm_melt<-na.omit(ecar_norm_melt)
dim(ecar_norm_melt)
ecar_norm_sum<-summarySE(ecar_norm_melt, measurevar="value", groupvars=c('group','time','cell_type','treatment'))
ecar_norm_sum$time<-as.numeric(levels(ecar_norm_sum$time))[ecar_norm_sum$time]

theme_set(theme_classic(base_size = 11,base_line_size = 0.05,base_rect_size = 0.05))

glyco_norm_tab1<-ggplot(glyco_norm_tab_grp, aes(x=treatment, y=basal_glyco_norm)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=3.5,binwidth = 0.01,aes(fill=cell_type),color='transparent')+
  stat_summary(fun=mean, geom="crossbar",size=0.15,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.3,)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='(mpH/min)/(ug/mL)')+ggtitle('Basal Glycolysis')+
  stat_pvalue_manual(test_basal_glyco_norm, label = "p.signif",size=7,face="bold")+facet_grid(.~ cell_type)+
  scale_y_continuous(limits = c(0,1),breaks=seq(0,1,by=0.2))+
  theme(legend.position="none",axis.title.y=element_text(size=20,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=16,face="bold"),plot.title = element_text(size= 20,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=20, color="black",face="bold"))

glyco_norm_tab2<-ggplot(glyco_norm_tab_grp, aes(x=treatment, y=max_glyco_norm)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=3.5,binwidth = 0.01,aes(fill=cell_type),color='transparent')+
  stat_summary(fun=mean, geom="crossbar",size=0.15,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.3,)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='(mpH/min)/(ug/mL)')+ggtitle('Maximum Glycolysis')+
  stat_pvalue_manual(test_max_glyco_norm, label = "p.signif",size=7,face="bold")+facet_grid(.~ cell_type)+
  scale_y_continuous(limits = c(0,1),breaks=seq(0,1,by=0.2))+
  theme(legend.position="none",axis.title.y=element_text(size=20,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=16,face="bold"),plot.title = element_text(size= 20,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=20, color="black",face="bold"))

glyco_norm_tab3<-ggplot(glyco_norm_tab_grp, aes(x=cell_type, y=spare_cap_glyco_norm_pct)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=3,binwidth =1,aes(fill=cell_type),color='transparent')+
  stat_summary(fun=mean, geom="crossbar",size=0.15,width=0.3, color="black")+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.15,size=0.3)+facet_grid(.~ treatment)+
  scale_fill_manual(values=mycolors1[1:2])+labs(y='Percentage (%)')+ggtitle('Spare Glycolytic Capacity')+
  stat_pvalue_manual(test_spare_cap_glyco_norm_pct, label = "p.signif",size=7,face="bold")+
  scale_y_continuous(limits = c(0,100),expand = c(0,0),breaks=seq(0,100,by=20))+
  theme(legend.position="none",axis.title.y=element_text(size=20,face="bold"),axis.title.x= element_blank()
        ,axis.text=element_text(size=16,face="bold"),plot.title = element_text(size=20,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=20, color="black",face="bold"))

glyco_norm_tab4<-ggplot(ecar_norm_sum, aes(x=time, y=value, group=group,color=group,linetype=group))+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=1,size=0.15)+geom_line(size=0.2)+geom_point(size=0.15)+
  scale_linetype_manual(values=c(rep(c('solid','dotted'),2)))+labs(x='Time (min)',y='(mpH/min)/(ug/mL)')+
  scale_y_continuous(limits = c(0,0.8),breaks=seq(0,0.8,by=0.2))+
  scale_color_manual(values=c(rep(mycolors1[1],2),rep(mycolors1[2],2)))+ggtitle('ECAR over time')+
  theme(legend.position="bottom",legend.title=element_blank(),legend.text=element_text(hjust=0.5,size=14)
        ,axis.title.y=element_text(size=16),axis.title.x=element_text(size=20),axis.text = element_text(size=16, face='bold')
        ,plot.title = element_text(size= 20,hjust = 0.5,face='bold'),legend.key.size = unit(2,"line",'point'))+
  guides(linetype = guide_legend(override.aes = list(size = 0.3)))


glyco_norm_grp2<-ggarrange(glyco_norm_tab1,NULL,glyco_norm_tab2,NULL,NULL,NULL,glyco_norm_tab3,NULL,glyco_norm_tab4,ncol=3,nrow=3
                          ,widths=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1),heights=c(1,0.05,1,0.05,0.05,0.05,1,0.05,1))+theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))

####Export Figures----
setwd('/results')#custom destination folder for file generation
ggexport(glyco_norm_grp2,filename='Fig1H.pdf',width=11.69, height=8.27, pointsize=8, res=250) #A4 landscape

setwd('/results')#custom destination folder for file generation
ggexport(resp_norm_grp,filename='Fig1I.pdf',width=11.69, height=8.27, pointsize=8, res=250) #A4 landscape

