#updated by E.S. de Camargo Magalhaes on 20/04/2024 using R v4.3.2 "Eye Holes"

####Cell cycle analysis----
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Rmisc))
suppressPackageStartupMessages(library(rstatix))

setwd('data/Datasets/Fig1G')#custom folder for importing files
cc_res<-as.data.frame(read_xlsx('analysis_cel_cycle.xlsx',sheet='cell_cycle'))

#Remove id values and separate df
cc_tab<-as.data.frame(cc_res[,-4])
cc_tab$phase<-factor(ifelse(cc_tab$phase=='G1','G0/G1',ifelse(cc_tab$phase=='G2','G2/M','S-phase')),levels=c('G0/G1','S-phase','G2/M'))
cc_tab$trt<-factor(ifelse(cc_tab$trt=='10uM','PFK15 10\u00b5M',ifelse(cc_tab$trt=='20uM','PFK15 20\u00b5M','DMSO')))
cc_tab$type<-factor(cc_tab$type)
cc_tab$type<-relevel(cc_tab$type,'H3.3K27M')

#Create color object
mycolors1<-c('red3','blue3','green4','darkmagenta','goldenrod4','darkorange','deeppink',
             'gray60','darkcyan','darkturquoise')

#Compare means1
test_cc1<-cc_tab%>% group_by(phase,type) %>% tukey_hsd(value~trt) %>%
  add_significance() %>% add_xy_position(x='phase')%>% mutate(y.position=y.position+8)

theme_set(theme_classic(base_size = 11,base_line_size = 0.05,base_rect_size = 0.05))

cc_plot<-ggplot(cc_tab, aes(x=phase, y=value, group=trt,color=trt))+ggtitle('Cell Cycle Analysis')+facet_grid(.~ type)+
  geom_bar(aes(phase, value, fill = trt),color='black',position = position_dodge(width=0.7), stat = "summary", fun="mean",width=0.5,size=0.25)+
  geom_text(size=7,fontface='bold',color='black',stat = "summary", fun = "mean",aes(y = stage(value, after_stat = y+8),label = round(after_stat(y),0)),position = position_dodge(width=0.7))+
  stat_summary(fun.data=mean_se,geom="errorbar",position = position_dodge(width=0.7), color="black", size=0.15,width=0.2)+
  stat_pvalue_manual(test_cc1, label = "p.adj.signif",size=7,face="bold")+
  geom_point(position = position_dodge(width=0.7), color='black',stat = "identity",size=0.25)+
  scale_fill_manual(values=mycolors1[3:5])+labs(y='Percentage(%)')+scale_y_continuous(limits = c(0,110),breaks=seq(0,100,by=20))+
  theme(legend.position="bottom",legend.title=element_blank(),legend.text=element_text(hjust=0.5,size=18,face='bold')
        ,axis.title.y=element_text(size=20,face="bold",margin=margin(t=0,r=15,b=0,l=0)),axis.title.x= element_blank()
        ,axis.text=element_text(size=18,face="bold"),plot.title = element_text(size=20,hjust = 0.5,face='bold')
        ,strip.text.x = element_text(size=20, color="black",face="bold"),legend.key.size = unit(2,"line"))

####Export Figures----
setwd('/results')#custom destination folder for file generation
ggexport(cc_plot,filename='Fig1G.pdf',width=21.21, height=15, pointsize=8, res=250) #A4 landscape

