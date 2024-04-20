#updated by E.S. de Camargo Magalhaes on 20/04/2024 using R v4.3.2 "Eye Holes"

####Trypan Blue 24h analysis----
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Rmisc))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(rstatix))
suppressPackageStartupMessages(library(drc))
suppressPackageStartupMessages(library(sandwich))
suppressPackageStartupMessages(library(lmtest))
suppressPackageStartupMessages(library(multcomp))
suppressPackageStartupMessages(library(ggbreak))
suppressPackageStartupMessages(library(BSDA))

setwd('data/Datasets/Fig1E')#custom folder for importing files
tpb24_res<-as.data.frame(read_xlsx('Trypan_blue_PFK15_24h_48h.xlsx',sheet='24H'))

#Create color object
mycolors1<-c('red3','blue3','green4','darkmagenta','goldenrod4','darkorange','deeppink',
             'gray60','darkcyan','darkturquoise')

#Extract values
tpb24_tab<-as.data.frame(tpb24_res[,c(1:3)])
colnames(tpb24_tab)<-c('pfk15','sample','value')
tpb24_tab$pfk15<-c(rep(0,6),rep(3,6),rep(10,6),rep(30,6),rep(100,6))
tpb24_tab$type<-factor(rep(c(rep('H3.3K27M',3),rep('H3.3K27M-KO',3)),5),levels=c('H3.3K27M','H3.3K27M-KO'))
tpb24_tab$log10_pfk15<-log10(tpb24_tab$pfk15)
tpb24_tab$log10_pfk15<-ifelse(tpb24_tab$log10_pfk15<0,0,tpb24_tab$log10_pfk15)
tpb24_tab<-as.data.frame(tpb24_tab[,-2])
tpb24_sum<-summarySE(tpb24_tab, measurevar="value",groupvars=c('type','pfk15','log10_pfk15'))
tpb24_btwt<-tpb24_tab[tpb24_tab$type=='H3.3K27M',]
tpb24_btko<-tpb24_tab[tpb24_tab$type=='H3.3K27M-KO',]

#Generate dose-response mod_tpb24els
#Global model
mod_tpb24<-drm(value~pfk15,type, data=tpb24_tab, fct = LL.4(names = c("Slope", "lower", "upper", "mod_tpb24_ed50")),pmodels=list(~type-1, ~1, ~1, ~type-1))
summary(mod_tpb24)
summary(glht(mod_tpb24))
mod_tpb24_ed50<-as.data.frame(ED(mod_tpb24,c(50),interval="delta"))
test24<-symnum(tsum.test(mean.x=mod_tpb24_ed50[1,1],s.x=mod_tpb24_ed50[1,2],n.x=3,mean.y=mod_tpb24_ed50[2,1],s.y=mod_tpb24_ed50[2,2],n.y=3)$p.value
               ,cutpoints = c(0, 0.001, 0.01, 0.05, 1),symbols = c("p<0.001", "p<0.01", "p<0.05", "ns "),abbr.colnames = F, na = 'N/A')[[1]]

#H3.3K27M mod_tpb24
mod_tpb24_btwt<-drm(value~pfk15,data=tpb24_btwt,fct=LL.4(names = c("Slope", "lower", "upper", "mod_tpb24_ed50")))
summary(mod_tpb24_btwt)
summary(glht(mod_tpb24_btwt))
mod_tpb24_ed50_btwt<-as.data.frame(ED(mod_tpb24_btwt,c(50),interval="delta"))

#H3.3K27M-KO mod_tpb24
mod_tpb24_btko<-drm(value~pfk15,data=tpb24_btko,fct=LL.4(names = c("Slope", "lower", "upper", "mod_tpb24_ed50")))
summary(mod_tpb24_btko)
summary(glht(mod_tpb24_btko))
mod_tpb24_ed50_btko<-as.data.frame(ED(mod_tpb24_btko,c(50),interval="delta"))

#Adjust data for plotting with ggplot2
mod_tpb242_btwt<-expand.grid(conc=exp(seq(log(0.5),log(100),length=100)))
pmod_tpb24_btwt<-predict(mod_tpb24_btwt,mod_tpb242_btwt,interval="confidence")
mod_tpb242_btwt$p<-pmod_tpb24_btwt[,1]
mod_tpb242_btwt$pmin<-pmod_tpb24_btwt[,2]
mod_tpb242_btwt$pmax<-pmod_tpb24_btwt[,3]
mod_tpb242_btwt$type<-factor(rep('H3.3K27M',100),levels=c('H3.3K27M','H3.3K27M-KO'))

mod_tpb242_btko<-expand.grid(conc=exp(seq(log(0.5),log(100),length=100)))
pmod_tpb24_btko<-predict(mod_tpb24_btko,mod_tpb242_btko,interval="confidence")
mod_tpb242_btko$p<-pmod_tpb24_btko[,1]
mod_tpb242_btko$pmin<-pmod_tpb24_btko[,2]
mod_tpb242_btko$pmax<-pmod_tpb24_btko[,3]
mod_tpb242_btko$type<-factor(rep('H3.3K27M-KO',100),levels=c('H3.3K27M','H3.3K27M-KO'))

#Recode variables for ggplot
mod_tpb242_btko$log10_conc<-log10(mod_tpb242_btko$conc)
mod_tpb242_btko$log10_conc<-ifelse(mod_tpb242_btko$log10_conc<0,0,mod_tpb242_btko$log10_conc)

mod_tpb242_btwt$log10_conc<-log10(mod_tpb242_btwt$conc)
mod_tpb242_btwt$log10_conc<-ifelse(mod_tpb242_btwt$log10_conc<0,0,mod_tpb242_btwt$log10_conc)

ed50_btko_tbp24<-round(log10(mod_tpb24_ed50$Estimate[2]),digits=2)
ed50_btwt_tbp24<-round(log10(mod_tpb24_ed50$Estimate[1]),digits=2)

#Plot models with ggplot2
theme_set(theme_classic(base_size = 11,base_line_size = 0.1,base_rect_size = 0.1))

tpb24_plot<-ggplot(tpb24_sum, aes(x=log10_pfk15,y=value)) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se,color=type), width=0.05,size=0.2)+
  scale_color_manual(values=rep(mycolors1[1:2],2))+annotate('text',x =1.6,y=75,label=test24, size=12)+
  geom_segment(data=mod_tpb24_ed50_btwt,aes(x=-Inf,y=50,xend=1.36,yend=50),linetype="dashed",color='black',size=0.15,linejoin="round")+
  geom_segment(data=mod_tpb24_ed50_btwt,aes(x=1.36,y=-Inf,xend=1.36,yend=50),linetype="dashed",color=mycolors1[1],size=0.15,linejoin="round")+
  geom_segment(data=mod_tpb24_ed50_btko,aes(x=1.30,y=-Inf,xend=1.30,yend=50),linetype="dashed",color=mycolors1[2],size=0.15,linejoin="round")+
  geom_line(data=mod_tpb242_btwt, aes(x=log10_conc, y=p,color=type),size=0.15)+ggtitle('Trypan Blue PFK15 24h')+
  geom_line(data=mod_tpb242_btko, aes(x=log10_conc, y=p,color=type),size=0.15)+geom_point(size=4,aes(color=type))+
  scale_x_continuous(limits=c(0,2),breaks=c(0,0.5,1,1.5,2.0))+
  annotate('text',x=ed50_btwt_tbp24+0.2,y=0,label=ed50_btwt_tbp24, size=12,color=mycolors1[1])+
  annotate('text',x=ed50_btko_tbp24-0.15,y=0,label=ed50_btko_tbp24, size=12,color=mycolors1[2])+
  xlab("Log10 [PFK15], \u00b5M") + ylab("Viability (%)")+scale_y_continuous(limits=c(0,110),breaks=seq(0,100,by=25))+
  theme(legend.position="bottom",legend.title=element_blank(),legend.text=element_text(hjust=0.5,size=32)
        ,axis.title.y=element_text(size=36),axis.title.x=element_text(size=38),axis.text = element_text(size=32)
        ,plot.title = element_text(size= 42,hjust = 0.5),legend.key.size = unit(5,"line",'point'))

print(tpb24_plot)


####Trypan Blue 48h analysis----
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Rmisc))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(rstatix))
suppressPackageStartupMessages(library(drc))
suppressPackageStartupMessages(library(sandwich))
suppressPackageStartupMessages(library(lmtest))
suppressPackageStartupMessages(library(multcomp))
suppressPackageStartupMessages(library(ggbreak))
suppressPackageStartupMessages(library(BSDA))

setwd('data/Datasets/Fig1E')#custom folder for importing files
tpb48_res<-as.data.frame(read_xlsx('Trypan_blue_PFK15_24h_48h.xlsx',sheet='48H'))

#Create color object
mycolors1<-c('red3','blue3','green4','darkmagenta','goldenrod4','darkorange','deeppink',
             'gray60','darkcyan','darkturquoise')

#Extract values
tpb48_tab<-as.data.frame(tpb48_res[,c(1:3)])
colnames(tpb48_tab)<-c('pfk15','sample','value')
tpb48_tab$pfk15<-c(rep(0,6),rep(3,6),rep(10,6),rep(30,6),rep(100,6))
tpb48_tab$type<-factor(rep(c(rep('H3.3K27M',3),rep('H3.3K27M-KO',3)),5),levels=c('H3.3K27M','H3.3K27M-KO'))
tpb48_tab$log10_pfk15<-log10(tpb48_tab$pfk15)
tpb48_tab$log10_pfk15<-ifelse(tpb48_tab$log10_pfk15<0,0,tpb48_tab$log10_pfk15)
tpb48_tab<-as.data.frame(tpb48_tab[,-2])
tpb48_sum<-summarySE(tpb48_tab, measurevar="value",groupvars=c('type','pfk15','log10_pfk15'))
tpb48_btwt<-tpb48_tab[tpb48_tab$type=='H3.3K27M',]
tpb48_btko<-tpb48_tab[tpb48_tab$type=='H3.3K27M-KO',]

#Generate dose-response mod_tpb48els
#Global model
mod_tpb48<-drm(value~pfk15,type, data=tpb48_tab, fct = LL.4(names = c("Slope", "lower", "upper", "mod_tpb48_ed50")),pmodels=list(~type-1, ~1, ~1, ~type-1))
summary(mod_tpb48)
summary(glht(mod_tpb48))
mod_tpb48_ed50<-as.data.frame(ED(mod_tpb48,c(50),interval="delta"))
test48<-symnum(tsum.test(mean.x=mod_tpb48_ed50[1,1],s.x=mod_tpb48_ed50[1,2],n.x=3,mean.y=mod_tpb48_ed50[2,1],s.y=mod_tpb48_ed50[2,2],n.y=3)$p.value
               ,cutpoints = c(0, 0.001, 0.01, 0.05, 1),symbols = c("p<0.001", "p<0.01", "p<0.05", "ns "),abbr.colnames = F, na = 'N/A')[[1]]

#BT-WT mod_tpb48el
mod_tpb48_btwt<-drm(value~pfk15,data=tpb48_btwt,fct=LL.4(names = c("Slope", "lower", "upper", "mod_tpb48_ed50")))
summary(mod_tpb48_btwt)
summary(glht(mod_tpb48_btwt))
mod_tpb48_ed50_btwt<-as.data.frame(ED(mod_tpb48_btwt,c(50),interval="delta"))

#H3.3K27M-KO mod_tpb48el
mod_tpb48_btko<-drm(value~pfk15,data=tpb48_btko,fct=LL.4(names = c("Slope", "lower", "upper", "mod_tpb48_ed50")))
summary(mod_tpb48_btko)
summary(glht(mod_tpb48_btko))
mod_tpb48_ed50_btko<-as.data.frame(ED(mod_tpb48_btko,c(50),interval="delta"))

#Adjust data for plotting with ggplot2
mod_tpb482_btwt<-expand.grid(conc=exp(seq(log(0.5),log(100),length=100)))
pmod_tpb48_btwt<-predict(mod_tpb48_btwt,mod_tpb482_btwt,interval="confidence")
mod_tpb482_btwt$p<-pmod_tpb48_btwt[,1]
mod_tpb482_btwt$pmin<-pmod_tpb48_btwt[,2]
mod_tpb482_btwt$pmax<-pmod_tpb48_btwt[,3]
mod_tpb482_btwt$type<-factor(rep('H3.3K27M',100),levels=c('H3.3K27M','H3.3K27M-KO'))

mod_tpb482_btko<-expand.grid(conc=exp(seq(log(0.5),log(100),length=100)))
pmod_tpb48_btko<-predict(mod_tpb48_btko,mod_tpb482_btko,interval="confidence")
mod_tpb482_btko$p<-pmod_tpb48_btko[,1]
mod_tpb482_btko$pmin<-pmod_tpb48_btko[,2]
mod_tpb482_btko$pmax<-pmod_tpb48_btko[,3]
mod_tpb482_btko$type<-factor(rep('H3.3K27M-KO',100),levels=c('H3.3K27M','H3.3K27M-KO'))

#Recode variables for ggplot
mod_tpb482_btko$log10_conc<-log10(mod_tpb482_btko$conc)
mod_tpb482_btko$log10_conc<-ifelse(mod_tpb482_btko$log10_conc<0,0,mod_tpb482_btko$log10_conc)

mod_tpb482_btwt$log10_conc<-log10(mod_tpb482_btwt$conc)
mod_tpb482_btwt$log10_conc<-ifelse(mod_tpb482_btwt$log10_conc<0,0,mod_tpb482_btwt$log10_conc)

ed50_btko_tbp48<-round(log10(mod_tpb48_ed50$Estimate[2]),digits=2)
ed50_btwt_tbp48<-round(log10(mod_tpb48_ed50$Estimate[1]),digits=2)

#Plot models with ggplot2
theme_set(theme_classic(base_size = 11,base_line_size = 0.1,base_rect_size = 0.1))

tpb48_plot<-ggplot(tpb48_sum, aes(x=log10_pfk15,y=value)) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se,color=type), width=0.05,size=0.2)+
  scale_color_manual(values=rep(mycolors1[1:2],2))+annotate('text',x =1.6,y=75,label=test48, size=12)+
  geom_segment(data=mod_tpb48_ed50_btwt,aes(x=-Inf,y=50,xend=log10(Estimate),yend=50),linetype="dashed",color='black',size=0.15,linejoin="round")+
  geom_segment(data=mod_tpb48_ed50_btwt,aes(x=log10(Estimate),y=-Inf,xend=log10(Estimate),yend=50),linetype="dashed",color=mycolors1[1],size=0.15,linejoin="round")+
  geom_segment(data=mod_tpb48_ed50_btko,aes(x=log10(Estimate),y=-Inf,xend=log10(Estimate),yend=50),linetype="dashed",color=mycolors1[2],size=0.15,linejoin="round")+
  geom_line(data=mod_tpb482_btwt, aes(x=log10_conc, y=p,color=type),size=0.15)+geom_point(size=4,aes(color=type))+
  geom_line(data=mod_tpb482_btko, aes(x=log10_conc, y=p,color=type),size=0.15)+ggtitle('Trypan Blue PFK15 48h')+
  scale_x_continuous(limits=c(0,2),breaks=c(0,0.5,1,1.5,2.0))+
  annotate('text',x=ed50_btwt_tbp48+0.175,y=0,label=ed50_btwt_tbp48, size=12,color=mycolors1[1])+
  annotate('text',x=ed50_btko_tbp48-0.175,y=0,label=ed50_btko_tbp48, size=12,color=mycolors1[2])+
  xlab("Log10 [PFK15], \u00b5M") + ylab("Viability (%)")+
  theme(legend.position="bottom",legend.title=element_blank(),legend.text=element_text(hjust=0.5,size=32)
        ,axis.title.y=element_text(size=36),axis.title.x=element_text(size=38),axis.text = element_text(size=32)
        ,plot.title = element_text(size= 42,hjust = 0.5),legend.key.size = unit(5,"line",'point'))

print(tpb48_plot)

####Export Figure----

setwd('/results')#custom destination folder for file generation
tpb_grp<-ggarrange(tpb24_plot,NULL,tpb48_plot,ncol=3,nrow=1,widths=c(1,0.05,1),heights=c(1,0.05,1))+theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))
ggsave(tpb_grp,filename="Fig1E.pdf", width = 21.21, height = 9.5, units = "in")



