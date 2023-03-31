library(ggplot2)
library(ggpubr)
library(gridExtra)

setwd("/Volumes/SizhuFiles/TMS+EEGstudy/Main/visualization/data")


## Oscillation: hc vs. patients ##
rm(list=ls(all=TRUE))
data <- read.csv("HCvsPats.csv",header = T)
Group <- c(rep(c(rep("HC",64),rep("MDD",53)),5))
Freq <- c(rep("Delta",117),rep("Theta",117),rep("Alpha",117),rep("Beta",117),rep("Gamma",117))
dist_cat_n <- c(rep(1,117),rep(2,117),rep(3,117),rep(4,117),rep(5,117))

col= 23
Value <- c(data[,col],data[,col+1],data[,col+2],data[,col+3],data[,col+4])
mydata <- data.frame(Group,Freq,Value,dist_cat_n)
mydata<-transform(mydata,scat_adj=ifelse(Group=="HC",-0.2,0.2))
mydata$Group <- factor(mydata$Group, levels=c('HC','MDD'))
mydata$Freq <- factor(mydata$Freq, levels=c('Delta','Theta','Alpha','Beta','Gamma'))
plot1<-ggplot(mydata,aes(Freq,Value,fill=factor(Group)))+
  geom_boxplot(outlier.size=0,aes(fill=factor(Group)),position = position_dodge(0.8),size=0.4)+
  geom_jitter(aes(dist_cat_n+scat_adj,Value,fill=factor(Group)),position=position_jitter(width=0.1,height=0),shape=21,size=1)+
  scale_fill_manual(values = c("white","white"))+
  scale_y_continuous(limits = c(0,0.25))+
  labs(title="lDLPFC", x="", y = "Reg-Amplitude")+
  stat_compare_means(aes(group=Group),method = 't.test',label.y=0.21,label='p.signif')+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5,face = "bold"),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10))


col= 28
Value <- c(data[,col],data[,col+1],data[,col+2],data[,col+3],data[,col+4])
mydata <- data.frame(Group,Freq,Value,dist_cat_n)
mydata<-transform(mydata,scat_adj=ifelse(Group=="HC",-0.2,0.2))
mydata$Group <- factor(mydata$Group, levels=c('HC','MDD'))
mydata$Freq <- factor(mydata$Freq, levels=c('Delta','Theta','Alpha','Beta','Gamma'))
plot2<-ggplot(mydata,aes(Freq,Value,fill=factor(Group)))+
  geom_boxplot(outlier.size=0,aes(fill=factor(Group)),position = position_dodge(0.8),size=0.4)+
  geom_jitter(aes(dist_cat_n+scat_adj,Value,fill=factor(Group)),position=position_jitter(width=0.1,height=0),shape=21,size=1)+
  scale_fill_manual(values = c("white","white"))+
  scale_y_continuous(limits = c(0,0.3))+
  labs(title="OFC", x="", y = "")+
  stat_compare_means(aes(group=Group),method = 't.test',label.y=0.25,label='p.signif')+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5,face = "bold"),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10))


col= 33
Value <- c(data[,col],data[,col+1],data[,col+2],data[,col+3],data[,col+4])
mydata <- data.frame(Group,Freq,Value,dist_cat_n)
mydata<-transform(mydata,scat_adj=ifelse(Group=="HC",-0.2,0.2))
mydata$Group <- factor(mydata$Group, levels=c('HC','MDD'))
mydata$Freq <- factor(mydata$Freq, levels=c('Delta','Theta','Alpha','Beta','Gamma'))
plot3<-ggplot(mydata,aes(Freq,Value,fill=factor(Group)))+
  geom_boxplot(outlier.size=0,aes(fill=factor(Group)),position = position_dodge(0.8),size=0.4)+
  geom_jitter(aes(dist_cat_n+scat_adj,Value,fill=factor(Group)),position=position_jitter(width=0.1,height=0),shape=21,size=1)+
  scale_fill_manual(values = c("white","white"))+
  scale_y_continuous(limits = c(0,0.3))+
  labs(title="HPC", x="", y = "")+
  stat_compare_means(aes(group=Group),method = 't.test',label.y=0.25,label='p.signif')+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5,face = "bold"),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10))

col= 38
Value <- c(data[,col],data[,col+1],data[,col+2],data[,col+3],data[,col+4])
mydata <- data.frame(Group,Freq,Value,dist_cat_n)
mydata<-transform(mydata,scat_adj=ifelse(Group=="HC",-0.2,0.2))
mydata$Group <- factor(mydata$Group, levels=c('HC','MDD'))
mydata$Freq <- factor(mydata$Freq, levels=c('Delta','Theta','Alpha','Beta','Gamma'))
plot4<-ggplot(mydata,aes(Freq,Value,fill=factor(Group)))+
  geom_boxplot(outlier.size=0,aes(fill=factor(Group)),position = position_dodge(0.8),size=0.4)+
  geom_jitter(aes(dist_cat_n+scat_adj,Value,fill=factor(Group)),position=position_jitter(width=0.1,height=0),shape=21,size=1)+
  scale_fill_manual(values = c("white","white"))+
  scale_y_continuous(limits = c(0,0.3))+
  labs(title="PPC", x="", y = "")+
  stat_compare_means(aes(group=Group),method = 't.test',label.y=0.25,label='p.signif')+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5,face = "bold"),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10))

grid.arrange(plot1, plot2, plot3,plot4,ncol = 4, nrow = 1, respect = T)



## Oscillation: pre vs. post ##
rm(list=ls(all=TRUE))
data <- read.csv("PrevsPost_Oscillation.csv",header = T)
data <- data[which(data$Treatment=='Active'),]
#data <- data[which(data$Treatment=='Sham'),]
nsub <- nrow(data)
Time <- c(rep(c(rep("Pre",nsub/2),rep("Post",nsub/2)),5))
Freq <- c(rep("Delta",nsub),rep("Theta",nsub),rep("Alpha",nsub),rep("Beta",nsub),rep("Gamma",nsub))
dist_cat_n <- c(rep(1,nsub),rep(2,nsub),rep(3,nsub),rep(4,nsub),rep(5,nsub))

col= 8
Value <- c(data[,col],data[,col+1],data[,col+2],data[,col+3],data[,col+4])
mydata <- data.frame(Time,Freq,Value,dist_cat_n)
mydata<-transform(mydata,scat_adj=ifelse(Time=="Pre",-0.2,0.2))
mydata$Time <- factor(mydata$Time, levels=c('Pre','Post'))
mydata$Freq <- factor(mydata$Freq, levels=c('Delta','Theta','Alpha','Beta','Gamma'))
plot1<-ggplot(mydata,aes(Freq,Value,fill=factor(Time)))+
  geom_boxplot(outlier.size=0,aes(fill=factor(Time)),position = position_dodge(0.8),size=0.4)+
  geom_jitter(aes(dist_cat_n+scat_adj,Value,fill=factor(Time)),position=position_jitter(width=0.1,height=0),shape=21,size=1)+
  scale_fill_manual(values = c('white','white'))+
  scale_y_continuous(limits = c(0,0.25))+
  labs(title="lDLPFC", x="", y = "Reg-Amplitude")+
  stat_compare_means(aes(group=Time),method = 't.test',label.y=0.18,label='p.signif',paired = TRUE,method.args = list(alternative = "greater"))+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5,face = "bold"),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10))


col= 13
Value <- c(data[,col],data[,col+1],data[,col+2],data[,col+3],data[,col+4])
mydata <- data.frame(Time,Freq,Value,dist_cat_n)
mydata<-transform(mydata,scat_adj=ifelse(Time=="Pre",-0.2,0.2))
mydata$Time <- factor(mydata$Time, levels=c('Pre','Post'))
mydata$Freq <- factor(mydata$Freq, levels=c('Delta','Theta','Alpha','Beta','Gamma'))
plot2<-ggplot(mydata,aes(Freq,Value,fill=factor(Time)))+
  geom_boxplot(outlier.size=0,aes(fill=factor(Time)),position = position_dodge(0.8),size=0.4)+
  geom_jitter(aes(dist_cat_n+scat_adj,Value,fill=factor(Time)),position=position_jitter(width=0.1,height=0),shape=21,size=1)+
  scale_fill_manual(values = c('white','white'))+
  scale_y_continuous(limits = c(0,0.25))+
  labs(title="OFC", x="", y = "")+
  stat_compare_means(aes(group=Time),method = 't.test',label.y=0.18,label='p.signif',paired = TRUE,method.args = list(alternative = "greater"))+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5,face = "bold"),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10))


col= 18
Value <- c(data[,col],data[,col+1],data[,col+2],data[,col+3],data[,col+4])
mydata <- data.frame(Time,Freq,Value,dist_cat_n)
mydata<-transform(mydata,scat_adj=ifelse(Time=="Pre",-0.2,0.2))
mydata$Time <- factor(mydata$Time, levels=c('Pre','Post'))
mydata$Freq <- factor(mydata$Freq, levels=c('Delta','Theta','Alpha','Beta','Gamma'))
plot3<-ggplot(mydata,aes(Freq,Value,fill=factor(Time)))+
  geom_boxplot(outlier.size=0,aes(fill=factor(Time)),position = position_dodge(0.8),size=0.4)+
  geom_jitter(aes(dist_cat_n+scat_adj,Value,fill=factor(Time)),position=position_jitter(width=0.1,height=0),shape=21,size=1)+
  scale_fill_manual(values = c('white','white'))+
  scale_y_continuous(limits = c(0,0.25))+
  labs(title="HPC", x="", y = "")+
  stat_compare_means(aes(group=Time),method = 't.test',label.y=0.18,label='p.signif',paired = TRUE,method.args = list(alternative = "greater"))+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5,face = "bold"),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10))

col= 23
Value <- c(data[,col],data[,col+1],data[,col+2],data[,col+3],data[,col+4])
mydata <- data.frame(Time,Freq,Value,dist_cat_n)
mydata<-transform(mydata,scat_adj=ifelse(Time=="Pre",-0.2,0.2))
mydata$Time <- factor(mydata$Time, levels=c('Pre','Post'))
mydata$Freq <- factor(mydata$Freq, levels=c('Delta','Theta','Alpha','Beta','Gamma'))
plot4<-ggplot(mydata,aes(Freq,Value,fill=factor(Time)))+
  geom_boxplot(outlier.size=0,aes(fill=factor(Time)),position = position_dodge(0.8),size=0.4)+
  geom_jitter(aes(dist_cat_n+scat_adj,Value,fill=factor(Time)),position=position_jitter(width=0.1,height=0),shape=21,size=1)+
  scale_fill_manual(values = c('white','white'))+
  scale_y_continuous(limits = c(0,0.25))+
  labs(title="PPC", x="", y = "")+
  stat_compare_means(aes(group=Time),method = 't.test',label.y=0.18,label='p.signif',paired = TRUE,method.args = list(alternative = "greater"))+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5,face = "bold"),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10))

grid.arrange(plot1, plot2, plot3,plot4,ncol = 4, nrow = 1, respect = T)

## plot correlations
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(ggExtra)

rm(list=ls(all=TRUE))
data_trt <- read.csv("/Volumes/SizhuFiles/TMS+EEGstudy/Source2_Data/PrevsPost_Oscillation.csv",header = T)

data_temp <- data_trt[which(data_trt$Treatment=='Active'),]
pre_active<-data_temp[which(data_temp$Time=='Pre'),]
post_active<-data_temp[which(data_temp$Time=='Post'),]

data_temp <- data_trt[which(data_trt$Treatment=='Sham'),]
pre_sham<-data_temp[which(data_temp$Time=='Pre'),]
post_sham<-data_temp[which(data_temp$Time=='Post'),]

ind <- 9
HDRS_Changes <- c(pre_active$HDRS - post_active$HDRS,pre_sham$HDRS - post_sham$HDRS)
lDLPFC_Changes <- c(pre_active[,c(ind)] - post_active[,c(ind)],pre_sham[,c(ind)] - post_sham[,c(ind)])
OFC_Changes <- c(pre_active[,c(ind+5)] - post_active[,c(ind+5)],pre_sham[,c(ind+5)] - post_sham[,c(ind+5)])
Hipp_Changes <- c(pre_active[,c(ind+10)] - post_active[,c(ind+10)],pre_sham[,c(ind+10)] - post_sham[,c(ind+10)])


data1 <- data.frame(Hipp_Changes,HDRS_Changes,Treatment=c(rep("Active",28),rep("Sham",25)))
my_data1 <- data1
data2 <- data.frame(lDLPFC_Changes,HDRS_Changes,Treatment=c(rep("Active",28),rep("Sham",25)))
my_data2 <- data2
data3 <- data.frame(OFC_Changes,Hipp_Changes,Treatment=c(rep("Active",28),rep("Sham",25)))
my_data3 <- data3

color5<-c('#F37E6F','#F6BE7A')
ggplot(my_data1, aes(Hipp_Changes, HDRS_Changes,colour=factor(Treatment)))+
  geom_point()+
  scale_color_manual(values=color5)+
  stat_smooth(method='lm',alpha=0.2)+stat_cor(method = "pearson",label.y =c(30,28))+
  theme_classic()+
  xlab("Changes in HPC_delta (pre-post)")+
  ylab("Changes in HDRS (pre-post)" )+
  theme(axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),legend.position = "none")


#### structural equation modeling ####

library(lavaan)
range = c(29:53)
lDLPFC<-lDLPFC_Changes[range] - mean(lDLPFC_Changes[range])
OFC<-OFC_Changes[range] - mean(OFC_Changes[range])
Hipp<-Hipp_Changes[range] - mean(Hipp_Changes[range])
my_data <- data.frame(lDLPFC,OFC,Hipp)

model0 <- 'Hipp ~  c0 * lDLPFC '

model1 <- 'OFC ~ a * lDLPFC
          Hipp ~  c * lDLPFC + b * OFC
          
          ie := a * b
          total := c + (a*b) 
          '
fit0 <- sem(model0, data = my_data)
summary(fit0)
standardizedsolution(fit0)



fit <- sem(model1, data = my_data)
summary(fit)
standardizedsolution(fit)


