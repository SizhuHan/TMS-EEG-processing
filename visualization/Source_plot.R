#install.packages("devtools")
library(devtools)
library(ggplot2)
library(reshape2) 
library(gridExtra)
library(ggpubr)
setwd("/Volumes/SizhuFiles/TMS+EEGstudy/Main/visualization/data")

#### plot LMFA-AUC/TEPs: bar plot ####

rm(list=ls(all=TRUE))
color1<-c("#222222","#E42600") # black=hc, red=patients
data <- read.csv("HCvsPats.csv",header = T)
Group <- factor(data[,c(2)])
ratio <- data[,c(8)]
data_temp <- data.frame(Group,ratio)
group=c('HC','Patients')


m1=mean(data_temp[which(data_temp$Group==group[1]),2],na.rm = TRUE)
m2=mean(data_temp[which(data_temp$Group==group[2]),2],na.rm = TRUE)
sd1=sd(data_temp[which(data_temp$Group==group[1]),2],na.rm = TRUE)/sqrt(64) #64 for hc,28 for active
sd2=sd(data_temp[which(data_temp$Group==group[2]),2],na.rm = TRUE)/sqrt(53) #53 for patients, 25 for sham

Ratio = c(m1,m2)
SD = c(sd1,sd2)
group=c('HC','MDD')
data_new = data.frame(group,Ratio,SD)
data_new$group = factor(data_new$group, levels=group)

p <- ggplot(data_new, aes(x=group, y=Ratio, fill=group)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Ratio-SD, ymax=Ratio+SD), width=.05,
                position=position_dodge(.9))

p + labs(title="", 
         x="", y = "P30")+
  scale_fill_manual(values=color1)+
  theme_classic()+theme(legend.position = "none")


#### hc vs. patients: SCD ####
library(ggpubr)
library(gridExtra)
rm(list=ls(all=TRUE))
color1<-c("#222222","#E42600") # black=hc, red=patients
data <- read.csv("HCvsPats.csv",header = T)
data$Group <- factor(data$Group,levels=c('HC','Patients'))
Group <- c(rep("HC",64),rep("MDD",53)) 


## hc vs. patients: SCD ##
col <- 15
Value <- c(data[,col])
mydata <- data.frame(Group,Value)
plot1<-ggplot(mydata,aes(Group,Value))+
  geom_violin(aes(fill=Group),trim=FALSE)+
  geom_boxplot(width=0.1)+
  scale_fill_manual(values = color1)+
  scale_y_continuous(limits = c(-0.25,0.5))+
  labs(title="lDLPFC", x="", y = "Reg-SCD")+
  stat_compare_means(comparisons = list(c('HC','MDD')),method = 't.test',label.y=0.45,label='p.signif')+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5,face = "bold",size=15),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10))

col <- 17
Value <- c(data[,col])
mydata <- data.frame(Group,Value)
plot2<-ggplot(mydata,aes(Group,Value))+
  geom_violin(aes(fill=Group),trim=FALSE)+
  geom_boxplot(width=0.1)+
  scale_fill_manual(values = color1)+
  scale_y_continuous(limits = c(-0.25,0.5))+
  labs(title="OFC", x="", y = "")+
  stat_compare_means(comparisons = list(c('HC','MDD')),method = 't.test',label.y=0.45,label='p.signif')+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5,face = "bold",size=15),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10))


col <- 19
Value <- c(data[,col])
mydata <- data.frame(Group,Value)
plot3<-ggplot(mydata,aes(Group,Value))+
  geom_violin(aes(fill=Group),trim=FALSE)+
  geom_boxplot(width=0.1)+
  scale_fill_manual(values = color1)+
  scale_y_continuous(limits = c(-0.25,0.75))+
  labs(title="HPC", x="", y = "")+
  stat_compare_means(comparisons = list(c('HC','MDD')),method = 't.test',label.y=0.7,label='p.signif')+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5,face = "bold",size=15),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10))

col <- 21
Value <- c(data[,col])
mydata <- data.frame(Group,Value)
plot4<-ggplot(mydata,aes(Group,Value))+
  geom_violin(aes(fill=Group),trim=FALSE)+
  geom_boxplot(width=0.1)+
  scale_fill_manual(values = color1)+
  scale_y_continuous(limits = c(-0.25,0.75))+
  labs(title="PPC", x="", y = "")+
  stat_compare_means(comparisons = list(c('HC','MDD')),method = 't.test',label.y=0.7,label='p.signif')+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5,face = "bold",size=15),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10))

grid.arrange(plot1, plot2, plot3,plot4, ncol = 4, nrow = 1, respect = T)

## hc vs. patients: SCS ## 

col <- 16
Value <- c(data[,col])
mydata <- data.frame(Group,Value)
plot5<-ggplot(mydata,aes(Group,Value))+
  geom_violin(aes(fill=Group),trim=FALSE)+
  geom_boxplot(width=0.1)+
  scale_fill_manual(values = color1)+
  scale_y_continuous(limits = c(-40,60))+
  labs(title="lDLPFC", x="", y = "SCS,mm")+
  stat_compare_means(comparisons = list(c('HC','MDD')),method = 'wilcox.test',label.y=50,label='p.signif')+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10))

col <- 18
Value <- c(data[,col])
mydata <- data.frame(Group,Value)
plot6<-ggplot(mydata,aes(Group,Value))+
  geom_violin(aes(fill=Group),trim=FALSE)+
  geom_boxplot(width=0.1)+
  scale_fill_manual(values = color1)+
  scale_y_continuous(limits = c(-60,160))+
  labs(title="OFC", x="", y = "")+
  stat_compare_means(comparisons = list(c('HC','MDD')),method = 'wilcox.test',label.y=140,label='p.signif')+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10))


col <- 20
Value <- c(data[,col])
mydata <- data.frame(Group,Value)
plot7<-ggplot(mydata,aes(Group,Value))+
  geom_violin(aes(fill=Group),trim=FALSE)+
  geom_boxplot(width=0.1)+
  scale_fill_manual(values = color1)+
  scale_y_continuous(limits = c(-60,160))+
  labs(title="HPC", x="", y = "")+
  stat_compare_means(comparisons = list(c('HC','MDD')),method = 'wilcox.test',label.y=140,label='p.signif')+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10))

col <- 22
Value <- c(data[,col])
mydata <- data.frame(Group,Value)
plot8<-ggplot(mydata,aes(Group,Value))+
  geom_violin(aes(fill=Group),trim=FALSE)+
  geom_boxplot(width=0.1)+
  scale_fill_manual(values = color1)+
  scale_y_continuous(limits = c(-60,160))+
  labs(title="PPC", x="", y = "")+
  stat_compare_means(comparisons = list(c('HC','MDD')),method = 'wilcox.test',label.y=140,label='p.signif')+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10))

grid.arrange(plot5, plot6, plot7, plot8, ncol = 4, nrow = 1, respect = T)


#### active vs. sham: LMFA/TEPs ####

library(ggpubr)
rm(list=ls(all=TRUE))
data <- read.csv("PrevsPost.csv",header = T)

data_temp <- data[which(data$Treatment=='Active'),]
pre_active<-data_temp[which(data_temp$Time=='Pre'),]
post_active<-data_temp[which(data_temp$Time=='Post'),]

data_temp <- data[which(data$Treatment=='Sham'),]
pre_sham<-data_temp[which(data_temp$Time=='Pre'),]
post_sham<-data_temp[which(data_temp$Time=='Post'),]

Treatment <- c(rep("Active",28),rep("Sham",25))

col <- 14
Value <- c(post_active[,col]- pre_active[,col],post_sham[,col]- pre_sham[,col])
mydata <- data.frame(Treatment,Value)
ggplot(mydata,aes(Treatment,Value))+
  geom_violin(aes(fill=Treatment),trim=FALSE)+
  geom_boxplot(width=0.15)+
  scale_fill_manual(values = c('#108F40','#848484'))+
  scale_y_continuous(limits = c(-0.5,0.75))+
  labs(title="", x="", y = "LMFA-AUC (post-pre)")+
  stat_compare_means(comparisons = list(c('Active','Sham')),method = 't.test',label.y=0.6,label='p.signif',paired = FALSE,method.args = list(alternative = "greater"))+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10))

#### SCD Plotting: active vs. sham ####

col <- 15
Value <- c(post_active[,col]- pre_active[,col],post_sham[,col]- pre_sham[,col])
mydata <- data.frame(Treatment,Value)
plot1<-ggplot(mydata,aes(Treatment,Value))+
  geom_violin(aes(fill=Treatment),trim=FALSE)+
  geom_boxplot(width=0.1)+
  scale_fill_manual(values = c('#108F40','#848484'))+
  scale_y_continuous(limits = c(-0.5,0.75))+
  labs(title="lDLPFC", x="", y = "Reg-SCD (post - pre)")+
  stat_compare_means(comparisons = list(c('Active','Sham')),method = 't.test',label.y=0.6,label='p.signif',paired = FALSE,method.args = list(alternative = "greater"))+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5,face = "bold",size=15),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10))

col <- 17
Value <- c(post_active[,col]- pre_active[,col],post_sham[,col]- pre_sham[,col])
mydata <- data.frame(Treatment,Value)
plot2<-ggplot(mydata,aes(Treatment,Value))+
  geom_violin(aes(fill=Treatment),trim=FALSE)+
  geom_boxplot(width=0.1)+
  scale_fill_manual(values = c('#108F40','#848484'))+
  scale_y_continuous(limits = c(-0.5,0.75))+
  labs(title="OFC", x="", y = "")+
  stat_compare_means(comparisons = list(c('Active','Sham')),method = 't.test',label.y=0.6,label='p.signif',paired = FALSE,method.args = list(alternative = "greater"))+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5,face = "bold",size=15),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10))


col <- 19
Value <- c(post_active[,col]- pre_active[,col],post_sham[,col]- pre_sham[,col])
mydata <- data.frame(Treatment,Value)
plot3<-ggplot(mydata,aes(Treatment,Value))+
  geom_violin(aes(fill=Treatment),trim=FALSE)+
  geom_boxplot(width=0.1)+
  scale_fill_manual(values = c('#108F40','#848484'))+
  scale_y_continuous(limits = c(-1,1))+
  labs(title="HPC", x="", y = "")+
  stat_compare_means(comparisons = list(c('Active','Sham')),method = 't.test',label.y = 0.80,label='p.signif',paired = FALSE,method.args = list(alternative = "greater"))+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5,face = "bold",size=15),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10))

col <-21
Value <- c(post_active[,col]- pre_active[,col],post_sham[,col]- pre_sham[,col])
mydata <- data.frame(Treatment,Value)
plot4<-ggplot(mydata,aes(Treatment,Value))+
  geom_violin(aes(fill=Treatment),trim=FALSE)+
  geom_boxplot(width=0.1)+
  scale_fill_manual(values = c('#108F40','#848484'))+
  scale_y_continuous(limits = c(-0.5,0.75))+
  labs(title="PPC", x="", y = "")+
  stat_compare_means(comparisons = list(c('Active','Sham')),method = 't.test',label.y=0.6,label='p.signif',paired = FALSE,method.args = list(alternative = "greater"))+
  theme_classic()+
  theme(legend.position = "none",
        plot.title=element_text(hjust=0.5,face = "bold",size=15),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10))

grid.arrange(plot1, plot2,  plot3, plot4, ncol = 4, nrow = 1, respect = T)

## SCD Plotting: pre vs. post ##

data <- read.csv("PrevsPost.csv",header = T)
data <- data[which(data$Treatment=='Active'),]
#data <- data[which(data$Treatment=='Sham'),]
nsub <- nrow(data)/2
Time <- c(rep("Pre",nsub),rep("Post",nsub))
pre_active<-data[which(data$Time=='Pre'),]
post_active<-data[which(data$Time=='Post'),]
paired = c(1:nsub,1:nsub)

col <- 15
Value <- c(pre_active[,col], post_active[,col])
mydata <- data.frame(Time,Value,paired)
mydata$Time <- factor(mydata$Time, levels=c('Pre','Post'))
plot1<-ggplot(mydata,aes(Time,Value,fill=Time))+
  geom_boxplot(width=0.3)+
  scale_fill_manual(values = c('white','white'))+
  geom_point(aes(fill=Time,group=paired),size=2,shape=21)+
  geom_line(aes(group = paired), size=0.5, color='gray', alpha=0.6)+
  scale_y_continuous(limits = c(-0.25,0.5))+
  labs(title="", x="", y = "Reg-SCD")+
  stat_compare_means(comparisons = list(c('Pre','Post')),method = 't.test',label.y=0.25,label='p.format',paired = TRUE,method.args = list(alternative = "less"))+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10))

col <- 17
Value <- c(pre_active[,col], post_active[,col])
mydata <- data.frame(Time,Value,paired)
mydata$Time <- factor(mydata$Time, levels=c('Pre','Post'))
plot2<-ggplot(mydata,aes(Time,Value,fill=Time))+
  geom_boxplot(width=0.3)+
  scale_fill_manual(values = c('white','white'))+
  geom_point(aes(fill=Time,group=paired),size=2,shape=21)+
  geom_line(aes(group = paired), size=0.5, color='gray', alpha=0.6)+
  scale_y_continuous(limits = c(-0.25,0.5))+
  labs(title="", x="", y = "")+
  stat_compare_means(comparisons = list(c('Pre','Post')),method = 't.test',label.y=0.25,label='p.format',paired = TRUE,method.args = list(alternative = "less"))+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10))



col <- 19
Value <- c(pre_active[,col], post_active[,col])
mydata <- data.frame(Time,Value,paired)
mydata$Time <- factor(mydata$Time, levels=c('Pre','Post'))
plot3<-ggplot(mydata,aes(Time,Value,fill=Time))+
  geom_boxplot(width=0.3)+
  scale_fill_manual(values = c('white','white'))+
  geom_point(aes(fill=Time,group=paired),size=2,shape=21)+
  geom_line(aes(group = paired), size=0.5, color='gray', alpha=0.6)+
  scale_y_continuous(limits = c(-0.25,0.5))+
  labs(title="", x="", y = "")+
  stat_compare_means(comparisons = list(c('Pre','Post')),method = 't.test',label.y = 0.25,label='p.format',paired = TRUE,method.args = list(alternative = "less"))+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10))

col <- 21
Value <- c(pre_active[,col], post_active[,col])
mydata <- data.frame(Time,Value,paired)
mydata$Time <- factor(mydata$Time, levels=c('Pre','Post'))
plot4<-ggplot(mydata,aes(Time,Value,fill=Time))+
  geom_boxplot(width=0.3)+
  scale_fill_manual(values = c('white','white'))+
  geom_point(aes(fill=Time,group=paired),size=2,shape=21)+
  geom_line(aes(group = paired), size=0.5, color='gray', alpha=0.6)+
  scale_y_continuous(limits = c(-0.25,0.5))+
  labs(title="", x="", y = "")+
  stat_compare_means(comparisons = list(c('Pre','Post')),method = 't.test',label.y=0.25,label='p.format',paired = TRUE,method.args = list(alternative = "less"))+
  theme_classic()+
  theme(legend.position = "none",
        plot.title=element_text(hjust=0.5),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10))

grid.arrange(plot1, plot2, plot3,plot4,ncol = 4, nrow = 1, respect = T)


## SCS Plotting: active vs. sham ##

col <- 16
Value <- c(post_active[,col]- pre_active[,col],post_sham[,col]- pre_sham[,col])
mydata <- data.frame(Treatment,Value)
plot5<-ggplot(mydata,aes(Treatment,Value))+
  geom_violin(aes(fill=Treatment),trim=FALSE)+
  geom_boxplot(width=0.1)+
  scale_fill_manual(values = c('#108F40','#848484'))+
  scale_y_continuous(limits = c(-60,80))+
  labs(title="", x="", y = "SCS (post-pre)")+
  stat_compare_means(comparisons = list(c('Active','Sham')),method = 't.test',label.y=75,label='p.signif')+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5),
        axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10))

col <- 18
Value <- c(post_active[,col]- pre_active[,col],post_sham[,col]- pre_sham[,col])
mydata <- data.frame(Treatment,Value)
plot6<-ggplot(mydata,aes(Treatment,Value))+
  geom_violin(aes(fill=Treatment),trim=FALSE)+
  geom_boxplot(width=0.1)+
  scale_fill_manual(values = c('#108F40','#848484'))+
  scale_y_continuous(limits = c(-160,160))+
  labs(title="", x="", y = "")+
  stat_compare_means(comparisons = list(c('Active','Sham')),method = 't.test',label.y=150,label='p.signif')+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5),
        axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10))


col <- 20
Value <- c(post_active[,col]- pre_active[,col],post_sham[,col]- pre_sham[,col])
mydata <- data.frame(Treatment,Value)
plot7<-ggplot(mydata,aes(Treatment,Value))+
  geom_violin(aes(fill=Treatment),trim=FALSE)+
  geom_boxplot(width=0.1)+
  scale_fill_manual(values = c('#108F40','#848484'))+
  scale_y_continuous(limits = c(-160,180))+
  labs(title="", x="", y = "")+
  stat_compare_means(comparisons = list(c('Active','Sham')),method = 't.test',label.y = 160,label='p.signif')+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5),
        axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10))

col <- 22
Value <- c(post_active[,col]- pre_active[,col],post_sham[,col]- pre_sham[,col])
mydata <- data.frame(Treatment,Value)
plot8<-ggplot(mydata,aes(Treatment,Value))+
  geom_violin(aes(fill=Treatment),trim=FALSE)+
  geom_boxplot(width=0.1)+
  scale_fill_manual(values = c('#108F40','#848484'))+
  scale_y_continuous(limits = c(-160,200))+
  labs(title="", x="", y = "")+
  stat_compare_means(comparisons = list(c('Active','Sham')),method = 't.test',label.y=180,label='p.signif')+
  theme_classic()+
  theme(legend.position = "none",
        plot.title=element_text(hjust=0.5),
        axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10))

grid.arrange(plot5, plot6,  plot7, plot8, ncol = 4, nrow = 1, respect = T)

## SCS Plotting: pre vs. post ##
library(gridExtra)
data <- read.csv("PrevsPost.csv",header = T)
data_temp <- data[which(data$Treatment=='Active'),]
#data_temp <- data[which(data$Treatment=='Sham'),]
nsub <- nrow(data_temp)/2
Time <- c(rep("Pre",nsub),rep("Post",nsub))
pre_active<-data_temp[which(data_temp$Time=='Pre'),]
post_active<-data_temp[which(data_temp$Time=='Post'),]
paired = c(1:nsub,1:nsub)

col <- 16
Value <- c(pre_active[,col], post_active[,col])
mydata <- data.frame(Time,Value)
mydata$Time <- factor(mydata$Time, levels=c('Pre','Post'))
plot5<-ggplot(mydata,aes(Time,Value,fill=Time))+
  geom_boxplot(width=0.3)+
  scale_fill_manual(values = c('white','white'))+
  geom_point(aes(fill=Time,group=paired),size=2,shape=21)+
  geom_line(aes(group = paired), size=0.5, color='gray', alpha=0.6)+
  scale_y_continuous(limits = c(-60,80))+
  labs(title="", x="", y = "SCS,mm")+
  stat_compare_means(comparisons = list(c('Pre','Post')),method = 't.test',label.y=75,label='p.signif')+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5),
        axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10))

col <- 18
Value <- c(pre_active[,col], post_active[,col])
mydata <- data.frame(Time,Value)
mydata$Time <- factor(mydata$Time, levels=c('Pre','Post'))
plot6<-ggplot(mydata,aes(Time,Value,fill=Time))+
  geom_boxplot(width=0.3)+
  scale_fill_manual(values = c('white','white'))+
  geom_point(aes(fill=Time,group=paired),size=2,shape=21)+
  geom_line(aes(group = paired), size=0.5, color='gray', alpha=0.6)+
  scale_y_continuous(limits = c(-160,160))+
  labs(title="", x="", y = "")+
  stat_compare_means(comparisons = list(c('Pre','Post')),method = 't.test',label.y=150,label='p.signif')+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5),
        axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10))

col <- 20
Value <- c(pre_active[,col], post_active[,col])
mydata <- data.frame(Time,Value)
mydata$Time <- factor(mydata$Time, levels=c('Pre','Post'))
plot7<-ggplot(mydata,aes(Time,Value,fill=Time))+
  geom_boxplot(width=0.3)+
  scale_fill_manual(values = c('white','white'))+
  geom_point(aes(fill=Time,group=paired),size=2,shape=21)+
  geom_line(aes(group = paired), size=0.5, color='gray', alpha=0.6)+
  scale_y_continuous(limits = c(-160,180))+
  labs(title="", x="", y = "")+
  stat_compare_means(comparisons = list(c('Pre','Post')),method = 't.test',label.y = 160,label='p.signif')+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5),
        axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10))

col <- 22
Value <- c(pre_active[,col], post_active[,col])
mydata <- data.frame(Time,Value)
mydata$Time <- factor(mydata$Time, levels=c('Pre','Post'))
plot8<-ggplot(mydata,aes(Time,Value,fill=Time))+
  geom_boxplot(width=0.3)+
  scale_fill_manual(values = c('white','white'))+
  geom_point(aes(fill=Time,group=paired),size=2,shape=21)+
  geom_line(aes(group = paired), size=0.5, color='gray', alpha=0.6)+
  scale_y_continuous(limits = c(-160,200))+
  labs(title="", x="", y = "")+
  stat_compare_means(comparisons = list(c('Pre','Post')),method = 't.test',label.y=180,label='p.signif')+
  theme_classic()+
  theme(legend.position = "none",
        plot.title=element_text(hjust=0.5),
        axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10))

grid.arrange(plot5, plot6,plot7, plot8, ncol = 4, nrow = 1, respect = T)



## two-sample t test ##
# hc vs. patients
data <- read.csv("HCvsPats.csv",header = T)
data$SuicideIdea <- factor(data$SuicideIdea,levels=c('HC','NO','YES'))
data$Course <- factor(data$Course,levels=c('HC','Short','Long'))

hc<-data[which(data$Group=='HC'),c(22)]
pat<-data[which(data$Group=='Patients'),c(22)]

mean_hc <- mean(hc)
sd_hc <- sd(hc)
mean_pat <- mean(pat)
sd_pat <- sd(pat)

t.test(hc,pat,paired=FALSE)
wilcox.test(hc,pat,paired=FALSE)

library(rcompanion)
wilcoxonR(x = hc,g = pat) # effect size r

library(effsize)
cliff.delta(d=hc,f=pat)
cohen.d( pat,hc)

rm(list=ls(all=TRUE))
data <- read.csv("PrevsPost.csv",header = T)
n <- 29
data_temp <- data[which(data$Treatment=='Active'),]
pre_active<-data_temp[which(data_temp$Time=='Pre'),]
post_active<-data_temp[which(data_temp$Time=='Post'),]
active_change = pre_active[,c(n)]-post_active[,c(n)]

mean_pre_active <- mean(pre_active[,c(n)])
sd_pre_active <- sd(pre_active[,c(n)])
mean_post_active <- mean(post_active[,c(n)])
sd_post_active <- sd(post_active[,c(n)])


data_temp <- data[which(data$Treatment=='Sham'),]
pre_sham<-data_temp[which(data_temp$Time=='Pre'),]
post_sham<-data_temp[which(data_temp$Time=='Post'),]
sham_change = pre_sham[,c(n)]-post_sham[,c(n)]

mean_pre_sham <- mean(pre_sham[,c(n)])
sd_pre_sham <- sd(pre_sham[,c(n)])
mean_post_sham <- mean(post_sham[,c(n)])
sd_post_sham <- sd(post_sham[,c(n)])

t.test(pre_active[,c(n)],pre_sham[,c(n)],paired=FALSE)
t.test(post_active[,c(n)],post_sham[,c(n)],paired=FALSE)

t.test(pre_active[,c(n)],post_active[,c(n)],paired=TRUE,alternative = 'less')
cohen.d(pre_active[,c(n)],post_active[,c(n)])

active_sd <- sd(active_change)
sham_sd <- sd(sham_change)
t.test(active_change,sham_change,paired=FALSE,alternative = 'less')
wilcox.test(active_change,sham_change,paired=FALSE)

cohen.d(active_change,sham_change)



#### plot baseline correltations ####
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(ggExtra)

rm(list=ls(all=TRUE))
color1<-c("#222222","#DC2421") # black=hc, red=patients
data <- read.csv("HCvsPats.csv",header = T)

data3 <- data[,c(2,5,14)] # lmfa & hdrs
data4 <- data[,c(2,6,14)] # lmfa & hama
data5 <- data[,c(2,7,14)] # lmfa& psqi


plots<- ggscatterhist(
  data3,  x ='LMFA.AUC', y = 'HDRS', 
  shape=21,color ="black",
  fill = "Group", # comment this line to plot all data in one color
  size =3, alpha = 0.8,
  palette = color1, # comment this line to plot all data in one color
  margin.plot = "density", #"histogram" or "density"
  margin.params = list(fill = "Group", color = "black", size = 0.2), # comment this line to plot all data in one color
  legend = "top",
  ggtheme = theme_minimal(),
  xlab = "LMFA-AUC",
  ylab = "HDRS",
  add = "reg.line",    # Add regression line
  conf.int = TRUE,     # Add confidence interval
  add.params = list(color = "black",  # using specific color (e.g., black) or group by subjects (e.g., Group)
                    fill = "lightgray"),
  )

plots$sp <- plots$sp+stat_cor(method = "pearson",label.x=0)  # uncomment this line to plot all data in one color
#plots$sp <- plots$sp+stat_cor(aes(color = Group)) # group by subjects

plots


#### plot treatment changes correlations ####
rm(list=ls(all=TRUE))
data_trt <- read.csv("PrevsPost.csv",header = T)
data_temp <- data_trt[which(data_trt$Treatment=='Active'),]
pre_active<-data_temp[which(data_temp$Time=='Pre'),]
post_active<-data_temp[which(data_temp$Time=='Post'),]

data_temp <- data_trt[which(data_trt$Treatment=='Sham'),]
pre_sham<-data_temp[which(data_temp$Time=='Pre'),]
post_sham<-data_temp[which(data_temp$Time=='Post'),]

LMFA.AUC_Changes <- c(pre_active$LMFA.AUC - post_active$LMFA.AUC,pre_sham$LMFA.AUC - post_sham$LMFA.AUC)

HDRS_Changes <- c(pre_active$HDRS - post_active$HDRS,pre_sham$HDRS - post_sham$HDRS)
HAMA_Changes <- c(pre_active$HAMA - post_active$HAMA,pre_sham$HAMA - post_sham$HAMA)
PSQI_Changes <- c(pre_active$PSQI - post_active$PSQI,pre_sham$PSQI - post_sham$PSQI)

lDLPFC_Changes <- c(post_active$lDLPFC.SCD - pre_active$lDLPFC.SCD,post_sham$lDLPFC.SCD - pre_sham$lDLPFC.SCD)
Hipp_Changes <- c(pre_active$Hipp.SCD - post_active$Hipp.SCD,pre_sham$Hipp.SCD - post_sham$Hipp.SCD)
OFC_Changes <- c(post_active$OFC.SCD - pre_active$OFC.SCD,post_sham$OFC.SCD - pre_sham$OFC.SCD)


data5 <- data.frame(LMFA.AUC_Changes,HDRS_Changes,Treatment=c(rep("Active",28),rep("Sham",25)))
data6 <- data.frame(Hipp_Changes,HDRS_Changes,Treatment=c(rep("Active",28),rep("Sham",25)))

data7 <- data.frame(LMFA.AUC_Changes,HAMA_Changes,Treatment=c(rep("Active",28),rep("Sham",25))) # hc & pats
data8 <- data.frame(Hipp_Changes,HAMA_Changes,Treatment=c(rep("Active",28),rep("Sham",25)))

data9 <- data.frame(LMFA.AUC_Changes,PSQI_Changes,Treatment=c(rep("Active",28),rep("Sham",25)))
data10 <- data.frame(Hipp_Changes,PSQI_Changes,Treatment=c(rep("Active",28),rep("Sham",25)))


###### significance plotting #####
color5<-c('#F37E6F','#F6BE7A')
mydata <- data5
ggplot(mydata, aes(LMFA.AUC_Changes, HDRS_Changes))+
  geom_point(aes(colour=factor(Treatment)))+
  scale_color_manual(values=color5)+
  stat_smooth(method='lm',color='black',alpha=0.2)+stat_cor(method = "pearson",label.y =-3)+
  theme_classic()+
  xlab("Changes in LMFA-AUC (pre-post)")+
  ylab("Changes in HDRS (pre-post)" )+
  theme(legend.position = "none")


mydata <- data6
ggplot(mydata, aes(Hipp_Changes, HDRS_Changes,colour=factor(Treatment)))+
  geom_point()+
  scale_color_manual(values=color5)+
  stat_smooth(method='lm',alpha=0.2)+stat_cor(method = "pearson",label.y =c(-1.5,-3.5))+
  theme_classic()+
  xlab("Changes in Hipp-SCD (pre-post)")+
  ylab("Changes in HDRS (pre-post)" )+
  theme(legend.position = "none")

## control analysis on HAMA & PSQI ##

mydata <- data7
k1<-ggplot(mydata, aes(LMFA.AUC_Changes, HAMA_Changes))+
  geom_point(aes(colour=factor(Treatment)))+
  scale_color_manual(values=color5)+
  stat_smooth(method='lm',color='black',alpha=0.2)+stat_cor(method = "pearson",label.y =-3)+
  theme_classic()+
  xlab("Changes in LMFA-AUC (pre-post)")+
  ylab("Changes in HAMA (pre-post)" )+
  theme(legend.position = "none")


k2<-ggplot(data8, aes(Hipp_Changes, HAMA_Changes,colour=factor(Treatment)))+
  geom_point()+
  scale_color_manual(values=color5)+
  stat_smooth(method='lm',alpha=0.2)+stat_cor(method = "pearson",label.y =c(-1.5,-3.5))+
  theme_classic()+
  xlab("Changes in Hipp-SCD (pre-post)")+
  ylab("Changes in HAMA (pre-post)" )+
  theme(legend.position = "none")

mydata <- data9
k3<-ggplot(mydata, aes(LMFA.AUC_Changes, PSQI_Changes))+
  geom_point(aes(colour=factor(Treatment)))+
  scale_color_manual(values=color5)+
  stat_smooth(method='lm',color='black',alpha=0.2)+stat_cor(method = "pearson",label.y =-3)+
  theme_classic()+
  xlab("Changes in LMFA-AUC (pre-post)")+
  ylab("Changes in PSQI (pre-post)" )+
  theme(legend.position = "none")

k4<-ggplot(data10, aes(Hipp_Changes, PSQI_Changes,colour=factor(Treatment)))+
  geom_point()+
  scale_color_manual(values=color5)+
  stat_smooth(method='lm',alpha=0.2)+stat_cor(method = "pearson",label.y =c(-1.5,-3.5))+
  theme_classic()+
  xlab("Changes in Hipp-SCD (pre-post)")+
  ylab("Changes in PSQI (pre-post)" )+
  theme(legend.position = "none")

grid.arrange(k1, k2, k3, k4, ncol = 2, nrow = 2, respect = T)

#### structural equation modeling ####
#install.packages("lavaan")
library(lavaan)
range = c(32:59)
lDLPFC<-lDLPFC_Changes[range] - mean(lDLPFC_Changes[range])
OFC<-OFC_Changes[range] - mean(OFC_Changes[range])
Hipp<-Hipp_Changes[range] - mean(Hipp_Changes[range])

my_data <- data.frame(lDLPFC,OFC,Hipp)

model0 <- 'Hipp ~ c * lDLPFC'

model1 <- 'OFC ~ a * lDLPFC
          Hipp ~  c * lDLPFC + b * OFC
          
          ie := a * b
          total := c + (a*b) 
          '

fit0 <- sem(model0, data = my_data)
summary(fit0)
standardizedsolution(fit0)

fit1 <- sem(model1, data = my_data)
summary(fit1)
standardizedsolution(fit1)


#### compute & plot ROC Curve ####
library(ROCR)
library(pROC)
library(GMCM)
library(car)
library(e1071)

rm(list=ls(all=TRUE))
data <- read.csv("HCvsPats.csv",header = T)
#data <- data[which(data$lDLPFCaddHipp > 0),]
#temp = data[,c(2,16)] # p60/n100
#temp = data[,c(2,16,18)] # p60/n100, lmfa,
temp = data[,c(2,16,18,21)] # p60/n100, lmfa,ofc_scd,
#temp = data[,c(2,16,18,21,19)] # p60/n100, lmfa,ofc_scd,dlpfc_scd
#temp = data[,c(2,16,18,21,19,27)] # p60/n100, lmfa,ofc_scd,dlpfc_scd,hpc_scd
#temp = data[,c(2,16,18,21,19,27,22)] # p60/n100, lmfa, ofc_scd,dlpfc_scd,hpc_scd,ofc_scs


# normalization
for (i in 2:length(temp)) {
   temp[,i] = (temp[,i] - mean(temp[,i]))/sd(temp[,i])
}

df = temp
nt=1000 # number of resample
x_all <- matrix(NaN,nrow = nt,ncol = 37) 
y_all <- matrix(NaN,nrow = nt,ncol = 37) 
auc_all <- matrix(NaN, nrow = nt, ncol =1)
rand_x_all <- matrix(NaN,nrow = nt,ncol = 37) 
rand_y_all <- matrix(NaN,nrow = nt,ncol = 37) 
rand_auc_all <- matrix(NaN, nrow = nt, ncol =1)
for (i in 1:nt) {
index = sample(1:nrow(df), size = .7 * nrow(df)) # resample
train = df[index,]
test = df[-index,]
rand_test = test

# shuffle test label
rand_test_ind <- sample(1:nrow(test))
Group <- test[rand_test_ind,c(1)]
rand_test$Group <- Group


model = glm(as.factor(Group)~.,family = binomial(link='logit'),data=train)
#model = svm(Group~.,data=train,kernel = "linear", cost = 10)
#vif(model) # variation inflation factors

#summary(model)
pred = predict(model,test,type="response")
modelroc <- roc(test$Group,pred)
#modelroc <- roc(test$Group,as.numeric(pred))


x_all[i,] = modelroc[["specificities"]]
y_all[i,] = modelroc[["sensitivities"]]
auc_all[i] = modelroc[["auc"]]


rand_pred = predict(model,rand_test,type="response")
rand_modelroc <- roc(rand_test$Group,rand_pred)
#rand_modelroc <- roc(rand_test$Group,as.numeric(rand_pred))

rand_x_all[i,] = rand_modelroc[["specificities"]]
rand_y_all[i,] = rand_modelroc[["sensitivities"]]
rand_auc_all[i] = rand_modelroc[["auc"]]

}

modelroc[["specificities"]] = colMeans(x_all) # average
modelroc[["sensitivities"]] = colMeans(y_all) # average

x_mean = 1-colMeans(x_all)
y_mean = colMeans(y_all)
y_sd = GMCM:::colSds(y_all)
auc_mean = colMeans(auc_all)

rand_modelroc[["specificities"]] = colMeans(rand_x_all) # average
rand_modelroc[["sensitivities"]] = colMeans(rand_y_all) # average

rand_x_mean = 1-colMeans(rand_x_all)
rand_y_mean = colMeans(rand_y_all)
rand_y_sd = GMCM:::colSds(rand_y_all)
rand_auc_mean = colMeans(rand_auc_all)


plot(NULL,xlim=c(0,1),ylim=c(0,1),xlab="1-Specificity",ylab="Sensitivity",axes=T,las=1)
# axis(side=1,at=seq(0,1,0.2),
#     labels=seq(0,1,0.2), cex.axis=1)
# axis(side=2,at=seq(0,1,0.2),
#      labels=seq(0,1,0.2), las=1,cex.axis=1)

lines(x=1-modelroc[["specificities"]],y=modelroc[["sensitivities"]],
     col = rgb(0.03,0.32,0.61,1),lwd=3)
polygon(c(x_mean, rev(x_mean)),
        c(y_mean-y_sd,rev(y_mean+y_sd)), col =rgb(0.03,0.32,0.61,0.3),border = NA)
text(0.7,0.4,paste("Real AUC = ", round(auc_mean,3)), cex=1.0, xlab=0.5,col=rgb(0.03,0.32,0.61,1))

lines(1-rand_modelroc[["specificities"]],rand_modelroc[["sensitivities"]], type ="l",
     col = rgb(0.8,0.8,0.8,1),lwd=3)
polygon(c(rand_x_mean, rev(rand_x_mean)),
        c(rand_y_mean-rand_y_sd,rev(rand_y_mean+rand_y_sd)), col =rgb(0.8,0.8,0.8,0.3),border = NA)
text(0.7,0.3,paste("Shuffled AUC = ", round(rand_auc_mean,3)), cex=1.0, xlab=0.5,col=rgb(0.8,0.8,0.8,1))
title(main="The winner model",cex.main = 2)


#### save & plotting ####
model1 <- c(auc_all,rand_auc_all)
model2 <- c(auc_all,rand_auc_all)
model3 <- c(auc_all,rand_auc_all)
model4 <- c(auc_all,rand_auc_all)
model5 <- c(auc_all,rand_auc_all)
model6 <- c(auc_all,rand_auc_all)
type <- rep(c(rep("Real",1000),rep("Shuffled",1000)),6)
model <- c(rep("model1",2000),rep("model2",2000),rep("model3",2000),rep("model4",2000),rep("model5",2000),rep("model6",2000))
value <- c(model1,model2,model3,model4,model5,model6)
AUC_data <- data.frame(type,model,value)
write.csv(AUC_data,"AUC_data.csv",row.names = TRUE)

AUC_data <-read.csv("AUC_data.csv")
real_auc <- AUC_data[which(AUC_data$type=="Real"),]
real_auc$model <- c(rep("Model1",1000),rep("Model2",1000),rep("Model3",1000),rep("Model4",1000),rep("Model5",1000),rep("Model6",1000))

ggplot(AUC_data,aes(x=model,y=value))+
  geom_violin(aes(fill=type),trim=FALSE, position = position_dodge(0.8))+
  geom_boxplot(aes(fill=type),width=0.1, position = position_dodge(0.8))+
  #scale_fill_discrete()+
  scale_fill_manual(values = c('#6baed6','white'))+
  scale_y_continuous(limits = c(0.25,1.2))+
  labs(title="", x="", y = "AUC")+
  stat_compare_means(aes(group=type),method = 't.test',label.y = 1.05,label='p.signif')+
  theme_classic()+
  theme(plot.title=element_text(hjust=0.5),
        axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10))


ggplot(real_auc,aes(x=model,y=value))+
  geom_violin(aes(fill=model),trim=FALSE)+
  geom_boxplot(width=0.1)+
  #scale_fill_brewer(palette = "Blues")+
  scale_fill_manual(values = c('#f7f7f7','#636363','#252525','#969696','#bdbdbd','#d9d9d9'))+
  scale_y_continuous(limits = c(0.25,1.3))+
  labs(title="Model comparison", x="", y = "AUC")+
  stat_compare_means(comparisons = list(c('Model1','Model2'),c('Model2','Model3'),c('Model1','Model3'),c('Model3','Model4'),c('Model3','Model5'),c('Model3','Model6')),
                     method = 't.test',label.y = c(1.0,1.03,1.1,1.03,1.1,1.17),label='p.signif',size=5)+
  theme_classic()+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5,face = "bold",size=25),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20,face = "bold"),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.line.x=element_line(size=1),
        axis.line.y=element_line(size=1))



