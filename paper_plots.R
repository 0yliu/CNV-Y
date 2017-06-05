#########################
## het density and phase concordance plot
rm(list = ls())
library(ggplot2)
library(ggExtra)

dat_pc_Cran = read.table("~/Dropbox/MD_Anderson/Projects/PaulScheet/CNV/Craniofacial_data/reRun_2016/update-V1/spCrossover/phaseConcordance_spCrossover_Cran.txt",header=TRUE,stringsAsFactors = FALSE)
dat_pc_AGP = read.table("~/Dropbox/MD_Anderson/Projects/PaulScheet/CNV/Autism_data/reRun_2016/spCrossover/phaseConcordance_spCrossover_AGP.txt",header=TRUE,stringsAsFactors=FALSE)

dat_hd_Cran = read.table("~/Dropbox/MD_Anderson/Projects/PaulScheet/CNV/Craniofacial_data/reRun_2016/update-V1/spCrossover/hetDensity_spCrossover_Cran.txt",header=TRUE,stringsAsFactors = FALSE)
dat_hd_AGP = read.table("~/Dropbox/MD_Anderson/Projects/PaulScheet/CNV/Autism_data/reRun_2016/spCrossover/hetDensity_spCrossover_AGP.txt",header=TRUE,stringsAsFactors=FALSE)

dat_Cran = read.table("~/Dropbox/MD_Anderson/Projects/PaulScheet/CNV/Craniofacial_data/reRun_2016/update-V1/spCrossover/outputDenovoList_Cran.txt",header=TRUE,stringsAsFactors = FALSE)
dat_AGP = read.table("~/Dropbox/MD_Anderson/Projects/PaulScheet/CNV/Autism_data/reRun_2016/spCrossover/outputDenovoList_AGP.txt",header=TRUE,stringsAsFactors=FALSE)

if (all(dat_pc_Cran$sample_name==dat_Cran$kid_name)){
  dat_Cran$phaseConcordance = dat_pc_Cran$phaseConcordance
  dat_Cran$lengthX = dat_pc_Cran$length.x.
}
if (all(dat_pc_AGP$sample_name==dat_AGP$kid_name)){
  dat_AGP$phaseConcordance = dat_pc_AGP$phaseConcordance
  dat_AGP$lengthX = dat_pc_AGP$length.x.
}

if (all(dat_hd_Cran$sample_name==dat_Cran$kid_name)){
  dat_Cran$expectedHD_gtype = dat_hd_Cran$expectedHD_gtype
  dat_Cran$hetD_gtype = dat_hd_Cran$hetD_gtype
  dat_Cran$num_Mkrs = dat_hd_Cran$num_Mkrs
}
if (all(dat_hd_AGP$sample_name==dat_AGP$kid_name)){
  dat_AGP$expectedHD_gtype = dat_hd_AGP$expectedHD_gtype
  dat_AGP$hetD_gtype = dat_hd_AGP$hetD_gtype
  dat_AGP$num_Mkrs = dat_hd_AGP$num_Mkrs
}
dat_Cran$sizeOfCNV = round(dat_Cran$sizeOfCNV/1000)

dat = rbind(dat_Cran,dat_AGP)
dat$dataSet = factor(c(rep('Craniofacial',dim(dat_Cran)[1]),rep('Autism',dim(dat_AGP)[1])), c('Craniofacial','Autism'))
dat$category = cut(dat$pu1_all, c(-0.01,0.1,0.89,1), labels=c("bi-allelic","undefined","tri-allelic"))

# het density plot
ggplot(dat,aes(x=pu1_all,y=expectedHD_gtype,shape=dataSet)) +
  geom_point(size=5, alpha=0.7, colour="#003300") +
  xlim(0,1) + 
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        strip.text=element_text(size=16,face="bold"),
        legend.position="bottom",
        legend.text=element_text(size=10),
        legend.title=element_text(size=10),
        legend.background = element_rect(colour = "black"),
        panel.margin=unit(3,"lines")) +
  labs(x="probability of tri-allelic dups", y="Normalized Het Density") 

# phase concordance plot
ggplot(dat,aes(x=pu1_all,y=phaseConcordance,shape=dataSet)) +
  geom_point(size=5, alpha=0.7, colour="#660066") +
  xlim(0,1) + 
  ylim(0,1) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        strip.text=element_text(size=16,face="bold"),
        legend.position="bottom",
        legend.text=element_text(size=10),
        legend.title=element_text(size=10),
        legend.background = element_rect(colour = "black"),
        panel.margin=unit(3,"lines")) +
  labs(x="probability of tri-allelic dups", y="Phase Concordance") 

# only consider two categories:
sub_dat = subset(dat, category=="bi-allelic" | category=="tri-allelic")
sub_dat$category = factor(sub_dat$category,levels=c("bi-allelic","tri-allelic"))

p_sc <- ggplot(sub_dat,aes(x=hetD_gtype,y=phaseConcordance,colour=category,shape=dataSet)) +
  geom_point(size=5, alpha=0.5) +
  xlim(0,1) + 
  ylim(0,1) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        strip.text=element_text(size=16,face="bold"),
        legend.position="bottom",
        legend.text=element_text(size=10),
        legend.title=element_text(size=10),
        legend.background = element_rect(colour = "black"),
        panel.margin=unit(2,"lines")) +
  labs(x="Het Density", y="Phase Concordance") +
  scale_colour_manual(values=c("blue", "red"),name="Duplication \n Configuration")
ggExtra::ggMarginal(p_sc, type="histogram")

# plot prediction error (with/without weights)
sub_dat$pu1[sub_dat$category=="bi-allelic"] <- 0
sub_dat$pu1[sub_dat$category=="tri-allelic"] <- 1
sub_dat$vari_pc = 1/sub_dat$lengthX * sub_dat$phaseConcordance * (1-sub_dat$phaseConcordance)
sub_dat$vari_hd = 1/sub_dat$num_Mkrs * sub_dat$hetD_gtype * (1-sub_dat$hetD_gtype)

# with weights:
PredictionError = numeric(dim(sub_dat)[1])
predictValue = numeric(dim(sub_dat)[1])
for (i in c(1:dim(sub_dat)[1])){
  Leave_one_dat = sub_dat[-i,]
  newdata = sub_dat[i,]
  fit1 <- glm(pu1~phaseConcordance+hetD_gtype, data=Leave_one_dat, family="binomial",
              weights=1/(vari_pc+vari_hd))
  predictValue[i] = predict(fit1, newdata, type="response")
  PredictionError[i] = newdata$pu1 - predictValue[i]
}
sub_dat$predictErr = PredictionError

a1 = 1*(predictValue>=0.9) + 1*(predictValue<=0.1) -1 + 1*(predictValue>=0.9)
# 1 -> tri-allelic, 0 -> bi-allelic, -1, unable
1- length(which(a1==-1))/length(a1)

print(sum(PredictionError^2))

reorderDat = sub_dat[with(sub_dat, order(sizeOfCNV)),]
p1_wgt <- ggplot(reorderDat, aes(x=c(1:dim(reorderDat)[1]),y=predictErr, colour=category))+
  geom_point(size=3,alpha=0.6) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        legend.position="bottom",
        legend.text=element_text(size=10),
        legend.title=element_text(size=10),
        legend.background = element_rect(colour = "black"),
        panel.margin=unit(2,"lines")) +
  labs(x="Index of ordered event size", y="prediction error") +
  scale_colour_manual(values=c("blue", "red"),name="Duplication \n Configuration")
ggExtra::ggMarginal(p1_wgt, margins="y", type="density")

# without weights:
PredictionError = numeric(dim(sub_dat)[1])
predictValue = numeric(dim(sub_dat)[1])
for (i in c(1:dim(sub_dat)[1])){
  Leave_one_dat = sub_dat[-i,]
  newdata = sub_dat[i,]
  fit1 <- glm(pu1~phaseConcordance+hetD_gtype, data=Leave_one_dat, family="binomial")
  predictValue[i] = predict(fit1, newdata, type="response")
  PredictionError[i] = newdata$pu1 - predictValue[i]
}
sub_dat$predictErr = PredictionError

a2 = 1*(predictValue>=0.9) + 1*(predictValue<=0.1) -1 + 1*(predictValue>=0.9)
# 1 -> tri-allelic, 0 -> bi-allelic, -1, unable
1- length(which(a2==-1))/length(a2)

print(sum(PredictionError^2))

reorderDat = sub_dat[with(sub_dat, order(sizeOfCNV)),]
p1 <- ggplot(reorderDat, aes(x=c(1:dim(reorderDat)[1]),y=predictErr, colour=category,shape=dataSet))+
  geom_point(size=3, alpha=0.6) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        legend.position="bottom",
        legend.text=element_text(size=10),
        legend.title=element_text(size=10),
        legend.background = element_rect(colour = "black"),
        panel.margin=unit(2,"lines")) +
  labs(x="Index of ordered event size",y="prediction error without weights") +
  scale_colour_manual(values=c("blue", "red"),name="Duplication \n Configuration")
ggExtra::ggMarginal(p1, margins="y", type="density")

####################################################
### Figure 1: histogram of the two datasets results
rm(list = ls())
setwd("~/Dropbox/MD_Anderson/Projects/PaulScheet/CNV/")
filename_Cran = "Craniofacial_data/reRun_2016/results/output_Craniofacial-data.results-new.txt"
data_Cran = read.table(filename_Cran,header=TRUE,stringsAsFactors=FALSE)
filename_AGP = "Autism_data/reRun_2016/cleanUpResults/output_AGP_combinedResults.txt"
data_AGP = read.table(filename_AGP,header=TRUE,stringsAsFactors=FALSE)
x1 = data_Cran$pu1_all
x2 = data_AGP$pu1_all

## Results histogram
category1 = rep('a',length(x1))
category1[x1<=0.1] = 'bi-allelic'
category1[x1>=0.9] = 'tri-allelic'
category1[intersect(which(x1>.1),which(x1<.9))] = 'undefined'
data_Cran$category = factor(category1,c("bi-allelic","undefined","tri-allelic"))

category2 = rep('a',length(x2))
category2[x2<=0.1] = 'bi-allelic'
category2[x2>=0.9] = 'tri-allelic'
category2[intersect(which(x2>.1),which(x2<.9))] = 'undefined'
data_AGP$category = factor(category2,c("bi-allelic","undefined","tri-allelic"))

DataSet = c(rep('Craniofacial',length(x1)),rep('Autism',length(x2)))
dat = rbind(data_Cran,data_AGP); dat$dataSet = factor(DataSet, c('Craniofacial','Autism'))

library(ggplot2)
ggplot(dat,aes(pu1_all,fill=category)) + 
  geom_histogram(breaks=seq(0,1,by=0.1)) +
  facet_wrap(~dataSet) +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        strip.text=element_text(size=16,face="bold"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.background = element_rect(colour = "black"),
        panel.margin=unit(2,"lines")) +
  labs(x="the Probability of \"tri-allelic\" duplications") +
  scale_fill_manual(values=c("blue", "darkgray", "red"),name="Duplication \n Configuration")

setEPS()
postscript("histResults_twoDatasets.eps",horizontal=FALSE, paper="special",height=10,width=20)
par(mfrow=c(1,2),mar=c(5,5,5,5))
#hist(x1,breaks=10,labels = TRUE,xlim=c(0,1),ylim=c(0,30),cex.lab=2,cex.main=3,cex.axis=2,col="black",
#     xlab=expression(p[u1]),ylab="counts",main="Craniofacial Data")
bp1 = hist(x1,breaks=10,labels=FALSE,xlim=c(0,1),ylim=c(0,30),cex.lab=2,cex.main=3,cex.axis=2,
           col="black",xlab=expression(p[u1]),ylab="counts",main="Craniofacial Data")
text(x=bp1$mids, y=bp1$counts, labels=bp1$counts, cex=2, pos=3)
#hist(x2,breaks=10,labels = TRUE,xlim=c(0,1),ylim=c(0,30),cex.lab=3,cex.main=3,cex.axis=2,col="black",
#     xlab=expression(p[u1]),ylab="counts",main="Autism Data")
bp2 = hist(x2,breaks=10,labels=FALSE,xlim=c(0,1),ylim=c(0,30),cex.lab=3,cex.main=3,cex.axis=2,
           col="black",xlab=expression(p[u1]),ylab="",main="Autism Data")
text(x=bp2$mids, y=bp2$counts, labels=bp2$counts, cex=2, pos=3)
dev.off()

## chi.square test for the counts in three groups
category1 = rep('a',length(x1))
category1[x1<=0.1] = 'bi-allelic'
category1[x1>=0.9] = 'tri-allelic'
category1[intersect(which(x1>.1),which(x1<.9))] = 'undefined'
data_Cran$category = factor(category1,c("bi-allelic","undefined","tri-allelic"))

category2 = rep('a',length(x2))
category2[x2<=0.1] = 'bi-allelic'
category2[x2>=0.9] = 'tri-allelic'
category2[intersect(which(x2>.1),which(x2<.9))] = 'undefined'
data_AGP$category = factor(category2,c("bi-allelic","undefined","tri-allelic"))

dat = data.frame(dataset=c(rep('Craniofacial',length(x1)),rep('Autism',length(x2))),category=c(category1,category2))
tbl = table(dat)
tbl = tbl[,-3]
(Xsq = chisq.test(tbl))
(fisher.test(tbl))
Xsq$observed   # observed counts (same as M)
Xsq$expected   # expected counts under the null
Xsq$residuals  # Pearson residuals
Xsq$stdres     # standardized residuals

## size in groups
library(ggplot2)
data_Cran$length=data_Cran$end_pos - data_Cran$start_pos+1
data_AGP$length=data_AGP$end_pos - data_AGP$start_pos+1
dat = rbind(data_Cran,data_AGP)
dat$dataSet = factor(c(rep("Craniofacial",58),rep("Autism",28)),c("Craniofacial","Autism"))

ggplot(dat,aes(category, length)) +
  geom_boxplot(aes(color=category),outlier.colour = "red") +
  coord_trans(y = "log10") +
  facet_wrap(~dataSet) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"),
        axis.title.x=element_blank(),
        strip.text=element_text(size=16,face="bold"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.background = element_rect(colour = "black"),
        plot.margin = unit(c(1,1,1,2), "line")) +
  scale_colour_manual(values=c("blue", "darkgray", "red"),name="Duplication \n Configuration")

x1 = data_Cran$length[which(data_Cran$category=="bi-allelic")];
x2 = data_Cran$length[which(data_Cran$category=="tri-allelic")];
x3 = data_AGP$length[which(data_AGP$category=="bi-allelic")];
x4 = data_AGP$length[which(data_AGP$category=="tri-allelic")];
wilcox.test(x1,x2)
wilcox.test(x2,x4)

## remove recombination ones:
new_d_Cran = data_Cran[-c(27,28,47),]
new_d_AGP = data_AGP[-c(15),]
x1 = new_d_Cran$length[which(new_d_Cran$category=="bi-allelic")];
x2 = new_d_Cran$length[which(new_d_Cran$category=="tri-allelic")];
x3 = new_d_AGP$length[which(new_d_AGP$category=="bi-allelic")];
x4 = new_d_AGP$length[which(new_d_AGP$category=="tri-allelic")];
wilcox.test(x1,x2,alternative="less")
wilcox.test(x2,x4,alternative="great")

dat = rbind(new_d_Cran,new_d_AGP)
dat$dataSet = factor(c(rep("Craniofacial",55),rep("Autism",27)),c("Craniofacial","Autism"))
ggplot(dat,aes(category, length)) +
  geom_boxplot(aes(color=category),outlier.colour = "red") +
  coord_trans(y = "log10") +
  facet_wrap(~dataSet) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"),
        axis.title.x=element_blank(),
        strip.text=element_text(size=16,face="bold"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.background = element_rect(colour = "black"),
        plot.margin = unit(c(1,1,1,2), "line")) +
  scale_colour_manual(values=c("blue", "darkgray", "red"),name="Duplication \n Configuration")

## baf.dev and lrr.dev
library(ggplot2)
require(gridExtra)

g1 <- ggplot(data_Cran,aes(x=baf.dev,y=lrr.dev,colour=category)) +
  geom_jitter(size=5,alpha=I(1/2)) +
  labs(title="Craniofacial Data - update list") +
  xlim(0.05,0.25) +
  ylim(0.1,0.45) +
  theme(legend.position="none")

g2 <- ggplot(data_AGP,aes(x=baf.dev,y=lrr.dev,colour=category)) +
  geom_jitter(size=5,alpha=I(1/2)) +
  labs(title="Autism Data - update list") +
  xlim(0.05,0.25) +
  ylim(0.1,0.45)

g <- grid.arrange(g1,g2,ncol=2,widths=c(1,1.3))
ggsave(g,file="devComp.pdf",width=10, height=5)


###########################
g1 <- ggplot(data_Cran,aes(x=baf.dev,y=lrr.dev)) +
  geom_jitter(size=5,alpha=I(1/2),colour="#009E73") +
  labs(title="Craniofacial Data") +
  xlim(0,0.25) +
  ylim(0.1,0.45) +
  theme(legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"),
        plot.title=element_text(size=20,face="bold"),
        plot.margin = unit(c(1,1,1,1), "line")) 

g2 <- ggplot(data_AGP,aes(x=baf.dev,y=lrr.dev)) +
  geom_jitter(size=5,alpha=I(1/2),colour="#009E73") +
  labs(title="Autism Data") +
  xlim(0,0.25) +
  ylim(0.1,0.45) +
  theme(legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"),
        plot.title=element_text(size=20,face="bold"),
        plot.margin = unit(c(1,1,1,1), "line"))
g <- grid.arrange(g1,g2,ncol=2,widths=c(1,1))
###########################
x1 = dat$lrr.dev
subData1 = subset(dat,dupCategories=="bi-allelic",select=c(baf.dev,lrr.dev))
subData3 = subset(dat,dupCategories=="undefinable",select=c(baf.dev,lrr.dev))
subData2 = subset(dat,dupCategories=="tri-allelic",select=c(baf.dev,lrr.dev))
mu1 = mean(subData1$lrr.dev)
mu2 = mean(subData2$lrr.dev)
mu3 = mean(subData3$lrr.dev)
x1[id1]=abs(subData1$lrr.dev-mu1)
x1[id2]=abs(subData2$lrr.dev-mu2)
x1[id3]=abs(subData3$lrr.dev-mu3)
dat$lrr.dev.dist = x1
dat$cat[id1]=rep(1,length(id1))
dat$cat[id2]=rep(2,length(id2))
dat$cat[id3]=rep(0,length(id3))

l = lm(cat~lrr.dev.dist,data=dat)

## dev.results plot
ggplot(dat,aes(x=baf.dev,y=lrr.dev.new,colour=dupCategories)) +
  geom_jitter(size=5,alpha=I(1/2)) +
  labs(title="Craniofacial Data - new list") #+
#  xlim(0,0.3) +
#  ylim(0,1)

#ggplot(dat,aes(x=log(size_kb),y=lrr.dev,colour=dupCategories)) +
#  geom_jitter(size=5,alpha=I(1/2)) 

# check het markers in duplication region from kids
rm(list=ls())
dat=read.table("hetMarkersList.txt",header=FALSE,stringsAsFactors=FALSE)
str(dat)
colnames(dat)=c("SNP.Name","Chr","MapInfo","BAF","LRR","GT1","GT2")
i1 = which(dat$GT1=="-")
i2 = which(dat$GT2=="-")
all(i1==i2)
x=dat$GT1
x[i1] = rep("missingGT",length(i1))
x[-i1] = rep("havingGT",length(x)-length(i1))
dat$Type = as.factor(x)

library(ggplot2)
ggplot(dat,aes(BAF,colour=Type))+
  geom_histogram(binwidth=0.01)+
  facet_wrap(~Type,ncol=1,scales="free_y") +
  labs(title="Craniofacial Data")+
  xlim(0,1)

# rm(list = ls())
# setwd("~/Dropbox/MD_Anderson/Projects/PaulScheet/CNV/Autism_data/reRun/")
# filename = "combinedResults.txt"
# data = read.table(filename,header=TRUE)
# all(data$kid_name==data$kid_name.1)
# plot(data$pu1_all,data$hetD_gtype)
# plot(data$pu1_all,data$hetD_BAF)
# plot(data$pu1_all,data$Expected_hetD)

rm(list = ls())
library(reshape)
library(ggplot2)
setwd("~/Dropbox/MD_Anderson/Projects/PaulScheet/CNV/Craniofacial_data/reRun_2016Jan/")
dat = read.table('overlapRate.txt',header=TRUE,stringsAsFactors=FALSE)
dat$reName = as.character(c(1:dim(dat)[1]))
for (i in 1:length(dat$reName)){
  dat$reName[i] = paste0(dat$reName[i],'-',dat$kidName[i],'-',dat$chr_num[i])
}
newData = dat[c('reName','prop_dist_1','prop_dist_2')]
colnames(newData) = c('reName','overlapRate_In_newList','overlapRate_In_oldList')
ggInput = melt(newData,id=c('reName'))
colnames(ggInput) = c('reName','group','rate')
g = ggplot(ggInput,aes(x=reName,y=rate,fill=group)) +
  geom_bar(stat='identity') +
  facet_grid(group~.) 
g

## lm model on size and pu1
filename = "~/Dropbox/MD_Anderson/Projects/PaulScheet/CNV/Craniofacial_data/reRun_2016Jan/output_Craniofacial-data.results-02092016.txt"
dat = read.table(filename,header=TRUE,stringsAsFactors = FALSE)
filename = "~/Dropbox/MD_Anderson/Projects/PaulScheet/CNV/Craniofacial_data/reRun_2016Jan/combinedResults.txt"
devDataCra = read.table(filename,header=TRUE,stringsAsFactors = FALSE)

## autism data, test dev
filename = "~/Dropbox/MD_Anderson/Projects/PaulScheet/CNV/Autism_data/reRun/dev.results.txt"
devDataAutism = read.table(filename,header=TRUE,stringsAsFactors = FALSE)

## Autism Data, check missing GT data
setwd("~/Dropbox/MD_Anderson/Projects/PaulScheet/CNV/Autism_data/reRun/")
filename = "missingGT-chr22.txt"
dat = read.table(filename,header=FALSE,stringsAsFactors=FALSE)
hist(as.numeric(dat$V5))

filename = "sampleBAF-check-chr22.txt"
dat = read.table(filename,header=FALSE,stringsAsFactors=FALSE)

####################################
## plot for chromosomes shape overlap segments
rm(list = ls())
setwd("~/Dropbox/MD_Anderson/Projects/PaulScheet/CNV/")

library(Gviz)
library(GenomicRanges)

chrId='22'
data = read.table("Craniofacial_data/reRun_2016/update-V1/input_list.txt",header=TRUE,stringsAsFactors=FALSE)
subsetData = subset(data,chr_num==chrId,select=c(kid_name,chr_num,start_pos,end_pos))
newdata <- subsetData[with(subsetData, order(start_pos,end_pos)),] 
colnames(newdata) = c('kid_name','chromosome','start','end')
newdata$width = newdata$end-newdata$start+1
View(newdata)

gen = "hg19"
chr = paste0('chr',chrId)
itrack <- IdeogramTrack(genome=gen,chromosome=chr)
#gtrack <- GenomeAxisTrack()
LL <- list(itrack)

resultsData = read.table('Craniofacial_data/reRun_2016/update-V1/results/output_Craniofacial-original.txt',header=TRUE,stringsAsFactors=FALSE)
l = dim(newdata)[1]
for (i in c(1:l)){
  kidName = newdata$kid_name[i]
  idx1 = which(resultsData$kid_name==kidName)
  if (length(idx1)>1){
    #print("larger than 1")
    idx2 = which(resultsData$start_pos==newdata$start[i])
    idx3 = which(resultsData$end_pos==newdata$end[i])
    idx = intersect(idx1,intersect(idx2,idx3))
  }else{
    idx = idx1
  }
  print(idx)
  pu1All = resultsData$pu1_all[idx]
  if (pu1All>=0.9){
    c = 'red'
  }else if (pu1All<=0.1){
    c = 'blue'
  }else{
    c = 'gray'
  }
  atrack <- AnnotationTrack(newdata[i,],genome=gen,chromosome=chr,
                            name=paste0("sample",as.character(i)),fill=c)
  LL <- c(LL,assign(paste0('atrack',as.character(i)), atrack))
}
class(LL)

chr22.end = 51304566
chr21.end = 48129895
chr20.end = 63025520
chr19.end = 59128983
chr18.end = 78077248
chr17.end = 81195210
chr16.end = 90354753
chr15.end = 102531392
chr14.end = 107349540
chr13.end = 115169878
chr12.end = 133851895
chr11.end = 135006516
chr10.end = 135534747
chr9.end = 141213431
chr8.end = 146364022
chr7.end = 159138663
chr6.end = 171115067
chr5.end =	180915260
chr4.end = 191154276
chr3.end = 198022430
chr2.end = 243199373
chr1.end = 249250621

png("Craniofacial_chr22_location.png",height=2,width=4,units='in',res=300)
plotTracks(LL,from=0,to=chr11.end,showId=FALSE,showBandId=TRUE,col=NULL)
dev.off()

###Get centromere information:
##centroData = read.table('Craniofacial_data/reRun_2016/general_inputfile/cytoBand.txt',
##                        header=FALSE,stringsAsFactors=FALSE)
centroData = data.frame(chr1.ctr=125000000,chr2.ctr=93300000,chr3.ctr=91000000,
                        chr4.ctr=50400000,chr5.ctr=48400000,chr6.ctr=61000000,
                        chr7.ctr=59900000,chr8.ctr=45600000,chr9.ctr=49000000,
                        chr10.ctr=40200000,chr11.ctr=53700000,chr12.ctr=35800000,
                        chr13.ctr=17900000,chr14.ctr=17600000,chr15.ctr=19000000,
                        chr16.ctr=36600000,chr17.ctr=24000000,chr18.ctr=17200000,
                        chr19.ctr=26500000,chr20.ctr=27500000,chr21.ctr=13200000,
                        chr22.ctr=14700000)
## centromere distance:
subData = resultsData[union(which(resultsData$pu1_all<=0.1),which(resultsData$pu1_all>=0.9)),]
subData$dupType = 1*(subData$pu1_all>=0.9) + 0*(subData$pu1_all<=0.1)
for (i in c(1:dim(subData)[1])){
  chrId = subData$chr_num[i]
  centroLoc = centroData[1,as.integer(chrId)]
  sDist = subData$start_pos[i]-centroLoc
  eDist = subData$end_pos[i]-centroLoc
  if (sDist<0 & eDist<0){
    #p-arm
    d = centroLoc-subData$end_pos[i]
    subData$d[i] = d
  }else if (sDist>0 & eDist>0){
    #q-arm
    d = sDist
    subData$d[i] = d
  }else{
    print(i)
  }
}
####Autism Data
## plot for chromosomes shape overlap segments
rm(list = ls())
setwd("~/Dropbox/MD_Anderson/Projects/PaulScheet/CNV/Autism_data/reRun_2016/")

library(Gviz)
library(GenomicRanges)

data = read.table("general_inputfile/input_list_updated.txt",header=TRUE,stringsAsFactors=FALSE)
chrId='8'
subsetData = subset(data,chr_num==chrId,select=c(kid_name,chr_num,start_pos,end_pos))
newdata <- subsetData[with(subsetData, order(start_pos,end_pos)),] 
colnames(newdata) = c('kid_name','chromosome','start','end')
newdata$width = newdata$end-newdata$start+1
View(newdata)

gen = "hg19"
chr = paste0('chr',chrId)
itrack <- IdeogramTrack(genome=gen,chromosome=chr)
#gtrack <- GenomeAxisTrack()
atrack1 <- AnnotationTrack(newdata[1,],genome=gen,chromosome=chr,name="sample1",fill='blue')
atrack2 <- AnnotationTrack(newdata[2,],genome=gen,chromosome=chr,name="sample2",fill='blue')
atrack3 <- AnnotationTrack(newdata[3,],genome=gen,chromosome=chr,name="sample3",fill='blue')
atrack4 <- AnnotationTrack(newdata[4,],genome=gen,chromosome=chr,name="sample4",fill='blue')
atrack5 <- AnnotationTrack(newdata[5,],genome=gen,chromosome=chr,name="sample5",fill='red')
atrack6 <- AnnotationTrack(newdata[6,],genome=gen,chromosome=chr,name="sample6",fill='blue')
atrack7 <- AnnotationTrack(newdata[9,],genome=gen,chromosome=chr,name="sample7",fill='gray')

png("Autism_chr2_location.png",height=1,width=4,units='in',res=300)
plotTracks(list(itrack,atrack1,atrack2,atrack3,atrack4,atrack5,atrack6,atrack7),from=0,to=190214555,
           showId=FALSE,showBandId=TRUE,col=NULL)
dev.off()

########################## model comparison:
data1 = read.table('Craniofacial_data/reRun_2016/update-V1/results/output_Craniofacial-original.txt',header=TRUE,stringsAsFactors=FALSE)
data2 = read.table('Craniofacial_data/reRun_2016/update-V1/results/output_Craniofacial-Logsmooth-v2.txt',header=TRUE,stringsAsFactors=FALSE)
dat = data.frame(x=c(1:dim(data1)[1]),original=data1$pu1_all,logSmooth=data2$pu1_all,size=data1$sizeOfCNV)
ggplot(dat, aes(x)) +  
  geom_point(aes(y=original), colour="red") +
  geom_point(aes(y=logSmooth),colour="blue")
