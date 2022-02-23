options(stringsAsFactors = FALSE)


survivalSigBreastLung<-read.table('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/NewRoute/ProbeCalculate/survival/BreastLungSigThreeDatabase/survivalSigBreastLung.txt',
                                  row.names = 1,header = TRUE,sep = '\t')
load('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/NewRoute/ResourceData/clinicalSurvivalUse.RData')
load('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/NewRoute/ResourceData/clinicalCaseControl.RData')
load('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/NewRoute/ResourceData/expTMMcasecontrol.RData')
load('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/NewRoute/ResourceData/expTMMsurvAllRaw.RData')
load('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/NewRoute/ResourceData/expTMMsurvEach.RData')

TMSprobeWeightedSum<-function(TMSgenetable,expdataset){
  
  ## up regulated gene weighted +1, down regulated gene weighted -1
  upprobe<-rownames(TMSgenetable[TMSgenetable$dysregulate %in% c('up'),])
  downprobe<-rownames(TMSgenetable[TMSgenetable$dysregulate %in% c('down'),])
  upprobeexp<-lapply(expdataset, function(x){x[upprobe,]})
  downprobeexp<-lapply(expdataset, function(x){x[downprobe,]})
  upTMScoreSum<-lapply(upprobeexp, function(x){apply(x, 2, sum)})
  downTMScoreSum<-lapply(downprobeexp, function(x){apply(x, 2, sum)})
  TMScoreWeightSum<-lapply(1:length(upTMScoreSum), function(x){
    upTMScoreSum[[x]]-downTMScoreSum[[x]]
  }) 
  names(TMScoreWeightSum)<-names(upTMScoreSum)
  return(TMScoreWeightSum)
}
TMScoreCaseControl<-TMSprobeWeightedSum(survivalSigBreastLung,expTMMcasecontrol)
TMScoreSurvAllRaw<-TMSprobeWeightedSum(survivalSigBreastLung,expTMMsurvAllRaw)
TMScoreSurvEach<-TMSprobeWeightedSum(survivalSigBreastLung,expTMMsurvEach)

TMScoreCaseControl1<-unlist(TMScoreCaseControl)
names(TMScoreCaseControl1)<-do.call(rbind,strsplit(names(TMScoreCaseControl1),'.',fixed = TRUE))[,2]

TMScoreSurvAllRaw1<-unlist(TMScoreSurvAllRaw)
names(TMScoreSurvAllRaw1)<-do.call(rbind,strsplit(names(TMScoreSurvAllRaw1),'.',fixed = TRUE))[,2]

TMScoreSurvEach1<-unlist(TMScoreSurvEach)
names(TMScoreSurvEach1)<-do.call(rbind,strsplit(names(TMScoreSurvEach1),'.',fixed = TRUE))[,2]

clinicalCaseControl<-cbind(clinicalCaseControl,TMScoreCaseControl=TMScoreCaseControl1[rownames(clinicalCaseControl)])

clinicalSurvivalAllRaw<-clinicalSurvivalUse[c(colnames(expTMMsurvAllRaw$lung),colnames(expTMMsurvAllRaw$breast),
                                              colnames(expTMMsurvAllRaw$colon),colnames(expTMMsurvAllRaw$kidney),
                                              colnames(expTMMsurvAllRaw$stomach)),]
clinicalSurvivalAllRaw<-cbind(clinicalSurvivalAllRaw,TMScoreSurvAllRaw=TMScoreSurvAllRaw1[rownames(clinicalSurvivalAllRaw)])

clinicalSurvivalUse<-cbind(clinicalSurvivalUse,TMScoreSurvEach=TMScoreSurvEach1[rownames(clinicalSurvivalUse)])

setwd('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/NewRoute/ProbeCalculate/survival/BreastLungSigThreeDatabase')
library(openxlsx)
write.xlsx(clinicalCaseControl,file = 'clinicalCaseControl.xlsx',colNames = TRUE, rowNames = TRUE)
write.xlsx(clinicalSurvivalAllRaw,file = 'clinicalSurvivalAllRaw.xlsx',colNames = TRUE, rowNames = TRUE)
write.xlsx(clinicalSurvivalUse,file = 'clinicalSurvivalUse.xlsx',colNames = TRUE, rowNames = TRUE)

############################################################
##########################################

library(openxlsx)
TMScoreCaseCtrl<-read.xlsx('E:/0/telomeres/Results/20210415NewRoute/ProbeCalculate/1survivalSigBreastLung_use/TMScore/clinicalCaseControl.xlsx',
                           sheet = 1,colNames = TRUE,rowNames = TRUE)
TMScoreSurvAllRaw<-read.xlsx('E:/0/telomeres/Results/20210415NewRoute/ProbeCalculate/1survivalSigBreastLung_use/TMScore/clinicalSurvivalAllRaw.xlsx',
                             sheet = 1,colNames = TRUE,rowNames = TRUE)
TMScoreSurvEach<-read.xlsx('E:/0/telomeres/Results/20210415NewRoute/ProbeCalculate/1survivalSigBreastLung_use/TMScore/clinicalSurvivalUse.xlsx',
                           sheet = 1,colNames = TRUE,rowNames = TRUE)

setwd('E:/0/telomeres/Results/20210415NewRoute/ProbeCalculate/1survivalSigBreastLung_use/TMScore/boxplot')


library(ggpubr)
p1<-ggviolin(data = subset(TMScoreCaseCtrl, tissue_type %in% c('tumor', 'normal')),  nrow=1,
             x = 'tissue_type', y = 'TMScoreCaseControl', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'tissue_type', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  stat_compare_means(comparisons = list(c('tumor','normal')),method='wilcox.test',
                     aes(group=tissue_type),label= 'p.signif',label.y = c(95))

p2<-ggviolin(data = subset(TMScoreSurvAllRaw, tissue_type %in% c('tumor', 'normal')), nrow=1, 
             x = 'tissue_type', y = 'TMScoreSurvAllRaw', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'tissue_type', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  stat_compare_means(comparisons = list(c('tumor','normal')),method='wilcox.test',
                     aes(group=tissue_type),label= 'p.signif',label.y = c(95))

ggarrange(p1,p2,nrow = 2,legend = 'top',common.legend = TRUE,font.label = list(size=6),heights = c(1,1))

ggviolin(data = subset(TMScoreCaseCtrl, stage_T %in% c('T1', 'T2', 'T3', 'T4')), 
         x = 'stage_T', y = 'TMScoreCaseControl', facet.by = 'cancer_type',
         palette = 'nc',
         fill = 'stage_T', add = c('boxplot',"jitter"),add.params = list(fill="white")) + 
  stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('T1','T2'),c('T2','T3'),c('T3','T4'),c('T1','T3'),c('T2','T4'),c('T1','T4')),
                     method='wilcox.test',label.y = c(90,90,90,95,97,100)) 

### stage
TMScoreCaseCtrl$stage<-factor(TMScoreCaseCtrl$stage,levels = c('I', 'II', 'III', 'IV'))
# TMScoreCaseCtrl$cancer_type<-factor(TMScoreCaseCtrl$cancer_type,levels = c())
p1<-ggviolin(data = subset(TMScoreCaseCtrl, stage %in% c('I', 'II', 'III', 'IV')), nrow=1,
             x = 'stage', y = 'TMScoreCaseControl', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'stage', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  stat_compare_means(label.y = 110)+
  stat_compare_means(comparisons = list(c('I','II'),c('II','III'),c('III','IV'),c('I','III'),c('II','IV'),c('I','IV')),
                     method='wilcox.test',label= 'p.signif',label.y = c(95,95,95,100,100,105)) # sig

TMScoreSurvAllRaw$stage[TMScoreSurvAllRaw$stage %in% c('1','1a','1b','IA','IB','I','Stage 1')]<-'I'
TMScoreSurvAllRaw$stage[TMScoreSurvAllRaw$stage %in% c('2','2b','2a','II','IIB','IIA','Stage 2')]<-'II'
TMScoreSurvAllRaw$stage[TMScoreSurvAllRaw$stage %in% c('3','3a','3b','IIIA','IIIB','IIIC','Stage 3','III')]<-'III'
TMScoreSurvAllRaw$stage[TMScoreSurvAllRaw$stage %in% c('4','IV','Stage 4')]<-'IV'
TMScoreSurvAllRaw$stage[TMScoreSurvAllRaw$stage %in% c('N/A','0','Unknown','Recurrence','Unkown')]<-NA

TMScoreSurvAllRaw$stage<-factor(TMScoreSurvAllRaw$stage,levels = c('I', 'II', 'III', 'IV'))
p2<-ggviolin(data = subset(TMScoreSurvAllRaw, stage %in% c('I', 'II', 'III', 'IV')), nrow=1,
             x = 'stage', y = 'TMScoreSurvAllRaw', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'stage', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  stat_compare_means(label.y = 110)+
  stat_compare_means(comparisons = list(c('I','II'),c('II','III'),c('III','IV'),c('I','III'),c('II','IV'),c('I','IV')),
                     method='wilcox.test',label= 'p.signif',label.y = c(95,95,95,100,100,105)) # sig 
ggarrange(p1,p2,nrow = 2,legend = 'top',common.legend = TRUE,font.label = list(size=6),heights = c(1,1))


TMScoreCaseCtrl$grade[TMScoreCaseCtrl$grade %in% c('1','I','G1')]<-'G1'
TMScoreCaseCtrl$grade[TMScoreCaseCtrl$grade %in% c('2','II','G2')]<-'G2'
TMScoreCaseCtrl$grade[TMScoreCaseCtrl$grade %in% c('3','III','G3')]<-'G3'
TMScoreCaseCtrl$grade[TMScoreCaseCtrl$grade %in% c('4')]<-'G4'
TMScoreCaseCtrl$grade[TMScoreCaseCtrl$grade %in% c('N/A')]<-NA

TMScoreCaseCtrl$grade<-factor(TMScoreCaseCtrl$grade,levels = c('G1','G2','G3','G4'))
p1<-ggviolin(data = subset(TMScoreCaseCtrl[TMScoreCaseCtrl$cancer_type=='breast',], grade %in% c('G1','G2', 'G3')), nrow=1, 
             x = 'grade', y = 'TMScoreCaseControl', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'grade', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('G1','G2'),c('G2','G3'),c('G1','G3')),
                     method='wilcox.test',label= 'p.signif',label.y = c(95,95,100)) # sig


TMScoreSurvAllRaw$grade[TMScoreSurvAllRaw$grade %in% c('1','G1','WELL DIFF','I')]<-'G1'
TMScoreSurvAllRaw$grade[TMScoreSurvAllRaw$grade %in% c('2','G2','MOD DIFF','II')]<-'G2'
TMScoreSurvAllRaw$grade[TMScoreSurvAllRaw$grade %in% c('3','G3','POORLY DIFF','III')]<-'G3'
TMScoreSurvAllRaw$grade[TMScoreSurvAllRaw$grade %in% c('4','IV')]<-'G4'
TMScoreSurvAllRaw$grade[TMScoreSurvAllRaw$grade %in% c('Unknown')]<-NA

TMScoreSurvAllRaw$grade<-factor(TMScoreSurvAllRaw$grade,levels = c('G1','G2','G3','G4'))
p2<-ggviolin(data = subset(TMScoreSurvAllRaw[TMScoreSurvAllRaw$cancer_type=='breast',], grade %in% c('G1','G2', 'G3')),  nrow=1, 
             x = 'grade', y = 'TMScoreSurvAllRaw', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'grade', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('G1','G2'),c('G2','G3'),c('G1','G3')),
                     method='wilcox.test',label= 'p.signif',label.y = c(95,95,100)) # sig 
p3<-ggarrange(p1,p2,nrow = 2,legend = 'top',common.legend = TRUE,font.label = list(size=6),heights = c(1,1))

### breast ER status
TMScoreCaseCtrl$ER_IHC_status[TMScoreCaseCtrl$ER_IHC_status %in% c('1','pos','POS','+')]<- 'pos'
TMScoreCaseCtrl$ER_IHC_status[TMScoreCaseCtrl$ER_IHC_status %in% c('0','neg','NEG','-')]<- 'neg'
TMScoreCaseCtrl$ER_IHC_status[TMScoreCaseCtrl$ER_IHC_status %in% c('N/A','EV')]<- NA

TMScoreCaseCtrl$ER_IHC_status<-factor(TMScoreCaseCtrl$ER_IHC_status,levels = c('pos','neg'))
p1<-ggviolin(data = subset(TMScoreCaseCtrl, ER_IHC_status %in% c('neg','pos')), 
             x = 'ER_IHC_status', y = 'TMScoreCaseControl', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'ER_IHC_status', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  # stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('neg','pos')),
                     method='wilcox.test',label= 'p.signif',label.y = c(100)) # sig 

TMScoreSurvAllRaw$ER_status[TMScoreSurvAllRaw$ER_status %in% c('1','positive','pos','Positive')]<- 'pos'
TMScoreSurvAllRaw$ER_status[TMScoreSurvAllRaw$ER_status %in% c('0','negative','neg','Negative')]<- 'neg'
TMScoreSurvAllRaw$ER_status[TMScoreSurvAllRaw$ER_status %in% c('EV')]<- NA

TMScoreSurvAllRaw$ER_status<-factor(TMScoreSurvAllRaw$ER_status,levels = c('pos','neg'))
p2<-ggviolin(data = subset(TMScoreSurvAllRaw, ER_status %in% c('neg','pos')), 
             x = 'ER_status', y = 'TMScoreSurvAllRaw', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'ER_status', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  # stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('neg','pos')),
                     method='wilcox.test',label= 'p.signif',label.y = c(100)) # sig
p4<-ggarrange(p1,p2,nrow = 2,legend = 'top',common.legend = TRUE,font.label = list(size=6),heights = c(1,1))
### breast PR status
TMScoreCaseCtrl$PR_IHC_status[TMScoreCaseCtrl$PR_IHC_status %in% c('1','pos','+')]<- 'pos'
TMScoreCaseCtrl$PR_IHC_status[TMScoreCaseCtrl$PR_IHC_status %in% c('0','neg','-')]<- 'neg'
TMScoreCaseCtrl$PR_IHC_status[TMScoreCaseCtrl$PR_IHC_status %in% c('N/A','EV')]<- NA

TMScoreCaseCtrl$PR_IHC_status<-factor(TMScoreCaseCtrl$PR_IHC_status,levels = c('pos','neg'))
p1<-ggviolin(data = subset(TMScoreCaseCtrl, PR_IHC_status %in% c('neg','pos')), 
             x = 'PR_IHC_status', y = 'TMScoreCaseControl', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'PR_IHC_status', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  # stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('neg','pos')),
                     method='wilcox.test',label= 'p.signif',label.y = c(100)) # sig


TMScoreSurvAllRaw$PgR_status[TMScoreSurvAllRaw$PgR_status %in% c('1','positive','pos','Positive')]<- 'pos'
TMScoreSurvAllRaw$PgR_status[TMScoreSurvAllRaw$PgR_status %in% c('0','negative','neg','Negative')]<- 'neg'
TMScoreSurvAllRaw$PgR_status[TMScoreSurvAllRaw$PgR_status %in% c('EV')]<- NA

TMScoreSurvAllRaw$PgR_status<-factor(TMScoreSurvAllRaw$PgR_status,levels = c('pos','neg'))
p2<-ggviolin(data = subset(TMScoreSurvAllRaw, PgR_status %in% c('neg','pos')), 
             x = 'PgR_status', y = 'TMScoreSurvAllRaw', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'PgR_status', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  # stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('neg','pos')),
                     method='wilcox.test',label= 'p.signif',label.y = c(100)) # sig 

p5<-ggarrange(p1,p2,nrow = 2,legend = 'top',common.legend = TRUE,font.label = list(size=6),heights = c(1,1))

p6<-ggarrange(p3,p4,p5,nrow = 1,ncol = 3,legend = 'top',common.legend = TRUE,font.label = list(size=6),widths = c(1,1,1))
### breast molecular subtype
TMScoreCaseCtrl$molecular_subtype[TMScoreCaseCtrl$molecular_subtype %in% c('Basal')]<- 'Basal'
TMScoreCaseCtrl$molecular_subtype[TMScoreCaseCtrl$molecular_subtype %in% c('ERBB2','HER2')]<- 'HER2'
TMScoreCaseCtrl$molecular_subtype[TMScoreCaseCtrl$molecular_subtype %in% c('Luminal A','LumA','LuminalA')]<- 'LumA'
TMScoreCaseCtrl$molecular_subtype[TMScoreCaseCtrl$molecular_subtype %in% c('Luminal B','LumB')]<- 'LumB'
TMScoreCaseCtrl$molecular_subtype[TMScoreCaseCtrl$molecular_subtype %in% c('Normal','normal')]<- 'Normal'
TMScoreCaseCtrl$molecular_subtype[TMScoreCaseCtrl$molecular_subtype %in% c('N/A')]<- NA

TMScoreCaseCtrl$molecular_subtype<-factor(TMScoreCaseCtrl$molecular_subtype,levels = c('LumA','Normal','LumB','HER2','Basal'))
p1<-ggviolin(data = subset(TMScoreCaseCtrl, molecular_subtype %in% c('LumA','Normal','LumB','HER2','Basal')), 
             x = 'molecular_subtype', y = 'TMScoreCaseControl', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'molecular_subtype', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1/2)+
  stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('LumA','Normal'),c('Normal','LumB'),c('LumB','HER2'),c('HER2','Basal'),c('LumA','LumB'),
                                        c('Normal','HER2'),c('LumB','Basal'),c('LumA','HER2'),c('Normal','Basal'),c('LumA','Basal')),
                     method='wilcox.test',label= 'p.signif',label.y = c(90,90,90,90,95,95,95,100,100,105)) # sig 


TMScoreSurvAllRaw$molecular_subtype[TMScoreSurvAllRaw$molecular_subtype %in% c('HER2','HER2-enriched','HER2-positive breast cancer','HER2 breast cancer tumor','ERBB2','Her2')]<-'HER2'
TMScoreSurvAllRaw$molecular_subtype[TMScoreSurvAllRaw$molecular_subtype %in% c('LumB','Luminal B','LuminalB')]<-'LumB'
TMScoreSurvAllRaw$molecular_subtype[TMScoreSurvAllRaw$molecular_subtype %in% c('Basal','Basal-like subtype','triple negative breast cancer (TNBC) patient','triple negative breast cancer','TNBC')]<-'Basal'
TMScoreSurvAllRaw$molecular_subtype[TMScoreSurvAllRaw$molecular_subtype %in% c('LumA','Luminal A','LuminalA')]<-'LumA'
TMScoreSurvAllRaw$molecular_subtype[TMScoreSurvAllRaw$molecular_subtype %in% c('Normal','Normal breast-like')]<-'Normal'
TMScoreSurvAllRaw$molecular_subtype[TMScoreSurvAllRaw$molecular_subtype %in% c('Uncla','--','Healthy')]<-NA
subdata<-subset(TMScoreSurvAllRaw, cancer_type %in% c('breast'))
subdata$molecular_subtype<-factor(subdata$molecular_subtype,levels = c('LumA','Normal','LumB','HER2','Basal'))
p2<-ggviolin(data = subset(subdata, molecular_subtype %in% c('LumA','Normal','LumB','HER2','Basal')), 
             x = 'molecular_subtype', y = 'TMScoreSurvAllRaw', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'molecular_subtype', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1/2)+
  stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('LumA','Normal'),c('Normal','LumB'),c('LumB','HER2'),c('HER2','Basal'),c('LumA','LumB'),
                                        c('Normal','HER2'),c('LumB','Basal'),c('LumA','HER2'),c('Normal','Basal'),c('LumA','Basal')),
                     method='wilcox.test',label= 'p.signif',label.y = c(90,90,90,90,95,95,95,100,100,105)) # sig 可用

p7<-ggarrange(p1,p2,nrow = 1,ncol = 2,legend = 'top',common.legend = TRUE,font.label = list(size=6),widths = c(1,1))
# p8<-ggarrange(p6,p7,nrow = 2,legend = 'top',common.legend = TRUE,font.label = list(size=6),heights =c(2,1))

############################################################



TMScoreCaseCtrl$stage_N_status[TMScoreCaseCtrl$stage_N_status %in% c('pos','1')]<-'pos'
TMScoreCaseCtrl$stage_N_status[TMScoreCaseCtrl$stage_N_status %in% c('neg','0')]<-'neg'
TMScoreCaseCtrl$stage_N_status<-factor(TMScoreCaseCtrl$stage_N_status,levels = c('neg','pos'))
ggviolin(data = subset(TMScoreCaseCtrl, stage_N_status %in% c('neg','pos')), 
         x = 'stage_N_status', y = 'TMScoreCaseControl', facet.by = 'cancer_type',
         palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
         fill = 'stage_N_status', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  # stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('neg','pos')),
                     method='wilcox.test',label= 'p.signif',label.y = c(100)) # ns

TMScoreCaseCtrl$P53_IHC_status[TMScoreCaseCtrl$P53_IHC_status %in% c('pos','+')]<-'pos'
TMScoreCaseCtrl$P53_IHC_status[TMScoreCaseCtrl$P53_IHC_status %in% c('neg','-')]<-'neg'
TMScoreCaseCtrl$P53_IHC_status<-factor(TMScoreCaseCtrl$P53_IHC_status,levels = c('neg','pos'))
ggviolin(data = subset(TMScoreCaseCtrl, P53_IHC_status %in% c('neg','pos')), 
         x = 'P53_IHC_status', y = 'TMScoreCaseControl', facet.by = 'cancer_type',
         palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
         fill = 'P53_IHC_status', add = c("jitter"),add.params = list(size=0.2)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  # stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('neg','pos')),
                     method='wilcox.test',label= 'p.signif',label.y = c(100)) # ns

TMScoreCaseCtrl$ERBB2_IHC_status[TMScoreCaseCtrl$ERBB2_IHC_status %in% c('pos','1')]<-'pos'
TMScoreCaseCtrl$ERBB2_IHC_status[TMScoreCaseCtrl$ERBB2_IHC_status %in% c('neg','0')]<-'neg'
TMScoreCaseCtrl$ERBB2_IHC_status<-factor(TMScoreCaseCtrl$ERBB2_IHC_status,levels = c('neg','pos'))
ggviolin(data = subset(TMScoreCaseCtrl, ERBB2_IHC_status %in% c('neg','pos')), 
         x = 'ERBB2_IHC_status', y = 'TMScoreCaseControl', facet.by = 'cancer_type',
         palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
         fill = 'ERBB2_IHC_status', add = c("jitter"),add.params = list(size=0.2)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  # stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('neg','pos')),
                     method='wilcox.test',label= 'p.signif',label.y = c(100)) # ns


breastCaseCtrl<-TMScoreCaseCtrl[TMScoreCaseCtrl$cancer_type=='breast',]
breastCaseCtrl$subtyping1[breastCaseCtrl$subtyping1 %in% c('IDC')]<-'IDC'
breastCaseCtrl$subtyping1[breastCaseCtrl$subtyping1 %in% c('LOC')]<-'ILC'
breastCaseCtrl$subtyping1[breastCaseCtrl$subtyping1 %in% c('MIX','other')]<-'other'
breastCaseCtrl$subtyping1<-factor(breastCaseCtrl$subtyping1,levels = c('normal','DCIS','ILC','IDC','other'))
ggviolin(data = subset(breastCaseCtrl, subtyping1 %in% c('normal','DCIS','ILC','IDC')), 
         x = 'subtyping1', y = 'TMScoreCaseControl', facet.by = 'cancer_type',
         palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
         fill = 'subtyping1', add = c("jitter"),add.params = list(size=0.2)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  # stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('normal','DCIS'),c('DCIS','ILC'),c('IDC','ILC'),c('normal','ILC'),c('DCIS','IDC'),c('normal','IDC')),
                     method='wilcox.test',label= 'p.signif',label.y = c(95,95,95,100,100,105)) 


breastSurv<-TMScoreSurvAllRaw[TMScoreSurvAllRaw$cancer_type=='breast',]
breastSurv$subtyping[breastSurv$subtyping %in% c(' ductal','DUC','IDC','ductal','invasive ductal carcinoma')]<-'IDC'
breastSurv$subtyping[breastSurv$subtyping %in% c(' lobular','ILC','lobular')]<-'ILC'
breastSurv$subtyping[breastSurv$subtyping %in% c(' squamous cell','other','TUB','MIX','MUC','medullary carcinoma',
                                                 'undifferentiated carcinoma infiltrating','metaplastic carcinoma',
                                                 'Trabecular','Solid','Classic','Alveolar','Mixed, non classic')]<-'other'

breastSurv$subtyping<-factor(breastSurv$subtyping,levels = c('ILC','IDC','other'))
ggviolin(data = subset(breastSurv, subtyping %in% c('IDC','ILC')), 
         x = 'subtyping', y = 'TMScoreSurvAllRaw', facet.by = 'cancer_type',
         palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
         fill = 'subtyping', add = c("jitter"),add.params = list(size=0.2)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  # stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('IDC','ILC')),
                     method='wilcox.test',label= 'p.signif',label.y = c(100)) # ***

TMScoreSurvAllRaw$HER2_status[TMScoreSurvAllRaw$HER2_status %in% c('0','normal','negative','neg','no','Negative')]<-'neg'
TMScoreSurvAllRaw$HER2_status[TMScoreSurvAllRaw$HER2_status %in% c('1','over-expression','HER2-positive breast cancer','HER2 breast cancer tumor','positive','pos','yes','Positive')]<-'pos'
TMScoreSurvAllRaw$HER2_status<-factor(TMScoreSurvAllRaw$HER2_status,levels = c('pos','neg'))
ggviolin(data = subset(TMScoreSurvAllRaw, HER2_status %in% c('pos','neg')), 
         x = 'HER2_status', y = 'TMScoreSurvAllRaw', facet.by = 'cancer_type',
         palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
         fill = 'HER2_status', add = c("jitter"),add.params = list(size=0.2)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  # stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('pos','neg')),
                     method='wilcox.test',label= 'p.signif',label.y = c(100)) # ns
## ****
TMScoreSurvAllRaw$P53_status[TMScoreSurvAllRaw$P53_status %in% c('0')]<-'neg'
TMScoreSurvAllRaw$P53_status[TMScoreSurvAllRaw$P53_status %in% c('1')]<-'pos'
TMScoreSurvAllRaw$P53_status<-factor(TMScoreSurvAllRaw$P53_status,levels = c('pos','neg'))
ggviolin(data = subset(TMScoreSurvAllRaw, P53_status %in% c('pos','neg')), 
         x = 'P53_status', y = 'TMScoreSurvAllRaw', facet.by = 'cancer_type',
         palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
         fill = 'P53_status', add = c("jitter"),add.params = list(size=0.2)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  # stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('pos','neg')),
                     method='wilcox.test',label= 'p.signif',label.y = c(100)) # ****

################ Lung cancer
## 1. subtyping
lungCaseCtrl<-TMScoreCaseCtrl[TMScoreCaseCtrl$cancer_type=='lung',]
lungCaseCtrl$subtyping1<-factor(lungCaseCtrl$subtyping1,levels = c('normal','LADC','LUSC','LCLC','SCLC'))
p1<-ggviolin(data = subset(lungCaseCtrl, subtyping1 %in% c('normal','LADC','LUSC','LCLC','SCLC')), 
             x = 'subtyping1', y = 'TMScoreCaseControl', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'subtyping1', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 0.5)+
  # stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('normal','LADC'),c('LADC','LUSC'),c('LUSC','LCLC'),c('SCLC','LCLC'),
                                        c('normal','LUSC'),c('LADC','LCLC'),c('LUSC','SCLC'),c('normal','LCLC'),
                                        c('LADC','SCLC'),c('normal','SCLC')),
                     method='wilcox.test',label= 'p.signif',label.y = c(105,105,105,105,110,110,110,115,115,120)) # sig 

lungSurv<-subset(TMScoreSurvAllRaw,cancer_type %in% c('lung'))
lungSurv$subtyping1[lungSurv$subtyping1 %in% c('Squamous','SCC','squamous','SQC','Squamous Cell Carcinoma of Lung')]<-'LUSC'
lungSurv$subtyping1[lungSurv$subtyping1 %in% c('Adenocarcinoma','ADC','adeno','LADC')]<-'LADC'
lungSurv$subtyping1[lungSurv$subtyping1 %in% c('LCC','large')]<-'LCLC'
lungSurv$subtyping1[lungSurv$subtyping1 %in% c('healthy','NTL')]<-'normal'
lungSurv$subtyping1[lungSurv$subtyping1 %in% c('LCNE','CARCI','BAS','Other')]<-NA

lungSurv$subtyping1<-factor(lungSurv$subtyping1,levels = c('normal','LADC','LUSC','LCLC'))
p2<-ggviolin(data = subset(lungSurv, subtyping1 %in% c('normal','LADC','LUSC','LCLC')), 
             x = 'subtyping1', y = 'TMScoreSurvAllRaw', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'subtyping1', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 0.5)+
  # stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('normal','LADC'),c('LADC','LUSC'),c('LUSC','LCLC'),
                                        c('normal','LUSC'),c('LADC','LCLC'),c('normal','LCLC')),
                     method='wilcox.test',label= 'p.signif',label.y = c(110,110,110,115,115,120)) # sig 
ggarrange(p1,p2,ncol = 2,widths = c(1,1),legend = 'top',common.legend = TRUE,font.label = list(size=6))
## 2. 
ggviolin(data = subset(lungCaseCtrl, stage_N %in% c("N0", "N1", "N2", "N3")), 
         x = 'stage_N', y = 'TMScoreCaseControl', facet.by = 'cancer_type',
         palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
         fill = 'stage_N', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 0.5)+
  # stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('N0','N1'),c('N1','N2'),c('N2','N3'),c('N0','N2'),
                                        c('N1','N3'),c('N0','N3')),
                     method='wilcox.test',label= 'p.signif',label.y = c(95,95,95,95,100,100,100,105,105,110)) # sig 
ggviolin(data = subset(lungSurv, KRAS_mutation_status %in% c("negative", "positive")), 
         x = 'KRAS_mutation_status', y = 'TMScoreSurvAllRaw', facet.by = 'cancer_type',
         palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
         fill = 'KRAS_mutation_status', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 0.5)+
  # stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('negative','positive')),
                     method='wilcox.test',label= 'p.signif',label.y = c(95,95,95,95,100,100,100,105,105,110)) # sig 

## Colon cancer
colonSurv<-subset(TMScoreSurvAllRaw,cancer_type %in% c('colon'))
colonSurv$MSI_status[colonSurv$MSI_status %in% c('MSI','MSI_H')]<-'MSI'
colonSurv$MSI_status[colonSurv$MSI_status %in% c('MSS')]<-'MSS'
colonSurv$MSI_status<-factor(colonSurv$MSI_status,levels = c('MSS','MSI'))
p1<-ggviolin(data = subset(colonSurv, MSI_status %in% c("MSI", "MSS")), 
             x = 'MSI_status', y = 'TMScoreSurvAllRaw', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'MSI_status', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  # stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('MSS','MSI')),
                     method='wilcox.test',label= 'p.signif',label.y = c(100)) # sig 



## Stomach cancer
stomachCaseCtrl<-TMScoreCaseCtrl[TMScoreCaseCtrl$cancer_type=='stomach',]
stomachCaseCtrl$MSI_feature[stomachCaseCtrl$MSI_feature %in% c('--')]<-NA
stomachCaseCtrl$MSI_feature<-factor(stomachCaseCtrl$MSI_feature,levels = c('MSS','MSI'))
ggviolin(data = subset(stomachCaseCtrl, MSI_feature %in% c("MSI", "MSS")), 
         x = 'MSI_feature', y = 'TMScoreCaseControl', facet.by = 'cancer_type',
         palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
         fill = 'MSI_feature', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  # stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('MSS','MSI')),
                     method='wilcox.test',label= 'p.signif',label.y = c(100)) # ns 

stomachSurv<-TMScoreSurvAllRaw[TMScoreSurvAllRaw$cancer_type=='stomach',]
stomachSurv$MSI_status<-factor(stomachSurv$MSI_status,levels = c('MSS','MSI'))
p2<-ggviolin(data = subset(stomachSurv, MSI_status %in% c("MSI", "MSS")), 
             x = 'MSI_status', y = 'TMScoreSurvAllRaw', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'MSI_status', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  # stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('MSS','MSI')),
                     method='wilcox.test',label= 'p.signif',label.y = c(100)) # sig
ggarrange(p1,p2,ncol = 2,widths = c(1,1),legend = 'top',common.legend = TRUE,font.label = list(size=6))


