options(stringsAsFactors = FALSE)


survivalSigBreastLung<-read.table('/data/survivalSigBreastLung.txt',
                                  row.names = 1,header = TRUE,sep = '\t')

load('/data/clinicalSurvival.RData')
load('/data/clinicalCaseControl.RData')
load('/data/expTMMcasecontrol.RData')
load('/data/expTMMsurv.RData')

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
TMScoreSurv<-TMSprobeWeightedSum(survivalSigBreastLung,expTMMsurv)

TMScoreCaseControl1<-unlist(TMScoreCaseControl)
names(TMScoreCaseControl1)<-do.call(rbind,strsplit(names(TMScoreCaseControl1),'.',fixed = TRUE))[,2]

TMScoreSurv1<-unlist(TMScoreSurv)
names(TMScoreSurv1)<-do.call(rbind,strsplit(names(TMScoreSurv1),'.',fixed = TRUE))[,2]


clinicalCaseControl<-cbind(clinicalCaseControl,TMScoreCaseControl=TMScoreCaseControl1[rownames(clinicalCaseControl)])

clinicalSurvival<-clinicalSurvivalUse[c(colnames(expTMMsurv$lung),colnames(expTMMsurv$breast), 
                                        colnames(expTMMsurv$colon),colnames(expTMMsurv$kidney),
                                        colnames(expTMMsurv$stomach)),]
clinicalSurvival<-cbind(clinicalSurvival,TMScoreSurv=TMScoreSurv1[rownames(clinicalSurvival)])


############################################################
##########################################

library(ggpubr)
p1<-ggviolin(data = subset(clinicalCaseControl, tissue_type %in% c('tumor', 'normal')),  nrow=1,
             x = 'tissue_type', y = 'TMScoreCaseControl', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'tissue_type', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  stat_compare_means(comparisons = list(c('tumor','normal')),method='wilcox.test',
                     aes(group=tissue_type),label= 'p.signif',label.y = c(95))

p2<-ggviolin(data = subset(clinicalSurvival, tissue_type %in% c('tumor', 'normal')), nrow=1, 
             x = 'tissue_type', y = 'TMScoreSurv', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'tissue_type', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  stat_compare_means(comparisons = list(c('tumor','normal')),method='wilcox.test',
                     aes(group=tissue_type),label= 'p.signif',label.y = c(95))

ggarrange(p1,p2,nrow = 2,legend = 'top',common.legend = TRUE,font.label = list(size=6),heights = c(1,1))


### stage
clinicalCaseControl$stage<-factor(clinicalCaseControl$stage,levels = c('I', 'II', 'III', 'IV'))
p1<-ggviolin(data = subset(clinicalCaseControl, stage %in% c('I', 'II', 'III', 'IV')), nrow=1,
             x = 'stage', y = 'TMScoreCaseControl', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'stage', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  stat_compare_means(label.y = 110)+
  stat_compare_means(comparisons = list(c('I','II'),c('II','III'),c('III','IV'),c('I','III'),c('II','IV'),c('I','IV')),
                     method='wilcox.test',label= 'p.signif',label.y = c(95,95,95,100,100,105)) # sig

clinicalSurvival$stage<-factor(clinicalSurvival$stage,levels = c('I', 'II', 'III', 'IV'))
p2<-ggviolin(data = subset(clinicalSurvival, stage %in% c('I', 'II', 'III', 'IV')), nrow=1,
             x = 'stage', y = 'TMScoreSurv', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'stage', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  stat_compare_means(label.y = 110)+
  stat_compare_means(comparisons = list(c('I','II'),c('II','III'),c('III','IV'),c('I','III'),c('II','IV'),c('I','IV')),
                     method='wilcox.test',label= 'p.signif',label.y = c(95,95,95,100,100,105)) # sig 
ggarrange(p1,p2,nrow = 2,legend = 'top',common.legend = TRUE,font.label = list(size=6),heights = c(1,1))


clinicalCaseControl$grade[clinicalCaseControl$grade %in% c('1','I','G1')]<-'G1'
clinicalCaseControl$grade[clinicalCaseControl$grade %in% c('2','II','G2')]<-'G2'
clinicalCaseControl$grade[clinicalCaseControl$grade %in% c('3','III','G3')]<-'G3'
clinicalCaseControl$grade[clinicalCaseControl$grade %in% c('4')]<-'G4'
clinicalCaseControl$grade[clinicalCaseControl$grade %in% c('N/A')]<-NA

clinicalCaseControl$grade<-factor(clinicalCaseControl$grade,levels = c('G1','G2','G3','G4'))
p1<-ggviolin(data = subset(clinicalCaseControl[clinicalCaseControl$cancer_type=='breast',], grade %in% c('G1','G2', 'G3')), nrow=1, 
             x = 'grade', y = 'TMScoreCaseControl', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'grade', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('G1','G2'),c('G2','G3'),c('G1','G3')),
                     method='wilcox.test',label= 'p.signif',label.y = c(95,95,100)) # sig


clinicalSurvival$grade<-factor(clinicalSurvival$grade,levels = c('G1','G2','G3','G4'))
p2<-ggviolin(data = subset(clinicalSurvival[clinicalSurvival$cancer_type=='breast',], grade %in% c('G1','G2', 'G3')),  nrow=1, 
             x = 'grade', y = 'TMScoreSurv', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'grade', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('G1','G2'),c('G2','G3'),c('G1','G3')),
                     method='wilcox.test',label= 'p.signif',label.y = c(95,95,100)) # sig 
p3<-ggarrange(p1,p2,nrow = 2,legend = 'top',common.legend = TRUE,font.label = list(size=6),heights = c(1,1))

### breast ER status
clinicalCaseControl$ER_IHC_status[clinicalCaseControl$ER_IHC_status %in% c('1','pos','POS','+')]<- 'pos'
clinicalCaseControl$ER_IHC_status[clinicalCaseControl$ER_IHC_status %in% c('0','neg','NEG','-')]<- 'neg'
clinicalCaseControl$ER_IHC_status[clinicalCaseControl$ER_IHC_status %in% c('N/A','EV')]<- NA

clinicalCaseControl$ER_IHC_status<-factor(clinicalCaseControl$ER_IHC_status,levels = c('pos','neg'))
p1<-ggviolin(data = subset(clinicalCaseControl, ER_IHC_status %in% c('neg','pos')), 
             x = 'ER_IHC_status', y = 'TMScoreCaseControl', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'ER_IHC_status', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  # stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('neg','pos')),
                     method='wilcox.test',label= 'p.signif',label.y = c(100)) # sig 


clinicalSurvival$ER_status<-factor(clinicalSurvival$ER_status,levels = c('pos','neg'))
p2<-ggviolin(data = subset(clinicalSurvival, ER_status %in% c('neg','pos')), 
             x = 'ER_status', y = 'TMScoreSurv', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'ER_status', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  # stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('neg','pos')),
                     method='wilcox.test',label= 'p.signif',label.y = c(100)) # sig
p4<-ggarrange(p1,p2,nrow = 2,legend = 'top',common.legend = TRUE,font.label = list(size=6),heights = c(1,1))
### breast PR status
clinicalCaseControl$PR_IHC_status[clinicalCaseControl$PR_IHC_status %in% c('1','pos','+')]<- 'pos'
clinicalCaseControl$PR_IHC_status[clinicalCaseControl$PR_IHC_status %in% c('0','neg','-')]<- 'neg'
clinicalCaseControl$PR_IHC_status[clinicalCaseControl$PR_IHC_status %in% c('N/A','EV')]<- NA

clinicalCaseControl$PR_IHC_status<-factor(clinicalCaseControl$PR_IHC_status,levels = c('pos','neg'))
p1<-ggviolin(data = subset(clinicalCaseControl, PR_IHC_status %in% c('neg','pos')), 
             x = 'PR_IHC_status', y = 'TMScoreCaseControl', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'PR_IHC_status', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  # stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('neg','pos')),
                     method='wilcox.test',label= 'p.signif',label.y = c(100)) # sig


clinicalSurvival$PgR_status<-factor(clinicalSurvival$PgR_status,levels = c('pos','neg'))
p2<-ggviolin(data = subset(clinicalSurvival, PgR_status %in% c('neg','pos')), 
             x = 'PgR_status', y = 'TMScoreSurv', facet.by = 'cancer_type',
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
clinicalCaseControl$molecular_subtype[clinicalCaseControl$molecular_subtype %in% c('Basal')]<- 'Basal'
clinicalCaseControl$molecular_subtype[clinicalCaseControl$molecular_subtype %in% c('ERBB2','HER2')]<- 'HER2'
clinicalCaseControl$molecular_subtype[clinicalCaseControl$molecular_subtype %in% c('Luminal A','LumA','LuminalA')]<- 'LumA'
clinicalCaseControl$molecular_subtype[clinicalCaseControl$molecular_subtype %in% c('Luminal B','LumB')]<- 'LumB'
clinicalCaseControl$molecular_subtype[clinicalCaseControl$molecular_subtype %in% c('Normal','normal')]<- 'Normal'
clinicalCaseControl$molecular_subtype[clinicalCaseControl$molecular_subtype %in% c('N/A')]<- NA

clinicalCaseControl$molecular_subtype<-factor(clinicalCaseControl$molecular_subtype,levels = c('LumA','Normal','LumB','HER2','Basal'))
p1<-ggviolin(data = subset(clinicalCaseControl, molecular_subtype %in% c('LumA','Normal','LumB','HER2','Basal')), 
             x = 'molecular_subtype', y = 'TMScoreCaseControl', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'molecular_subtype', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1/2)+
  stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('LumA','Normal'),c('Normal','LumB'),c('LumB','HER2'),c('HER2','Basal'),c('LumA','LumB'),
                                        c('Normal','HER2'),c('LumB','Basal'),c('LumA','HER2'),c('Normal','Basal'),c('LumA','Basal')),
                     method='wilcox.test',label= 'p.signif',label.y = c(90,90,90,90,95,95,95,100,100,105)) # sig 


subdata<-subset(clinicalSurvival, cancer_type %in% c('breast'))
subdata$molecular_subtype<-factor(subdata$molecular_subtype,levels = c('LumA','Normal','LumB','HER2','Basal'))
p2<-ggviolin(data = subset(subdata, molecular_subtype %in% c('LumA','Normal','LumB','HER2','Basal')), 
             x = 'molecular_subtype', y = 'TMScoreSurv', facet.by = 'cancer_type',
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



################ Lung cancer
## 1. subtyping
lungCaseCtrl<-clinicalCaseControl[clinicalCaseControl$cancer_type=='lung',]
lungCaseCtrl$subtyping<-factor(lungCaseCtrl$subtyping,levels = c('normal','LADC','LUSC','LCLC','SCLC'))
p1<-ggviolin(data = subset(lungCaseCtrl, subtyping %in% c('normal','LADC','LUSC','LCLC','SCLC')), 
             x = 'subtyping', y = 'TMScoreCaseControl', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'subtyping', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 0.5)+
  # stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('normal','LADC'),c('LADC','LUSC'),c('LUSC','LCLC'),c('SCLC','LCLC'),
                                        c('normal','LUSC'),c('LADC','LCLC'),c('LUSC','SCLC'),c('normal','LCLC'),
                                        c('LADC','SCLC'),c('normal','SCLC')),
                     method='wilcox.test',label= 'p.signif',label.y = c(105,105,105,105,110,110,110,115,115,120)) # sig 

lungSurv<-subset(clinicalSurvival,cancer_type %in% c('lung'))
lungSurv$subtyping1[lungSurv$subtyping1 %in% c('Squamous','SCC','squamous','SQC','Squamous Cell Carcinoma of Lung')]<-'LUSC'
lungSurv$subtyping1[lungSurv$subtyping1 %in% c('Adenocarcinoma','ADC','adeno','LADC')]<-'LADC'
lungSurv$subtyping1[lungSurv$subtyping1 %in% c('LCC','large')]<-'LCLC'
lungSurv$subtyping1[lungSurv$subtyping1 %in% c('healthy','NTL')]<-'normal'
lungSurv$subtyping1[lungSurv$subtyping1 %in% c('LCNE','CARCI','BAS','Other')]<-NA

lungSurv$subtyping1<-factor(lungSurv$subtyping1,levels = c('normal','LADC','LUSC','LCLC'))
p2<-ggviolin(data = subset(lungSurv, subtyping1 %in% c('normal','LADC','LUSC','LCLC')), 
             x = 'subtyping1', y = 'TMScoreSurv', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'subtyping1', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 0.5)+
  # stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('normal','LADC'),c('LADC','LUSC'),c('LUSC','LCLC'),
                                        c('normal','LUSC'),c('LADC','LCLC'),c('normal','LCLC')),
                     method='wilcox.test',label= 'p.signif',label.y = c(110,110,110,115,115,120)) # sig 
ggarrange(p1,p2,ncol = 2,widths = c(1,1),legend = 'top',common.legend = TRUE,font.label = list(size=6))


## Colon cancer
colonSurv<-subset(clinicalSurvival,cancer_type %in% c('colon'))
colonSurv$MSI_status[colonSurv$MSI_status %in% c('MSI','MSI_H')]<-'MSI'
colonSurv$MSI_status[colonSurv$MSI_status %in% c('MSS')]<-'MSS'
colonSurv$MSI_status<-factor(colonSurv$MSI_status,levels = c('MSS','MSI'))
p1<-ggviolin(data = subset(colonSurv, MSI_status %in% c("MSI", "MSS")), 
             x = 'MSI_status', y = 'TMScoreSurv', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'MSI_status', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  # stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('MSS','MSI')),
                     method='wilcox.test',label= 'p.signif',label.y = c(100)) # sig 



## Stomach cancer
stomachSurv<-clinicalSurvival[clinicalSurvival$cancer_type=='stomach',]
stomachSurv$MSI_status<-factor(stomachSurv$MSI_status,levels = c('MSS','MSI'))
p2<-ggviolin(data = subset(stomachSurv, MSI_status %in% c("MSI", "MSS")), 
             x = 'MSI_status', y = 'TMScoreSurv', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'MSI_status', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  # stat_compare_means(label.y = 105)+
  stat_compare_means(comparisons = list(c('MSS','MSI')),
                     method='wilcox.test',label= 'p.signif',label.y = c(100)) # sig
ggarrange(p1,p2,ncol = 2,widths = c(1,1),legend = 'top',common.legend = TRUE,font.label = list(size=6))


