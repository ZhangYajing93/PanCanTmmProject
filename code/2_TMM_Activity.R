######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## TMM activity by GSVA
options(stringsAsFactors = FALSE)
load('/data/probeTMMOrder.RData')
# load('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/OfficialVersion/ResourceData/expTMMcasecontrol.RData')
load('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/RawData/CaseControl/NormalizedAndClinicalDataToUse/CombatRemoveBatchEffect/expdataCombat.RData')
load('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/OfficialVersion/ResourceData/clinicalCaseControl.RData')


library(GSVA)

allTMMWayGsva<-lapply(expCombat, function(x){
  gsva(expr = x, gset.idx.list = split(rownames(probeTMMOrder),as.factor(probeTMMOrder$TMMway)), method = "gsva")
})
allpathWayGsva<-lapply(expCombat, function(x){
  gsva(expr = x, gset.idx.list = split(rownames(probeTMMOrder),as.factor(probeTMMOrder$pathway)), method = "gsva")
})

difTMMgsva<-lapply(expCombat, function(x){
  gsva(expr = x, gset.idx.list = list(rownames(probeTMMOrder[probeTMMOrder$difexpCancers %in% c(4:7),])), method = "gsva")
})
difTMMwayGsva<-lapply(expCombat, function(x){
  gsva(expr = x, gset.idx.list = split(rownames(probeTMMOrder[probeTMMOrder$difexpCancers %in% c(4:7),]),
                                       as.factor(probeTMMOrder[probeTMMOrder$difexpCancers %in% c(4:7),'TMMway'])), 
       method = "gsva")
})

allGsva<-lapply(expCombat, function(x){
  gsva(expr = x, gset.idx.list = list(rownames(probeTMMOrder)), method = "gsva")
})

setwd('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/OfficialVersion/difTMM/gsva')


load('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/OfficialVersion/ResourceData/tissueTypeCaseCtrl.RData')
allTMMWayGsva<-lapply(allTMMWayGsva, function(x){x[,intersect(rownames(tissueTypeCaseCtrl),colnames(x))]})
allpathWayGsva<-lapply(allpathWayGsva, function(x){x[,intersect(rownames(tissueTypeCaseCtrl),colnames(x))]})
difTMMgsva<-lapply(difTMMgsva, function(x){x[,intersect(rownames(tissueTypeCaseCtrl),colnames(x))]})
difTMMwayGsva<-lapply(difTMMwayGsva, function(x){x[,intersect(rownames(tissueTypeCaseCtrl),colnames(x))]})
allGsva<-lapply(allGsva, function(x){x[,intersect(rownames(tissueTypeCaseCtrl),colnames(x))]})
## change result to 0-1
allTMMWayGsva<-lapply(allTMMWayGsva, function(x){
  t(apply(x, 1, function(z){
    s=z-min(z)
    s=s/max(s)
  }))
})
allpathWayGsva<-lapply(allpathWayGsva, function(x){
  t(apply(x, 1, function(z){
    s=z-min(z)
    s=s/max(s)
  }))
})
difTMMgsva<-lapply(difTMMgsva, function(z){
  s=z-min(z)
  s=s/max(s)
})
difTMMwayGsva<-lapply(difTMMwayGsva, function(x){
  t(apply(x, 1, function(z){
    s=z-min(z)
    s=s/max(s)
  }))
})

allGsva<-lapply(allGsva, function(z){
  s=z-min(z)
  s=s/max(s)
})


allTMMWayGsva1<-cbind(allTMMWayGsva$lung,allTMMWayGsva$breast,allTMMWayGsva$colon,
                      allTMMWayGsva$liver,allTMMWayGsva$kidney,allTMMWayGsva$thyroid,
                      allTMMWayGsva$stomach)
allTMMWayGsva1<-allTMMWayGsva1[,rownames(tissueTypeCaseCtrl)]
rownames(allTMMWayGsva1)<-paste0(rownames(allTMMWayGsva1),'_allTMMWay')

allpathWayGsva1<-cbind(allpathWayGsva$lung,allpathWayGsva$breast,allpathWayGsva$colon,
                       allpathWayGsva$liver,allpathWayGsva$kidney,allpathWayGsva$thyroid,
                       allpathWayGsva$stomach)
allpathWayGsva1<-allpathWayGsva1[,rownames(tissueTypeCaseCtrl)]
rownames(allpathWayGsva1)<-paste0(rownames(allpathWayGsva1),'_allpathWay')

difTMMgsva1<-unlist(difTMMgsva)
names(difTMMgsva1)<-do.call(rbind,strsplit(names(difTMMgsva1),'.',fixed = TRUE))[,2]
difTMMgsva1<-difTMMgsva1[rownames(tissueTypeCaseCtrl)]

difTMMwayGsva1<-cbind(difTMMwayGsva$lung,difTMMwayGsva$breast,difTMMwayGsva$colon,
                      difTMMwayGsva$liver,difTMMwayGsva$kidney,difTMMwayGsva$thyroid,
                      difTMMwayGsva$stomach)
difTMMwayGsva1<-difTMMwayGsva1[,rownames(tissueTypeCaseCtrl)]
rownames(difTMMwayGsva1)<-paste0(rownames(difTMMwayGsva1),'_difTMMway')

allGsva1<-unlist(allGsva)
names(allGsva1)<-do.call(rbind,strsplit(names(allGsva1),'.',fixed = TRUE))[,2]
allGsva1<-allGsva1[rownames(tissueTypeCaseCtrl)]

tissueTypeCaseCtrl<-cbind(tissueTypeCaseCtrl,t(allTMMWayGsva1),t(allpathWayGsva1),difTMM=difTMMgsva1,t(difTMMwayGsva1),allTMM=allGsva1)

setwd('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/OfficialVersion/difTMM/gsva')
write.table(tissueTypeCaseCtrl,file = 'TMMgsvaAllGene.txt',row.names = TRUE,col.names = TRUE,sep = '\t')
save.image(file = 'TMMgsvaAllGene.RData')

#####################
TMMgsvaAllGene<-read.table('/data/TMMgsvaAllGene.txt',
                           header = TRUE,row.names = 1,sep = '\t')
library(ggpubr)
p<-ggviolin(data = subset(TMMgsvaAllGene, tissue_type2 %in% c('tumor', 'normal')), nrow=1,
            x = 'tissue_type2', y = 'allTMM', facet.by = 'cancer_type',
            palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
            fill = 'tissue_type2', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  stat_compare_means(comparisons = list(c('tumor','normal')),label= 'p.signif',method='wilcox.test') # sig
p1<-ggviolin(data = subset(TMMgsvaAllGene, tissue_type2 %in% c('tumor', 'normal')), nrow=1,
             x = 'tissue_type2', y = 'Telomerase_allTMMWay', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'tissue_type2', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  stat_compare_means(comparisons = list(c('tumor','normal')),label= 'p.signif',method='wilcox.test') # sig
p2<-ggviolin(data = subset(TMMgsvaAllGene, tissue_type2 %in% c('tumor', 'normal')), 
             x = 'tissue_type2', y = 'ALT_allTMMWay', facet.by = 'cancer_type',nrow=1,
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'tissue_type2', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  stat_compare_means(comparisons = list(c('tumor','normal')),label= 'p.signif',method='wilcox.test') # sig
ggarrange(p,p1,p2,nrow = 3,legend = 'top',common.legend = TRUE,font.label = list(size=6),heights = c(1,1,1))



p<-ggviolin(data = subset(TMMgsvaAllGene[TMMgsvaAllGene$cancer_type %in% c('colon','liver','thyroid'),], tissue_type2 %in% c('tumor','adenoma', 'normal')), 
            nrow=1,order = c('normal','adenoma','tumor'),
            x = 'tissue_type2', y = 'allTMM', facet.by = 'cancer_type',
            palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
            fill = 'tissue_type2', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  stat_compare_means(label.y= 1.2)+
  stat_compare_means(comparisons = list(c('normal','adenoma'),c('adenoma','tumor'),c('normal','tumor')),label= 'p.signif',
                     method='wilcox.test',label.y = c(1.05,1.1,1.15))
p1<-ggviolin(data = subset(TMMgsvaAllGene[TMMgsvaAllGene$cancer_type %in% c('colon','liver','thyroid'),], tissue_type2 %in% c('tumor','adenoma', 'normal')), 
             nrow=1,order = c('normal','adenoma','tumor'),
             x = 'tissue_type2', y = 'Telomerase_allTMMWay', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'tissue_type2', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  stat_compare_means(label.y= 1.2)+
  stat_compare_means(comparisons = list(c('normal','adenoma'),c('adenoma','tumor'),c('normal','tumor')),label= 'p.signif',
                     method='wilcox.test',label.y = c(1.05,1.1,1.15))
p2<-ggviolin(data = subset(TMMgsvaAllGene[TMMgsvaAllGene$cancer_type %in% c('colon','liver','thyroid'),], tissue_type2 %in% c('tumor','adenoma', 'normal')), 
             nrow=1,order = c('normal','adenoma','tumor'),
             x = 'tissue_type2', y = 'ALT_allTMMWay', facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),font.label = list(size=6),
             fill = 'tissue_type2', add = c("jitter"),add.params = list(size=0.2,alpha=0.5)) + 
  theme_classic()+
  theme(aspect.ratio = 1)+
  stat_compare_means(label.y= 1.2)+
  stat_compare_means(comparisons = list(c('normal','adenoma'),c('adenoma','tumor'),c('normal','tumor')),label= 'p.signif',
                     method='wilcox.test',label.y = c(1.05,1.1,1.15))
ggarrange(p,p1,p2,nrow = 3,legend = 'top',common.legend = TRUE,font.label = list(size=6),heights = c(1,1,1))


allTelomeraseplot<-TMMgsvaAllGene[,c('cancer_type','tissue_type2','TERT_pathway_allpathWay','TERC_DKC1_pathway_allpathWay')]
colnames(allTelomeraseplot)<-c("cancer_type",'tissue_type2')
allTelomeraseplot<-rbind(allTelomeraseplot[,c(1,2,3)],allTelomeraseplot[,c(1,2,4)])
allTelomeraseplot<-cbind(allTelomeraseplot,TMMway=c(rep('TERT',nrow(TMMgsvaAllGene)),rep('TERC',nrow(TMMgsvaAllGene))))
colnames(allTelomeraseplot)<-c('cancer_type','tissue_type2','TMMgsvaAllGene','TMMWay')

library(ggpubr)
p1<-ggviolin(data = subset(allTelomeraseplot, tissue_type2 %in% c('tumor', 'normal')), 
             x = 'TMMWay', y = 'TMMgsvaAllGene', color = 'tissue_type2',facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),fill = 'tissue_type2',width = 0.9,nrow=1, orientation = "horiz")+
  theme_classic()+
  theme(aspect.ratio = 1)+
  stat_compare_means(aes(group=tissue_type2),label= 'p.signif')

allALTplot<-TMMgsvaAllGene[,c('cancer_type','tissue_type2','HR_allpathWay','Chromatin_Decompaction_pathway_allpathWay','PML_allpathWay','TERRA.Telomere_instability_allpathWay')]
colnames(allALTplot)<-c("cancer_type",'tissue_type2')
allALTplot<-rbind(allALTplot[,c(1,2,3)],allALTplot[,c(1,2,4)],allALTplot[,c(1,2,5)],allALTplot[,c(1,2,6)])
allALTplot<-cbind(allALTplot,TMMway=c(rep('HR',nrow(TMMgsvaAllGene)),rep('ChDe',nrow(TMMgsvaAllGene)),rep('PML',nrow(TMMgsvaAllGene)),rep('TERRA',nrow(TMMgsvaAllGene))))
colnames(allALTplot)<-c('cancer_type','tissue_type2','TMMgsvaAllGene','TMMWay')

library(ggpubr)
p2<-ggviolin(data = subset(allALTplot, tissue_type2 %in% c('tumor', 'normal')), 
             x = 'TMMWay', y = 'TMMgsvaAllGene', color = 'tissue_type2',facet.by = 'cancer_type',
             palette = 'nc',draw_quantiles = c(0.25,0.5,0.75),fill = 'tissue_type2',width = 0.9, orientation = "horiz",nrow=1)+
  theme_classic()+
  theme(aspect.ratio = 2)+
  stat_compare_means(aes(group=tissue_type2),label= 'p.signif')
ggarrange(p1,p2,nrow = 2,legend = 'top',common.legend = TRUE,font.label = list(size=6),heights = c(1,1))


