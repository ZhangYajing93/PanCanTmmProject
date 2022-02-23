##‘ @ Choose probes corresponding to TMM-related genes

load('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/RawData/CaseControl/Results/TMMgenelist/TMMgenelist.RData')
load('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/RawData/CaseControl/NormalizedAndClinicalDataToUse/CombatRemoveBatchEffect/expdataCombat.RData')
load('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/RawData/SurvivalDataset/NormalizedDataSurvival/expdataSurvivalCombatRaw.RData')
setwd('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/RawData/CaseControl/NormalizedAndClinicalDataToUse/CombatRemoveBatchEffect')
gpl570<-read.csv("/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/RawData/CaseControl/NormalizedAndClinicalDataToUse/CombatRemoveBatchEffect/GPL570-55999.txt",
                 comment.char = "#",sep = "\t",header = T,stringsAsFactors = F,row.names = 1)
gpl570<-gpl570[,c(10,11)]
colnames(gpl570)<-c('SYMBOL','ENTREZID')

probeTMM<-gpl570[!is.na(match(gpl570$ENTREZID,TMMset$ENTREZID)),] #587 probes of 218 unique genes

probeTMMOrder<-lapply(1:nrow(TMMgene), function(x){
  probeTMM[probeTMM$ENTREZID==TMMgene$ENTREZID[x],]
})
probeTMMOrder<-do.call(rbind,probeTMMOrder)

# TMMset[is.na(match(TMMset$ENTREZID,unique(probeTMM$ENTREZID))),]

expTMM<-lapply(expCombat, function(x){
  x[rownames(probeTMM),]
})

expTMMsurvAllRaw<-lapply(norExpDataSurvivalCombatRaw, function(x){
  x[rownames(probeTMM),]
})
setwd('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/NewRoute/ResourceData')
save(expTMMsurvAllRaw,file = 'expTMMsurvAllRaw.RData')


## Calculate differential expression probes
## use limma R package
library(limma)
GroupDesign<-function(tissuetype){
  group<-factor(tissuetype,levels = unique(tissuetype))
  groupdesign<-model.matrix(~0+group)
  colnames(groupdesign)<-unique(tissuetype)
  return(groupdesign)
}

groupdesign<-lapply(names(expTMM), function(x){
  GroupDesign(tissueType[tissueType$cancer_type==x,'tissue_type1'])
})
names(groupdesign)<-names(expTMM)

## change low/high-grade dysplastic in liver cancer
colnames(groupdesign$liver)<-c('tumor','normal','cirrhosis',
                               'low_grade_dysplastic','high_grade_dysplastic')



DifExp<-function(eset,designgroup,tissuetype){
  fit<-lmFit(eset,designgroup)
  if(tissuetype%in%c('lung','breast','kidney','stomach')){
    cont.matrix<-makeContrasts(TvsN=tumor-normal,levels = designgroup)
    fit2<-contrasts.fit(fit,cont.matrix)
    fit2<-eBayes(fit2)
  }else if(tissuetype=='colon'){
    cont.matrix<-makeContrasts(colitis-normal,adenoma-colitis,tumor-adenoma,
                               adenoma-normal,tumor-colitis,
                               tumor-normal,levels = designgroup)
    fit2<-contrasts.fit(fit,cont.matrix)
    fit2<-eBayes(fit2)
  }else if(tissuetype=='liver'){
    cont.matrix<-makeContrasts(cirrhosis-normal,low_grade_dysplastic-cirrhosis,high_grade_dysplastic-low_grade_dysplastic,tumor-high_grade_dysplastic,
                               low_grade_dysplastic-normal,high_grade_dysplastic-cirrhosis,tumor-low_grade_dysplastic,
                               high_grade_dysplastic-normal,tumor-cirrhosis,
                               tumor-normal,levels = designgroup)
    fit2<-contrasts.fit(fit,cont.matrix)
    fit2<-eBayes(fit2)
  }else{
    cont.matrix<-makeContrasts(adenoma-normal,tumor-adenoma,
                               tumor-normal,levels = designgroup)
    fit2<-contrasts.fit(fit,cont.matrix)
    fit2<-eBayes(fit2)
  }
  
}

difexpprobe<-lapply(names(groupdesign), function(x){
  DifExp(expTMM[[x]],groupdesign[[x]],x)
})
names(difexpprobe)<-names(groupdesign) 
# topTable(difexpprobe$lung,coef = 1,number = 10)


difexpprobe<-lapply(difexpprobe, function(x){
  topTable(x,coef = ncol(x[['coefficients']]),number = 587,adjust.method = 'fdr')
})

difexpprobeOrder<-lapply(difexpprobe, function(x){x[rownames(probeTMMOrder),]})

## logFC 
difexpprobeOrderFC<-lapply(difexpprobeOrder, function(x){x[,'logFC']})
difexpprobeOrderFC<-do.call(cbind,difexpprobeOrderFC)
rownames(difexpprobeOrderFC)<-rownames(probeTMMOrder)
difexpprobeOrderFC<-t(difexpprobeOrderFC)
## difexpprobeOrder
difexpprobeOrderP<-lapply(difexpprobeOrder, function(x){x[,'adj.P.Val']})
difexpprobeOrderP<-do.call(cbind,difexpprobeOrderP)
rownames(difexpprobeOrderP)<-rownames(probeTMMOrder)
difexpprobeOrderP<-t(difexpprobeOrderP)


difexpprobeOrderFCplot<-difexpprobeOrderFC
difexpprobeOrderFCplot[which(difexpprobeOrderFCplot>log2(1.5))]<-1
difexpprobeOrderFCplot[which(difexpprobeOrderFCplot<log2(1/1.5))]<-(-1)
difexpprobeOrderFCplot[which(abs(difexpprobeOrderFCplot)<log2(1.5))]<-0


difexpprobeOrderPplot<-matrix(ifelse(difexpprobeOrderP < 0.05, "*", ""), nrow(difexpprobeOrderP))
colnames(difexpprobeOrderPplot)<-colnames(difexpprobeOrderP)
rownames(difexpprobeOrderPplot)<-rownames(difexpprobeOrderP)


probeTMMOrderDif<-probeTMMOrder[probeTMMOrder$difexpCancers %in% c(seq(1,7)),]
library(pheatmap)

pheatmap(difexpprobeOrderFCplot[,rownames(probeTMMOrderDif[probeTMMOrderDif$TMMway=='Telomerase',])],
         cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = TRUE,cellwidth = 6,
         cellheight = 6,annotation_col = probeTMMOrderDif[probeTMMOrderDif$TMMway=='Telomerase',c(4,3)],
         display_numbers = difexpprobeOrderPplot[,rownames(probeTMMOrderDif[probeTMMOrderDif$TMMway=='Telomerase',])],
         labels_col = probeTMMOrderDif[probeTMMOrderDif$TMMway=='Telomerase','SYMBOL'],
         angle_col = 90,fontsize_col = 6,filename = 'difTelomeraseProbe.pdf')

pheatmap(difexpprobeOrderFCplot[,rownames(probeTMMOrderDif[probeTMMOrderDif$TMMway=='ALT',])],
         cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = TRUE,cellwidth = 6,
         cellheight = 6,annotation_col = probeTMMOrderDif[probeTMMOrderDif$TMMway=='ALT',c(4,3)],
         display_numbers = difexpprobeOrderPplot[,rownames(probeTMMOrderDif[probeTMMOrderDif$TMMway=='ALT',])],
         labels_col = probeTMMOrderDif[probeTMMOrderDif$TMMway=='ALT','SYMBOL'],
         angle_col = 90,fontsize_col = 6,filename = 'difALTProbe.pdf')
dev.off()




library(ggpubr)
probeTMMOrderDifBar<-cbind(probe=rownames(probeTMMOrderDif),probeTMMOrderDif)
setwd('E:/0/telomeres/Results/OfficialVersion/TMMgene/pheatmap/barplot')
pdf(file = 'difTelomeraseProbe.pdf')
ggbarplot(probeTMMOrderDifBar[probeTMMOrderDifBar$TMMway=='Telomerase',],x='probe',y='difexpCancers',
          width = 1,fill = 'black',alpha=0.8,color = 'white')
dev.off()
pdf(file = 'difALTProbe.pdf')
ggbarplot(probeTMMOrderDifBar[probeTMMOrderDifBar$TMMway=='ALT',],x='probe',y='difexpCancers',
          width = 1,fill = 'black',alpha=0.8,color = 'white')
dev.off()

probeTMMOrderDif<-probeTMMOrder[probeTMMOrder$difexpCancers %in% c(seq(3,7)),]
library(pheatmap)
setwd('E:/0/telomeres/Results/OfficialVersion/TMMgene/pheatmap')
pheatmap(difexpprobeOrderFCplot[,rownames(probeTMMOrderDif)],
         cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = TRUE,cellwidth = 6,
         cellheight = 6,annotation_col = probeTMMOrderDif[,c(4,3)],
         display_numbers = difexpprobeOrderPplot[,rownames(probeTMMOrderDif)],
         labels_col = probeTMMOrderDif[,'SYMBOL'],
         angle_col = 90,fontsize_col = 6,filename = 'difTelomeraseProbe3up.pdf')

library(ggpubr)
probeTMMOrderDifBar<-cbind(probe=rownames(probeTMMOrderDif),probeTMMOrderDif)
setwd('E:/0/telomeres/Results/OfficialVersion/TMMgene/pheatmap/barplot')
pdf(file = 'difTelomeraseProbe3up.pdf')
ggbarplot(probeTMMOrderDifBar,x='probe',y='difexpCancers',
          width = 1,fill = 'black',alpha=0.8,color = 'white')
dev.off()


