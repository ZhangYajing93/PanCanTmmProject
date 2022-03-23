## load data
load('/data/TMMgenes.RData')
load('/data/probeTMM.RData')
load('/data/expTMMcasecontrol.RData')
load('/data/expTMMsurv.RData')


## Calculate differential expression probes
## use limma R package
library(limma)
GroupDesign<-function(tissuetype){
  group<-factor(tissuetype,levels = unique(tissuetype))
  groupdesign<-model.matrix(~0+group)
  colnames(groupdesign)<-unique(tissuetype)
  return(groupdesign)
}

groupdesign<-lapply(names(expTMMcasecontrol), function(x){
  GroupDesign(tissueTypeCaseCtrl[tissueTypeCaseCtrl$cancer_type==x,'tissue_type1'])
})
names(groupdesign)<-names(expTMMcasecontrol)

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
  DifExp(expTMMcasecontrol[[x]],groupdesign[[x]],x)
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
pheatmap(difexpprobeOrderFCplot[,rownames(probeTMMOrderDif)],
         cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = TRUE,cellwidth = 6,
         cellheight = 6,annotation_col = probeTMMOrderDif[,c(4,3)],
         display_numbers = difexpprobeOrderPplot[,rownames(probeTMMOrderDif)],
         labels_col = probeTMMOrderDif[,'SYMBOL'],
         angle_col = 90,fontsize_col = 6,filename = 'difTelomeraseProbe3up.pdf')

library(ggpubr)

pdf(file = 'difTelomeraseProbe3upBar.pdf')
ggbarplot(probeTMMOrderDifBar,x='probe',y='difexpCancers',
          width = 1,fill = 'black',alpha=0.8,color = 'white')
dev.off()


