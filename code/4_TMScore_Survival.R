## TMScore Survival Analysis
survivalSigBreastLung<-read.table('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/NewRoute/ProbeCalculate/survival/BreastLungSigThreeDatabase/survivalSigBreastLung.txt',
                                  row.names = 1,header = TRUE,sep = '\t')

load('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/NewRoute/ResourceData/expTMMsurvAllRaw.RData')
load('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/NewRoute/ResourceData/clinicalSurvivalUse.RData')
clinicalSurvivalUse<-clinicalSurvivalUse[c(colnames(expTMMsurvAllRaw$lung),colnames(expTMMsurvAllRaw$breast),
                                           colnames(expTMMsurvAllRaw$colon),colnames(expTMMsurvAllRaw$kidney),
                                           colnames(expTMMsurvAllRaw$stomach)),]

TMSprobeWeightedSum<-function(TMSgenetable,expdataset){
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
TMScore<-TMSprobeWeightedSum(survivalSigBreastLung,expTMMsurvAllRaw)


SurvDataThreeGroupAll<-function(clinical.data,surv.type,sample.score){
  ## surv.type: 'OS','PFS','RFS','DFS','DMFS','MFS','DSS'
  ## sample.score:TMScore,list 
  clisur<-clinical.data[clinical.data$tissue_type %in% c('tumor'),]
  clidata<-switch(surv.type,
                  OS=subset(clisur,!(is.na(OS_status)|is.na(OS_month)))[,c('cancer_type','GEO_number','geo_accession','OS_status','OS_month')],
                  PFS=subset(clisur,!(is.na(PFS_status)|is.na(PFS_month)))[,c('cancer_type','GEO_number','geo_accession','PFS_status','PFS_month')],
                  RFS=subset(clisur,!(is.na(RFS_status)|is.na(RFS_month)))[,c('cancer_type','GEO_number','geo_accession','RFS_status','RFS_month')],
                  DFS=subset(clisur,!(is.na(DFS_status)|is.na(DFS_month)))[,c('cancer_type','GEO_number','geo_accession','DFS_status','DFS_month')],
                  DMFS=subset(clisur,!(is.na(DMFS_status)|is.na(DMFS_month)))[,c('cancer_type','GEO_number','geo_accession','DMFS_status','DMFS_month')],
                  MFS=subset(clisur,!(is.na(MFS_status)|is.na(MFS_month)))[,c('cancer_type','GEO_number','geo_accession','MFS_status','MFS_month')],
                  DSS=subset(clisur,!(is.na(DSS_status)|is.na(DSS_month)))[,c('cancer_type','GEO_number','geo_accession','DSS_status','DSS_month')])
  sam.score<-lapply(sample.score, function(x){x[intersect(rownames(clidata),names(x))]})
  sam.score<-sam.score[lapply(sam.score,length)!=0]
  tertiles<-lapply(sam.score, function(x){quantile(x,c(0,1/3,2/3,1))})
  Label<-lapply(1:length(sam.score), function(x){
    cut(sam.score[[x]],tertiles[[x]],labels = FALSE,include.lowest = TRUE)
  })
  names(Label)<-names(sam.score)
  for (i in 1:length(sam.score)) {
    names(Label[[i]])<-names(sam.score[[i]])
  }
  Label<-unlist(Label)
  a<-names(Label)
  b<-strsplit(a,'.',fixed = TRUE)
  c<-do.call(rbind,b)
  d<-c[,2]
  names(Label)<-d
  common <- intersect(rownames(clidata), names(Label))
  clidata<-clidata[common,]
  Label<-Label[common]
  clidata<-cbind(clidata,Label)
  colnames(clidata)<-c('cancer_type','GEO_number','Patient_ID','event','time','sample.label')
  clidata$event<-as.numeric(clidata$event)
  clidata$time<-as.numeric(clidata$time)
  clidata$sample.label<-as.factor(clidata$sample.label)
  clidata<-split(clidata,as.factor(clidata$cancer_type))
  clidata<-lapply(clidata, function(x){
    split(x,as.factor(x$cancer_type))
  })
  return(clidata)
}

surcliOSTertilesWeight<-SurvDataThreeGroupAll(clinicalSurvivalUse,'OS',TMScore)
surcliPFSTertilesWeight<-SurvDataThreeGroupAll(clinicalSurvivalUse,'PFS',TMScore)
surcliRFSTertilesWeight<-SurvDataThreeGroupAll(clinicalSurvivalUse,'RFS',TMScore)
surcliDFSTertilesWeight<-SurvDataThreeGroupAll(clinicalSurvivalUse,'DFS',TMScore)
surcliDMFSTertilesWeight<-SurvDataThreeGroupAll(clinicalSurvivalUse,'DMFS',TMScore)
surcliMFSTertilesWeight<-SurvDataThreeGroupAll(clinicalSurvivalUse,'MFS',TMScore)
surcliDSSTertilesWeight<-SurvDataThreeGroupAll(clinicalSurvivalUse,'DSS',TMScore)

source(file = '/pub5/xiaoyun/Jobs/J22/ZYJ2020/Rscript/km_survival_plot.R')

setwd('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/NewRoute/ProbeCalculate/survival/BreastLungSigThreeDatabase/survival/survAllRaw')

for (i in 1:length(surcliOSTertilesWeight)) {
  pdf(file = paste0('OS',names(surcliOSTertilesWeight)[i],'.pdf'))
  lapply(1:length(surcliOSTertilesWeight[[i]]), function(x){
    print(plot.surv(surcliOSTertilesWeight[[i]][[x]], upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
                    surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
                    conf.int = TRUE, main = paste0(names(surcliOSTertilesWeight)[i],names(surcliOSTertilesWeight[[i]])[[x]]), ylab = "OS"))
  })
  dev.off()
}

for (i in 1:length(surcliPFSTertilesWeight)) {
  pdf(file = paste0('PFS',names(surcliPFSTertilesWeight)[i],'.pdf'))
  lapply(1:length(surcliPFSTertilesWeight[[i]]), function(x){
    print(plot.surv(surcliPFSTertilesWeight[[i]][[x]], upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
                    surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
                    conf.int = TRUE, main = paste0(names(surcliPFSTertilesWeight)[i],names(surcliPFSTertilesWeight[[i]])[[x]]), ylab = "PFS"))
  })
  dev.off()
}

for (i in 1:length(surcliRFSTertilesWeight)) {
  pdf(file = paste0('RFS',names(surcliRFSTertilesWeight)[i],'.pdf'))
  lapply(1:length(surcliRFSTertilesWeight[[i]]), function(x){
    print(plot.surv(surcliRFSTertilesWeight[[i]][[x]], upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
                    surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
                    conf.int = TRUE, main = paste0(names(surcliRFSTertilesWeight)[i],names(surcliRFSTertilesWeight[[i]])[[x]]), ylab = "RFS"))
  })
  dev.off()
}

for (i in 1:length(surcliDFSTertilesWeight)) {
  pdf(file = paste0('DFS',names(surcliDFSTertilesWeight)[i],'.pdf'))
  lapply(1:length(surcliDFSTertilesWeight[[i]]), function(x){
    print(plot.surv(surcliDFSTertilesWeight[[i]][[x]], upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
                    surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
                    conf.int = TRUE, main = paste0(names(surcliDFSTertilesWeight)[i],names(surcliDFSTertilesWeight[[i]])[[x]]), ylab = "DFS"))
  })
  dev.off()
}

for (i in 1:length(surcliDMFSTertilesWeight)) {
  pdf(file = paste0('DMFS',names(surcliDMFSTertilesWeight)[i],'.pdf'))
  lapply(1:length(surcliDMFSTertilesWeight[[i]]), function(x){
    print(plot.surv(surcliDMFSTertilesWeight[[i]][[x]], upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
                    surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
                    conf.int = TRUE, main = paste0(names(surcliDMFSTertilesWeight)[i],names(surcliDMFSTertilesWeight[[i]])[[x]]), ylab = "DMFS"))
  })
  dev.off()
}

for (i in 1:length(surcliMFSTertilesWeight)) {
  pdf(file = paste0('MFS',names(surcliMFSTertilesWeight)[i],'.pdf'))
  lapply(1:length(surcliMFSTertilesWeight[[i]]), function(x){
    print(plot.surv(surcliMFSTertilesWeight[[i]][[x]], upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
                    surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
                    conf.int = TRUE, main = paste0(names(surcliMFSTertilesWeight)[i],names(surcliMFSTertilesWeight[[i]])[[x]]), ylab = "MFS"))
  })
  dev.off()
}

for (i in 1:length(surcliDSSTertilesWeight)) {
  pdf(file = paste0('DSS',names(surcliDSSTertilesWeight)[i],'.pdf'))
  lapply(1:length(surcliDSSTertilesWeight[[i]]), function(x){
    print(plot.surv(surcliDSSTertilesWeight[[i]][[x]], upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
                    surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
                    conf.int = TRUE, main = paste0(names(surcliDSSTertilesWeight)[i],names(surcliDSSTertilesWeight[[i]])[[x]]), ylab = "DSS"))
  })
  dev.off()
}





############################################ Cox regression analysis#######################################################
################ breast cancer ################ 
source('/pub5/xiaoyun/Jobs/J22/EvoClass2.0/Section3/RScripts/Cox.function.R')
clinicalSurvivalUse<-clinicalSurvivalUse[c(colnames(expTMMsurvAllRaw$lung),colnames(expTMMsurvAllRaw$breast),
                                           colnames(expTMMsurvAllRaw$colon),colnames(expTMMsurvAllRaw$kidney),
                                           colnames(expTMMsurvAllRaw$stomach)),]
clinicalSurvivalUse$HER2_status[clinicalSurvivalUse$HER2_status %in% c('1','over-expression','HER2-positive breast cancer','HER2 breast cancer tumor','positive','pos','yes','Positive')]<-'pos'
clinicalSurvivalUse$HER2_status[clinicalSurvivalUse$HER2_status %in% c('0','normal','negative','neg','no','Negative')]<-'neg'

clinicalSurvivalUse$ER_status[clinicalSurvivalUse$ER_status %in% c('1','positive','pos','Positive')]<- 'pos'
clinicalSurvivalUse$ER_status[clinicalSurvivalUse$ER_status %in% c('0','negative','neg','Negative')]<- 'neg'
clinicalSurvivalUse$ER_status[clinicalSurvivalUse$ER_status %in% c('EV')]<- NA

clinicalSurvivalUse$PgR_status[clinicalSurvivalUse$PgR_status %in% c('1','positive','pos','Positive')]<- 'pos'
clinicalSurvivalUse$PgR_status[clinicalSurvivalUse$PgR_status %in% c('0','negative','neg','Negative')]<- 'neg'
clinicalSurvivalUse$PgR_status[clinicalSurvivalUse$PgR_status %in% c('EV')]<- NA

clinicalSurvivalUse$node_status[clinicalSurvivalUse$node_status %in% c('1','6','2','11','8','10','7','20','4','5','3','14','22',
                                                                       '16','15','17','13','9','18','28','12','positive','pos','20')]<-'pos'
clinicalSurvivalUse$node_status[clinicalSurvivalUse$node_status %in% c('0','negative')]<-'neg'

clinicalSurvivalUse$age<-as.numeric(clinicalSurvivalUse$age)
clinicalSurvivalUse$OS_status<-as.numeric(clinicalSurvivalUse$OS_status)
clinicalSurvivalUse$OS_month<-as.numeric(clinicalSurvivalUse$OS_month)
clinicalSurvivalUse$DFS_status<-as.numeric(clinicalSurvivalUse$DFS_status)
clinicalSurvivalUse$DFS_month<-as.numeric(clinicalSurvivalUse$DFS_month)
clinicalSurvivalUse$RFS_status<-as.numeric(clinicalSurvivalUse$RFS_status)
clinicalSurvivalUse$RFS_month<-as.numeric(clinicalSurvivalUse$RFS_month)
clinicalSurvivalUse$DMFS_status<-as.numeric(clinicalSurvivalUse$DMFS_status)
clinicalSurvivalUse$DMFS_month<-as.numeric(clinicalSurvivalUse$DMFS_month)


clinicalSurvivalUse$molecular_subtype[clinicalSurvivalUse$molecular_subtype %in% c('HER2','HER2-enriched','HER2-positive breast cancer','HER2 breast cancer tumor','ERBB2','Her2')]<-'HER2'
clinicalSurvivalUse$molecular_subtype[clinicalSurvivalUse$molecular_subtype %in% c('LumB','Luminal B','LuminalB')]<-'LumB'
clinicalSurvivalUse$molecular_subtype[clinicalSurvivalUse$molecular_subtype %in% c('Basal','Basal-like subtype','triple negative breast cancer (TNBC) patient','triple negative breast cancer','TNBC')]<-'Basal'
clinicalSurvivalUse$molecular_subtype[clinicalSurvivalUse$molecular_subtype %in% c('LumA','Luminal A','LuminalA')]<-'LumA'
clinicalSurvivalUse$molecular_subtype[clinicalSurvivalUse$molecular_subtype %in% c('Normal','Normal breast-like')]<-'Normal'
clinicalSurvivalUse$molecular_subtype[clinicalSurvivalUse$molecular_subtype %in% c('Uncla','--','Healthy')]<-NA

clinicalSurvivalUse$grade[clinicalSurvivalUse$grade %in% c('1','G1','WELL DIFF','I')]<-'G1'
clinicalSurvivalUse$grade[clinicalSurvivalUse$grade %in% c('2','G2','MOD DIFF','II')]<-'G2'
clinicalSurvivalUse$grade[clinicalSurvivalUse$grade %in% c('3','G3','POORLY DIFF','III')]<-'G3'
clinicalSurvivalUse$grade[clinicalSurvivalUse$grade %in% c('4','IV')]<-'G4'
clinicalSurvivalUse$grade[clinicalSurvivalUse$grade %in% c('Unknown')]<-NA

clinicalSurvivalUse$stage[clinicalSurvivalUse$stage %in% c('1','1a','1b','IA','IB','I','Stage 1')]<-'I'
clinicalSurvivalUse$stage[clinicalSurvivalUse$stage %in% c('2','2b','2a','II','IIB','IIA','Stage 2')]<-'II'
clinicalSurvivalUse$stage[clinicalSurvivalUse$stage %in% c('3','3a','3b','IIIA','IIIB','IIIC','Stage 3','III')]<-'III'
clinicalSurvivalUse$stage[clinicalSurvivalUse$stage %in% c('4','IV','Stage 4')]<-'IV'
clinicalSurvivalUse$stage[clinicalSurvivalUse$stage %in% c('N/A','0','Unknown','Recurrence','Unkown')]<-NA

clinicalSurvivalUse$stage_T[clinicalSurvivalUse$stage_T %in% c('T1','T1b','T1a','T1c','1','1c','pT1','pT1b','pT1a')]<-'T1'
clinicalSurvivalUse$stage_T[clinicalSurvivalUse$stage_T %in% c('T2','T2a','T2b','2','pT2')]<-'T2'
clinicalSurvivalUse$stage_T[clinicalSurvivalUse$stage_T %in% c('T3','3','pT3','pT3b','pT3a')]<-'T3'
clinicalSurvivalUse$stage_T[clinicalSurvivalUse$stage_T %in% c('T4','4','pT4')]<-'T4'
clinicalSurvivalUse$stage_T[clinicalSurvivalUse$stage_T %in% c('T0','Tis','0')]<-'T0'
clinicalSurvivalUse$stage_T[clinicalSurvivalUse$stage_T %in% c('NTL','TX','N/A','pTX')]<-NA

clinicalSurvivalUse$stage_N[clinicalSurvivalUse$stage_N %in% c('N1','N1a','N1mi','N1b','N1c','N1biii','N1bi','1','pN1')]<-'N1'
clinicalSurvivalUse$stage_N[clinicalSurvivalUse$stage_N %in% c('N2','N2a','2','pN2')]<-'N2'
clinicalSurvivalUse$stage_N[clinicalSurvivalUse$stage_N %in% c('N3','N3b','N3a','3')]<-'N3'
clinicalSurvivalUse$stage_N[clinicalSurvivalUse$stage_N %in% c('N0','N0 (i-)','N0 (i+)','0','pN0')]<-'N0'
clinicalSurvivalUse$stage_N[clinicalSurvivalUse$stage_N %in% c('NTL','NX','N/A','X','N+','pNX')]<-NA

clinicalSurvivalUse$stage_M[clinicalSurvivalUse$stage_M %in% c('M0','0')]<-'M0'
clinicalSurvivalUse$stage_M[clinicalSurvivalUse$stage_M %in% c('M1','1')]<-'M1'
clinicalSurvivalUse$stage_M[clinicalSurvivalUse$stage_M %in% c('MX','NTL','X','N/A')]<-NA

clinicalSurvivalUse$MMR_status[clinicalSurvivalUse$MMR_status %in% c('N/A')]<-NA
clinicalSurvivalUse$KRAS_mutation_status[clinicalSurvivalUse$KRAS_mutation_status %in% c('N/A')]<-NA
clinicalSurvivalUse$P53_status[clinicalSurvivalUse$P53_status %in% c('N/A')]<-NA
clinicalSurvivalUse$BRAF_mutation[clinicalSurvivalUse$BRAF_mutation %in% c('N/A')]<-NA
clinicalSurvivalUse$CIN_status[clinicalSurvivalUse$CIN_status %in% c('N/A')]<-NA
clinicalSurvivalUse$CIMP_status[clinicalSurvivalUse$CIMP_status %in% c('N/A')]<-NA

clinicalSurvivalUse$smoking_status[clinicalSurvivalUse$smoking_status %in% c('1','Ever-smoker')]<-'1'
clinicalSurvivalUse$smoking_status[clinicalSurvivalUse$smoking_status %in% c('0','Never-smoker')]<-'0'

clinicalSurvivalUse[clinicalSurvivalUse$cancer_type %in% c('lung'),'subtyping1'][clinicalSurvivalUse[clinicalSurvivalUse$cancer_type %in% c('lung'),'subtyping1'] %in% c('Squamous','SCC','squamous','SQC','Squamous Cell Carcinoma of Lung')]<-'LUSC'
clinicalSurvivalUse[clinicalSurvivalUse$cancer_type %in% c('lung'),'subtyping1'][clinicalSurvivalUse[clinicalSurvivalUse$cancer_type %in% c('lung'),'subtyping1'] %in% c('Adenocarcinoma','ADC','adeno','LADC')]<-'LADC'
clinicalSurvivalUse[clinicalSurvivalUse$cancer_type %in% c('lung'),'subtyping1'][clinicalSurvivalUse[clinicalSurvivalUse$cancer_type %in% c('lung'),'subtyping1'] %in% c('LCC','large')]<-'LCLC'
clinicalSurvivalUse[clinicalSurvivalUse$cancer_type %in% c('lung'),'subtyping1'][clinicalSurvivalUse[clinicalSurvivalUse$cancer_type %in% c('lung'),'subtyping1'] %in% c('LCNE','CARCI','BAS','Other')]<-NA
clinicalSurvivalUse[clinicalSurvivalUse$cancer_type %in% c('lung'),'subtyping1'][clinicalSurvivalUse[clinicalSurvivalUse$cancer_type %in% c('lung'),'subtyping1'] %in% c('healthy','NTL')]<-'normal'


load('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/NewRoute/ProbeCalculate/survival/BreastLungSigThreeDatabase/survival/x.RData')
clinicalSurvivalUse[clinicalSurvivalUse$GEO_number=='GSE39582','MSI_status']<-x$MSI_status
clinicalSurvivalUse$MSI_status[clinicalSurvivalUse$MSI_status %in% c('N/A')]<-NA
clinicalSurvivalUse$MSI_status[clinicalSurvivalUse$MSI_status %in% c('MSI','MSI_H')]<-'MSI'
clinicalSurvivalUse$MSI_status[clinicalSurvivalUse$MSI_status %in% c('MSS')]<-'MSS'

setwd('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/NewRoute/ResourceData')
save(clinicalSurvivalUse,file = 'clinicalSurvivalUse.RData')

CancerSurv<-function(clinical.data,surv.type){
  clisur<-clinical.data[clinical.data$tissue_type %in% c('tumor'),]
  clidata<-switch(surv.type,
                  OS=subset(clisur,!(is.na(OS_status)|is.na(OS_month))),
                  PFS=subset(clisur,!(is.na(PFS_status)|is.na(PFS_month))),
                  RFS=subset(clisur,!(is.na(RFS_status)|is.na(RFS_month))),
                  DFS=subset(clisur,!(is.na(DFS_status)|is.na(DFS_month))),
                  DMFS=subset(clisur,!(is.na(DMFS_status)|is.na(DMFS_month))),
                  MFS=subset(clisur,!(is.na(MFS_status)|is.na(MFS_month))),
                  DSS=subset(clisur,!(is.na(DSS_status)|is.na(DSS_month))))
}

cancerSurv<-split(clinicalSurvivalUse,factor(clinicalSurvivalUse$cancer_type))
cancerSurv<-lapply(names(cancerSurv), function(x){cbind(cancerSurv[[x]],TMScore=TMScore[[x]][rownames(cancerSurv[[x]])])})
names(cancerSurv)<-names(split(clinicalSurvivalUse,factor(clinicalSurvivalUse$cancer_type)))



############### breast cancer OS Cox analysis
cancerSurvOS<-lapply(cancerSurv, function(x){
  CancerSurv(x,'OS')
})

setwd('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/NewRoute/ProbeCalculate/survival/BreastLungSigThreeDatabase/survival/survAllRaw/Cox')
save(cancerSurvOS,file = 'cancerSurvOS.RData')

############### breast cancer OS Cox analysis
#' @GSE20711 age grade molecular_subtype HER2_status ER_status node_status tumor_size TMScore
#' @GSE48390 molecular_subtype HER2_status ER_status
#' @GSE20685 age stage_T stage_N stage_M
#' @GSE135565 age stage
#' @GSE103091 age
#' @GSE65216 age molecular_subtype post_menopausal_status node_status tumor_size height_cm weight_kg 
#' @GSE88770 grade 'HER2_status','ER_status', PgR_status Ki67_status node_status


breastSurv<-cancerSurvOS$breast[,c('cancer_type','GEO_number','geo_accession','OS_status','OS_month',
                                   'age','stage','grade','molecular_subtype','HER2_status','ER_status',
                                   'PgR_status','node_status','tumor_size','stage_T','stage_N','stage_M',
                                   'post_menopausal_status','TMScore')]
breastSurv$tumor_size[c(12:14,44,74,83)]<-c('1.7','2.8','2.6','1.2','2.1','4.5')
breastSurv$tumor_size<-as.numeric(breastSurv$tumor_size)
breastSurv$tumor_size[688:817]<-breastSurv$tumor_size[688:817]/10
breastSurv$tumor_size1<-breastSurv$tumor_size
breastSurv$tumor_size1[breastSurv$tumor_size1 <=2]<-'1'
breastSurv$tumor_size1<-as.numeric(breastSurv$tumor_size1)
breastSurv$tumor_size1[breastSurv$tumor_size1 >2]<-'2'

## age stage_T stage_N stage_M 327 
#' @ GSE20685
a<-breastSurv[,c(1:5,6,15:17,19)]
b<-na.omit(a)
Cox.function(b$OS_month,b$OS_status,b,clinical.variate=c(6:ncol(b)))
c<-cut(b$TMScore,quantile(b$TMScore,c(0,1/3,2/3,1)),labels = FALSE,include.lowest = TRUE)
b$TMScore<-c
b$TMScore<-as.character(b$TMScore)
d<-Cox.function(b$OS_month,b$OS_status,b,clinical.variate=c(6:ncol(b)))

Forestplot <- function(cox.result){
  require(forestplot)
  
  # 数据重组织
  cox.result <- subset(cox.result, variate %in% c("agediagnosis") |cox.result[, 2] != " ")
  rownames(cox.result) <- c(cox.result[1, 1], cox.result[-c(1), 2])
  cox.result <- cox.result[c('2 VS 1', '3 VS 1', 'agediagnosis','T2 VS T1', 'T3 VS T1',
                             'T4 VS T1','N1 VS N0','N2 VS N0','N3 VS N0','M1 VS M0'), 
                           c('multiv HR (95% CI for HR)',  'multiv p value')]
  
  # 构建需要输入的参数信息
  varibles <- c("Breast cancer OS varibles", 'TMScore','intermediate vs. low','high vs. low',
                "Age at dignosis", 'Stage T',"T2 vs. T1", "T3 vs. T1", "T4 vs. T1",'Stage N', 
                "N1 vs. N0", "N2 vs. N0", "N3 vs. N0",'Stage M','M1 vs. M0') 
  
  info.index <- c(3,4,5,7:9,11:13,15)
  # 用于绘图的HR信息
  HR <- rep(NA, length(varibles))
  multiCox.HR <- as.character(cox.result[, 1])
  HR[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(x," \\(")[[1]][1]))
  lowerCI <- rep(NA, length(varibles))
  lowerCI[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(strsplit(x," \\(")[[1]][2], "-")[[1]][1]))
  upperCI <- rep(NA, length(varibles))
  upperCI[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(strsplit(x,"-")[[1]][2], "\\)")[[1]][1]))
  # 文字信息
  HR.info <- rep(NA, length(varibles))
  HR.info[info.index] <- paste(round(HR[info.index],2)," (", round(lowerCI[info.index],2), "-", round(upperCI[info.index],2),")",sep="")
  multiCox.P <- as.character(cox.result[,2])
  P.value <- rep(NA, length(varibles))
  P.value[info.index] <- multiCox.P
  tabletext <- cbind(varibles, HR.info, P.value)
  colnames(tabletext) <- c("Varible", "HR(95%CI)", "P")
  
  # 绘制森林图
  
  
  forest.plot <- forestplot(tabletext[,1:3],
                            HR, lowerCI, upperCI, 
                            txt_gp = fpTxtGp(ticks = gpar(cex = 0.3), xlab = gpar(cex = 0.5), cex = 0.5),
                            graph.pos = 2, # 森林图为位置
                            zero = 1, # 设置无效线
                            boxsize = .1, # 规定方块的大小一致
                            col = fpColors(box = "royalblue",line = "black"), 
                            # vertices = TRUE, # 给线增加端点
                            clip = c(.1, 8), # 超出范围的区间增加箭头
                            xlab ="HR",
                            line.margin = unit(5,"mm"),
                            lineheight = unit(6,"mm"),
                            xticks = (c(0, 0.5, 1:8)),
                            graphwidth = unit(30,"mm"))
  
  return(forest.plot)
}

library(ggpubr)
ggsave(plot = Forestplot(d), filename='breastOS.pdf')

############### breast cancer DFS Cox analysis
cancerSurvDFS<-lapply(cancerSurv, function(x){
  CancerSurv(x,'DFS')
})

setwd('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/NewRoute/ProbeCalculate/survival/BreastLungSigThreeDatabase/survival/survAllRaw/Cox')
save(cancerSurvDFS,file = 'cancerSurvDFS.RData')

############### breast cancer DFS Cox analysis
#' @GSE61304 age stage stage_T stage_N grade subtyping 'ER_status', PgR_status
#' @GSE31448 age stage_T subtyping molecular_subtype 'HER2_status','ER_status', PgR_status Ki67_status P53_status SBR_grade node_status
#' @GSE21653 age stage_T stage_N subtyping 'HER2_status','ER_status', PgR_status Ki67_status P53_status SBR_grade

breastSurv<-cancerSurvDFS$breast[,c('cancer_type','GEO_number','geo_accession','DFS_status','DFS_month',
                                    'age','stage','stage_T','stage_N','grade','subtyping','molecular_subtype','HER2_status','ER_status',
                                    'PgR_status','Ki67_status','P53_status','SBR_grade','node_status','TMScore')]

breastSurv$subtyping[breastSurv$subtyping %in% c(' ductal','DUC','IDC')]<-'IDC'
breastSurv$subtyping[breastSurv$subtyping %in% c(' lobular','ILC')]<-'ILC'
breastSurv$subtyping[breastSurv$subtyping %in% c(' squamous cell','other','TUB','MIX','MUC')]<-NA
breastSurv$Ki67_status[breastSurv$Ki67_status %in% c('0')]<-'0'
breastSurv$Ki67_status[breastSurv$Ki67_status %in% c('1','pos')]<-'1'

## age stage_T subtyping ER_status PgR_status 478 
#' @： "GSE61304" "GSE31448" "GSE21653"
a<-breastSurv[,c(1:5,6,8,11,14,15,20)]
b<-na.omit(a)
Cox.function(b$DFS_month,b$DFS_status,b,clinical.variate=c(6:ncol(b)))
c<-cut(b$TMScore,quantile(b$TMScore,c(0,1/3,2/3,1)),labels = FALSE,include.lowest = TRUE)
b$TMScore<-c
b$TMScore<-as.character(b$TMScore)
d<-Cox.function(b$DFS_month,b$DFS_status,b,clinical.variate=c(6:ncol(b)))

Forestplot <- function(cox.result){
  require(forestplot)
  
  # 有重复名称，更改
  cox.result[3,2]<-'ER_status pos vs neg'
  cox.result[5,2]<-'PgR_status pos vs neg'
  
  # 数据重组织
  cox.result <- subset(cox.result, variate %in% c("agediagnosis") |cox.result[, 2] != " ")
  rownames(cox.result) <- c(cox.result[1, 1], cox.result[-c(1), 2])
  cox.result <- cox.result[c('2 VS 1', '3 VS 1', 'agediagnosis','T2 VS T1', 'T3 VS T1',
                             'ER_status pos vs neg','PgR_status pos vs neg','ILC VS IDC'),
                           c('multiv HR (95% CI for HR)','multiv p value')]
  # 构建需要输入的参数信息
  varibles <- c("Breast cancer DFS varibles", 'TMScore','intermediate vs. low','high vs. low',
                "Age at dignosis",'Stage T',"T2 vs. T1","T3 vs. T1", 'ER status','positive vs. negative',
                'PR status','positive vs. negative','subtyping','ILC vs. IDC') 
  
  info.index <- c(3,4,5,7,8,10,12,14)
  # 用于绘图的HR信息
  HR <- rep(NA, length(varibles))
  multiCox.HR <- as.character(cox.result[, 1])
  HR[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(x," \\(")[[1]][1]))
  lowerCI <- rep(NA, length(varibles))
  lowerCI[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(strsplit(x," \\(")[[1]][2], "-")[[1]][1]))
  upperCI <- rep(NA, length(varibles))
  upperCI[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(strsplit(x,"-")[[1]][2], "\\)")[[1]][1]))
  # 文字信息
  HR.info <- rep(NA, length(varibles))
  HR.info[info.index] <- paste(round(HR[info.index],2)," (", round(lowerCI[info.index],2), "-", round(upperCI[info.index],2),")",sep="")
  multiCox.P <- as.character(cox.result[,2])
  P.value <- rep(NA, length(varibles))
  P.value[info.index] <- multiCox.P
  tabletext <- cbind(varibles, HR.info, P.value)
  colnames(tabletext) <- c("Varible", "HR(95%CI)", "P")
  
  # 绘制森林图
  
  
  forest.plot <- forestplot(tabletext[,1:3],
                            HR, lowerCI, upperCI, 
                            txt_gp = fpTxtGp(ticks = gpar(cex = 0.3), xlab = gpar(cex = 0.5), cex = 0.5),
                            graph.pos = 2, # 森林图为位置
                            zero = 1, # 设置无效线
                            boxsize = .1, # 规定方块的大小一致
                            col = fpColors(box = "royalblue",line = "black"), 
                            # vertices = TRUE, # 给线增加端点
                            clip = c(.1, 8), # 超出范围的区间增加箭头
                            xlab ="HR",
                            line.margin = unit(5,"mm"),
                            lineheight = unit(6,"mm"),
                            xticks = (c(0, 0.5, 1:8)),
                            graphwidth = unit(30,"mm"))
  
  return(forest.plot)
}

library(ggpubr)
ggsave(plot = Forestplot(d), filename='breastDFS.pdf')


############### breast cancer RFS Cox analysis
cancerSurvRFS<-lapply(cancerSurv, function(x){
  CancerSurv(x,'RFS')
})

setwd('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/NewRoute/ProbeCalculate/survival/BreastLungSigThreeDatabase/survival/survAllRaw/Cox')
save(cancerSurvRFS,file = 'cancerSurvRFS.RData')

############### breast cancer RFS Cox analysis
#' @GSE16391 age grade HER2_status','ER_status', PgR_status post_menopausal_status node_status tumor_size
#' @GSE9195  age grade ER_status', PgR_status node_status tumor_size
#' @GSE20711 age grade molecular_subtype HER2_status ER_status node_status tumor_size TMScore
#' @GSE20685 age stage_T stage_N stage_M 
#' @GSE71258 HER2_status','ER_status', PgR_status dtc_IHC_positive post_menopausal_status
#' @GSE6532  age grade ER_status', PgR_status node_status tumor_size
#' @GSE65216 age molecular_subtype post_menopausal_status node_status tumor_size height_cm weight_kg 
#' @GSE88770 grade 'HER2_status','ER_status', PgR_status Ki67_status node_status


breastSurv<-cancerSurvRFS$breast[,c('cancer_type','GEO_number','geo_accession','RFS_status','RFS_month',
                                    'age','stage_T','stage_N','stage_M','grade','molecular_subtype',
                                    'HER2_status','ER_status','PgR_status','Ki67_status','dtc_IHC_positive',
                                    'post_menopausal_status','node_status','tumor_size','height_cm','weight_kg','TMScore')]
breastSurv$post_menopausal_status[breastSurv$post_menopausal_status %in% c('before_chemotherapy','Post','post')]<-'post'
breastSurv$post_menopausal_status[breastSurv$post_menopausal_status %in% c('after chemotherapy','Pre','pre')]<-'pre'

breastSurv$tumor_size[c(137:139,169,199,208)]<-c('1.7','2.8','2.6','1.2','2.1','4.5')
breastSurv$tumor_size<-as.numeric(breastSurv$tumor_size)
breastSurv$tumor_size[708:837]<-breastSurv$tumor_size[708:837]/10

breastSurv$tumor_size[49:213][breastSurv$tumor_size[49:213] <= 2]<-'1'
breastSurv$tumor_size<-as.numeric(breastSurv$tumor_size)
breastSurv$tumor_size[621:837][breastSurv$tumor_size[621:837] <= 2]<-'1'
breastSurv$tumor_size<-as.numeric(breastSurv$tumor_size)
breastSurv$tumor_size[49:213][breastSurv$tumor_size[49:213] > 2]<-'2'
breastSurv$tumor_size<-as.numeric(breastSurv$tumor_size)
breastSurv$tumor_size[621:837][breastSurv$tumor_size[621:837] > 2]<-'2'


## grade ER_status PgR_status node_status 288 
#' @ "GSE16391" "GSE9195"  "GSE6532"  "GSE88770"
a<-breastSurv[,c(1:5,10,13,14,18,22)]
b<-na.omit(a)
Cox.function(b$RFS_month,b$RFS_status,b,clinical.variate=c(6:ncol(b)))
c<-cut(b$TMScore,quantile(b$TMScore,c(0,1/3,2/3,1)),labels = FALSE,include.lowest = TRUE)
b$TMScore<-c
b$TMScore<-as.character(b$TMScore)
d<-Cox.function(b$RFS_month,b$RFS_status,b,clinical.variate=c(6:ncol(b)))

Forestplot <- function(cox.result){
  require(forestplot)
  
  # 有重复名称，更改
  cox.result[1,2]<-'ER_status pos vs neg'
  cox.result[3,2]<-'ER_status pos vs neg'
  cox.result[8,2]<-'node_status pos vs neg'
  cox.result[10,2]<-'PgR_status pos vs neg'
  cox.result<-cox.result[-c(3),]
  
  # 数据重组织
  cox.result <- subset(cox.result, cox.result[, 2] != " ")
  rownames(cox.result) <- cox.result[, 2]
  cox.result <- cox.result[c('2 VS 1', '3 VS 1', 'G2 VS G1','G3 VS G1',
                             'ER_status pos vs neg','PgR_status pos vs neg','node_status pos vs neg'),
                           c('multiv HR (95% CI for HR)','multiv p value')]
  # 构建需要输入的参数信息
  varibles <- c("Breast cancer RFS varibles", 'TMScore','intermediate vs. low','high vs. low',
                "Grade",'G2 vs. G1',"G3 vs. G1",'ER status','positive vs. negative',
                'PR status','positive vs. negative','Lymph node involvement','positive vs. negative') 
  
  info.index <- c(3,4,6,7,9,11,13)
  # 用于绘图的HR信息
  HR <- rep(NA, length(varibles))
  multiCox.HR <- as.character(cox.result[, 1])
  HR[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(x," \\(")[[1]][1]))
  lowerCI <- rep(NA, length(varibles))
  lowerCI[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(strsplit(x," \\(")[[1]][2], "-")[[1]][1]))
  upperCI <- rep(NA, length(varibles))
  upperCI[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(strsplit(x,"-")[[1]][2], "\\)")[[1]][1]))
  # 文字信息
  HR.info <- rep(NA, length(varibles))
  HR.info[info.index] <- paste(round(HR[info.index],2)," (", round(lowerCI[info.index],2), "-", round(upperCI[info.index],2),")",sep="")
  multiCox.P <- as.character(cox.result[,2])
  P.value <- rep(NA, length(varibles))
  P.value[info.index] <- multiCox.P
  tabletext <- cbind(varibles, HR.info, P.value)
  colnames(tabletext) <- c("Varible", "HR(95%CI)", "P")
  
  # 绘制森林图
  
  
  forest.plot <- forestplot(tabletext[,1:3],
                            HR, lowerCI, upperCI, 
                            txt_gp = fpTxtGp(ticks = gpar(cex = 0.3), xlab = gpar(cex = 0.5), cex = 0.5),
                            graph.pos = 2, # 森林图为位置
                            zero = 1, # 设置无效线
                            boxsize = .1, # 规定方块的大小一致
                            col = fpColors(box = "royalblue",line = "black"), 
                            # vertices = TRUE, # 给线增加端点
                            clip = c(.1, 8), # 超出范围的区间增加箭头
                            xlab ="HR",
                            line.margin = unit(5,"mm"),
                            lineheight = unit(6,"mm"),
                            xticks = (c(0, 0.5, 1:8)),
                            graphwidth = unit(30,"mm"))
  
  return(forest.plot)
}

library(ggpubr)
ggsave(plot = Forestplot(d), filename='breastRFS.pdf')


############### breast cancer DMFS Cox analysis
cancerSurvDMFS<-lapply(cancerSurv, function(x){
  CancerSurv(x,'DMFS')
})

setwd('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/NewRoute/ProbeCalculate/survival/BreastLungSigThreeDatabase/survival/survAllRaw/Cox')
save(cancerSurvDMFS,file = 'cancerSurvDMFS.RData')

############### breast cancer DMFS Cox analysis
#' @GSE61304 age stage stage_T stage_N grade subtyping 'ER_status', PgR_status
#' @GSE9195  age grade ER_status', PgR_status node_status tumor_size
#' @GSE58984 age molecular_subtype HER2_status','ER_status', PgR_status node_status tumor_size
#' @GSE20685 age stage_T stage_N stage_M
#' @GSE6532  age grade ER_status', PgR_status node_status tumor_size

breastSurv<-cancerSurvDMFS$breast[,c('cancer_type','GEO_number','geo_accession','DMFS_status','DMFS_month',
                                     'age','stage','stage_T','stage_N','stage_M','grade','subtyping','molecular_subtype',
                                     'HER2_status','ER_status','PgR_status','node_status','tumor_size','TMScore')]

breastSurv$tumor_size<-as.numeric(breastSurv$tumor_size)
breastSurv$tumor_size[135:228]<-breastSurv$tumor_size[135:228]/10
breastSurv$tumor_size1<-breastSurv$tumor_size
breastSurv$tumor_size1[breastSurv$tumor_size1 <= 2]<-'1'
breastSurv$tumor_size1<-as.numeric(breastSurv$tumor_size1)
breastSurv$tumor_size1[breastSurv$tumor_size1 > 2]<-'2'

## age stage_T stage_N 382 
#' "GSE61304" "GSE20685"
a<-breastSurv[,c(1:5,6,8,9,19)]
b<-na.omit(a)
Cox.function(b$DMFS_month,b$DMFS_status,b,clinical.variate=c(6:ncol(b)))
c<-cut(b$TMScore,quantile(b$TMScore,c(0,1/3,2/3,1)),labels = FALSE,include.lowest = TRUE)
b$TMScore<-c
b$TMScore<-as.character(b$TMScore)
d<-Cox.function(b$DMFS_month,b$DMFS_status,b,clinical.variate=c(6:ncol(b)))

Forestplot <- function(cox.result){
  require(forestplot)
  
  # 数据重组织
  cox.result <- subset(cox.result, variate %in% c("agediagnosis") |cox.result[, 2] != " ")
  rownames(cox.result) <- c(cox.result[1, 1], cox.result[-c(1), 2])
  cox.result <- cox.result[c('2 VS 1', '3 VS 1', 'agediagnosis','T2 VS T1', 'T3 VS T1',
                             'T4 VS T1','N1 VS N0','N2 VS N0','N3 VS N0'), 
                           c('multiv HR (95% CI for HR)',  'multiv p value')]
  
  # 构建需要输入的参数信息
  varibles <- c("Breast cancer DMFS varibles", 'TMScore','intermediate vs. low','high vs. low',
                "Age at dignosis", 'Stage T',"T2 vs. T1", "T3 vs. T1", "T4 vs. T1",'Stage N', 
                "N1 vs. N0", "N2 vs. N0", "N3 vs. N0") 
  
  info.index <- c(3,4,5,7:9,11:13)
  # 用于绘图的HR信息
  HR <- rep(NA, length(varibles))
  multiCox.HR <- as.character(cox.result[, 1])
  HR[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(x," \\(")[[1]][1]))
  lowerCI <- rep(NA, length(varibles))
  lowerCI[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(strsplit(x," \\(")[[1]][2], "-")[[1]][1]))
  upperCI <- rep(NA, length(varibles))
  upperCI[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(strsplit(x,"-")[[1]][2], "\\)")[[1]][1]))
  # 文字信息
  HR.info <- rep(NA, length(varibles))
  HR.info[info.index] <- paste(round(HR[info.index],2)," (", round(lowerCI[info.index],2), "-", round(upperCI[info.index],2),")",sep="")
  multiCox.P <- as.character(cox.result[,2])
  P.value <- rep(NA, length(varibles))
  P.value[info.index] <- multiCox.P
  tabletext <- cbind(varibles, HR.info, P.value)
  colnames(tabletext) <- c("Varible", "HR(95%CI)", "P")
  
  # 绘制森林图
  
  
  forest.plot <- forestplot(tabletext[,1:3],
                            HR, lowerCI, upperCI, 
                            txt_gp = fpTxtGp(ticks = gpar(cex = 0.3), xlab = gpar(cex = 0.5), cex = 0.5),
                            graph.pos = 2, # 森林图为位置
                            zero = 1, # 设置无效线
                            boxsize = .1, # 规定方块的大小一致
                            col = fpColors(box = "royalblue",line = "black"), 
                            # vertices = TRUE, # 给线增加端点
                            clip = c(.1, 8), # 超出范围的区间增加箭头
                            xlab ="HR",
                            line.margin = unit(5,"mm"),
                            lineheight = unit(6,"mm"),
                            xticks = (c(0, 0.5, 1:8)),
                            graphwidth = unit(30,"mm"))
  
  return(forest.plot)
}

library(ggpubr)
ggsave(plot = Forestplot(d), filename='breastDMFS.pdf')


################ lung cancer ################ 
############### lung cancer OS Cox analysis
#' @GSE29013 age gender stage subtyping1 TMScore
#' @GSE19188 gender subtyping1
#' @GSE37745 age gender stage subtyping1 who_performance_status
#' @GSE31210 age gender stage subtyping1 smoking_status ALK_fusion KRAS_mutation_status EGFR_mutation_status MYC MYC_copy
#' @GSE30219 age gender stage_T stage_N stage_M subtyping1
#' @GSE157010 age gender stage_T subtyping1


lungSurv<-cancerSurvOS$lung[,c('cancer_type','GEO_number','geo_accession','OS_status','OS_month',
                               'age','gender','stage','stage_T', 'stage_N', 'stage_M','subtyping1','smoking_status',
                               'ALK_fusion','KRAS_mutation_status','EGFR_mutation_status','MYC',
                               'MYC_copy','who_performance_status','TMScore')]
lungSurv$subtyping1[lungSurv$subtyping1 %in% c('Squamous','SCC','squamous','SQC','Squamous Cell Carcinoma of Lung')]<-'LUSC'
lungSurv$subtyping1[lungSurv$subtyping1 %in% c('Adenocarcinoma','ADC','adeno','LADC')]<-'LADC'
lungSurv$subtyping1[lungSurv$subtyping1 %in% c('LCC','large','LCNE','CARCI','BAS','Other')]<-NA
lungSurv$MYC[lungSurv$MYC %in% c('nd')]<-NA

## age gender stage subtyping1 smoking_status 281
#' "GSE29013" "GSE37745" "GSE31210"
a<-lungSurv[,c(1:5,6,7,8,12,13,20)]
b<-na.omit(a)
Cox.function(b$OS_month,b$OS_status,b,clinical.variate=c(6:ncol(b)))
c<-cut(b$TMScore,quantile(b$TMScore,c(0,1/3,2/3,1)),labels = FALSE,include.lowest = TRUE)
b$TMScore<-c
b$TMScore<-as.character(b$TMScore)
d<-Cox.function(b$OS_month,b$OS_status,b,clinical.variate=c(6:ncol(b)))

Forestplot <- function(cox.result){
  require(forestplot)
  
  # 数据重组织
  cox.result <- subset(cox.result, variate %in% c("agediagnosis") |cox.result[, 2] != " ")
  rownames(cox.result) <- c(cox.result[1, 1], cox.result[-c(1), 2])
  cox.result <- cox.result[c('2 VS 1', '3 VS 1', 'agediagnosis','M VS F','II VS I', 'III VS I',
                             '1 VS 0','LUSC VS LADC'), 
                           c('multiv HR (95% CI for HR)',  'multiv p value')]
  
  # 构建需要输入的参数信息
  varibles <- c("Lung cancer OS varibles", 'TMScore','intermediate vs. low','high vs. low',
                "Age at dignosis",'Gender','male vs. female', 'Stage',"II vs. I","III vs. I", 
                "Smoking", "yes vs. no",'Subtyping', "LUSC vs. LADC") 
  
  info.index <- c(3,4,5,7,9:10,12,14)
  # 用于绘图的HR信息
  HR <- rep(NA, length(varibles))
  multiCox.HR <- as.character(cox.result[, 1])
  HR[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(x," \\(")[[1]][1]))
  lowerCI <- rep(NA, length(varibles))
  lowerCI[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(strsplit(x," \\(")[[1]][2], "-")[[1]][1]))
  upperCI <- rep(NA, length(varibles))
  upperCI[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(strsplit(x,"-")[[1]][2], "\\)")[[1]][1]))
  # 文字信息
  HR.info <- rep(NA, length(varibles))
  HR.info[info.index] <- paste(round(HR[info.index],2)," (", round(lowerCI[info.index],2), "-", round(upperCI[info.index],2),")",sep="")
  multiCox.P <- as.character(cox.result[,2])
  P.value <- rep(NA, length(varibles))
  P.value[info.index] <- multiCox.P
  tabletext <- cbind(varibles, HR.info, P.value)
  colnames(tabletext) <- c("Varible", "HR(95%CI)", "P")
  
  # 绘制森林图
  
  
  forest.plot <- forestplot(tabletext[,1:3],
                            HR, lowerCI, upperCI, 
                            txt_gp = fpTxtGp(ticks = gpar(cex = 0.3), xlab = gpar(cex = 0.5), cex = 0.5),
                            graph.pos = 2, # 森林图为位置
                            zero = 1, # 设置无效线
                            boxsize = .1, # 规定方块的大小一致
                            col = fpColors(box = "royalblue",line = "black"), 
                            # vertices = TRUE, # 给线增加端点
                            clip = c(.1, 8), # 超出范围的区间增加箭头
                            xlab ="HR",
                            line.margin = unit(5,"mm"),
                            lineheight = unit(6,"mm"),
                            xticks = (c(0, 0.5, 1:8)),
                            graphwidth = unit(30,"mm"))
  
  return(forest.plot)
}

library(ggpubr)
ggsave(plot = Forestplot(d), filename='lungOS.pdf')


############### lung cancer RFS Cox analysis
#' @GSE37745 age gender stage subtyping1 who_performance_status
#' @GSE31210 age gender stage subtyping1 smoking_status ALK_fusion KRAS_mutation_status EGFR_mutation_status MYC MYC_copy
#' @GSE30219 age gender stage_T stage_N stage_M subtyping1

lungSurv<-cancerSurvRFS$lung[,c('cancer_type','GEO_number','geo_accession','RFS_status','RFS_month',
                                'age','gender','stage','stage_T', 'stage_N', 'stage_M','subtyping1','smoking_status',
                                'ALK_fusion','KRAS_mutation_status','EGFR_mutation_status','MYC',
                                'MYC_copy','who_performance_status','TMScore')]

lungSurv$subtyping1[lungSurv$subtyping1 %in% c('Squamous','SCC','squamous','SQC','Squamous Cell Carcinoma of Lung')]<-'LUSC'
lungSurv$subtyping1[lungSurv$subtyping1 %in% c('Adenocarcinoma','ADC','adeno','LADC')]<-'LADC'
lungSurv$subtyping1[lungSurv$subtyping1 %in% c('LCC','large','LCNE','CARCI','BAS','Other')]<-NA
lungSurv$MYC[lungSurv$MYC %in% c('nd')]<-NA


## age gender stage subtyping1 309 
#' "GSE37745" "GSE31210"
a<-lungSurv[,c(1:5,6:8,12,20)]
b<-na.omit(a)
Cox.function(b$RFS_month,b$RFS_status,b,clinical.variate=c(6:ncol(b)))
c<-cut(b$TMScore,quantile(b$TMScore,c(0,1/3,2/3,1)),labels = FALSE,include.lowest = TRUE)
b$TMScore<-c
b$TMScore<-as.character(b$TMScore)
d<-Cox.function(b$RFS_month,b$RFS_status,b,clinical.variate=c(6:ncol(b)))

Forestplot <- function(cox.result){
  require(forestplot)
  
  # 数据重组织
  cox.result <- subset(cox.result, variate %in% c("agediagnosis") |cox.result[, 2] != " ")
  rownames(cox.result) <- c(cox.result[1, 1], cox.result[-c(1), 2])
  cox.result <- cox.result[c('2 VS 1', '3 VS 1', 'agediagnosis','M VS F', 'II VS I',
                             'III VS I','IV VS I','LUSC VS LADC'),
                           c('multiv HR (95% CI for HR)','multiv p value')]
  # 构建需要输入的参数信息
  varibles <- c("Lung cancer RFS varibles", 'TMScore','intermediate vs. low','high vs. low',
                "Age at dignosis",'Gender',"male vs. female","Stage", 'II vs. I','III vs. I',
                'IV vs. I','Subtyping','LUSC VS LADC') 
  
  info.index <- c(3,4,5,7,9:11,13)
  # 用于绘图的HR信息
  HR <- rep(NA, length(varibles))
  multiCox.HR <- as.character(cox.result[, 1])
  HR[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(x," \\(")[[1]][1]))
  lowerCI <- rep(NA, length(varibles))
  lowerCI[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(strsplit(x," \\(")[[1]][2], "-")[[1]][1]))
  upperCI <- rep(NA, length(varibles))
  upperCI[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(strsplit(x,"-")[[1]][2], "\\)")[[1]][1]))
  # 文字信息
  HR.info <- rep(NA, length(varibles))
  HR.info[info.index] <- paste(round(HR[info.index],2)," (", round(lowerCI[info.index],2), "-", round(upperCI[info.index],2),")",sep="")
  multiCox.P <- as.character(cox.result[,2])
  P.value <- rep(NA, length(varibles))
  P.value[info.index] <- multiCox.P
  tabletext <- cbind(varibles, HR.info, P.value)
  colnames(tabletext) <- c("Varible", "HR(95%CI)", "P")
  
  # 绘制森林图
  
  
  forest.plot <- forestplot(tabletext[,1:3],
                            HR, lowerCI, upperCI, 
                            txt_gp = fpTxtGp(ticks = gpar(cex = 0.3), xlab = gpar(cex = 0.5), cex = 0.5),
                            graph.pos = 2, # 森林图为位置
                            zero = 1, # 设置无效线
                            boxsize = .1, # 规定方块的大小一致
                            col = fpColors(box = "royalblue",line = "black"), 
                            # vertices = TRUE, # 给线增加端点
                            clip = c(.1, 8), # 超出范围的区间增加箭头
                            xlab ="HR",
                            line.margin = unit(5,"mm"),
                            lineheight = unit(6,"mm"),
                            xticks = (c(0, 0.5, 1:8)),
                            graphwidth = unit(30,"mm"))
  
  return(forest.plot)
}

library(ggpubr)
ggsave(plot = Forestplot(d), filename='lungRFS.pdf')

## age gender stage smoking_status ALK_fusion KRAS_mutation_status EGFR_mutation_status MYC 224 LADC
#'"GSE31210"
a<-lungSurv[,c(1:5,6:8,13,14,15,16,17,20)]
b<-na.omit(a)
Cox.function(b$RFS_month,b$RFS_status,b,clinical.variate=c(6:ncol(b)))
c<-cut(b$TMScore,quantile(b$TMScore,c(0,1/3,2/3,1)),labels = FALSE,include.lowest = TRUE)
b$TMScore<-c
b$TMScore<-as.character(b$TMScore)
d<-Cox.function(b$RFS_month,b$RFS_status,b,clinical.variate=c(6:ncol(b)))

Forestplot <- function(cox.result){
  require(forestplot)
  
  # 有重复名称，更改
  cox.result[3,2]<-'ALK_fusion WT VS M'
  cox.result[5,2]<-'EGFR_mutation WT VS M'
  cox.result[9,2]<-'KRAS_mutation WT VS M'
  
  # 数据重组织
  cox.result <- subset(cox.result, variate %in% c("agediagnosis") |cox.result[, 2] != " ")
  rownames(cox.result) <- c(cox.result[1, 1], cox.result[-c(1), 2])
  cox.result <- cox.result[c('2 VS 1', '3 VS 1', 'agediagnosis','M VS F', 'II VS I','1 VS 0',
                             'EGFR_mutation WT VS M','KRAS_mutation WT VS M','ALK_fusion WT VS M',
                             'Low VS High'),
                           c('multiv HR (95% CI for HR)','multiv p value')]
  # 构建需要输入的参数信息
  varibles <- c("LADC RFS varibles", 'TMScore','intermediate vs. low','high vs. low',
                "Age at dignosis",'Gender',"male vs. female","Stage", 'II vs. I','Smoking','yes vs. no',
                'EGFR mutation','mutation vs. wild type',
                'KRAS mutation','wild type vs. mutation','ALK fusion','yes vs. no',
                'MYC expression','low vs. high') 
  
  info.index <- c(3,4,5,7,9,11,13,15,17,19)
  # 用于绘图的HR信息
  HR <- rep(NA, length(varibles))
  multiCox.HR <- as.character(cox.result[, 1])
  HR[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(x," \\(")[[1]][1]))
  lowerCI <- rep(NA, length(varibles))
  lowerCI[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(strsplit(x," \\(")[[1]][2], "-")[[1]][1]))
  upperCI <- rep(NA, length(varibles))
  upperCI[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(strsplit(x,"-")[[1]][2], "\\)")[[1]][1]))
  # 文字信息
  HR.info <- rep(NA, length(varibles))
  HR.info[info.index] <- paste(round(HR[info.index],2)," (", round(lowerCI[info.index],2), "-", round(upperCI[info.index],2),")",sep="")
  multiCox.P <- as.character(cox.result[,2])
  P.value <- rep(NA, length(varibles))
  P.value[info.index] <- multiCox.P
  tabletext <- cbind(varibles, HR.info, P.value)
  colnames(tabletext) <- c("Varible", "HR(95%CI)", "P")
  
  # 绘制森林图
  
  
  forest.plot <- forestplot(tabletext[,1:3],
                            HR, lowerCI, upperCI, 
                            txt_gp = fpTxtGp(ticks = gpar(cex = 0.3), xlab = gpar(cex = 0.5), cex = 0.5),
                            graph.pos = 2, # 森林图为位置
                            zero = 1, # 设置无效线
                            boxsize = .1, # 规定方块的大小一致
                            col = fpColors(box = "royalblue",line = "black"), 
                            # vertices = TRUE, # 给线增加端点
                            clip = c(.1, 8), # 超出范围的区间增加箭头
                            xlab ="HR",
                            line.margin = unit(5,"mm"),
                            lineheight = unit(6,"mm"),
                            xticks = (c(0, 0.5, 1:8)),
                            graphwidth = unit(30,"mm"))
  
  return(forest.plot)
}

library(ggpubr)
ggsave(plot = Forestplot(d), filename='LADCRFS.pdf')


################ colon cancer ################ 
############### colon cancer OS Cox analysis
#' @GSE29621 gender stage stage_T stage_N stage_M grade TMScore
#' @GSE17538 age gender stage grade
#' @GSE39582 age gender stage stage_T stage_N stage_M MMR_status KRAS_mutation_status P53_status BRAF_mutation CIN_status
#' @GSE72970 age gender stage_T stage_N
#' @GSE39084 age gender stage stage_T stage_N stage_M MSI_status KRAS_mutation_status P53_status BRAF_mutation PIK3CA_mutation_status CIMP_status lynch_syndrom

colonSurv<-cancerSurvOS$colon[,c('cancer_type','GEO_number','geo_accession','OS_status','OS_month',
                                 'age','gender','stage','stage_T', 'stage_N', 'stage_M','grade','MMR_status',
                                 'MSI_status','KRAS_mutation_status','P53_status','BRAF_mutation',
                                 'PIK3CA_mutation_status','CIN_status','CIMP_status','lynch_syndrom','TMScore')]

## age gender stage_T KRAS_mutation_status 585 
#'  "GSE39582" "GSE39084"
a<-colonSurv[,c(1:5,6,7,9,15,22)]
b<-na.omit(a)
Cox.function(b$OS_month,b$OS_status,b,clinical.variate=c(6:ncol(b)))
c<-cut(b$TMScore,quantile(b$TMScore,c(0,1/3,2/3,1)),labels = FALSE,include.lowest = TRUE)
b$TMScore<-c
b$TMScore<-as.character(b$TMScore)
d<-Cox.function(b$OS_month,b$OS_status,b,clinical.variate=c(6:ncol(b)))

Forestplot <- function(cox.result){
  require(forestplot)
  
  # 数据重组织
  cox.result <- subset(cox.result, variate %in% c("agediagnosis") |cox.result[, 2] != " ")
  rownames(cox.result) <- c(cox.result[1, 1], cox.result[-c(1), 2])
  cox.result <- cox.result[c('2 VS 1', '3 VS 1', 'agediagnosis','M VS F','T2 VS T1',
                             'T3 VS T1','T4 VS T1','WT VS M'), 
                           c('multiv HR (95% CI for HR)',  'multiv p value')]
  
  # 构建需要输入的参数信息
  varibles <- c("Colon cancer OS varibles", 'TMScore','intermediate vs. low','high vs. low',
                "Age at dignosis",'Gender','male vs. female', 'Stage T',"T2 vs. T1", 
                "T3 vs. T1", "T4 vs. T1",'KRAS mutation', "wild type vs. mutation") 
  
  info.index <- c(3,4,5,7,9:11,13)
  # 用于绘图的HR信息
  HR <- rep(NA, length(varibles))
  multiCox.HR <- as.character(cox.result[, 1])
  HR[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(x," \\(")[[1]][1]))
  lowerCI <- rep(NA, length(varibles))
  lowerCI[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(strsplit(x," \\(")[[1]][2], "-")[[1]][1]))
  upperCI <- rep(NA, length(varibles))
  upperCI[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(strsplit(x,"-")[[1]][2], "\\)")[[1]][1]))
  # 文字信息
  HR.info <- rep(NA, length(varibles))
  HR.info[info.index] <- paste(round(HR[info.index],2)," (", round(lowerCI[info.index],2), "-", round(upperCI[info.index],2),")",sep="")
  multiCox.P <- as.character(cox.result[,2])
  P.value <- rep(NA, length(varibles))
  P.value[info.index] <- multiCox.P
  tabletext <- cbind(varibles, HR.info, P.value)
  colnames(tabletext) <- c("Varible", "HR(95%CI)", "P")
  
  # 绘制森林图
  
  
  forest.plot <- forestplot(tabletext[,1:3],
                            HR, lowerCI, upperCI, 
                            txt_gp = fpTxtGp(ticks = gpar(cex = 0.3), xlab = gpar(cex = 0.5), cex = 0.5),
                            graph.pos = 2, # 森林图为位置
                            zero = 1, # 设置无效线
                            boxsize = .1, # 规定方块的大小一致
                            col = fpColors(box = "royalblue",line = "black"), 
                            # vertices = TRUE, # 给线增加端点
                            clip = c(.1, 8), # 超出范围的区间增加箭头
                            xlab ="HR",
                            line.margin = unit(5,"mm"),
                            lineheight = unit(6,"mm"),
                            xticks = (c(0, 0.5, 1:8)),
                            graphwidth = unit(30,"mm"))
  
  return(forest.plot)
}

library(ggpubr)
ggsave(plot = Forestplot(d), filename='colonOS.pdf')

############### colon cancer RFS Cox analysis
#' @GSE39582 age gender stage stage_T stage_N stage_M MMR_status KRAS_mutation_status P53_status BRAF_mutation CIN_status

colonSurv<-cancerSurvRFS$colon[,c('cancer_type','GEO_number','geo_accession','RFS_status','RFS_month',
                                  'age','gender','stage','stage_T', 'stage_N', 'stage_M','MMR_status',
                                  'KRAS_mutation_status','P53_status','BRAF_mutation','CIN_status','TMScore')]

## age gender stage MMR_status KRAS_mutation_status P53_status BRAF_mutation 303 
#'"GSE39582"
a<-colonSurv[,c(1:5,6:8,12,13,14,15,17)]
b<-na.omit(a)
Cox.function(b$RFS_month,b$RFS_status,b,clinical.variate=c(6:ncol(b)))
c<-cut(b$TMScore,quantile(b$TMScore,c(0,1/3,2/3,1)),labels = FALSE,include.lowest = TRUE)
b$TMScore<-c
b$TMScore<-as.character(b$TMScore)
d<-Cox.function(b$RFS_month,b$RFS_status,b,clinical.variate=c(6:ncol(b)))

Forestplot <- function(cox.result){
  require(forestplot)
  
  # 有重复名称，更改
  cox.result[3,2]<-'BRAF_mutation WT VS M'
  cox.result[7,2]<-'KRAS_mutation WT VS M'
  cox.result[11,2]<-'P53_status WT VS M'
  
  # 数据重组织
  cox.result <- subset(cox.result, variate %in% c("agediagnosis") |cox.result[, 2] != " ")
  rownames(cox.result) <- c(cox.result[1, 1], cox.result[-c(1), 2])
  cox.result <- cox.result[c('2 VS 1', '3 VS 1', 'agediagnosis','M VS F', 'II VS I',
                             'III VS I','IV VS I','KRAS_mutation WT VS M','P53_status WT VS M',
                             'BRAF_mutation WT VS M'),
                           c('multiv HR (95% CI for HR)','multiv p value')]
  # 构建需要输入的参数信息
  varibles <- c("Colon cancer RFS varibles", 'TMScore','intermediate vs. low','high vs. low',
                "Age at dignosis",'Gender',"male vs. female","Stage", 'II vs. I','III vs. I',
                'IV vs. I','KRAS mutation','wild type vs. mutation','P53 status','wild type vs. mutation',
                'BRAF status','wild type vs. mutation','MMR status','pMMR vs. dMMR') 
  
  info.index <- c(3,4,5,7,9:11,13,15,17,19)
  # 用于绘图的HR信息
  HR <- rep(NA, length(varibles))
  multiCox.HR <- as.character(cox.result[, 1])
  HR[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(x," \\(")[[1]][1]))
  lowerCI <- rep(NA, length(varibles))
  lowerCI[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(strsplit(x," \\(")[[1]][2], "-")[[1]][1]))
  upperCI <- rep(NA, length(varibles))
  upperCI[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(strsplit(x,"-")[[1]][2], "\\)")[[1]][1]))
  # 文字信息
  HR.info <- rep(NA, length(varibles))
  HR.info[info.index] <- paste(round(HR[info.index],2)," (", round(lowerCI[info.index],2), "-", round(upperCI[info.index],2),")",sep="")
  multiCox.P <- as.character(cox.result[,2])
  P.value <- rep(NA, length(varibles))
  P.value[info.index] <- multiCox.P
  tabletext <- cbind(varibles, HR.info, P.value)
  colnames(tabletext) <- c("Varible", "HR(95%CI)", "P")
  
  # 绘制森林图
  
  
  forest.plot <- forestplot(tabletext[,1:3],
                            HR, lowerCI, upperCI, 
                            txt_gp = fpTxtGp(ticks = gpar(cex = 0.3), xlab = gpar(cex = 0.5), cex = 0.5),
                            graph.pos = 2, # 森林图为位置
                            zero = 1, # 设置无效线
                            boxsize = .1, # 规定方块的大小一致
                            col = fpColors(box = "royalblue",line = "black"), 
                            # vertices = TRUE, # 给线增加端点
                            clip = c(.1, 8), # 超出范围的区间增加箭头
                            xlab ="HR",
                            line.margin = unit(5,"mm"),
                            lineheight = unit(6,"mm"),
                            xticks = (c(0, 0.5, 1:8)),
                            graphwidth = unit(30,"mm"))
  
  return(forest.plot)
}

library(ggpubr)
ggsave(plot = Forestplot(d), filename='colonRFS.pdf')



################ stomach cancer ################ 
############### stomach cancer OS Cox analysis
#' @GSE57303 age gender stage stage_T stage_N stage_M grade subtyping TMScore
#' @GSE62254 age gender stage stage_T stage_N stage_M subtyping MSI_status EBV_status P53_status MLH1_status
#' @GSE15460 age gender stage subtyping

stomachSurv<-cancerSurvOS$stomach[,c('cancer_type','GEO_number','geo_accession','OS_status','OS_month',
                                     'age','gender','stage','stage_T', 'stage_N', 'stage_M','grade',
                                     'subtyping','MSI_status','EBV_status','P53_status','MLH1_status','TMScore')]

stomachSurv$subtyping[stomachSurv$subtyping %in% c('diffuse type','diffuse','Diffuse')]<-'diffuse'
stomachSurv$subtyping[stomachSurv$subtyping %in% c('intestinal type','intestinal','Intestinal')]<-'intestinal'
stomachSurv$subtyping[stomachSurv$subtyping %in% c('mixed type','mixed','Mixed')]<-'mixed'
stomachSurv$subtyping[stomachSurv$subtyping %in% c('unknown')]<-NA
# stomachSurv$subtyping[stomachSurv$subtyping %in% c('mixed')]<-NA


stomachSurv$MLH1_status[stomachSurv$MLH1_status %in% c('negative','Negative')]<-'neg'
stomachSurv$MLH1_status[stomachSurv$MLH1_status %in% c('positive','Positive','positive; MSH2 mutation (+)')]<-'pos'
stomachSurv$MLH1_status[stomachSurv$MLH1_status %in% c('partial loss')]<-NA

## age gender stage stage_T stage_N stage_M subtyping 368 
#'  "GSE57303" "GSE62254"
a<-stomachSurv[,c(1:5,6,7,8,9:11,13,18)]
b<-na.omit(a)
Cox.function(b$OS_month,b$OS_status,b,clinical.variate=c(6:ncol(b)))
c<-cut(b$TMScore,quantile(b$TMScore,c(0,1/3,2/3,1)),labels = FALSE,include.lowest = TRUE)
b$TMScore<-c
b$TMScore<-as.character(b$TMScore)
d<-Cox.function(b$OS_month,b$OS_status,b,clinical.variate=c(6:ncol(b)))

Forestplot <- function(cox.result){
  require(forestplot)
  
  # 数据重组织
  cox.result <- subset(cox.result, variate %in% c("agediagnosis") |cox.result[, 2] != " ")
  rownames(cox.result) <- c(cox.result[1, 1], cox.result[-c(1), 2])
  cox.result <- cox.result[c('2 VS 1', '3 VS 1', 'agediagnosis','M VS F','II VS I',
                             'III VS I','IV VS I','T3 VS T2', 'T4 VS T2','N1 VS N0',
                             'N2 VS N0','N3 VS N0','M1 VS M0','intestinal VS diffuse',
                             'mixed VS diffuse'),c('multiv HR (95% CI for HR)','multiv p value')]
  # 构建需要输入的参数信息
  varibles <- c("Stomach cancer OS varibles", 'TMScore','intermediate vs. low','high vs. low',
                "Age at dignosis",'Gender','male vs. female', 'Stage','II vs. I','III vs. I',
                'IV vs. I','Stage T',"T3 vs. T2","T4 vs. T2", 'Stage N','N1 vs. N0','N2 vs. N0',
                'N3 vs. N0','Stage M','M1 vs. M0','subtyping','intestinal vs. diffuse','mixed vs. diffuse') 
  
  info.index <- c(3,4,5,7,9:11,13,14,16:18,20,22,23)
  # 用于绘图的HR信息
  HR <- rep(NA, length(varibles))
  multiCox.HR <- as.character(cox.result[, 1])
  HR[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(x," \\(")[[1]][1]))
  lowerCI <- rep(NA, length(varibles))
  lowerCI[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(strsplit(x," \\(")[[1]][2], "-")[[1]][1]))
  upperCI <- rep(NA, length(varibles))
  upperCI[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(strsplit(x,"-")[[1]][2], "\\)")[[1]][1]))
  # 文字信息
  HR.info <- rep(NA, length(varibles))
  HR.info[info.index] <- paste(round(HR[info.index],2)," (", round(lowerCI[info.index],2), "-", round(upperCI[info.index],2),")",sep="")
  multiCox.P <- as.character(cox.result[,2])
  P.value <- rep(NA, length(varibles))
  P.value[info.index] <- multiCox.P
  tabletext <- cbind(varibles, HR.info, P.value)
  colnames(tabletext) <- c("Varible", "HR(95%CI)", "P")
  
  # 绘制森林图
  
  
  forest.plot <- forestplot(tabletext[,1:3],
                            HR, lowerCI, upperCI, 
                            txt_gp = fpTxtGp(ticks = gpar(cex = 0.3), xlab = gpar(cex = 0.5), cex = 0.5),
                            graph.pos = 2, # 森林图为位置
                            zero = 1, # 设置无效线
                            boxsize = .1, # 规定方块的大小一致
                            col = fpColors(box = "royalblue",line = "black"), 
                            # vertices = TRUE, # 给线增加端点
                            clip = c(.1, 8), # 超出范围的区间增加箭头
                            xlab ="HR",
                            line.margin = unit(5,"mm"),
                            lineheight = unit(6,"mm"),
                            xticks = (c(0, 0.5, 1:8)),
                            graphwidth = unit(30,"mm"))
  
  return(forest.plot)
}

library(ggpubr)
ggsave(plot = Forestplot(d), filename='stomachOS.pdf')
