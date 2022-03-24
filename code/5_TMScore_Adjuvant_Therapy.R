survivalSigBreastLung<-read.table('/data/survivalSigBreastLung.txt',
                                  row.names = 1,header = TRUE,sep = '\t')
load('/data/clinicalSurvival.RData')
load('/data/expTMMsurv.RData')

clinicalSurvivalUse<-clinicalSurvivalUse[c(colnames(expTMMsurv$lung),colnames(expTMMsurv$breast),
                                           colnames(expTMMsurv$colon),colnames(expTMMsurv$kidney),
                                           colnames(expTMMsurv$stomach)),]

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
TMScore<-TMSprobeWeightedSum(survivalSigBreastLung,expTMMsurv)


SurvDataThreeGroupAll<-function(clinical.data,surv.type,sample.score){
  clisur<-clinical.data[clinical.data$tissue_type %in% c('tumor'),]
  # clisur<-clisur[clisur$GEO_number %in% c('GSE29013','GSE37745'),]
  # clisur<-clisur[clisur$adjuvant_chemo %in% c('Yes','yes'),]
  # clisur<-clisur[clisur$adjuvant_chemo %in% c('No','no'),]
  # clisur<-clisur[clisur$adjuvant_chemo %in% c('No','no','Yes','yes'),]
  # clisur<-clisur[clisur$GEO_number %in% c('GSE9195','GSE58984','GSE16391'),]
  # clisur<-clisur[clisur$GEO_number %in% c('GSE20685','GSE65216'),]
  # clisur<-clisur[clisur$treatment %in% c('yes'),]
  # data1<-clisur[clisur$GEO_number %in% c('GSE103091'),]
  # data1<-data1[data1$adjuvant_chemo %in% c('1'),]
  # clisur<-rbind(data1,clisur)
  # print(dim(clisur))
  # clisur<-clisur[clisur$GEO_number %in% c('GSE92921','GSE17538','GSE33113','GSE38832','GSE39084'),]
  # clisur<-clisur[clisur$GEO_number %in% c('GSE39582'),]
  # clisur<-clisur[clisur$adjuvant_chemo %in% c('N'),]
  # clisur<-clisur[clisur$adjuvant_chemo %in% c('Y'),]
  # clisur<-clisur[clisur$GEO_number %in% c('GSE14333'),]
  # clisur<-clisur[clisur$adjuvant_chemo %in% c('N'),]
  # clisur<-clisur[clisur$adjuvant_chemo %in% c('Y'),]
  # clisur<-clisur[clisur$Radiotherapy %in% c('N'),]
  # clisur<-clisur[clisur$GEO_number %in% c('GSE29621','GSE39582','GSE72970','GSE143985'),]
  # clisur<-clisur[clisur$adjuvant_chemo %in% c('N'),]
  # clisur<-clisur[clisur$adjuvant_chemo %in% c('Y','1'),]
  # print(unique(clisur$adjuvant_chemo))
  # clisur<-clisur[clisur$GEO_number %in% c('GSE29621','GSE39582','GSE143985'),]
  # 
  # clisur<-clisur[clisur$adjuvant_chemo %in% c('Y','1','N'),]
  # clisur$adjuvant_chemo[clisur$adjuvant_chemo %in% c('Y','1')]<-'Y'
  # clisur$adjuvant_chemo[clisur$adjuvant_chemo %in% c('N')]<-'N'
  
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

source(file = '/code/km_survival_plot.R')

pdf(file = paste0('OS','lungadjtherapy','.pdf'))
plot.surv(surcliOSTertilesWeight$lung$lung, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'lungadjtherapy', ylab = "OS")
dev.off()

pdf(file = paste0('RFS','lungadjtherapy','.pdf'))
plot.surv(surcliRFSTertilesWeight$lung$lung, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'lungadjtherapy', ylab = "RFS")
dev.off()

pdf(file = paste0('PFS','lungadjtherapy','.pdf'))
plot.surv(surcliPFSTertilesWeight$lung$lung, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'lungadjtherapy', ylab = "PFS")
dev.off()

pdf(file = paste0('OS','lungNoadjtherapy','.pdf'))
plot.surv(surcliOSTertilesWeight$lung$lung, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'lungNoadjtherapy', ylab = "OS")
dev.off()

pdf(file = paste0('RFS','lungNoadjtherapy','.pdf'))
plot.surv(surcliRFSTertilesWeight$lung$lung, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'lungNoadjtherapy', ylab = "OS")
dev.off()


pdf(file = paste0('OS','colonGSE39582','.pdf')) #
plot.surv(surcliOSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonGSE39582', ylab = "OS")
dev.off()

pdf(file = paste0('RFS','colonGSE39582','.pdf')) #
plot.surv(surcliRFSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonGSE39582', ylab = "RFS")
dev.off()

pdf(file = paste0('OS','colonNoTherGSE39582','.pdf')) 
plot.surv(surcliOSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonNoTherGSE39582', ylab = "OS")
dev.off()

pdf(file = paste0('RFS','colonNoTherGSE39582','.pdf')) 
plot.surv(surcliRFSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonNoTherGSE39582', ylab = "RFS")
dev.off()


pdf(file = paste0('OS','colonTherGSE39582','.pdf')) 
plot.surv(surcliOSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonTherGSE39582', ylab = "OS")
dev.off()

pdf(file = paste0('RFS','colonTherGSE39582','.pdf')) 
plot.surv(surcliRFSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonTherGSE39582', ylab = "RFS")
dev.off()


pdf(file = paste0('OS','colonchemodataall','.pdf')) 
plot.surv(surcliOSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonchemodataall', ylab = "OS")
dev.off()


pdf(file = paste0('RFS','colonchemodataall','.pdf')) 
plot.surv(surcliRFSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonchemodataall', ylab = "RFS")
dev.off()



pdf(file = paste0('OS','colonchemodataN','.pdf')) 
plot.surv(surcliOSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonchemodataN', ylab = "OS")
dev.off()

pdf(file = paste0('RFS','colonchemodataN','.pdf')) 
plot.surv(surcliRFSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonchemodataN', ylab = "RFS")
dev.off()




pdf(file = paste0('OS','colonchemodataY','.pdf')) 
plot.surv(surcliOSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonchemodataY', ylab = "OS")
dev.off()


pdf(file = paste0('RFS','colonchemodataY','.pdf')) 
plot.surv(surcliRFSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonchemodataY', ylab = "RFS")
dev.off()

