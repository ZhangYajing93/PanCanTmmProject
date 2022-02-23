
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
  clisur<-clinical.data[clinical.data$tissue_type %in% c('tumor'),]
  # clisur<-clisur[clisur$molecular_subtype %in% c('Basal'),]
  # clisur<-clisur[clisur$subtyping1 %in% c('LADC'),]
  # clisur<-clisur[clisur$subtyping1 %in% c('LUSC'),]
  # clisur<-clisur[clisur$stage %in% c('I','II'),]
  # clisur<-clisur[clisur$stage %in% c('III','IV'),]
  # clisur<-clisur[clisur$ER_status %in% c('neg'),]
  clisur<-clisur[clisur$GEO_number %in% c('GSE29013','GSE37745'),]
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
  # clidata<-switch(surv.type,
  #                 OS=subset(clisur,!(is.na(OS_status)|is.na(OS_month)))[,c('cancer_type','GEO_number','geo_accession','OS_status','OS_month','adjuvant_chemo')],
  #                 PFS=subset(clisur,!(is.na(PFS_status)|is.na(PFS_month)))[,c('cancer_type','GEO_number','geo_accession','PFS_status','PFS_month','adjuvant_chemo')],
  #                 RFS=subset(clisur,!(is.na(RFS_status)|is.na(RFS_month)))[,c('cancer_type','GEO_number','geo_accession','RFS_status','RFS_month','adjuvant_chemo')],
  #                 DFS=subset(clisur,!(is.na(DFS_status)|is.na(DFS_month)))[,c('cancer_type','GEO_number','geo_accession','DFS_status','DFS_month','adjuvant_chemo')],
  #                 DMFS=subset(clisur,!(is.na(DMFS_status)|is.na(DMFS_month)))[,c('cancer_type','GEO_number','geo_accession','DMFS_status','DMFS_month','adjuvant_chemo')],
  #                 MFS=subset(clisur,!(is.na(MFS_status)|is.na(MFS_month)))[,c('cancer_type','GEO_number','geo_accession','MFS_status','MFS_month','adjuvant_chemo')],
  #                 DSS=subset(clisur,!(is.na(DSS_status)|is.na(DSS_month)))[,c('cancer_type','GEO_number','geo_accession','DSS_status','DSS_month','adjuvant_chemo')])
  # colnames(clidata)<-c('cancer_type','GEO_number','Patient_ID','event','time','sample.label')
  
  
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

setwd('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/NewRoute/ProbeCalculate/survival/BreastLungSigThreeDatabase/survival/survAllRaw/subtype')

pdf(file = paste0('OS','BreastBasal','.pdf'))
plot.surv(surcliOSTertilesWeight$breast$breast, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'Breast Cancer Basal', ylab = "OS")
dev.off()

pdf(file = paste0('RFS','BreastBasal','.pdf'))
plot.surv(surcliRFSTertilesWeight$breast$breast, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'Breast Cancer Basal', ylab = "RFS")
dev.off()

pdf(file = paste0('DFS','BreastBasal','.pdf'))
plot.surv(surcliDFSTertilesWeight$breast$breast, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'Breast Cancer Basal', ylab = "DFS")
dev.off()

pdf(file = paste0('MFS','BreastBasal','.pdf'))
plot.surv(surcliMFSTertilesWeight$breast$breast, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'Breast Cancer Basal', ylab = "MFS")
dev.off()

pdf(file = paste0('OS','LADC','.pdf')) 
plot.surv(surcliOSTertilesWeight$lung$lung, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'LADC', ylab = "OS")
dev.off()

pdf(file = paste0('RFS','LADC','.pdf'))
plot.surv(surcliRFSTertilesWeight$lung$lung, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'LADC', ylab = "RFS")
dev.off()

pdf(file = paste0('PFS','LADC','.pdf'))
plot.surv(surcliPFSTertilesWeight$lung$lung, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'LADC', ylab = "PFS")
dev.off()

pdf(file = paste0('OS','LUSC','.pdf')) 
plot.surv(surcliOSTertilesWeight$lung$lung, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'LUSC', ylab = "OS")
dev.off()

pdf(file = paste0('RFS','LUSC','.pdf'))
plot.surv(surcliRFSTertilesWeight$lung$lung, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'LUSC', ylab = "RFS")
dev.off()

pdf(file = paste0('PFS','LUSC','.pdf'))
plot.surv(surcliPFSTertilesWeight$lung$lung, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'LUSC', ylab = "PFS")
dev.off()


setwd('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/NewRoute/ProbeCalculate/survival/BreastLungSigThreeDatabase/survival/survAllRaw/subtype/stage')
for (i in 1:length(surcliOSTertilesWeight)) {
  pdf(file = paste0('OS_stageIIIandIV',names(surcliOSTertilesWeight)[i],'.pdf'))
  lapply(1:length(surcliOSTertilesWeight[[i]]), function(x){
    print(plot.surv(surcliOSTertilesWeight[[i]][[x]], upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
                    surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
                    conf.int = TRUE, main = paste0(names(surcliOSTertilesWeight)[i],names(surcliOSTertilesWeight[[i]])[[x]]), ylab = "OS"))
  })
  dev.off()
}

for (i in 1:length(surcliPFSTertilesWeight)) {
  pdf(file = paste0('PFS_stageIIIandIV',names(surcliPFSTertilesWeight)[i],'.pdf'))
  lapply(1:length(surcliPFSTertilesWeight[[i]]), function(x){
    print(plot.surv(surcliPFSTertilesWeight[[i]][[x]], upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
                    surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
                    conf.int = TRUE, main = paste0(names(surcliPFSTertilesWeight)[i],names(surcliPFSTertilesWeight[[i]])[[x]]), ylab = "PFS"))
  })
  dev.off()
}

for (i in 1:length(surcliRFSTertilesWeight)) {
  pdf(file = paste0('RFS_stageIIIandIV',names(surcliRFSTertilesWeight)[i],'.pdf'))
  lapply(1:length(surcliRFSTertilesWeight[[i]]), function(x){
    print(plot.surv(surcliRFSTertilesWeight[[i]][[x]], upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
                    surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
                    conf.int = TRUE, main = paste0(names(surcliRFSTertilesWeight)[i],names(surcliRFSTertilesWeight[[i]])[[x]]), ylab = "RFS"))
  })
  dev.off()
}

for (i in 1:length(surcliDFSTertilesWeight)) {
  pdf(file = paste0('DFS_stageIIIandIV',names(surcliDFSTertilesWeight)[i],'.pdf'))
  lapply(1:length(surcliDFSTertilesWeight[[i]]), function(x){
    print(plot.surv(surcliDFSTertilesWeight[[i]][[x]], upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
                    surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
                    conf.int = TRUE, main = paste0(names(surcliDFSTertilesWeight)[i],names(surcliDFSTertilesWeight[[i]])[[x]]), ylab = "DFS"))
  })
  dev.off()
}

for (i in 1:length(surcliDMFSTertilesWeight)) {
  pdf(file = paste0('DMFS_stageIIIandIV',names(surcliDMFSTertilesWeight)[i],'.pdf'))
  lapply(1:length(surcliDMFSTertilesWeight[[i]]), function(x){
    print(plot.surv(surcliDMFSTertilesWeight[[i]][[x]], upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
                    surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
                    conf.int = TRUE, main = paste0(names(surcliDMFSTertilesWeight)[i],names(surcliDMFSTertilesWeight[[i]])[[x]]), ylab = "DMFS"))
  })
  dev.off()
}


setwd('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/NewRoute/ProbeCalculate/survival/BreastLungSigThreeDatabase/survival/survAllRaw/subtype')

pdf(file = paste0('OS','BreastER-','.pdf'))
plot.surv(surcliOSTertilesWeight$breast$breast, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'BreastER-', ylab = "OS")
dev.off()

pdf(file = paste0('RFS','BreastER-','.pdf'))
plot.surv(surcliRFSTertilesWeight$breast$breast, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'BreastER-', ylab = "RFS")
dev.off()

pdf(file = paste0('DFS','BreastER-','.pdf'))
plot.surv(surcliDFSTertilesWeight$breast$breast, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'BreastER-', ylab = "DFS")
dev.off()

pdf(file = paste0('MFS','BreastER-','.pdf'))
plot.surv(surcliMFSTertilesWeight$breast$breast, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'BreastER-', ylab = "MFS")
dev.off()

pdf(file = paste0('DMFS','BreastER-','.pdf'))
plot.surv(surcliDMFSTertilesWeight$breast$breast, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'BreastER-', ylab = "DMFS")
dev.off()

setwd('/pub5/xiaoyun/Jobs/J22/ZYJ2020/Telomere/NewRoute/ProbeCalculate/survival/BreastLungSigThreeDatabase/survival/survAllRaw/therapy')
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


lungtherapy<-clinicalSurvivalUse[clinicalSurvivalUse$GEO_number %in% c('GSE29013','GSE37745'),]
lungtherapy<-cbind(lungtherapy,TMScore=TMScore$lung[rownames(lungtherapy)])
Clidata<-function(clisur,surv.type){
  switch(surv.type,
         OS=subset(clisur,!(is.na(OS_status)|is.na(OS_month)))[,c('cancer_type','GEO_number','geo_accession','OS_status','OS_month','TMScore','adjuvant_chemo')],
         PFS=subset(clisur,!(is.na(PFS_status)|is.na(PFS_month)))[,c('cancer_type','GEO_number','geo_accession','PFS_status','PFS_month','TMScore','adjuvant_chemo')],
         RFS=subset(clisur,!(is.na(RFS_status)|is.na(RFS_month)))[,c('cancer_type','GEO_number','geo_accession','RFS_status','RFS_month','TMScore','adjuvant_chemo')],
         DFS=subset(clisur,!(is.na(DFS_status)|is.na(DFS_month)))[,c('cancer_type','GEO_number','geo_accession','DFS_status','DFS_month','TMScore','adjuvant_chemo')],
         DMFS=subset(clisur,!(is.na(DMFS_status)|is.na(DMFS_month)))[,c('cancer_type','GEO_number','geo_accession','DMFS_status','DMFS_month','TMScore','adjuvant_chemo')],
         MFS=subset(clisur,!(is.na(MFS_status)|is.na(MFS_month)))[,c('cancer_type','GEO_number','geo_accession','MFS_status','MFS_month','TMScore','adjuvant_chemo')],
         DSS=subset(clisur,!(is.na(DSS_status)|is.na(DSS_month)))[,c('cancer_type','GEO_number','geo_accession','DSS_status','DSS_month','TMScore','adjuvant_chemo')])
  
}
clidata<-Clidata(lungtherapy,'OS')
clidata<-clidata[clidata$adjuvant_chemo %in% c('No','no','Yes','yes'),]
clidata$adjuvant_chemo[clidata$adjuvant_chemo %in% c('No','no')]<-0
clidata$adjuvant_chemo[clidata$adjuvant_chemo %in% c('Yes','yes')]<-1
clidata<-clidata[clidata$TMScore<quantile(clidata$TMScore,c(0,1/3,2/3,1))[3],]
colnames(clidata)<-c(colnames(clidata)[1:3],'event','time','TMScore','sample.label')

pdf(file = paste0('OS','lungHighTMScore','.pdf'))
plot.surv(clidata, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'lungHighTMScore', ylab = "OS")
dev.off() # high/low TMScore ns





pdf(file = paste0('RFS','breastchemo','.pdf')) 
plot.surv(surcliRFSTertilesWeight$breast$breast, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'breastchemo', ylab = "RFS")
dev.off()

pdf(file = paste0('DMFS','breastchemo','.pdf'))
plot.surv(surcliDMFSTertilesWeight$breast$breast, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'breastchemo', ylab = "DMFS")
dev.off()


pdf(file = paste0('OS','breastadjchemo','.pdf'))
plot.surv(surcliOSTertilesWeight$breast$breast, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'breastadjchemo', ylab = "OS")
dev.off()

pdf(file = paste0('RFS','breastadjchemo','.pdf'))
plot.surv(surcliRFSTertilesWeight$breast$breast, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'breastadjchemo', ylab = "RFS")
dev.off()

pdf(file = paste0('DMFS','breastadjchemo','.pdf')) 
plot.surv(surcliDMFSTertilesWeight$breast$breast, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'breastadjchemo', ylab = "DMFS")
dev.off()

pdf(file = paste0('MFS','breastadjchemo','.pdf')) #
plot.surv(surcliMFSTertilesWeight$breast$breast, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'breastadjchemo', ylab = "MFS")
dev.off()

pdf(file = paste0('OS','colonNoTher','.pdf')) #
plot.surv(surcliOSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonNoTher', ylab = "OS")
dev.off()

pdf(file = paste0('DFS','colonNoTher','.pdf')) #
plot.surv(surcliDFSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonNoTher', ylab = "DFS")
dev.off()

pdf(file = paste0('DSS','colonNoTher','.pdf')) #
plot.surv(surcliDSSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonNoTher', ylab = "DSS")
dev.off()

pdf(file = paste0('PFS','colonNoTher','.pdf')) #
plot.surv(surcliPFSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonNoTher', ylab = "PFS")
dev.off()

pdf(file = paste0('OS1','colonGSE39582','.pdf')) #
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

pdf(file = paste0('DFS','colonNoTherGSE14333','.pdf')) 
plot.surv(surcliDFSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonNoTherGSE14333', ylab = "DFS")
dev.off()

pdf(file = paste0('DFS','colonTherGSE14333','.pdf')) 
plot.surv(surcliDFSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonTherGSE14333', ylab = "DFS")
dev.off()

pdf(file = paste0('DFS','colonNoRadioTherGSE14333','.pdf')) 
plot.surv(surcliDFSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonNoRadioTherGSE14333', ylab = "DFS")
dev.off()

pdf(file = paste0('DFS','colonGSE14333','.pdf')) 
plot.surv(surcliDFSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonGSE14333', ylab = "DFS")
dev.off()


pdf(file = paste0('OS','colonchemodataall','.pdf')) 
plot.surv(surcliOSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonchemodataall', ylab = "OS")
dev.off()

pdf(file = paste0('PFS','colonchemodataall','.pdf')) 
plot.surv(surcliPFSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonchemodataall', ylab = "PFS")
dev.off()

pdf(file = paste0('DFS','colonchemodataall','.pdf')) 
plot.surv(surcliDFSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonchemodataall', ylab = "DFS")
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

pdf(file = paste0('PFS','colonchemodataN','.pdf')) 
plot.surv(surcliPFSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonchemodataN', ylab = "PFS")
dev.off()

pdf(file = paste0('DFS','colonchemodataN','.pdf')) 
plot.surv(surcliDFSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonchemodataN', ylab = "DFS")
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

pdf(file = paste0('PFS','colonchemodataY','.pdf')) 
plot.surv(surcliPFSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonchemodataY', ylab = "PFS")
dev.off()

pdf(file = paste0('DFS','colonchemodataY','.pdf')) #
plot.surv(surcliDFSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonchemodataY', ylab = "DFS")
dev.off()

pdf(file = paste0('RFS','colonchemodataY','.pdf')) 
plot.surv(surcliRFSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonchemodataY', ylab = "RFS")
dev.off()


pdf(file = paste0('OS','colonchemodataYvsN','.pdf')) 
plot.surv(surcliOSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonchemodataYvsN', ylab = "OS")
dev.off()

pdf(file = paste0('PFS','colonchemodataYvsN','.pdf')) 
plot.surv(surcliPFSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonchemodataYvsN', ylab = "PFS")
dev.off()

pdf(file = paste0('DFS','colonchemodataYvsN','.pdf')) #
plot.surv(surcliDFSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonchemodataYvsN', ylab = "DFS")
dev.off()

pdf(file = paste0('RFS','colonchemodataYvsN','.pdf')) 
plot.surv(surcliRFSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonchemodataYvsN', ylab = "RFS")
dev.off()

pdf(file = paste0('OS','colonGSE39582chemodataYvsN','.pdf')) 

plot.surv(surcliOSTertilesWeight$colon$colon, upper.time = NULL, xscale = 1, xlab = "Time (Month)", median.time = FALSE, 
          surv.median.line = "none", HR = TRUE, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, main = 'colonGSE39582chemodataYvsN', ylab = "OS")
dev.off()