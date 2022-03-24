
#函数三 plot.surv
#' @description 绘制相应的KM生存曲线，给出对应的risk.table，当绘制图像仅包括两组时，在图中标注出HR及置信区间，每组样本对应的中位生存时间
#'   在副标题中给出
#' @param clinical.data: 具有分类标签的临床数据信息，至少包含四列：Patient_ID表示临床样本ID；event表示样本事件结局(数值型或逻辑型)；
#'   time表示样本事件的时间（数值型）；最后一列是通过上面的分类函数产生的分类标签sample.label，必须是factor
#' @param upper.time: 数值型，生存时间上限，必须与clinical.data中的单位一致，默认为NULL；如果使用，超过年限的样本将被去掉
#' @param xscale: 字符型，允许的选项包括"d_m"，"d_y"，"m_d"，"m_y"，"y_d"和"y_m"，其中d =天，m =月和y =年。 
#'   例如，xscale =“d_m”会将x轴单位从几天转换为几个月
#' @param xlab: 字符向量，x轴的标签
#' @param median.time: 逻辑向量，是否在图中给出中位生存时间及其95%置信区间，其单位会随着x轴单位的改变而改变，默认为TRUE
#' @param surv.median.line: 在中位生存期时绘制水平/垂直线的字符向量，允许的值包括c（"none"，"hv"，"h"，"v"）中的一个，
#'   v：垂直，h：水平，默认不给出（none）
#' @param HR: 逻辑值，是否在图中标注出HR信息。注意：当且仅当样本分为两组时可以使用 
#' @param risk.table: 逻辑向量，是否绘制risk.table，默认为TRUE
#' @param conf.int: 逻辑向量，是否画出置信区间，默认为FALSE
#' @param pval: 逻辑向量，是否给出log rank p值，默认为TRUE
#' @param ylab: 字符向量，y轴的标签
#' @param main: 主标题的名字

plot.surv <- function(clinical.data, upper.time = NULL, xscale = 1, xlab = "Time", median.time = TRUE, 
                      surv.median.line = "none", HR = FALSE, risk.table = TRUE, pval = TRUE, 
                      conf.int = FALSE, main = NULL, ylab = "Survival probability") {
  
  #载入相关R包
  require(survival)
  require(survminer)
  require(RColorBrewer)
  require(gridExtra)
  
	#确定事件类型和时间的单位
  # survival.event <- survival.event[1];
  # unit.xlabel <- unit.xlabel[1];
  
  #如果设置upper.time，则去除生存时间超过upper.time的样本
  if (!is.null(upper.time)) clinical.data <- clinical.data[clinical.data$time <= upper.time,]
  
  # #选择日期格式 
  # xSL <- data.frame(xScale=c(1,7,30,365.25),xLab=c("Days","Weeks","Months","Years"), stringsAsFactors=FALSE)
  # switch(unit.xlabel, year={xScale <- 365.25;}, month={xScale <- 30;}, week={xScale <- 7;}, day={xScale <- 1})
  # xLab <- xSL[which(xSL[,1]==xScale),2];
  
  #构造颜色
	if (!is.factor(clinical.data$sample.label)) 
    clinical.data$sample.label <- as.factor(clinical.data$sample.label)
  
  t.name <- levels(clinical.data$sample.label)

	if (length(t.name) > 6) stop("样本分组>6，超过函数接受范围")
    colors <- c("#808080","#EA4335","#4285F4","#FBBC05","#34A853","#000000") # 顺序：灰，红，蓝，黄，绿，黑
    t.col <- colors[1:length(t.name)]
    
  # 构造生存对象
  km.curves <- survfit(Surv(time, event)~sample.label, data=clinical.data)
    
    # 计算HR值和95%CI
	if (length(t.name) == 2) {
	  if (HR) {
        cox.obj <- coxph(Surv(time, event)~sample.label, data=clinical.data)
        tmp <- summary(cox.obj)
        HRs <- round(tmp$coefficients[ ,2], digits = 2)
        HR.confint.lower <- round(tmp$conf.int[,"lower .95"], 2)
        HR.confint.upper <- round(tmp$conf.int[,"upper .95"], 2)
        HRs <- paste0(HRs, " (", HR.confint.lower, "-", HR.confint.upper, ")")	  
	  }
  }

  # 构造生存图像中图例显示文字
  legend.content <- substr(names(km.curves$strata), start = 14, stop = 1000)
  
	# x轴刻度单位转换
  if (is.numeric(xscale) | (xscale %in% c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"))) {
    xscale = xscale
  } else {
    stop('xscale should be numeric or one of c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m").')
  }
    
	# 隐函数：转换生存时间单位
  .format_xticklabels <- function(labels, xscale){
    # 1 year = 365.25 days
    # 1 month = 365.25/12 = 30.4375 days
    if (is.numeric(xscale)) xtrans <- 1/xscale
    else
      xtrans <- switch(xscale,
                       d_m = 12/365.25,
                       d_y = 1/365.25,
                       m_d = 365.25/12,
                       m_y = 1/12,
                       y_d = 365.25,
                       y_m = 12,
                       1
      )
    round(labels*xtrans,2)
  }

  # 在图中添加中位生存时间及其95%CI,放在副标题位置
  subtitle <- NULL
  if (median.time) {
    if (is.numeric(xscale)) {
      median.km.obj = km.curves
    } else if (xscale %in% c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m")) {
      clinical.data$time <- .format_xticklabels(labels = clinical.data$time, xscale = xscale)
      median.km.obj <- survfit(Surv(time, event)~sample.label, data=clinical.data)
    }
    survival.time.info <- NULL
    survival.time.info <- rbind(survival.time.info, summary(median.km.obj)$table)
    median.survival <- round(survival.time.info[!duplicated(survival.time.info[,7:9]),7:9], digits = 2) # 注意：这里取得的置信区间上界可能为NA
    if (length(levels(clinical.data$sample.label)) == 1) {
      tmp1 <- levels(clinical.data$sample.label)
    } else {
      tmp1 <- do.call(rbind,strsplit(rownames(summary(median.km.obj)$table), split = "="))[,2]
    }
    tmp2 <- paste(median.survival[,1], "(", median.survival[,2], "-", median.survival[,3], ")")
    subtitle <- paste(tmp1, tmp2, sep = ":", collapse = "\n")
  }
    
  # ggsurvplot绘制生存图像
	ggsurv <- ggsurvplot(km.curves,               # survfit object with calculated statistics.
						 data = clinical.data,             # data used to fit survival curves.
						 palette = t.col,
						 
						 #图的主题构架
						 risk.table = risk.table,       # show risk table.
						 pval = pval,             # show p-value of log-rank test.
						 surv.median.line = surv.median.line,  # add the median survival pointer.
						 title = main,     #主标题名字
						 subtitle = subtitle, #副标题
						 font.main = 15,       #主标题字体大小              
						 xlab = xlab,   # customize X axis label.
						 ylab = ylab,   # customize Y axis label
						 xscale = xscale,


						 #图例设置
						 legend.title = "", #图例标题，一般不用，设为空
						 legend.labs = legend.content, #图例文字描述
						 legend = c(0.8,0.9), #图例的位置，取值在【0,1】之间
						 font.legend = 9,     #图例字体大小
						 
						 #risk table设置
						 tables.theme = theme_cleantable(),#table主题
						 risk.table.title = "No. at risk:",#table标题
						 risk.table.y.text.col = T, # 使用颜色代替Y轴文字
						 risk.table.y.text = FALSE, # Y轴不使用文字注释
						 tables.height = 0.15,      # table的高度
						 risk.table.fontsize = 3    # risk table内文字的大小
						)
	# HR标注的位置限定
  if (length(t.name) == 2) {
	  if (HR) 
      ggsurv$plot <- ggsurv$plot + ggplot2::annotate("text", x = max(km.curves$time)/13, 
                                                     y = 0.15, size = 5, label = paste("HR=", HRs))
  } 
	#图的标题居中
  ggsurv$plot <- ggsurv$plot + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(size = 10), 
								                     plot.margin = unit(c(5.5, 5.5, 5.5, 50), "points"))
  #表的标题
  ggsurv$table <- ggsurv$table + theme(plot.title = element_text(hjust = -0.04),
                                       plot.margin = unit(c(5.5, 5.5, 5.5, 50), "points"))

  # 判断分类的类数，如果只有两类，就不必计算两两之间的log rank p值
  if(length(t.name) > 2) {
  # 计算pairwise的log rank的p值
  res <- pairwise_survdiff(Surv(time, event)~sample.label, data=clinical.data);
  pairwise.pvalue <- round(res$p.value, digits = 4);
  pairwise.pvalue[which(pairwise.pvalue < 0.0001)] <- "<0.0001";
  pairwise.pvalue[is.na(pairwise.pvalue)] <- "-"
  
  # 添加表格
  tt <- ttheme_minimal(core = list(fg_params = list(col = "black"),bg_params = list(fill = NA, col = "black")),
                       colhead = list(fg_params = list(col = NA),bg_params = list(fill = t.col, col = "black")),
                       rowhead = list(fg_params = list(col = NA, hjust = 1),bg_params = list(fill = c("white",t.col[-1]), col = "black"))
                      )
  pairwise.table <- tableGrob(pairwise.pvalue, theme = tt)
  ggsurv <- ggarrange(ggarrange(ggsurv$plot, ggsurv$table, nrow=2, heights=c(2,0.5)),
                      pairwise.table, nrow=2, heights = c(2,0.5),
                      labels = c("","p from pairwise comparisons"),
                      hjust = 0, font.label = list(size = 15, face = "plain"))
  }

  ggsurv	
}

