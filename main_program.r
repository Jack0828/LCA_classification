# #####1. 加载包
library(poLCA)
library(tidyLPA)
library(survival)
library(dplyr)
library(purrr)
library(cmprsk)
library(ggplot2)
library(survminer)
library(adjustedCurves)
library(tidyr)
library(broom)
# 2.1设置工作路径
# setwd("C:\\....")


#######################第一部分：seer数据库的相关分析#############################

# 3.2读取文件，设置用于LCA的formula，处理各变量为因子类型，准备部分变量供存放结果。
seerdata <- read.csv("LCASeer.csv")
# 移除 Patient_ID 重复两次或以上的所有行
seerdata_clean <- seerdata %>%
  group_by(Patient_ID) %>%
  filter(n() == 1) %>%
  ungroup()
seerdata_clean1 <- seerdata %>%
  group_by(Patient_ID) %>%
  filter(n() == 1) %>%
  ungroup()

seerdata_clean <- seerdata_clean[, c("Race", "Age_class", "Pathology", "Stage")]
seerdata_clean$Race <- as.integer(as.factor(seerdata_clean$Race))
seerdata_clean$Age_class <- as.integer(as.factor(seerdata_clean$Age_class))
seerdata_clean$Pathology <- as.integer(as.factor(seerdata_clean$Pathology))
seerdata_clean$Stage <- as.integer(as.factor(seerdata_clean$Stage))
formula_s <- as.formula(paste("cbind( Race, Age_class, Pathology, Stage) ~ 1"))
results_s1 <- results_s <- vector("list", 6) # 为4个模型创建一个长度为4的列表
lmr_p_values_s <- min_averages_s <- numeric(6)  
# 3.3开始进行LCA循环
set.seed(123)
i <- 1
j <- 1
for (i in 1:6) {
  lca_result <- poLCA(formula_s, data = seerdata_clean, nclass = i,maxiter = 10000, calc.se = FALSE)
  colname <- paste("LCA", i, "class", sep = "") 
  seerdata_clean[[colname]] <- lca_result$predclass
  results_s1 <- lca_result
  results_s[[i]] <- list(AIC = lca_result$aic, BIC = lca_result$bic, SABIC = lca_result$bic - log(nrow(seerdata_clean)) * (i - 1) * length(lca_result$probs[[1]]), Entropy = poLCA.entropy(lca_result), Min_pro = min(lca_result$P), llik= lca_result$llik,npar=lca_result$npar,P=lca_result$P)
  if(i > 1) {
    combined_result <- as.data.frame(cbind(lca_result$posterior, predclass = lca_result$predclass))
    class_counts <- table(combined_result$predclass)
    results_s[[i]]$class_counts <- class_counts
    datasets <- split(combined_result, combined_result$predclass)
    class_averages <- numeric(length(datasets))
    for (j in 1:length(datasets)) {
      class_averages[j] <- mean(datasets[[j]][, j])
    }
    min_averages_s[i] <- min(class_averages)
    n_obs <- nrow(seerdata_clean)
    lmr_result <- calc_lrt(n_obs, results_s[[i-1]]$llik, results_s[[i-1]]$npar, i-1, lca_result$llik,lca_result$npar, i)
    lmr_p_values_s[i] <- lmr_result[4]
  } else {
    min_averages_s[i] <- 1
    lmr_p_values_s[i] <- NA
  }
}
results_df_s <- map_dfr(lapply(results_s, function(x) {x$class_counts <- NULL
return(x)}), ~as.data.frame(t(unlist(.x))))
results_df_s$lmr_p_values <- lmr_p_values_s
results_df_s$min_averages <- min_averages_s
results_df_s

#3.4选定了最佳类别之后，再单独跑一次LCA以提取分类的信息。
LCASEER4class <- poLCA(formula_s, data = seerdata_clean, nclass = 4, maxiter = 15000,calc.se = FALSE)
seerdata_clean1$LCA4class <- LCASEER4class$predclass

seerCSSdata <- seerdata_clean1
reference_seer <- read.csv("reference.csv")
seerCSSdata <- merge(seerCSSdata, reference_seer, by = "Patient_ID")
reference_seer_after <- read.csv("reference_seer_after.csv")
seerCSSdata <- merge(seerCSSdata, reference_seer_after, by = "Patient_ID")

# 首先要把表格中的除了生存时间之外的变量转换为因子变量
seerCSSdata$CSS <- as.factor(seerCSSdata$CSS)
seerCSSdata$LCA4class <- factor(seerCSSdata$LCA4class, levels = sort(unique(seerCSSdata$LCA4class)))
seerCSSdata$Race <- as.factor(seerCSSdata$Race.x)
seerCSSdata$Age_class <- as.factor(seerCSSdata$Age_class)
seerCSSdata$Stage <- as.factor(seerCSSdata$Stage)
# seerCSSdata$LCA4class <- relevel(seerCSSdata$LCA4class, ref = "2")
fit_seer <- survfit(Surv(OS_time, OS_status) ~ LCA4class, data = seerCSSdata)
ggsurvplot(fit_seer, 
           data = seerCSSdata,
           pval = TRUE,                      
           risk.table = TRUE,
           risk.table.y.text = FALSE,
           strata = FALSE,
           conf.int = FALSE,                
           palette = c("#ffc000","#f84d4c","#8f7ef1","#00b0f1"), # 自定义颜色
           title = "Survival curves for LCA 4 class",
           xlab = "Time",
           ylab = "Survival Probability")


CSS01 <- subset(seerCSSdata, CSS == 1 | CSS == 0)
CSS02 <- subset(seerCSSdata, CSS == 2 | CSS == 0)
fitCSS01 <- coxph(Surv(OS_time, OS_status) ~ LCA4class, data = CSS01, robust = TRUE)
fitCSS02 <- coxph(Surv(OS_time, OS_status) ~ LCA4class, data = CSS02, robust = TRUE)
summary(fitCSS01)
summary(fitCSS02)

seerCSSdata$LCA4class <- as.character(seerCSSdata$LCA4class)
seerCSSdata$LCA4class[seerCSSdata$LCA4class == "1"] <- "a"
seerCSSdata$LCA4class[seerCSSdata$LCA4class == "2"] <- "b"
seerCSSdata$LCA4class[seerCSSdata$LCA4class == "3"] <- "c"
seerCSSdata$LCA4class[seerCSSdata$LCA4class == "4"] <- "d"

seerCSSdata$LCA4class[seerCSSdata$LCA4class == "a"] <- "3"
seerCSSdata$LCA4class[seerCSSdata$LCA4class == "b"] <- "1"
seerCSSdata$LCA4class[seerCSSdata$LCA4class == "c"] <- "4"
seerCSSdata$LCA4class[seerCSSdata$LCA4class == "d"] <- "2"


seerCSSdata$LCA4class <- as.integer(seerCSSdata$LCA4class)
#绘制桑基图数据
sankydata <- data.frame(LCA4class = c(1, 2, 3, 4), CSS0 = c(0, 0, 0, 0), CSS1 = c(0, 0, 0, 0), CSS2 = c(0, 0, 0, 0))
for (i in 1:4) {
  sankydata$CSS0[i] <- sum(seerCSSdata$LCA4class == i & seerCSSdata$CSS == 0)
  sankydata$CSS1[i] <- sum(seerCSSdata$LCA4class == i & seerCSSdata$CSS == 1)
  sankydata$CSS2[i] <- sum(seerCSSdata$LCA4class == i & seerCSSdata$CSS == 2)
}
#绘制森林图数据
histo_age <- data.frame(LCA4class = c(1, 2, 3, 4), Age1 = c(0, 0, 0, 0), Age2 = c(0, 0, 0, 0), Age3 = c(0, 0, 0, 0), Age4 = c(0, 0, 0, 0))
for (i in 1:4) {
  histo_age$Age1[i] <- sum(seerCSSdata$LCA4class == i & seerCSSdata$Age_class == 1)
  histo_age$Age2[i] <- sum(seerCSSdata$LCA4class == i & seerCSSdata$Age_class == 2)
  histo_age$Age3[i] <- sum(seerCSSdata$LCA4class == i & seerCSSdata$Age_class == 3)
  histo_age$Age4[i] <- sum(seerCSSdata$LCA4class == i & seerCSSdata$Age_class == 4)
}

histo_stage <- data.frame(LCA4class = c(1, 2, 3, 4), Stage1 = c(0, 0, 0, 0), Stage2 = c(0, 0, 0, 0), Stage3 = c(0, 0, 0, 0), Stage4 = c(0, 0, 0, 0))
for (i in 1:4) {
  histo_stage$Stage1[i] <- sum(seerCSSdata$LCA4class == i & seerCSSdata$Stage == 1)
  histo_stage$Stage2[i] <- sum(seerCSSdata$LCA4class == i & seerCSSdata$Stage == 2)
  histo_stage$Stage3[i] <- sum(seerCSSdata$LCA4class == i & seerCSSdata$Stage == 3)
  histo_stage$Stage4[i] <- sum(seerCSSdata$LCA4class == i & seerCSSdata$Stage == 4)
}

histo_Pathology <- data.frame(LCA4class = c(1, 2, 3, 4), Pathology1 = c(0, 0, 0, 0), Pathology2 = c(0, 0, 0, 0), Pathology3 = c(0, 0, 0, 0), Pathology4 = c(0, 0, 0, 0), Pathology5 = c(0, 0, 0, 0))
for (i in 1:4) {
  histo_Pathology$Pathology1[i] <- sum(seerCSSdata$LCA4class == i & seerCSSdata$Pathology == 1)
  histo_Pathology$Pathology2[i] <- sum(seerCSSdata$LCA4class == i & seerCSSdata$Pathology == 2)
  histo_Pathology$Pathology3[i] <- sum(seerCSSdata$LCA4class == i & seerCSSdata$Pathology == 3)
  histo_Pathology$Pathology4[i] <- sum(seerCSSdata$LCA4class == i & seerCSSdata$Pathology == 4)
  histo_Pathology$Pathology5[i] <- sum(seerCSSdata$LCA4class == i & seerCSSdata$Pathology == 5)
}

histo_Race <- data.frame(LCA4class = c(1, 2, 3, 4), Race1 = c(0, 0, 0, 0), Race2 = c(0, 0, 0, 0), Race3 = c(0, 0, 0, 0), Race4 = c(0, 0, 0, 0))
for (i in 1:4) {
  histo_Race$Race1[i] <- sum(seerCSSdata$LCA4class == i & seerCSSdata$Race == 1)
  histo_Race$Race2[i] <- sum(seerCSSdata$LCA4class == i & seerCSSdata$Race == 2)
  histo_Race$Race3[i] <- sum(seerCSSdata$LCA4class == i & seerCSSdata$Race == 3)
  histo_Race$Race4[i] <- sum(seerCSSdata$LCA4class == i & seerCSSdata$Race == 4)
}

seerCSSdata$LCA4class <- as.character(seerCSSdata$LCA4class)
seerCSSdata$LCA4class[seerCSSdata$LCA4class == "3"] <- "c"
seerCSSdata$LCA4class[seerCSSdata$LCA4class == "4"] <- "d"
seerCSSdata$LCA4class[seerCSSdata$LCA4class == "c"] <- "4"
seerCSSdata$LCA4class[seerCSSdata$LCA4class == "d"] <- "3"
seerCSSdata$LCA4class <- as.integer(seerCSSdata$LCA4class)

#接下来，我需要把seerCSSdata的RT列中的0和2统一替换为0，然后把1,3,4,5，统一替换为1.
seerCSSdata$RT[seerCSSdata$RT == 0] <- 0
seerCSSdata$RT[seerCSSdata$RT == 2] <- 0
seerCSSdata$RT[seerCSSdata$RT == 1] <- 1
seerCSSdata$RT[seerCSSdata$RT == 3] <- 1
seerCSSdata$RT[seerCSSdata$RT == 4] <- 1
seerCSSdata$RT[seerCSSdata$RT == 5] <- 1
seerCSSdata$RT <- as.factor(seerCSSdata$RT)

seerstage1 <- seerCSSdata[seerCSSdata$Stage == 1,]
seerstage2 <- seerCSSdata[seerCSSdata$Stage == 2,]
seerstage3 <- seerCSSdata[seerCSSdata$Stage == 3,]
seerstage4 <- seerCSSdata[seerCSSdata$Stage == 4,]

cox_seer_stage1 <- coxph(Surv(OS_time, OS_status) ~ RT+LCA4class, data = seerstage1)
cox_seer_stage2 <- coxph(Surv(OS_time, OS_status) ~ RT+LCA4class, data = seerstage2)
cox_seer_stage3 <- coxph(Surv(OS_time, OS_status) ~ RT+LCA4class, data = seerstage3)
cox_seer_stage4 <- coxph(Surv(OS_time, OS_status) ~ RT+LCA4class, data = seerstage4)
summary(cox_seer_stage1)
summary(cox_seer_stage2)
summary(cox_seer_stage3)
summary(cox_seer_stage4)

seerCSSdata$RT <- as.factor(seerCSSdata$RT)
seerCSSdata$RT <- relevel(seerCSSdata$RT, ref = "0")

seer_cox_results <- data.frame(Stage=integer(), Age_class=integer(), LCA4class=integer(), RT0_count=integer(), RT1_count=integer(),TotalCount=integer(), HR=character(), CI=character(), pval=character())

for(stage in unique(seerCSSdata$Stage)) {
  for(age in unique(seerCSSdata$Age_class)) {
    for(lca in unique(seerCSSdata$LCA4class)) {
      subset_data <- seerCSSdata[seerCSSdata$Stage == stage & seerCSSdata$Age_class == age & seerCSSdata$LCA4class == lca, ]
      RT0_count <- sum(subset_data$RT == "0")
      RT1_count <- sum(subset_data$RT == "1")
      TotalCount <- nrow(subset_data)
      if(TotalCount > 0){
        cox_model <- coxph(Surv(OS_time, OS_status) ~ RT, data = subset_data)
        summary_cox <- summary(cox_model)
        
        HR <- round(summary_cox$conf.int[,"exp(coef)"], 2)
        CI_lower <- round(summary_cox$conf.int[,"lower .95"], 2)
        CI_upper <- round(summary_cox$conf.int[,"upper .95"], 2)
        pval <- round(summary_cox$coefficients[,"Pr(>|z|)"], 3)
      } else {
        HR <- CI_lower <- CI_upper <- pval <- "—"
      }
      seer_cox_results <- rbind(seer_cox_results, data.frame(Stage=stage, Age_class=age, LCA4class=lca, RT0_count=RT0_count, RT1_count=RT1_count,TotalCount=TotalCount, HR=HR, CI=paste(CI_lower, CI_upper, sep="-"), pval=pval))
    }
  }
}


perform_cox <- function(data, variable) {
  results <- data.frame(Variable=character(), HR=numeric(), CI=character(), pval=numeric(), stringsAsFactors=FALSE)
  
  for(val in unique(data[[variable]])) {
    subset_data <- data[data[[variable]] == val,]
    if(nrow(subset_data) > 0) {
      cox_model <- coxph(Surv(OS_time, OS_status) ~ LCA4class, data = subset_data)
      summary_cox <- summary(cox_model)
      if (summary(cox_model)$n > 0) {
        HR <- round(summary_cox$conf.int[,"exp(coef)"], 2)
        CI <- paste(round(summary_cox$conf.int[,"lower .95"], 2), 
                    round(summary_cox$conf.int[,"upper .95"], 2), sep="-")
        pval <- round(summary_cox$coefficients[,"Pr(>|z|)"], 3)
      } else {HR <- CI <- pval <- NA}} else {HR <- CI <- pval <- NA}
    results <- rbind(results, data.frame(Variable=as.character(val), HR=HR, CI=CI, pval=pval))
  }
  return(results)
}

cox_stage_result <- perform_cox(seerCSSdata, "Stage")
cox_race_result <- perform_cox(seerCSSdata, "Race")
cox_age_result <- perform_cox(seerCSSdata, "Age_class")
cox_pathology_result <- perform_cox(seerCSSdata, "Pathology")
# cox_lca_result <- perform_cox(seerCSSdata, "LCA4class")
cox_stage_result
cox_age_result
cox_pathology_result
# cox_lca_result

###########绘制CIF图###############



cif_data <- data.frame()
for(class in unique(seerCSSdata$LCA4class)) {
  sub_data <- seerCSSdata[seerCSSdata$LCA4class == class, ]
  cif_result <- cuminc(ftime = sub_data$OS_time, fstatus = sub_data$CSS)
  css_death <- data.frame(
    time = cif_result[[2]]$time,
    est = cif_result[[2]]$est,
    type = "CSS Death",
    class = class
  )
  non_css_death <- data.frame(
    time = cif_result[[1]]$time,
    est = cif_result[[1]]$est,
    type = "Non-CSS Death",
    class = class
  )
  cif_data <- rbind(cif_data, css_death, non_css_death)
}
cif_data$class <- factor(cif_data$class, levels = sort(unique(seerCSSdata$LCA4class)))
colors <- c("#FFC000","#F84D4D", "#8E7EF0",  "#00B0F0")
ggplot(cif_data, aes(x = time, y = est, color = class, linetype = type)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors) +
  labs(x = "Time", y = "Cumulative Incidence", title = "CIF for CSS Death by LCA4class") +
  theme_minimal() +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())  # 移除图例标题



###### 第二部分：TCGA#########
tcgadata <- read.csv("LCATCGA.csv")####可以在这段文本之后统一进行替换从而得到统一的变量名赋值。
#如果想导入txt文件，
tcgadata$Race <- as.factor(tcgadata$Race)
tcgadata$Age_class <- as.factor(tcgadata$Age_class)
tcgadata$Pathology <- as.factor(tcgadata$Pathology)
tcgadata$Stage <- as.factor(tcgadata$Stage)
formula_t <- as.formula(paste("cbind(Race, Age_class, Pathology, Stage) ~ 1"))# 创建formula用于LCA的公式
results_t <- vector("list", 6) # 为6个模型创建一个长度为6的列表
lmr_p_values_t <- min_averages_t <- numeric(6)   # 存储每个模型的最小平均值# 存储LMR检验的p值
set.seed(123)
# 2.3开始进行LCA循环，并展示结果
i <- 1
j <- 1
for (i in 1:6) {
  lca_result <- poLCA(formula_t, data = tcgadata, nclass = i, maxiter = 15000, nrep = 30, calc.se = FALSE, verbose =FALSE,)
  results_t[[i]] <- list(AIC = lca_result$aic, BIC = lca_result$bic, SABIC = lca_result$bic - log(nrow(tcgadata)) * (i - 1) * length(lca_result$probs[[1]]), Entropy = poLCA.entropy(lca_result), Min_pro =round(min(lca_result$P),3),llik= lca_result$llik, npar=lca_result$npar,P=round(lca_result$P,3))# 对1到6个潜在类别进行LCA
  if(i > 1) {
    combined_result <- as.data.frame(cbind(lca_result$posterior, predclass = lca_result$predclass))# 合并posterior和predclass，并转换为数据框
    class_counts <- table(combined_result$predclass)# 新增 - 统计每个类别的数量
    results_t[[i]]$class_counts <- class_counts
    datasets <- split(combined_result, combined_result$predclass)# 根据predclass分割数据
    class_averages <- numeric(length(datasets))# 初始化向量存储每个类均值
    for (j in 1:length(datasets)) {
      class_averages[j] <- mean(datasets[[j]][, j])
    }# 计算每个类别的平均值
    min_averages_t[i] <- round(min(class_averages),3)# 找到最小的平均值并存储
    n_obs <- nrow(tcgadata)# 计算LMR检验的p值
    lmr_result <- calc_lrt(n_obs, results_t[[i-1]]$llik, results_t[[i-1]]$npar, i-1, lca_result$llik,lca_result$npar, i)
    lmr_p_values_t[i] <- lmr_result[4] # LMR检验的p值是结果的第四个元素
  } else {
    min_averages_t[i] <- 1
    lmr_p_values_t[i] <- NA
  }
}
# list(min_averages_t = min_averages_t, lmr_p_values_t = lmr_p_values_t)#展示结果
#####汇总每个输出表得到模型拟合的参数比较数据框，直接可以赋值粘贴到表格中。
results_df_t <- map_dfr(lapply(results_t, function(x) {x$class_counts <- NULL
return(x)}), ~as.data.frame(t(unlist(.x))))
results_df_t$llik <- results_df_t$npar <- NULL
results_df_t$lmr_p_values <- lmr_p_values_t
results_df_t$min_averages <- min_averages_t
results_df_t
#####
reference_tcga_after <- read.csv("reference_tcga_after.csv")
tcgadata <- merge(tcgadata, reference_tcga_after, by = "Seq")
#####
# 2.4根据参数选定了最佳类别之后，再单独跑一次LCA以提取分类的信息。
LCATCGA <- poLCA(formula_t, data = tcgadata, nclass = 4, maxiter = 10000, nrep=30,calc.se = FALSE)
tcgadata$LCA4class <- LCATCGA$predclass


# 5.1准备进行生存分析和cox回归，比较原始PAM50模型在LCA分类结果的权重调整前后的性能
tcgadata$OS_status <- as.numeric(tcgadata$OS_status) 
tcgadata$PAM50 <- as.factor(tcgadata$PAM50)
tcgadata$LCA4class <- as.factor(tcgadata$LCA4class)
tcgadata$PAM50 <- factor(tcgadata$PAM50)
tcgadata$PAM50 <- relevel(tcgadata$PAM50, ref = "Normal")
fit_total <- survfit(Surv(OS_time, OS_status) ~ PAM50, data = tcgadata)

########二、PAM50原始模型的生存分析和C指数、鲁棒性等信息###########

fit_PAM50 <- coxph(Surv(OS_time, OS_status) ~ PAM50, data = tcgadata,robust = TRUE)
summary(fit_PAM50)

########二、经【边际方法】调整后的性能信息###########
tcgadata$LCA4class <- as.factor(tcgadata$LCA4class)
fit_LCA4class <- coxph(Surv(OS_time, OS_status) ~ LCA4class, data = tcgadata,robust = TRUE)# 运行Cox模型，得到C指数和鲁棒性评分(比较PAM50得到PAM50>LCA4class)
summary(fit_LCA4class)
tcgadata$LCA4class <- factor(tcgadata$LCA4class, levels = sort(unique(tcgadata$LCA4class)))# 按照数字大小对因子的水平进行排序
hazard_ratios_LCA4class <- exp(coef(fit_LCA4class))# 计算每个LCA4class的风险比，并为参考类别（类别"1"）添加风险比1
hazard_ratios_LCA4class <- c(1, hazard_ratios_LCA4class)  # 类别“1”作为参考类别
weights_LCA4class <- sapply(tcgadata$LCA4class, function(x) 1 / hazard_ratios_LCA4class[as.integer(x)])# 设置权重（这里以1/风险比为权重）
fit_PAM50_after <- coxph(Surv(OS_time, OS_status) ~ PAM50, data = tcgadata, weights = weights_LCA4class, robust = TRUE)# 使用权重拟合鲁棒性Cox模型
tcgadata$weights_LCA4class <- weights_LCA4class# 将权重添加到数据框中
summary(fit_PAM50_after)# 查看模型摘要，包括鲁棒性得分(使用LCA4class分类结合边际方法对PAM50分型进行优化，总体的结果比较明显)
###############二、子分类的两两组合进行cox##########################
# 定义一个函数来执行两两组合的Cox分析
perform_cox_analysis_1 <- function(class1, class2, data) {
  subset_data <- data[data$PAM50 %in% c(class1, class2), ]
  cox_model <- coxph(Surv(OS_time, OS_status) ~ PAM50, data = subset_data, robust = TRUE)
  return(summary(cox_model))
}
pam50_classes <- unique(tcgadata$PAM50)
cox_results_1 <- cox_results_2 <-list()
combinations <- combn(pam50_classes, 2, simplify = FALSE) 
cox_results_df_1 <- cox_results_df_2 <- data.frame(C_Index = numeric(), Robustness = numeric(), P_Value = numeric(), row.names = character())
for (comb in combinations) {
  cox_results_1[[paste(comb[1], "vs", comb[2])]] <- perform_cox_analysis_1(comb[1], comb[2], tcgadata)
}

for (comb_name in names(cox_results_1)) {
  c_index <- cox_results_1[[comb_name]][["concordance"]][["C"]]
  robustness_score <- cox_results_1[[comb_name]][["robscore"]][["test"]]
  p_value <- cox_results_1[[comb_name]][["robscore"]][["pvalue"]]
  cox_results_df_1[comb_name, ] <- c(C_Index = round(c_index,3), Robustness = round(robustness_score,3), P_Value = round(p_value,3))}
cox_results_df_1
perform_cox_analysis_2 <- function(class1, class2, data) {
  subset_data <- data[data$PAM50 %in% c(class1, class2), ]
  cox_model <- coxph(Surv(OS_time, OS_status) ~ PAM50, data = subset_data, weights = subset_data$weights_LCA4class, robust = TRUE)
  return(summary(cox_model))
}
pam50_classes <- unique(tcgadata$PAM50)
cox_results_2 <- list()
for (comb in combinations) {
  cox_results_2[[paste(comb[1], "vs", comb[2])]] <- perform_cox_analysis_2(comb[1], comb[2], tcgadata)
}
for (comb_name in names(cox_results_2)) {
  c_index <- cox_results_2[[comb_name]][["concordance"]][["C"]]
  robustness_score <- cox_results_2[[comb_name]][["robscore"]][["test"]]
  p_value <- cox_results_2[[comb_name]][["robscore"]][["pvalue"]]
  cox_results_df_2[comb_name, ] <- c(C_Index = round(c_index,3), Robustness = round(robustness_score,3), P_Value = round(p_value,3))}
cox_results_df_1
cox_results_df_2

#######################二、使用ggsurvplot绘制生存曲线###############################
fit_total
C1_1 <- coxph(Surv(OS_time, OS_status) ~ PAM50, data = tcgadata)
S1_1 <- summary(C1_1)
df1_1 <- paste0(round(S1_1$conf.int,2)[,1],"(",round(S1_1$conf.int,2)[,3],"-",round(S1_1$conf.int,2)[,4],")")
ggsurvplot(fit_total, 
           data = tcgadata,
           pval = TRUE,                      
           risk.table = TRUE,               
           risk.table.title="",
           risk.table.y.text = FALSE,           
           conf.int = FALSE,                
           palette = c("#B6EE97","#68A7F8","#EB7D62","#7788E6","#F2D19E"), 
           legend.labs = c("Basal","Her2","LumA","LumB","Normal"),
           title = "",
           xlab = "Time(year)",
           ylab = "Survival Probability",
           xlim = c(0, 3650),               
           ylim = c(0, 1),                  
           break.time.by = 365
)          

# 使用ggsurvplot绘制基于LCA4class分类的生存曲线
fit_total_after <- survfit(Surv(OS_time, OS_status) ~ PAM50, data = tcgadata,weights = weights_LCA4class)
C1_2 <- coxph(Surv(OS_time, OS_status) ~ PAM50, data = tcgadata,weights = weights_LCA4class)
S1_2 <- summary(C1_2)
df1_2 <- paste0(round(S1_2$conf.int,2)[,1],"(",round(S1_2$conf.int,2)[,3],"-",round(S1_2$conf.int,2)[,4],")")
ggsurvplot(fit_total_after, 
           data = tcgadata,
           pval = TRUE,                      
           risk.table = TRUE,               
           risk.table.title="",
           risk.table.y.text = FALSE,           
           conf.int = FALSE,                
           palette = c("#B6EE97","#68A7F8","#EB7D62","#7788E6","#F2D19E"), 
           legend.labs = c("Basal","Her2","LumA","LumB","Normal"),
           title = "",
           xlab = "Time(year)",
           ylab = "Survival Probability",
           xlim = c(0, 3650),               
           ylim = c(0, 1),                  
           break.time.by = 365
)          


#####PFS生存曲线#####
tcgadata$PFS_status <- as.numeric(tcgadata$PFS_status)

fit_PFS <- survfit(Surv(Followup_time, PFS_status) ~ PAM50, data = tcgadata)
C2_1 <- coxph(Surv(Followup_time, PFS_status) ~ PAM50, data = tcgadata)
S2_1 <- summary(C2_1)
df2_1 <- paste0(round(S2_1$conf.int,2)[,1],"(",round(S2_1$conf.int,2)[,3],"-",round(S2_1$conf.int,2)[,4],")")
ggsurvplot(fit_PFS, 
           data = tcgadata,
           pval = TRUE,                      
           risk.table = TRUE,               
           risk.table.title="",
           risk.table.y.text = FALSE,           
           conf.int = FALSE,                
           palette = c("#B6EE97","#68A7F8","#EB7D62","#7788E6","#F2D19E"), 
           legend.labs = c("Basal","Her2","LumA","LumB","Normal"),
           title = "",
           xlab = "Time(year)",
           ylab = "Survival Probability",
           xlim = c(0, 3650),               
           ylim = c(0, 1),                  
           break.time.by = 365
)          

fit_PFS_after <- survfit(Surv(Followup_time, PFS_status) ~ PAM50, data = tcgadata,weights = weights_LCA4class)
C2_2 <- coxph(Surv(Followup_time, PFS_status) ~ PAM50, data = tcgadata,weights = weights_LCA4class)
S2_2 <- summary(C2_2)
df2_2 <- paste0(round(S2_2$conf.int,2)[,1],"(",round(S2_2$conf.int,2)[,3],"-",round(S2_2$conf.int,2)[,4],")")
ggsurvplot(fit_PFS_after, 
           data = tcgadata,
           pval = TRUE,                      
           risk.table = TRUE,               
           risk.table.title="",
           risk.table.y.text = FALSE,           
           conf.int = FALSE,                
           palette = c("#B6EE97","#68A7F8","#EB7D62","#7788E6","#F2D19E"), 
           legend.labs = c("Basal","Her2","LumA","LumB","Normal"),
           title = "",
           xlab = "Time(year)",
           ylab = "Survival Probability",
           xlim = c(0, 3650),               
           ylim = c(0, 1),                  
           break.time.by = 365
)          

#####DFS生存曲线#####
tcgadata$DFS_status <- as.numeric(tcgadata$DFS_status) 
fit_DFS <- survfit(Surv(Followup_time, DFS_status) ~ PAM50, data = tcgadata)
C3_1 <- coxph(Surv(Followup_time, DFS_status) ~ PAM50, data = tcgadata)
S3_1 <- summary(C3_1)
df3_1 <- paste0(round(S3_1$conf.int,2)[,1],"(",round(S3_1$conf.int,2)[,3],"-",round(S3_1$conf.int,2)[,4],")")
ggsurvplot(fit_DFS, 
           data = tcgadata,
           pval = TRUE,                      
           risk.table = TRUE,               
           risk.table.title="",
           risk.table.y.text = FALSE,           
           conf.int = FALSE,                
           palette = c("#B6EE97","#68A7F8","#EB7D62","#7788E6","#F2D19E"), 
           legend.labs = c("Basal","Her2","LumA","LumB","Normal"),
           title = "",
           xlab = "Time(year)",
           ylab = "Survival Probability",
           xlim = c(0, 3650),               
           ylim = c(0, 1),                  
           break.time.by = 365
)          

fit_DFS_after <- survfit(Surv(Followup_time, DFS_status) ~ PAM50, data = tcgadata,weights = weights_LCA4class)
C3_2 <- coxph(Surv(Followup_time, DFS_status) ~ PAM50, data = tcgadata,weights = weights_LCA4class)
S3_2 <- summary(C3_2)
df3_2 <- paste0(round(S3_2$conf.int,2)[,1],"(",round(S3_2$conf.int,2)[,3],"-",round(S3_2$conf.int,2)[,4],")")
ggsurvplot(fit_DFS_after, 
           data = tcgadata,
           pval = TRUE,                      
           risk.table = TRUE,               
           risk.table.title="",
           risk.table.y.text = FALSE,           
           conf.int = FALSE,                
           palette = c("#B6EE97","#68A7F8","#EB7D62","#7788E6","#F2D19E"), 
           legend.labs = c("Basal","Her2","LumA","LumB","Normal"),
           title = "",
           xlab = "Time(year)",
           ylab = "Survival Probability",
           xlim = c(0, 3650),               
           ylim = c(0, 1),                  
           break.time.by = 365
)          

#####LC生存曲线#####
tcgadata$LC_status <- as.numeric(tcgadata$LC_status) 
fit_LC <- survfit(Surv(Followup_time, LC_status) ~ PAM50, data = tcgadata)
C4_1 <- coxph(Surv(Followup_time, LC_status) ~ PAM50, data = tcgadata)
S4_1 <- summary(C4_1)
df4_1 <- paste0(round(S4_1$conf.int,2)[,1],"(",round(S4_1$conf.int,2)[,3],"-",round(S4_1$conf.int,2)[,4],")")
ggsurvplot(fit_LC, 
           data = tcgadata,
           pval = TRUE,                      
           risk.table = TRUE,               
           risk.table.title="",
           risk.table.y.text = FALSE,           
           conf.int = FALSE,                
           palette = c("#B6EE97","#68A7F8","#EB7D62","#7788E6","#F2D19E"), 
           legend.labs = c("Basal","Her2","LumA","LumB","Normal"),
           title = "",
           xlab = "Time(year)",
           ylab = "Survival Probability",
           xlim = c(0, 3650),               
           ylim = c(0, 1),                  
           break.time.by = 365
)          


fit_LC_after <- survfit(Surv(Followup_time, LC_status) ~ PAM50, data = tcgadata,weights = weights_LCA4class)
C4_2 <- coxph(Surv(Followup_time, LC_status) ~ PAM50, data = tcgadata,weights = weights_LCA4class)
S4_2 <- summary(C4_2)
df4_2 <- paste0(round(S4_2$conf.int,2)[,1],"(",round(S4_2$conf.int,2)[,3],"-",round(S4_2$conf.int,2)[,4],")")
ggsurvplot(fit_LC_after, 
           data = tcgadata,
           pval = TRUE,                      
           risk.table = TRUE,               
           risk.table.title="",
           risk.table.y.text = FALSE,           
           conf.int = FALSE,                
           palette = c("#B6EE97","#68A7F8","#EB7D62","#7788E6","#F2D19E"), 
           legend.labs = c("Basal","Her2","LumA","LumB","Normal"),
           title = "",
           xlab = "Time(year)",
           ylab = "Survival Probability",
           xlim = c(0, 3650),               
           ylim = c(0, 1),                  
           break.time.by = 365
)          

#####DM生存曲线#####
tcgadata$DM_status <- as.numeric(tcgadata$DM_status) 
fit_DM <- survfit(Surv(Followup_time, DM_status) ~ PAM50, data = tcgadata)
C5_1 <- coxph(Surv(Followup_time, DM_status) ~ PAM50, data = tcgadata)
S5_1 <- summary(C5_1)
df5_1 <- paste0(round(S5_1$conf.int,2)[,1],"(",round(S5_1$conf.int,2)[,3],"-",round(S5_1$conf.int,2)[,4],")")
ggsurvplot(fit_DM, 
           data = tcgadata,
           pval = TRUE,                      
           risk.table = TRUE,               
           risk.table.title="",
           risk.table.y.text = FALSE,           
           conf.int = FALSE,                
           palette = c("#B6EE97","#68A7F8","#EB7D62","#7788E6","#F2D19E"), 
           legend.labs = c("Basal","Her2","LumA","LumB","Normal"),
           title = "",
           xlab = "Time(year)",
           ylab = "Survival Probability",
           xlim = c(0, 3650),               
           ylim = c(0, 1),                  
           break.time.by = 365
)          

fit_DM_after <- survfit(Surv(Followup_time, DM_status) ~ PAM50, data = tcgadata,weights = weights_LCA4class)
C5_2 <- coxph(Surv(Followup_time, DM_status) ~ PAM50, data = tcgadata,weights = weights_LCA4class)
S5_2 <- summary(C5_2)
df5_2 <- paste0(round(S5_2$conf.int,2)[,1],"(",round(S5_2$conf.int,2)[,3],"-",round(S5_2$conf.int,2)[,4],")")
ggsurvplot(fit_DM_after, 
           data = tcgadata,
           pval = TRUE,                      
           risk.table = TRUE,               
           risk.table.title="",
           risk.table.y.text = FALSE,           
           conf.int = FALSE,                
           palette = c("#B6EE97","#68A7F8","#EB7D62","#7788E6","#F2D19E"), 
           legend.labs = c("Basal","Her2","LumA","LumB","Normal"),
           title = "",
           xlab = "Time(year)",
           ylab = "Survival Probability",
           xlim = c(0, 3650),               
           ylim = c(0, 1),                  
           break.time.by = 365
)          
#我需要把df1_1,df6_2,df2_1,df2_2,df3_1,df3_2,df4_1,df4_2,df5_1,df5_2这十个变量的值写入一个表格中，每一个变量都是一列，这十个变量都是1*1的矩阵，
df5x4 <- data.frame(cbind(df1_1,df1_2,df2_1,df2_2,df4_1,df4_2,df5_1,df5_2))
df5x4
pvalue <- data.frame(cbind(S1_1$robscore['pvalue'],S1_2$robscore['pvalue'],S2_1$robscore['pvalue'],S2_2$robscore['pvalue'],S4_1$robscore['pvalue'],S4_2$robscore['pvalue'],S5_1$robscore['pvalue'],S5_2$robscore['pvalue']))





##################第三部分：10Gene################################
# 4.2读取文件，设置用于LCA的formula，处理各变量为因子类型，准备部分变量供存放结果。
Tengene <- c("BAG1", "BCL2", "BIRC5", "CCNB1", "ERBB2", "ESR1", "GRB7", "MKI67", "MYBL2", "PGR")
formula_g <- as.formula(paste("cbind(", paste(Tengene, collapse = ", "), ") ~ 1"))
results_g <- vector("list", 6) 
lmr_p_values_g <- min_averages_g <- numeric(6)
set.seed(123)
# 4.3开始进行LCA循环，并展示结果
i <- 1
j <- 1
for (i in 1:6) {
  lca_result <- poLCA(formula_g, data = tcgadata, nclass = i, maxiter = 10000,nrep =30, calc.se = FALSE,verbose= FALSE)
  results_g[[i]] <- list(AIC = lca_result$aic, BIC = lca_result$bic, SABIC = lca_result$bic - log(nrow(tcgadata)) * (i - 1) * length(lca_result$probs[[1]]), Entropy = poLCA.entropy(lca_result), Min_pro =round(min(lca_result$P) ,3), llik= lca_result$llik,npar=lca_result$npar)# 对1到6个潜在类别进行LCA
  if(i > 1) {
    combined_result <- as.data.frame(cbind(lca_result$posterior, predclass = lca_result$predclass))# 合并posterior和predclass，并转换为数据框
    class_counts <- table(combined_result$predclass)# 新增 - 统计每个类别的数量
    results_g[[i]]$class_counts <- class_counts
    datasets <- split(combined_result, combined_result$predclass)# 根据predclass分割数据
    class_averages <- numeric(length(datasets))# 初始化向量存储每个类别的平均值
    for (j in 1:length(datasets)) {
      class_averages[j] <- mean(datasets[[j]][, j])
    }# 计算每个类别的平均值
    min_averages_g[i] <- round(min(class_averages),3)# 找到最小的平均值并存储
    n_obs <- nrow(tcgadata)# 计算LMR检验的p值
    lmr_result <- calc_lrt(n_obs, results_g[[i-1]]$llik, results_g[[i-1]]$npar, i-1, lca_result$llik,lca_result$npar, i)
    lmr_p_values_g[i] <- round(lmr_result[4],3)# LMR检验的p值是结果的第四个元素
  } else {
    min_averages_g[i] <- 1
    lmr_p_values_g[i] <- NA
  }
}
results_df_g <- map_dfr(lapply(results_g, function(x) {x$class_counts <- NULL
return(x)}), ~as.data.frame(t(unlist(.x))))
results_df_g$llik <- results_df_g$npar <- NULL
results_df_g$lmr_p_values <- lmr_p_values_g
results_df_g$min_averages <- min_averages_g
results_df_g

#4.4选定了最佳类别之后，再单独跑一次LCA以提取分类的信息。
LCAgene3class <- poLCA(formula_g, data = tcgadata, nclass = 3, maxiter = 10000, nrep=30, calc.se = TRUE)
tcgadata$LCA3class <- LCAgene3class$predclass
#####三、经【边际方法】调整后的性能信息###################
tcgadata$LCA3class <- as.factor(tcgadata$LCA3class)
tcgadata$LCA3class <- factor(tcgadata$LCA3class, levels = sort(unique(tcgadata$LCA3class)))
fit_LCA3class <- coxph(Surv(OS_time, OS_status) ~ LCA3class, data = tcgadata,robust = TRUE)
fit_LCA3class_after <- coxph(Surv(OS_time, OS_status) ~ LCA3class, data = tcgadata, weights = weights_LCA4class, robust = TRUE)
summary(fit_LCA3class)
summary(fit_LCA3class_after)

#分别把每个类别作为参考组
levels(tcgadata$LCA3class)
# tcgadata$LCA3class <- relevel(tcgadata$LCA3class, ref = "2")
fit_LCA3class <- coxph(Surv(OS_time, OS_status) ~ LCA3class, data = tcgadata,robust = TRUE)# 运行Cox模型
fit_LCA3class_after <- coxph(Surv(OS_time, OS_status) ~ LCA3class, data = tcgadata, weights = weights_LCA4class, robust = TRUE)# 使用权重拟合鲁棒性Cox模型
summary(fit_LCA3class)# 查看模型摘要，包括鲁棒性得分
summary(fit_LCA3class_after)

fit_LCA3class_total <- survfit(Surv(OS_time, OS_status) ~ LCA3class, data = tcgadata,robust =TRUE)
fit_LCA3class_total_after <- survfit(Surv(OS_time, OS_status) ~ LCA3class, data = tcgadata, weights = weights_LCA4class, robust = TRUE)
# 使用ggsurvplot绘制生存曲线
ggsurvplot(fit_LCA3class_total, 
           data = tcgadata,
           pval = TRUE,                      
           risk.table = TRUE,               
           conf.int = FALSE,                
           palette = c("#F7A882", "#EB7D62", "#68A7F8"), # 自定义颜色
           title = "Survival curves for NewClass Categories",
           xlab = "Time",
           ylab = "Survival Probability")

ggsurvplot(fit_LCA3class_total_after, 
           data = tcgadata,
           pval = FALSE,                      
           risk.table = TRUE,               
           conf.int = FALSE,                
           palette = c("#F7A882", "#EB7D62", "#68A7F8"), # 自定义颜色，根据分类数量调整
           title = "Survival curves for Newclass-asjusted Categories",
           xlab = "Time",
           ylab = "Survival Probability")


#执行两两组合的Cox分析（3是未加权重，4是增加了权重）
cox_results_df_3 <- cox_results_df_4 <- data.frame(C_Index = numeric(), Robustness = numeric(), P_Value = numeric(), row.names = character())

# 定义一个函数来执行两两组合的Cox分析
perform_cox_analysis_3 <- function(class1, class2, data) {
  subset_data <- data[data$LCA3class %in% c(class1, class2), ]
  cox_model <- coxph(Surv(OS_time, OS_status) ~ LCA3class, data = subset_data, robust = TRUE)
  return(summary(cox_model))
}
class_LCA3 <- unique(tcgadata$LCA3class)# 获取LCAclass的不同类别并执行所有可能的两两组合的Cox分析
cox_results_3 <- list()
combinations_3 <- combn(class_LCA3, 2, simplify = FALSE) # 生成所有可能的两两组合
for (comb in combinations_3) {
  cox_results_3[[paste(comb[1], "vs", comb[2])]] <- perform_cox_analysis_3(comb[1], comb[2], tcgadata)
}

# 遍历cox_results列表并提取所需信息
for (comb_name in names(cox_results_3)) {
  c_index <- cox_results_3[[comb_name]][["concordance"]][["C"]]
  robustness_score <- cox_results_3[[comb_name]][["robscore"]][["test"]]
  p_value <- cox_results_3[[comb_name]][["robscore"]][["pvalue"]]
  cox_results_df_3[comb_name, ] <- c(C_Index = c_index, Robustness = robustness_score, P_Value = p_value)
}

# 定义一个函数来执行两两组合的Cox分析
perform_cox_analysis_4 <- function(class1, class2, data) {
  subset_data <- data[data$LCA3class %in% c(class1, class2), ]
  cox_model <- coxph(Surv(OS_time, OS_status) ~ LCA3class, data = subset_data, weights = subset_data$weights_LCA4class, robust = TRUE)
  return(summary(cox_model))
}
class_LCA3 <- unique(tcgadata$LCA3class)
cox_results_4 <- list()
combinations_4 <- combn(class_LCA3, 2, simplify = FALSE) # 生成所有可能的两两组合
for (comb in combinations_4) {
  cox_results_4[[paste(comb[1], "vs", comb[2])]] <- perform_cox_analysis_4(comb[1], comb[2], tcgadata)
}
# 遍历cox_results列表并提取所需信息
for (comb_name in names(cox_results_4)) {
  c_index <- cox_results_4[[comb_name]][["concordance"]][["C"]]
  robustness_score <- cox_results_4[[comb_name]][["robscore"]][["test"]]
  p_value <- cox_results_4[[comb_name]][["robscore"]][["pvalue"]]
  cox_results_df_4[comb_name, ] <- c(C_Index = c_index, Robustness = robustness_score, P_Value = p_value)
}

# 查看结果
cox_results_df_3
cox_results_df_4

#######################三、使用ggsurvplot绘制生存曲线###############################
fit_gene <- survfit(Surv(OS_time, OS_status) ~ LCA3class, data = tcgadata,robust = TRUE)
C6_1 <- coxph(Surv(OS_time, OS_status) ~ LCA3class, data = tcgadata)
S6_1 <- summary(C6_1)
df6_1 <- paste0(round(S6_1$conf.int,2)[,1],"(",round(S6_1$conf.int,2)[,3],"-",round(S6_1$conf.int,2)[,4],")")
ggsurvplot(fit_gene, 
           data = tcgadata,
           pval = TRUE,                      
           risk.table = TRUE,               
           risk.table.title="",
           risk.table.y.text = FALSE,           
           conf.int = FALSE,                
           palette = c("#B6EE97","#68A7F8","#EB7D62","#7788E6","#F2D19E"), 
           legend.labs = c(1,2,3),
           title = "",
           xlab = "Time(year)",
           ylab = "Survival Probability",
           xlim = c(0, 3650),               
           ylim = c(0, 1),                  
           break.time.by = 365
)          

# 使用ggsurvplot绘制基于LCA4class分类的生存曲线
fit_gene_after <- survfit(Surv(OS_time, OS_status) ~ LCA3class, data = tcgadata,weights = weights_LCA4class)
C6_2 <- coxph(Surv(OS_time, OS_status) ~ LCA3class, data = tcgadata,weights = weights_LCA4class)
S6_2 <- summary(C6_2)
df6_2 <- paste0(round(S6_2$conf.int,2)[,1],"(",round(S6_2$conf.int,2)[,3],"-",round(S6_2$conf.int,2)[,4],")")
#在ggsurvplot中添加图例需要使用legend.labs参数，但是只有这个参数是不够的，还需要使用palette参数来指定颜色，否则会报错
ggsurvplot(fit_gene_after, 
           data = tcgadata,
           pval = FALSE,                      
           risk.table = TRUE,               
           risk.table.title="",
           risk.table.y.text = FALSE,           
           conf.int = FALSE,                
           palette = c("#B6EE97","#68A7F8","#EB7D62","#7788E6","#F2D19E"), 
           legend.labs = c(1,2,3),
           title = "",
           xlab = "Time(year)",
           ylab = "Survival Probability",
           xlim = c(0, 3650),               
           ylim = c(0, 1),                  
           break.time.by = 365
) 




# 
# seerCSSdata$RT <- as.factor(seerCSSdata$RT)
# seerCSSdata$RT <- relevel(seerCSSdata$RT, ref = "0")

HR6_1 <- round(S6_1$conf.int[,"exp(coef)"], 2)
CI_lower6_1 <- round(S6_1$conf.int[,"lower .95"], 2)
CI_upper6_1 <- round(S6_1$conf.int[,"upper .95"], 2)
pval6_1 <- round(S6_1$coefficients[,"Pr(>|z|)"], 3)
results6_1<- data.frame(HR=HR6_1, CI=paste(CI_lower6_1, CI_upper6_1, sep="-"), pval=pval6_1)
results6_1
HR6_2 <- round(S6_2$conf.int[,"exp(coef)"], 2)
CI_lower6_2 <- round(S6_2$conf.int[,"lower .95"], 2)
CI_upper6_2 <- round(S6_2$conf.int[,"upper .95"], 2)
pval6_2 <- round(S6_2$coefficients[,"Pr(>|z|)"], 3)
results6_2<- data.frame(HR=HR6_2, CI=paste(CI_lower6_2, CI_upper6_2, sep="-"), pval=pval6_2)
results6_2


tcgaclass1 <- tcgadata[tcgadata$LCA3class == 1,]
tcgaclass2 <- tcgadata[tcgadata$LCA3class == 2,]
tcgaclass3 <- tcgadata[tcgadata$LCA3class == 3,]

c_t1 <- coxph(Surv(OS_time, OS_status) ~ OS_status, data = tcgaclass1)
c_t2 <- coxph(Surv(OS_time, OS_status) ~ OS_status, data = tcgaclass2)
c_t3 <- coxph(Surv(OS_time, OS_status) ~ OS_status, data = tcgaclass3)

fit_DM_after <- survfit(Surv(Followup_time, DM_status) ~ PAM50, data = tcgadata,weights = weights_LCA4class)
C5_2 <- coxph(Surv(Followup_time, DM_status) ~ PAM50, data = tcgadata,weights = weights_LCA4class)
S5_2 <- summary(C5_2)
df5_2 <- paste0(round(S5_2$conf.int,2)[,1],"(",round(S5_2$conf.int,2)[,3],"-",round(S5_2$conf.int,2)[,4],")")







reference_tcga <- read.csv("reference_tcga.csv")
tcgadata <- merge(tcgadata, reference_tcga, by = "Seq")


tcgadata$LCA3class <- factor(tcgadata$LCA3class, levels = sort(unique(tcgadata$LCA3class)))
tcgadata$Race <- as.factor(tcgadata$Race)
tcgadata$Age_class <- as.factor(tcgadata$Age_class)
tcgadata$Stage <- as.factor(tcgadata$Stage)
tcgadata$LCA3class <- as.integer(tcgadata$LCA3class)
tcgadata$Mode_of_therapy <- factor(tcgadata$Mode_of_therapy, levels = sort(unique(tcgadata$Mode_of_therapy)))

#绘制森林图数据
forest_age <- data.frame(LCA3class = c(1, 2, 3), Age1 = c(0, 0, 0), Age2 = c(0, 0, 0), Age3 = c(0, 0, 0), Age4 = c(0, 0, 0))
for (i in 1:3) {
  forest_age$Age1[i] <- sum(tcgadata$LCA3class == i & tcgadata$Age_class == 1)
  forest_age$Age2[i] <- sum(tcgadata$LCA3class == i & tcgadata$Age_class == 2)
  forest_age$Age3[i] <- sum(tcgadata$LCA3class == i & tcgadata$Age_class == 3)
  forest_age$Age4[i] <- sum(tcgadata$LCA3class == i & tcgadata$Age_class == 4)
}

forest_stage <- data.frame(LCA3class = c(1, 2, 3), Stage1 = c(0, 0, 0), Stage2 = c(0, 0, 0), Stage3 = c(0, 0, 0), Stage4 = c(0, 0, 0))
for (i in 1:3) {
  forest_stage$Stage1[i] <- sum(tcgadata$LCA3class == i & tcgadata$Stage == 1)
  forest_stage$Stage2[i] <- sum(tcgadata$LCA3class == i & tcgadata$Stage == 2)
  forest_stage$Stage3[i] <- sum(tcgadata$LCA3class == i & tcgadata$Stage == 3)
  forest_stage$Stage4[i] <- sum(tcgadata$LCA3class == i & tcgadata$Stage == 4)
}

forest_Pathology <- data.frame(LCA3class = c(1, 2, 3), Pathology1 = c(0, 0, 0), Pathology2 = c(0, 0, 0), Pathology3 = c(0, 0, 0), Pathology4 = c(0, 0, 0), Pathology5 = c(0, 0, 0))
for (i in 1:3) {
  forest_Pathology$Pathology1[i] <- sum(tcgadata$LCA3class == i & tcgadata$Pathology == 1)
  forest_Pathology$Pathology2[i] <- sum(tcgadata$LCA3class == i & tcgadata$Pathology == 2)
  forest_Pathology$Pathology3[i] <- sum(tcgadata$LCA3class == i & tcgadata$Pathology == 3)
  forest_Pathology$Pathology4[i] <- sum(tcgadata$LCA3class == i & tcgadata$Pathology == 4)
  forest_Pathology$Pathology5[i] <- sum(tcgadata$LCA3class == i & tcgadata$Pathology == 5)
}

forest_Race <- data.frame(LCA3class = c(1, 2, 3), Race1 = c(0, 0, 0), Race2 = c(0, 0, 0), Race3 = c(0, 0, 0), Race4 = c(0, 0, 0))
for (i in 1:3) {
  forest_Race$Race1[i] <- sum(tcgadata$LCA3class == i & tcgadata$Race == 1)
  forest_Race$Race2[i] <- sum(tcgadata$LCA3class == i & tcgadata$Race == 2)
  forest_Race$Race3[i] <- sum(tcgadata$LCA3class == i & tcgadata$Race == 3)
  forest_Race$Race4[i] <- sum(tcgadata$LCA3class == i & tcgadata$Race == 4)
}

# tcgadata$LCA3class <- as.character(tcgadata$LCA3class)
# tcgadata$LCA3class[tcgadata$LCA3class == "3"] <- "c"
# tcgadata$LCA3class[tcgadata$LCA3class == "4"] <- "d"
# tcgadata$LCA3class[tcgadata$LCA3class == "c"] <- "4"
# tcgadata$LCA3class[tcgadata$LCA3class == "d"] <- "3"
# tcgadata$LCA3class <- as.integer(tcgadata$LCA3class)
# 
# #接下来，我需要把tcgadata的RT列中的0和2统一替换为0，然后把1,3,4,5，统一替换为1.
# tcgadata$RT[tcgadata$RT == 0] <- 0
# tcgadata$RT[tcgadata$RT == 2] <- 0
# tcgadata$RT[tcgadata$RT == 1] <- 1
# tcgadata$RT[tcgadata$RT == 3] <- 1
# tcgadata$RT[tcgadata$RT == 4] <- 1
# tcgadata$RT[tcgadata$RT == 5] <- 1
# tcgadata$RT <- as.factor(tcgadata$RT)
# 
# seerstage1 <- tcgadata[tcgadata$Stage == 1,]
# seerstage2 <- tcgadata[tcgadata$Stage == 2,]
# seerstage3 <- tcgadata[tcgadata$Stage == 3,]
# seerstage4 <- tcgadata[tcgadata$Stage == 4,]
# 
# cox_seer_stage1 <- coxph(Surv(OS_time, OS_status) ~ RT+LCA3class, data = seerstage1)
# cox_seer_stage2 <- coxph(Surv(OS_time, OS_status) ~ RT+LCA3class, data = seerstage2)
# cox_seer_stage3 <- coxph(Surv(OS_time, OS_status) ~ RT+LCA3class, data = seerstage3)
# cox_seer_stage4 <- coxph(Surv(OS_time, OS_status) ~ RT+LCA3class, data = seerstage4)
# summary(cox_seer_stage1)
# summary(cox_seer_stage2)
# summary(cox_seer_stage3)
# summary(cox_seer_stage4)

#我想计算每个Stage下、每个LCA3class下、每个Mode_of_therapy分类的数量，统计在一个表格中count_MT里面，前三列分别按照Stage、LCA4class、Mode_of_therapy分类，后三列分别是Mode_of_therapy为0、1、2、3、4的数量。
library(dplyr)

count_MT <- tcgadata %>%
  group_by(Stage, LCA3class, Mode_of_therapy) %>%
  summarise(
    Count_MT0 = sum(Mode_of_therapy == 0, na.rm = TRUE),
    Count_MT1 = sum(Mode_of_therapy == 1, na.rm = TRUE),
    Count_MT2 = sum(Mode_of_therapy == 2, na.rm = TRUE),
    Count_MT3 = sum(Mode_of_therapy == 3, na.rm = TRUE),
    Count_MT4 = sum(Mode_of_therapy == 4, na.rm = TRUE)
  ) %>%
  ungroup()


tcgadata$Mode_of_therapy <- relevel(tcgadata$Mode_of_therapy, ref = "0")

tcga_cox_results <- data.frame(Stage=integer(), Age_class=integer(), LCA3class=integer(), RT0_count=integer(), RT1_count=integer(),TotalCount=integer(), HR=character(), CI=character(), pval=character())

for(stage in unique(tcgadata$Stage)) {
  for(lca in unique(tcgadata$LCA3class)) {
      subset_data <- tcgadata[tcgadata$Stage == stage & tcgadata$LCA3class == lca, ]
      M0_count <- sum(subset_data$Mode_of_therapy == "0")
      M1_count <- sum(subset_data$Mode_of_therapy == "1")
      M2_count <- sum(subset_data$Mode_of_therapy == "2")
      M3_count <- sum(subset_data$Mode_of_therapy == "3")
      M4_count <- sum(subset_data$Mode_of_therapy == "4")
      TotalCount <- nrow(subset_data)
      if(TotalCount > 0){
        cox_model <- coxph(Surv(OS_time, OS_status) ~ Mode_of_therapy, data = subset_data)
        summary_cox <- summary(cox_model)
        
        HR <- round(summary_cox$conf.int[,"exp(coef)"], 2)
        CI_lower <- round(summary_cox$conf.int[,"lower .95"], 2)
        CI_upper <- round(summary_cox$conf.int[,"upper .95"], 2)
        pval <- round(summary_cox$coefficients[,"Pr(>|z|)"], 3)
      } else {
        HR <- CI_lower <- CI_upper <- pval <- "—"
      }
      tcga_cox_results <- rbind(tcga_cox_results, data.frame(Stage=stage, LCA3class=lca, M0=M0_count, M1=M1_count,M2=M2_count,M3=M3_count,M4=M4_count,TotalCount=TotalCount, HR=HR, CI=paste(CI_lower, CI_upper, sep="-"), pval=pval))
  }
}

tcgadata$Radiation_Therapy <- factor(tcgadata$Radiation_Therapy, levels = sort(unique(tcgadata$Radiation_Therapy)))
# tcgadata$Radiation_Therapy <- relevel(tcgadata$Radiation_Therapy, ref = "0")

tcgaRT_cox_results <- data.frame(Stage=integer(),  LCA3class=integer(), RT0_count=integer(), RT1_count=integer(), HR=character(), CI=character(), pval=character())

for(stage in unique(tcgadata$Stage)) {
  for(lca in unique(tcgadata$LCA3class)) {
    subset_data <- tcgadata[tcgadata$Stage == stage & tcgadata$LCA3class == lca, ]
    RT0_count <- sum(subset_data$Radiation_Therapy == "0")
    RT1_count <- sum(subset_data$Radiation_Therapy == "1")
    if(TotalCount > 0){
      cox_model <- coxph(Surv(OS_time, OS_status) ~ Radiation_Therapy, data = subset_data)
      summary_cox <- summary(cox_model)
      
      HR <- round(summary_cox$conf.int[,"exp(coef)"], 2)
      CI_lower <- round(summary_cox$conf.int[,"lower .95"], 2)
      CI_upper <- round(summary_cox$conf.int[,"upper .95"], 2)
      pval <- round(summary_cox$coefficients[,"Pr(>|z|)"], 3)
    } else {
      HR <- CI_lower <- CI_upper <- pval <- "—"
    }
    tcgaRT_cox_results <- rbind(tcgaRT_cox_results, data.frame(Stage=stage, LCA3class=lca, RT0=RT0_count, RT1=RT1_count, HR=HR, CI=paste(CI_lower, CI_upper, sep="-"), pval=pval))
  }
}




tcgadata$Pharmaceutical_Therapy <- factor(tcgadata$Pharmaceutical_Therapy, levels = sort(unique(tcgadata$Pharmaceutical_Therapy)))
# tcgadata$Pharmaceutical_Therapy <- relevel(tcgadata$Pharmaceutical_Therapy, ref = "0")

tcgaCH_cox_results <- data.frame(Stage=integer(),  LCA3class=integer(), CH0_count=integer(), CH1_count=integer(), HR=character(), CI=character(), pval=character())

for(stage in unique(tcgadata$Stage)) {
  for(lca in unique(tcgadata$LCA3class)) {
    subset_data <- tcgadata[tcgadata$Stage == stage & tcgadata$LCA3class == lca, ]
    CH0_count <- sum(subset_data$Pharmaceutical_Therapy == "0")
    CH1_count <- sum(subset_data$Pharmaceutical_Therapy == "1")
    if(TotalCount > 0){
      cox_model <- coxph(Surv(OS_time, OS_status) ~ Pharmaceutical_Therapy, data = subset_data)
      summary_cox <- summary(cox_model)
      

      # 创建一个空列表来存储模型
      cox_tengene <- list()
      
      for (i in 1:3) {
        for (j in Tengene) {
          formula <- as.formula(paste("Surv(OS_time, OS_status) ~", j))
          cox_model <- coxph(formula, data = subset(tcgadata, LCA3class == i))
          cox_tengene[[paste("LCA3class", i, "Gene", j)]] <- cox_model
        }
      }
      
      # 初始化森林图数据的数据框
      forest_plot_data_tengene <- data.frame()
      
      # 遍历每个 Cox 模型
      for(model_name in names(cox_tengene)) {
        model <- cox_tengene[[model_name]]
        if (!is.null(model)) {
          summary_cox <- summary(model)
          
          # 提取风险比 (HR) 和置信区间
          HR <- round(summary_cox$conf.int[,"exp(coef)"], 2)
          CI_lower <- round(summary_cox$conf.int[,"lower .95"], 2)
          CI_upper <- round(summary_cox$conf.int[,"upper .95"], 2)
          pval <- round(summary_cox$coefficients[,"Pr(>|z|)"], 3)
          
          # 将结果添加到数据框
          forest_plot_data_tengene <- rbind(forest_plot_data_tengene, data.frame(
            model = model_name,
            HR = HR,
            lower_ci = CI_lower,
            upper_ci = CI_upper,
            p_value = pval
          ))
        }
      }
      
      # 绘制森林图
      ggplot(subset(forest_plot_data_tengene,class == 1), aes(x = HR, y = model1)) +
        geom_point() +
        geom_errorbarh(aes(xmin = pmax(lower_ci, 0.1), xmax = pmin(upper_ci, 10)), height = 0.3) +
        geom_segment(data = subset(forest_plot_data_tengene, lower_ci < 0.1),
                     aes(x = 0.1, xend = 0.099, y = model1, yend = model1),
                     arrow = arrow(type = "closed", length = unit(0.1, "inches"))) +
        geom_segment(data = subset(forest_plot_data_tengene, upper_ci > 10),
                     aes(x = 10, xend = 11, y = model1, yend = model1),
                     arrow = arrow(type = "closed", length = unit(0.1, "inches"))) +
        scale_x_log10(limits = c(0.1, 10)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        theme_minimal() +
        theme(panel.grid.major = element_blank(),  # 移除主要网格线
              panel.grid.minor = element_blank()) + # 移除次要网格线
        xlab("Hazard Ratio (HR)") +
        ylab("Model Identifier")      
      

      
      ggplot(subset(forest_plot_data_tengene,class == 2), aes(x = HR, y = model1)) +
        geom_point() +
        geom_errorbarh(aes(xmin = pmax(lower_ci, 0.1), xmax = pmin(upper_ci, 10)), height = 0.3) +
        geom_segment(data = subset(forest_plot_data_tengene, lower_ci < 0.1),
                     aes(x = 0.1, xend = 0.099, y = model1, yend = model1),
                     arrow = arrow(type = "closed", length = unit(0.1, "inches"))) +
        geom_segment(data = subset(forest_plot_data_tengene, upper_ci > 10),
                     aes(x = 10, xend = 11, y = model1, yend = model1),
                     arrow = arrow(type = "closed", length = unit(0.1, "inches"))) +
        scale_x_log10(limits = c(0.1, 10)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        theme_minimal() +
        theme(panel.grid.major = element_blank(),  # 移除主要网格线
              panel.grid.minor = element_blank()) + # 移除次要网格线
        xlab("Hazard Ratio (HR)") +
        ylab("Model Identifier")      

      
      
      ggplot(subset(forest_plot_data_tengene,class == 3), aes(x = HR, y = model1)) +
        geom_point() +
        geom_errorbarh(aes(xmin = pmax(lower_ci, 0.1), xmax = pmin(upper_ci, 10)), height = 0.3) +
        geom_segment(data = subset(forest_plot_data_tengene, lower_ci < 0.1),
                     aes(x = 0.1, xend = 0.099, y = model1, yend = model1),
                     arrow = arrow(type = "closed", length = unit(0.1, "inches"))) +
        geom_segment(data = subset(forest_plot_data_tengene, upper_ci > 10),
                     aes(x = 10, xend = 11, y = model1, yend = model1),
                     arrow = arrow(type = "closed", length = unit(0.1, "inches"))) +
        scale_x_log10(limits = c(0.1, 10)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        theme_minimal() +
        theme(panel.grid.major = element_blank(),  # 移除主要网格线
              panel.grid.minor = element_blank()) + # 移除次要网格线
        xlab("Hazard Ratio (HR)") +
        ylab("Model Identifier")      

      



      
      
      
      
      
      
      
      
      
      HR <- round(summary_cox$conf.int[,"exp(coef)"], 2)
      CI_lower <- round(summary_cox$conf.int[,"lower .95"], 2)
      CI_upper <- round(summary_cox$conf.int[,"upper .95"], 2)
      pval <- round(summary_cox$coefficients[,"Pr(>|z|)"], 3)
    } else {
      HR <- CI_lower <- CI_upper <- pval <- "—"
    }
    tcgaCH_cox_results <- rbind(tcgaCH_cox_results, data.frame(Stage=stage, LCA3class=lca, CH0=CH0_count, CH1=CH1_count, HR=HR, CI=paste(CI_lower, CI_upper, sep="-"), pval=pval))
  }
}


for(i in c("Age_class","Pathology","Stage","Marry","Chemotherapy","Laterality","Primary_Site","T","N","M","RT")){
  print(table(seerCSSdata[,i],seerCSSdata$CSS))
}

for(i in c("Age_class","Pathology","Stage","Marry","Chemotherapy","Laterality","Primary_Site","T","N","M","RT")){
  print(table(seerCSSdata[,i],seerCSSdata$LCA4class))
}

for(i in c("Age_class","Pathology","Stage","Race","Pharmaceutical_Therapy","Radiation_Therapy","Mode_of_therapy","T","N","M")){
  print(table(tcgadata[,i],tcgadata$PAM50))
}

for(i in c("Age_class","Pathology","Race","Stage")){
  print(table(tcgadata[,i],tcgadata$LCA4class))
}

for(i in Tengene){# Tengene <- c("BAG1", "BCL2", "BIRC5", "CCNB1", "ERBB2", "ESR1", "GRB7", "MKI67", "MYBL2", "PGR")
  print(table(tcgadata[,i],tcgadata$LCA3class))
}


for(i in c("Age_class","Race","Stage","Pathology")){
  print(table(seerCSSdata[,i],seerCSSdata$LCA4class))
}

#绘制基因三分类的KM曲线

# 分割数据集
tcgadata_LCA3class1 <- subset(tcgadata, LCA3class == 1)
tcgadata_LCA3class2 <- subset(tcgadata, LCA3class == 2)
tcgadata_LCA3class3 <- subset(tcgadata, LCA3class == 3)

# 定义颜色
palette1 <- c("#039278", "#05AF90", "#69C6AD", "#AEDBCD")
palette2 <- c("#467FB8", "#5598DC", "#84B4EE", "#B8D0F4")
palette3 <- c("#C33B3B", "#E84747", "#F97D7D", "#FBB5B5")

# 绘制生存曲线图
# LCA3class = 1
fit_LCA3class1 <- survfit(Surv(OS_time, OS_status) ~ LCA4class, data = tcgadata_LCA3class1)
ggsurvplot(fit_LCA3class1, data = tcgadata_LCA3class1, pval = TRUE, risk.table = FALSE, conf.int = FALSE, 
           palette = palette1, title = "Survival curves for LCA4class with LCA3class = 1", 
           xlab = "Time", ylab = "Survival Probability", xlim = c(0, 1825), ylim = c(0, 1))

# LCA3class = 2
fit_LCA3class2 <- survfit(Surv(OS_time, OS_status) ~ LCA4class, data = tcgadata_LCA3class2)
ggsurvplot(fit_LCA3class2, data = tcgadata_LCA3class2, pval = TRUE, risk.table = FALSE, conf.int = FALSE, 
           palette = palette2, title = "Survival curves for LCA4class with LCA3class = 2", 
           xlab = "Time", ylab = "Survival Probability", xlim = c(0, 1825), ylim = c(0, 1))

# LCA3class = 3
fit_LCA3class3 <- survfit(Surv(OS_time, OS_status) ~ LCA4class, data = tcgadata_LCA3class3)
ggsurvplot(fit_LCA3class3, data = tcgadata_LCA3class3, pval = TRUE, risk.table = FALSE, conf.int = FALSE, 
           palette = palette3, title = "Survival curves for LCA4class with LCA3class = 3", 
           xlab = "Time", ylab = "Survival Probability", xlim = c(0, 1825), ylim = c(0, 1))



#绘制基因四分类的KM曲线
# 分割数据集
tcgadata_LCA4class1 <- subset(tcgadata, LCA4class == 1)
tcgadata_LCA4class2 <- subset(tcgadata, LCA4class == 2)
tcgadata_LCA4class3 <- subset(tcgadata, LCA4class == 3)
tcgadata_LCA4class4 <- subset(tcgadata, LCA4class == 4)

# 定义颜色
palette_LCA4class1 <- c("#039278", "#467FB8", "#C33B3B")
palette_LCA4class2 <- c("#05AF90", "#5598DC", "#E84747")
palette_LCA4class3 <- c("#69C6AD", "#84B4EE", "#F97D7D")
palette_LCA4class4 <- c("#AEDBCD", "#B8D0F4", "#FBB5B5")

# 绘制生存曲线图
# LCA4class = 1
fit_LCA4class1 <- survfit(Surv(OS_time, OS_status) ~ LCA3class, data = tcgadata_LCA4class1)
ggsurvplot(fit_LCA4class1, data = tcgadata_LCA4class1, pval = TRUE, risk.table = FALSE, conf.int = FALSE, 
           palette = palette_LCA4class1, title = "Survival curves for LCA3class with LCA4class = 1", 
           xlab = "Time", ylab = "Survival Probability", xlim = c(0, 1825), ylim = c(0, 1))

# LCA4class = 2
fit_LCA4class2 <- survfit(Surv(OS_time, OS_status) ~ LCA3class, data = tcgadata_LCA4class2)
ggsurvplot(fit_LCA4class2, data = tcgadata_LCA4class2, pval = TRUE, risk.table = FALSE, conf.int = FALSE, 
           palette = palette_LCA4class2, title = "Survival curves for LCA3class with LCA4class = 2", 
           xlab = "Time", ylab = "Survival Probability", xlim = c(0, 1825), ylim = c(0, 1))

# LCA4class = 3
fit_LCA4class3 <- survfit(Surv(OS_time, OS_status) ~ LCA3class, data = tcgadata_LCA4class3)
ggsurvplot(fit_LCA4class3, data = tcgadata_LCA4class3, pval = TRUE, risk.table = FALSE, conf.int = FALSE, 
           palette = palette_LCA4class3, title = "Survival curves for LCA3class with LCA4class = 3", 
           xlab = "Time", ylab = "Survival Probability", xlim = c(0, 1825), ylim = c(0, 1))

# LCA4class = 4
fit_LCA4class4 <- survfit(Surv(OS_time, OS_status) ~ LCA3class, data = tcgadata_LCA4class4)
ggsurvplot(fit_LCA4class4, data = tcgadata_LCA4class4, pval = TRUE, risk.table = FALSE, conf.int = FALSE, 
           palette = palette_LCA4class4, title = "Survival curves for LCA3class with LCA4class = 4", 
           xlab = "Time", ylab = "Survival Probability", xlim = c(0, 1825), ylim = c(0, 1))

# 提取生存曲线数据
extract_survival_data <- function(fit) {
  summary_fit <- summary(fit)
  data.frame(
    Time = summary_fit$time,
    SurvivalProbability = summary_fit$surv,
    RiskNumber = summary_fit$n.risk,
    EventsNumber = summary_fit$n.event,
    ConfIntLower = summary_fit$lower,
    ConfIntUpper = summary_fit$upper
  )
}

# 使用函数提取数据
survival_data_LCA3class1 <- extract_survival_data(fit_LCA3class1)
survival_data_LCA3class2 <- extract_survival_data(fit_LCA3class2)
survival_data_LCA3class3 <- extract_survival_data(fit_LCA3class3)
survival_data_LCA4class1 <- extract_survival_data(fit_LCA4class1)
survival_data_LCA4class2 <- extract_survival_data(fit_LCA4class2)
survival_data_LCA4class3 <- extract_survival_data(fit_LCA4class3)
survival_data_LCA4class4 <- extract_survival_data(fit_LCA4class4)
# 进行统计检验（如果需要）
logrank_test3_1 <- survdiff(Surv(OS_time, OS_status) ~ LCA4class, data = tcgadata_LCA3class1)
logrank_test3_2 <- survdiff(Surv(OS_time, OS_status) ~ LCA4class, data = tcgadata_LCA3class2)
logrank_test3_3 <- survdiff(Surv(OS_time, OS_status) ~ LCA4class, data = tcgadata_LCA3class3)
logrank_test4_1 <- survdiff(Surv(OS_time, OS_status) ~ LCA3class, data = tcgadata_LCA4class1)
logrank_test4_2 <- survdiff(Surv(OS_time, OS_status) ~ LCA3class, data = tcgadata_LCA4class2)
logrank_test4_3 <- survdiff(Surv(OS_time, OS_status) ~ LCA3class, data = tcgadata_LCA4class3)
logrank_test4_4 <- survdiff(Surv(OS_time, OS_status) ~ LCA3class, data = tcgadata_LCA4class4)

# 查看提取的数据
print(survival_data_LCA3class1)
print(survival_data_LCA3class2)
print(survival_data_LCA3class3)
print(survival_data_LCA4class1)
print(survival_data_LCA4class2)
print(survival_data_LCA4class3)
print(survival_data_LCA4class4)
print(logrank_test3_1)
print(logrank_test3_2)
print(logrank_test3_3)
print(logrank_test4_1)
print(logrank_test4_2)
print(logrank_test4_3)
print(logrank_test4_4)


#输出关键结果

results_df_s
head(seerdata_clean)
head(seerdata_clean1)
head(seerCSSdata)
sankydata
summary(fitCSS01)
summary(fitCSS02)
seer_cox_results
cox_age_result
cox_pathology_result
cox_race_result
cox_stage_result
results_df_t
summary(weights_LCA4class)
head(tcgadata)
summary(fit_PAM50)
summary(fit_PAM50_after)
summary(fit_PFS)
summary(fit_PFS_after)
summary(fit_LC)
summary(fit_LC_after)
summary(fit_DM)
summary(fit_DM_after)
cox_results_df
cox_results_df_2
results_df_g
Tengene
fit_gene
fit_gene_after
cox_results_df_3
cox_results_df_4
summary(fit_LCA3class)
summary(fit_LCA3class_after)
summary(fit_LCA3class_total)
summary(fit_LCA3class_total_after)
summary(fit_LCA4class)


#我需要从TCGAdata中提取PAM50=normal的数据出来，然后再进行生存分析，这些数据放在normal_case中
normal_case <- tcgadata[tcgadata$PAM50 == "Normal", ]
