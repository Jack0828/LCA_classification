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
# setwd("C:\\....")

#######################Part.1 #############################
seerdata <- read.csv("LCASeer.csv")
# Remove any rows where Patient_ID is repeated twice or more
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
results_s1 <- results_s <- vector("list", 6) 
lmr_p_values_s <- min_averages_s <- numeric(6)  
# 1.2 Start LCA
set.seed(123)
i <- 1
j <- 1
for (i in 1:6) {
  lca_result <- poLCA(formula_s, data = seerdata_clean, nclass = i,maxiter = 15000, calc.se = FALSE)
  colname <- paste("LCAs", i, "class", sep = "") 
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

#1.3
LCASEER4class <- poLCA(formula_s, data = seerdata_clean, nclass = 4, maxiter = 15000,calc.se = FALSE)
seerdata_clean1$LCA4class <- LCASEER4class$predclass

seerCSSdata <- seerdata_clean1
reference_seer <- read.csv("reference.csv")
seerCSSdata <- merge(seerCSSdata, reference_seer, by = "Patient_ID")
reference_seer_after <- read.csv("reference_seer_after.csv")
seerCSSdata <- merge(seerCSSdata, reference_seer_after, by = "Patient_ID")


seerCSSdata$CSS <- as.factor(seerCSSdata$CSS)
seerCSSdata$LCA4class <- factor(seerCSSdata$LCA4class, levels = sort(unique(seerCSSdata$LCA4class)))
seerCSSdata$Race <- as.factor(seerCSSdata$Race.x)
seerCSSdata$Age_class <- as.factor(seerCSSdata$Age_class)
seerCSSdata$Stage <- as.factor(seerCSSdata$Stage)
fit_seer <- survfit(Surv(OS_time, OS_status) ~ LCA4class, data = seerCSSdata)
ggsurvplot(fit_seer, 
           data = seerCSSdata,
           pval = TRUE,                      
           risk.table = TRUE,
           risk.table.y.text = FALSE,
           strata = FALSE,
           conf.int = FALSE,                
           palette = c("#ffc000","#f84d4c","#8f7ef1","#00b0f1"), 
           title = "Survival curves for LCA 4 class",
           xlab = "Time",
           ylab = "Survival Probability")


CSS01 <- subset(seerCSSdata, CSS == 1 | CSS == 0)
CSS02 <- subset(seerCSSdata, CSS == 2 | CSS == 0)
fitCSS01 <- coxph(Surv(OS_time, OS_status) ~ LCA4class, data = CSS01, robust = TRUE)
fitCSS02 <- coxph(Surv(OS_time, OS_status) ~ LCA4class, data = CSS02, robust = TRUE)

# Sankey plot
sankydata <- data.frame(LCA4class = c(1, 2, 3, 4), CSS0 = c(0, 0, 0, 0), CSS1 = c(0, 0, 0, 0), CSS2 = c(0, 0, 0, 0))
for (i in 1:4) {
  sankydata$CSS0[i] <- sum(seerCSSdata$LCA4class == i & seerCSSdata$CSS == 0)
  sankydata$CSS1[i] <- sum(seerCSSdata$LCA4class == i & seerCSSdata$CSS == 1)
  sankydata$CSS2[i] <- sum(seerCSSdata$LCA4class == i & seerCSSdata$CSS == 2)
}

###seerCSSdata$RT standardization
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

###########  CIF plot  ###############
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
        panel.grid.minor = element_blank())  

###### Part.2 TCGA Database #########
tcgadata <- read.csv("LCATCGA.csv")
tcgadata$Race <- as.factor(tcgadata$Race)
tcgadata$Age_class <- as.factor(tcgadata$Age_class)
tcgadata$Pathology <- as.factor(tcgadata$Pathology)
tcgadata$Stage <- as.factor(tcgadata$Stage)
formula_t <- as.formula(paste("cbind(Race, Age_class, Pathology, Stage) ~ 1"))
results_t <- vector("list", 6)
lmr_p_values_t <- min_averages_t <- numeric(6) 
set.seed(123)
# 2.3 start LCA
i <- 1
j <- 1
for (i in 1:6) {
  lca_result <- poLCA(formula_t, data = tcgadata, nclass = i, maxiter = 15000, nrep = 30, calc.se = FALSE, verbose =FALSE,)
  results_t[[i]] <- list(AIC = lca_result$aic, BIC = lca_result$bic, SABIC = lca_result$bic - log(nrow(tcgadata)) * (i - 1) * length(lca_result$probs[[1]]), Entropy = poLCA.entropy(lca_result), Min_pro =round(min(lca_result$P),3),llik= lca_result$llik, npar=lca_result$npar,P=round(lca_result$P,3))
  if(i > 1) {
    combined_result <- as.data.frame(cbind(lca_result$posterior, predclass = lca_result$predclass))
    class_counts <- table(combined_result$predclass)
    results_t[[i]]$class_counts <- class_counts
    datasets <- split(combined_result, combined_result$predclass)
    class_averages <- numeric(length(datasets))
    for (j in 1:length(datasets)) {
      class_averages[j] <- mean(datasets[[j]][, j])
    }
    min_averages_t[i] <- round(min(class_averages),3)
    n_obs <- nrow(tcgadata)
    lmr_result <- calc_lrt(n_obs, results_t[[i-1]]$llik, results_t[[i-1]]$npar, i-1, lca_result$llik,lca_result$npar, i)
    lmr_p_values_t[i] <- lmr_result[4]
  } else {
    min_averages_t[i] <- 1
    lmr_p_values_t[i] <- NA
  }
}
results_df_t <- map_dfr(lapply(results_t, function(x) {x$class_counts <- NULL
return(x)}), ~as.data.frame(t(unlist(.x))))

reference_tcga_after <- read.csv("reference_tcga_after.csv")
tcgadata <- merge(tcgadata, reference_tcga_after, by = "Seq")
LCATCGA <- poLCA(formula_t, data = tcgadata, nclass = 4, maxiter = 10000, nrep=30,calc.se = FALSE)
tcgadata$LCA4class <- LCATCGA$predclass


# 5.1 start cox and survival analysis
tcgadata$OS_status <- as.numeric(tcgadata$OS_status) 
tcgadata$PAM50 <- as.factor(tcgadata$PAM50)
tcgadata$LCA4class <- as.factor(tcgadata$LCA4class)
tcgadata$PAM50 <- factor(tcgadata$PAM50)
tcgadata$PAM50 <- relevel(tcgadata$PAM50, ref = "Normal")
fit_total <- survfit(Surv(OS_time, OS_status) ~ PAM50, data = tcgadata)
fit_PAM50 <- coxph(Surv(OS_time, OS_status) ~ PAM50, data = tcgadata,robust = TRUE)
tcgadata$LCA4class <- as.factor(tcgadata$LCA4class)
fit_LCA4class <- coxph(Surv(OS_time, OS_status) ~ LCA4class, data = tcgadata,robust = TRUE)
tcgadata$LCA4class <- factor(tcgadata$LCA4class, levels = sort(unique(tcgadata$LCA4class)))

hazard_ratios_LCA4class <- exp(coef(fit_LCA4class))
hazard_ratios_LCA4class <- c(1, hazard_ratios_LCA4class) 
weights_LCA4class <- sapply(tcgadata$LCA4class, function(x) 1 / hazard_ratios_LCA4class[as.integer(x)])
                            
fit_PAM50_after <- coxph(Surv(OS_time, OS_status) ~ PAM50, data = tcgadata, weights = weights_LCA4class, robust = TRUE)
tcgadata$weights_LCA4class <- weights_LCA4class

###############Part.2.2 coxph model adjustment and comparision in different PAM50 subtype ##########################

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

#######################Part.2.2 Survival curves with different events ###############################
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
##### PFS #####
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
##### DFS #####
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
##### LC #####
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
##### DM ####
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
df_2rd <- data.frame(cbind(df1_1,df1_2,df2_1,df2_2,df4_1,df4_2,df5_1,df5_2))
pvalue_2rd <- data.frame(cbind(S1_1$robscore['pvalue'],S1_2$robscore['pvalue'],S2_1$robscore['pvalue'],S2_2$robscore['pvalue'],S4_1$robscore['pvalue'],S4_2$robscore['pvalue'],S5_1$robscore['pvalue'],S5_2$robscore['pvalue']))
                 
#data statistics
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

for(i in c("Age_class","Race","Stage","Pathology")){
  print(table(seerCSSdata[,i],seerCSSdata$LCA4class))
}
