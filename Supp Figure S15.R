##Script for Figure S15.

library(survival)
library(survminer)
library(dplyr)
library(readr)
library(pheatmap)

clinical_data <- read.delim("clinical.cases_selection.2019-05-01/clinical.tsv")

clinical_data <- clinical_data[,c("submitter_id", "age_at_diagnosis", "vital_status", "days_to_death", "days_to_last_follow_up", "tumor_stage")]

clinical_data$Patient <- substr(clinical_data$submitter_id, 9, 12)
clinical_data$submitter_id <- NULL


clinical_data$Censored <- ifelse(clinical_data$days_to_death != "--", 1, 0)
clinical_data$Survival_time <- clinical_data$days_to_death
clinical_data$Survival_time <- as.character(clinical_data$Survival_time)
clinical_data$Survival_time[clinical_data$Survival_time == "--"] <- as.character(clinical_data$days_to_last_follow_up[clinical_data$Survival_time == "--"])
clinical_data <- clinical_data[clinical_data$Survival_time != "--",]


load("ssGSEA_module_results.Rdata")

tumour <- "LGG"
local_ssGSEA <- ssGSEA_module_results[["tumour"]][[tumour]]
local_ssGSEA <- local_ssGSEA[rownames(local_ssGSEA) != "grey",]
local_ssGSEA_q <- t(apply(local_ssGSEA, 1, function(x){
  q <- quantile(x, c(1/2))
  y <- rep(NA, length(x))
  y[x < q[1]] <- "Bellow"
  y[x >= q[1]] <- "Above"
  return(y)
}))
colnames(local_ssGSEA_q) <- colnames(local_ssGSEA)

local_clinical_data <- clinical_data[clinical_data$Patient %in% colnames(local_ssGSEA_q),]

module <- "purple"

local_ssGSEA_q2 <- local_ssGSEA_q[rownames(local_ssGSEA_q) == module,]
Above <- names(local_ssGSEA_q2)[local_ssGSEA_q2 == "Above"]
Bellow <- names(local_ssGSEA_q2)[local_ssGSEA_q2 == "Bellow"]
Above <- Above[!is.na(Above)]
Bellow <- Bellow[!is.na(Bellow)]
local_clinical_data2 <- local_clinical_data
local_clinical_data2$Module_activity <- ifelse(local_clinical_data2$Patient %in% Above, "High",
                                               ifelse(local_clinical_data2$Patient %in% Bellow, "Low", NA))
local_clinical_data2 <- local_clinical_data2[!is.na(local_clinical_data2$Module_activity),]

surv_object <- Surv(time = as.numeric(local_clinical_data2$Survival_time), event = local_clinical_data2$Censored)
fit1 <- survfit(surv_object ~ Module_activity, data = local_clinical_data2)

pdf("Figure_S15.pdf", height=5, width=5)
g <- ggsurvplot(fit1, data = local_clinical_data2, pval = TRUE, 
                title=paste(tumour, module, sep=" - "))
print(g)

dev.off()