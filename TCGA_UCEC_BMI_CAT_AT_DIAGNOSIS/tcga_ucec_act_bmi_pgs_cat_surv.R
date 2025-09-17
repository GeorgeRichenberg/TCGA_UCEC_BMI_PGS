### Set working directory and load required packages:
setwd()
library(tidyverse)
library(survival)
library(survminer)

### Load required data:
clinical1 <- read_delim("TCGA.UCEC.sampleMap_UCEC_clinicalMatrix.txt", delim = "\t")
pgs <- read_delim("ucec_bmi_female_primary_std.txt", delim = "\t")
surv1 <- read_delim("UCEC_survival.txt", delim = "\t")

### Formatting clinical1:
colnames(clinical1)
which(duplicated(clinical1$`_PATIENT`))
clinical2 <- clinical1 %>% 
  distinct(`_PATIENT`, .keep_all = TRUE)
clinical3 <- clinical2 %>% 
  select(`_PATIENT`, CDE_ID_3226963, age_at_initial_pathologic_diagnosis, clinical_stage, height, weight)
clinical3$clinical_stage <- as_factor(clinical3$clinical_stage)
levels(clinical3$clinical_stage)
clinical3 <- clinical3 %>% 
  mutate(stage = recode(clinical_stage, "Stage I"="1", 
                        "Stage II"="2", 
                        "Stage III"="3", 
                        "Stage IV"="4", 
                        "Stage IA"="1", 
                        "Stage IB"="1", 
                        "Stage IC"="1", 
                        "Stage IIA"="2", 
                        "Stage IIB"="2",
                        "Stage IIC"="2",
                        "Stage IIIA"="3", 
                        "Stage IIIB"="3", 
                        "Stage IIIC"="3", 
                        "Stage IIIC1"="3", 
                        "Stage IIIC2"="3", 
                        "Stage IVA"="4",
                        "Stage IVB"="4"))
clinical3 <- clinical3 %>% 
  rename(tcgaid = `_PATIENT`, 
         age = age_at_initial_pathologic_diagnosis,
         ms_status = CDE_ID_3226963)
clinical3 <- clinical3 %>% 
  select(tcgaid, age, ms_status, stage, height, weight)

### Calculating BMI:
clinical3 <- clinical3 %>%
  mutate(height_m =  height / 100)
clinical3 <- clinical3 %>%
  mutate(height_squared =  height_m * height_m)
clinical3 <- clinical3 %>%
  mutate(bmi = weight/height_squared)

### Removing cases with BMI<15 and BMI>50:
clinical3 <- clinical3 %>%
  filter(bmi>15) %>%
  filter(bmi<50) %>%
  select(-c(height, weight, height_m, height_squared))
min(clinical3$bmi)
max(clinical3$bmi)

### Removing rows with missing and indeterminate microsatellite status:
dim(clinical3)
sum(is.na(clinical3$ms_status))
clinical4 <- clinical3 %>%
  drop_na(ms_status)
dim(clinical4)
clinical4$ms_status <- as_factor(clinical4$ms_status)
levels(clinical4$ms_status)
which(clinical4$ms_status == "Indeterminate")
clinical5 <- clinical4 %>%
  filter(ms_status != "Indeterminate")
dim(clinical5)

### Re-defining microsatellite status:
clinical5$ms_status <- as.character(clinical5$ms_status)
clinical6 <- clinical5 %>%
  mutate(ms_status = replace(ms_status, ms_status == "MSS", 0)) %>%
  mutate(ms_status = replace(ms_status, ms_status == "MSI-L", 1)) %>%
  mutate(ms_status = replace(ms_status, ms_status == "MSI-H", 1))
clinical6$ms_status <- as.numeric(clinical6$ms_status)

### Categorising BMI into five groups:
### (i) Underweight (BMI < 18.5) ~ 1
### (ii) Healthy weight (BMI 18.5 - 24.9) ~ 2
### (iii) Overweight (BMI 25.0 - 29.9) ~ 3
### (iv) Obese (BMI 30.0 - 39.9) ~ 4
### (v) Severely obese (BMI >= 40.0) ~ 5
clinical6 <- clinical6 %>%
  mutate(bmi = round(bmi, 1)) %>%
  mutate(bmi_cat = case_when(
    bmi < 18.5 ~ 1,
    between(bmi, 18.5, 24.9) ~ 1,
    between(bmi, 25.0, 29.9) ~ 2,
    between(bmi, 30.0, 39.9) ~ 3,
    bmi >= 40.0 ~ 4))
sum(is.na(clinical6$bmi_cat))
clinical6 <- clinical6 %>%
  select(-bmi)

### Format survival data:
dim(surv1)
surv1 <- surv1 %>% 
  arrange(sample)
surv2 <- surv1 %>%
  distinct(`_PATIENT`, .keep_all = TRUE)
dim(surv2)
surv3 <- surv2 %>%
  filter(is.na(Redaction))
dim(surv3)
surv3 <- surv3 %>%
  select(-c(sample, Redaction))
surv3 <- surv3 %>%
  rename(tcgaid = `_PATIENT`)

### Merge pgs, clinical3, ancestry2 and mutsig4 for multivariable analysis:
merge_surv_covars <- pgs %>% 
  inner_join(clinical6, by = "tcgaid") %>% 
  inner_join(surv3, by = "tcgaid")
merge_surv_covars <- merge_surv_covars %>%
  select(-pgs)
dim(merge_surv_covars)
merge_surv_covars %>%
  group_by(bmi_cat) %>%
  count()

### Run survival analysis:
os_results <- coxph(Surv(OS.time, OS) ~ bmi_cat + age + ms_status + stage, 
                    data = merge_surv_covars)
os_results2 <- summary(os_results)
os <- as.data.frame(os_results2$coefficients)
os$lower_95_ci <- os$'exp(coef)' - 1.96 * os$'se(coef)'
os$upper_95_ci <- os$'exp(coef)' + 1.96 * os$'se(coef)'
os <- os[, c(1:3, 6:7, 4:5)]
os <- os %>%
  rownames_to_column(var = "Variable") %>%
  dplyr::mutate("Survival endpoint" = "Overall survival") %>%
  select("Survival endpoint", everything())
os

dss_results <- coxph(Surv(DSS.time, DSS) ~ bmi_cat + age + ms_status + stage, 
                     data = merge_surv_covars)
dss_results2 <- summary(dss_results)
dss <- as.data.frame(dss_results2$coefficients)
dss$lower_95_ci <- dss$'exp(coef)' - 1.96 * dss$'se(coef)'
dss$upper_95_ci <- dss$'exp(coef)' + 1.96 * dss$'se(coef)'
dss <- dss[, c(1:3, 6:7, 4:5)]
dss <- dss %>%
  rownames_to_column(var = "Variable") %>%
  dplyr::mutate("Survival endpoint" = "Disease-specific survival") %>%
  select("Survival endpoint", everything())
dss

dfi_results <- coxph(Surv(DFI.time, DFI) ~ bmi_cat + age + ms_status + stage, 
                     data = merge_surv_covars)
dfi_results2 <- summary(dfi_results)
dfi <- as.data.frame(dfi_results2$coefficients)
dfi$lower_95_ci <- dfi$'exp(coef)' - 1.96 * dfi$'se(coef)'
dfi$upper_95_ci <- dfi$'exp(coef)' + 1.96 * dfi$'se(coef)'
dfi <- dfi[, c(1:3, 6:7, 4:5)]
dfi <- dfi %>%
  rownames_to_column(var = "Variable") %>%
  dplyr::mutate("Survival endpoint" = "Disease-free interval") %>%
  select("Survival endpoint", everything())
dfi

pfi_results <- coxph(Surv(PFI.time, PFI) ~ bmi_cat + age + ms_status + stage, 
                     data = merge_surv_covars)
pfi_results2 <- summary(pfi_results)
pfi <- as.data.frame(pfi_results2$coefficients)
pfi$lower_95_ci <- pfi$'exp(coef)' - 1.96 * pfi$'se(coef)'
pfi$upper_95_ci <- pfi$'exp(coef)' + 1.96 * pfi$'se(coef)'
pfi <- pfi[, c(1:3, 6:7, 4:5)]
pfi <- pfi %>%
  rownames_to_column(var = "Variable") %>%
  dplyr::mutate("Survival endpoint" = "Progression-free interval") %>%
  select("Survival endpoint", everything())
pfi

### Combine:
survival <- rbind(os, dss, dfi, pfi)
survival
write_csv(survival, file = "bmi_four_cat_survival.csv")

### Kaplan-Meier plot by BMI cat:
merge_surv_covars2 <- merge_surv_covars %>%
  mutate(bmi_cat = case_when(
    bmi_cat == 1 ~ "Underweight + Healthy weight",
    bmi_cat == 2 ~ "Overweight",
    bmi_cat == 3 ~ "Obese",
    bmi_cat == 4 ~ "Severely obese"))
merge_surv_covars2$bmi_cat <- factor(merge_surv_covars2$bmi_cat,
                                     levels = c("Underweight + Healthy weight",
                                                "Overweight",
                                                "Obese",
                                                "Severely obese"))
fit1 <- survfit(Surv(OS.time, OS) ~ bmi_cat, data = merge_surv_covars2)
kmplot <- ggsurvplot(fit1,
                     conf.int = FALSE,
                     risk.table = TRUE,
                     pval = TRUE,
                     tables.theme = theme_cleantable(),
                     legend = "bottom",
                     ggtheme = theme_minimal())
kmplot
pdf("act_bmi_pgs_cat_OS_curve.pdf")
print(kmplot, newpage = FALSE)
dev.off()