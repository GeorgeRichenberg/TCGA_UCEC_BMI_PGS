### Load required packages:
library(tidyverse)
library(vroom)
library(broom)
library(DescTools)
library(pscl)
library(countreg)

### Load required data:
clinical1 <- read_delim("TCGA.UCEC.sampleMap_UCEC_clinicalMatrix.txt", delim = "\t")
mutsig1 <- vroom("TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv", delim = ",")
ancestry1 <- read_delim("WashU_PCA_ethnicity_assigned.tsv", delim = "\t")
pgs <- read_delim("ucec_bmi_female_primary_std.txt", delim = "\t")

### Formatting clinical1:
colnames(clinical1)
which(duplicated(clinical1$`_PATIENT`))
clinical2 <- clinical1 %>% 
  distinct(`_PATIENT`, .keep_all = TRUE)
clinical3 <- clinical2 %>% 
  dplyr::select(`_PATIENT`, CDE_ID_3226963, age_at_initial_pathologic_diagnosis, clinical_stage)
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
  dplyr::select(tcgaid, age, ms_status, stage)

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

### Format mutsig1:
mutsig2 <- mutsig1 %>% 
  filter(`Cancer Types` == "Uterus-AdenoCa")
mutsig3 <- mutsig2 %>%
  dplyr::select(-`Cancer Types`)
mutsig3 <- mutsig3 %>%
  rename(tcgaid = `Sample Names`)
mutsig3 <- mutsig3 %>%
  arrange(tcgaid)
mutsig3 <- mutsig3 %>%
  separate(tcgaid, into = c("tcgaid1","tcgaid2","tcgaid3","primary"), sep = "-")
levels(as_factor(mutsig3$primary))
mutsig3 <- mutsig3 %>%
  filter(str_detect(primary, "06A") != TRUE)
mutsig3 <- mutsig3 %>%
  dplyr::select(-primary)
mutsig3 <- mutsig3 %>%
  unite(tcgaid, tcgaid1, tcgaid2, tcgaid3, sep = "-")
which(duplicated(mutsig3$tcgaid))

### Restricting analysis to mutational signatures observed in at least 10% of endometrial samples:
zero_index_mutsig3 <- mutsig3 %>%
  summarise_all(list(~ sum(.==0)))
proportion_zeros <- colMeans(mutsig3 == 0)
proportion_zeros
which(zero_index_mutsig3<(0.90*(nrow(mutsig3))))
mutsig4 <- mutsig3 %>%
  dplyr::select(c(1, 2, 4, 15, 16, 19, 21, 47, 51))

### Format ancestry1:
which(duplicated(ancestry1$Case))
ancestry1 <- ancestry1 %>%
  arrange(Sample)
ancestry2 <- ancestry1 %>%
  distinct(Case, .keep_all = TRUE)
ancestry2 <- ancestry2 %>%
  dplyr::select(Case, PC1:PC10)
ancestry2 <- ancestry2 %>%
  rename(tcgaid = Case)

### Merge data for multivariable analysis:
merge_mutsig_covars <- pgs %>% 
  inner_join(clinical6, by = "tcgaid") %>% 
  inner_join(ancestry2, by = "tcgaid") %>% 
  inner_join(mutsig4, by = "tcgaid")
dim(merge_mutsig_covars)

### Run ZINB analysis:
sbs2_zinb <- zeroinfl(SBS2 ~ pgs + age + ms_status + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Accuracy,
                      data = merge_mutsig_covars, 
                      dist = "negbin")
summary_sbs2_zinb <- summary(sbs2_zinb)
sbs2_zinb_results <- coef(summary_sbs2_zinb)
sbs2_zinb_results_count <- sbs2_zinb_results$count
sbs2_zinb_results_zero <- sbs2_zinb_results$zero
sbs2_results <- rbind(sbs2_zinb_results_count, sbs2_zinb_results_zero)
write.csv(sbs2_results, file = "sbs2_zinb.csv")

sbs10a_zinb <- zeroinfl(SBS10a ~ pgs + age + ms_status + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Accuracy,
                        data = merge_mutsig_covars, 
                        dist = "negbin")
summary_sbs10a_zinb <- summary(sbs10a_zinb)
sbs10a_zinb_results <- coef(summary_sbs10a_zinb)
sbs10a_zinb_results_count <- sbs10a_zinb_results$count
sbs10a_zinb_results_zero <- sbs10a_zinb_results$zero
sbs10a_results <- rbind(sbs10a_zinb_results_count, sbs10a_zinb_results_zero)
write.csv(sbs10a_results, file = "sbs10a_zinb.csv")

sbs10b_zinb <- zeroinfl(SBS10b ~ pgs + age + ms_status + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Accuracy,
                        data = merge_mutsig_covars, 
                        dist = "negbin")
summary_sbs10b_zinb <- summary(sbs10b_zinb)
sbs10b_zinb_results <- coef(summary_sbs10b_zinb)
sbs10b_zinb_results_count <- sbs10b_zinb_results$count
sbs10b_zinb_results_zero <- sbs10b_zinb_results$zero
sbs10b_results <- rbind(sbs10b_zinb_results_count, sbs10b_zinb_results_zero)
write.csv(sbs10b_results, file = "sbs10b_zinb.csv")

sbs13_zinb <- zeroinfl(SBS13 ~ pgs + age + ms_status + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Accuracy,
                       data = merge_mutsig_covars, 
                       dist = "negbin")
summary_sbs13_zinb <- summary(sbs13_zinb)
sbs13_zinb_results <- coef(summary_sbs13_zinb)
sbs13_zinb_results_count <- sbs13_zinb_results$count
sbs13_zinb_results_zero <- sbs13_zinb_results$zero
sbs13_results <- rbind(sbs13_zinb_results_count, sbs13_zinb_results_zero)
write.csv(sbs13_results, file = "sbs13_zinb.csv")

sbs15_zinb <- zeroinfl(SBS15 ~ pgs + age + ms_status + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Accuracy,
                       data = merge_mutsig_covars, 
                       dist = "negbin")
summary_sbs15_zinb <- summary(sbs15_zinb)
sbs15_zinb_results <- coef(summary_sbs15_zinb)
sbs15_zinb_results_count <- sbs15_zinb_results$count
sbs15_zinb_results_zero <- sbs15_zinb_results$zero
sbs15_results <- rbind(sbs15_zinb_results_count, sbs15_zinb_results_zero)
write.csv(sbs15_results, file = "sbs15_zinb.csv")

sbs40_zinb <- zeroinfl(SBS40 ~ pgs + age + ms_status + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Accuracy,
                       data = merge_mutsig_covars, 
                       dist = "negbin")
summary_sbs40_zinb <- summary(sbs40_zinb)
sbs40_zinb_results <- coef(summary_sbs40_zinb)
sbs40_zinb_results_count <- sbs40_zinb_results$count
sbs40_zinb_results_zero <- sbs40_zinb_results$zero
sbs40_results <- rbind(sbs40_zinb_results_count, sbs40_zinb_results_zero)
write.csv(sbs40_results, file = "sbs40_zinb.csv")

sbs44_zinb <- zeroinfl(SBS44 ~ pgs + age + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Accuracy,
                       data = merge_mutsig_covars, 
                       dist = "negbin") 
summary_sbs44_zinb <- summary(sbs44_zinb)
sbs44_zinb_results <- coef(summary_sbs44_zinb)
sbs44_zinb_results_count <- sbs44_zinb_results$count
sbs44_zinb_results_zero <- sbs44_zinb_results$zero
sbs44_results <- rbind(sbs44_zinb_results_count, sbs44_zinb_results_zero)
write.csv(sbs44_results, file = "sbs44_zinb.csv")

sbs_all <- rbind(sbs2_results,  
                 sbs10a_results, 
                 sbs10b_results, 
                 sbs13_results, 
                 sbs15_results, 
                 sbs40_results,
                 sbs44_results)
write.csv(sbs_all, file = "sbs_zinb.csv")

### MSS/MSI separate analyses:
### Selecting MSS cases:
clinical7 <- clinical5 %>%
  filter(ms_status == "MSS")
dim(clinical7)
clinical7$ms_status <- as.character(clinical7$ms_status)
clinical7 <- clinical7 %>%
  mutate(ms_status = replace(ms_status, ms_status == "MSS", 1))

### Merging data for MSS cases:
merge_mutsig_covars_mss <- pgs %>% 
  inner_join(clinical7, by = "tcgaid") %>% 
  inner_join(ancestry2, by = "tcgaid") %>% 
  inner_join(mutsig4, by = "tcgaid")
dim(merge_mutsig_covars_mss)

### Selecting MSI-H + MSI-L cases:
clinical8 <- clinical5 %>%
  filter(ms_status != "MSS")
dim(clinical8)
clinical8$ms_status <- as.character(clinical8$ms_status)
clinical8 <- clinical8 %>%
  mutate(ms_status = replace(ms_status, ms_status == "MSI-H", 1)) %>%
  mutate(ms_status = replace(ms_status, ms_status == "MSI-L", 1))

### Merging data for MSI-H and MSI-L cases:
merge_mutsig_covars_msi <- pgs %>%
  inner_join(clinical8, by = "tcgaid") %>% 
  inner_join(ancestry2, by = "tcgaid") %>% 
  inner_join(mutsig4, by = "tcgaid")
dim(merge_mutsig_covars_msi)

### Run ZINB analysis for MSS cases:
sbs2_zinb_mss <- zeroinfl(SBS2 ~ pgs + age + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Accuracy,
                          data = merge_mutsig_covars_mss, 
                          dist = "negbin")
summary_sbs2_zinb_mss <- summary(sbs2_zinb_mss)
sbs2_zinb_results_mss <- coef(summary_sbs2_zinb_mss)
sbs2_zinb_results_count_mss <- sbs2_zinb_results_mss$count
sbs2_zinb_results_zero_mss <- sbs2_zinb_results_mss$zero
sbs2_results_mss <- rbind(sbs2_zinb_results_count_mss, sbs2_zinb_results_zero_mss)
write.csv(sbs2_results_mss, file = "sbs2_zinb_mss.csv")

sbs10a_zinb_mss <- zeroinfl(SBS10a ~ pgs + age + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Accuracy,
                            data = merge_mutsig_covars_mss, 
                            dist = "negbin")
summary_sbs10a_zinb_mss <- summary(sbs10a_zinb_mss)
sbs10a_zinb_results_mss <- coef(summary_sbs10a_zinb_mss)
sbs10a_zinb_results_count_mss <- sbs10a_zinb_results_mss$count
sbs10a_zinb_results_zero_mss <- sbs10a_zinb_results_mss$zero
sbs10a_results_mss <- rbind(sbs10a_zinb_results_count_mss, sbs10a_zinb_results_zero_mss)
write.csv(sbs10a_results_mss, file = "sbs10a_zinb_mss.csv")

sbs10b_zinb_mss <- zeroinfl(SBS10b ~ pgs + age + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Accuracy,
                            data = merge_mutsig_covars_mss, 
                            dist = "negbin")
summary_sbs10b_zinb_mss <- summary(sbs10b_zinb_mss)
sbs10b_zinb_results_mss <- coef(summary_sbs10b_zinb_mss)
sbs10b_zinb_results_count_mss <- sbs10b_zinb_results_mss$count
sbs10b_zinb_results_zero_mss <- sbs10b_zinb_results_mss$zero
sbs10b_results_mss <- rbind(sbs10b_zinb_results_count_mss, sbs10b_zinb_results_zero_mss)
write.csv(sbs10b_results_mss, file = "sbs10b_zinb_mss.csv")

sbs13_zinb_mss <- zeroinfl(SBS13 ~ pgs + age + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Accuracy,
                           data = merge_mutsig_covars_mss, 
                           dist = "negbin")
summary_sbs13_zinb_mss <- summary(sbs13_zinb_mss)
sbs13_zinb_results_mss <- coef(summary_sbs13_zinb_mss)
sbs13_zinb_results_count_mss <- sbs13_zinb_results_mss$count
sbs13_zinb_results_zero_mss <- sbs13_zinb_results_mss$zero
sbs13_results_mss <- rbind(sbs13_zinb_results_count_mss, sbs13_zinb_results_zero_mss)
write.csv(sbs13_results_mss, file = "sbs13_zinb_mss.csv")

sbs15_zinb_mss <- zeroinfl(SBS15 ~ pgs + age + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Accuracy,
                           data = merge_mutsig_covars_mss, 
                           dist = "negbin")
summary_sbs15_zinb_mss <- summary(sbs15_zinb_mss)
sbs15_zinb_results_mss <- coef(summary_sbs15_zinb_mss)
sbs15_zinb_results_count_mss <- sbs15_zinb_results_mss$count
sbs15_zinb_results_zero_mss <- sbs15_zinb_results_mss$zero
sbs15_results_mss <- rbind(sbs15_zinb_results_count_mss, sbs15_zinb_results_zero_mss)
write.csv(sbs15_results_mss, file = "sbs15_zinb_mss.csv")

sbs40_zinb_mss <- zeroinfl(SBS40 ~ pgs + age + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Accuracy,
                           data = merge_mutsig_covars_mss, 
                           dist = "negbin")
summary_sbs40_zinb_mss <- summary(sbs40_zinb_mss)
sbs40_zinb_results_mss <- coef(summary_sbs40_zinb_mss)
sbs40_zinb_results_count_mss <- sbs40_zinb_results_mss$count
sbs40_zinb_results_zero_mss <- sbs40_zinb_results_mss$zero
sbs40_results_mss <- rbind(sbs40_zinb_results_count_mss, sbs40_zinb_results_zero_mss)
write.csv(sbs40_results_mss, file = "sbs40_zinb_mss.csv")

sbs44_zinb_mss <- zeroinfl(SBS44 ~ pgs + age + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Accuracy,
                           data = merge_mutsig_covars_mss, 
                           dist = "negbin")
summary_sbs44_zinb_mss <- summary(sbs44_zinb_mss)
sbs44_zinb_results_mss <- coef(summary_sbs44_zinb_mss)
sbs44_zinb_results_count_mss <- sbs44_zinb_results_mss$count
sbs44_zinb_results_zero_mss <- sbs44_zinb_results_mss$zero
sbs44_results_mss <- rbind(sbs44_zinb_results_count_mss, sbs44_zinb_results_zero_mss)
write.csv(sbs44_results_mss, file = "sbs44_zinb_mss.csv")

sbs_all_mss <- rbind(sbs2_results_mss,  
                     sbs10a_results_mss, 
                     sbs10b_results_mss, 
                     sbs13_results_mss, 
                     sbs15_results_mss, 
                     sbs40_results_mss,
                     sbs44_results_mss)
write.csv(sbs_all_mss, file = "sbs_zinb_mss.csv")

### Run ZINB analysis for MSI cases:
sbs2_zinb_msi <- zeroinfl(SBS2 ~ pgs + age + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Accuracy,
                          data = merge_mutsig_covars_msi, 
                          dist = "negbin")
summary_sbs2_zinb_msi <- summary(sbs2_zinb_msi)
sbs2_zinb_results_msi <- coef(summary_sbs2_zinb_msi)
sbs2_zinb_results_count_msi <- sbs2_zinb_results_msi$count
sbs2_zinb_results_zero_msi <- sbs2_zinb_results_msi$zero
sbs2_results_msi <- rbind(sbs2_zinb_results_count_msi, sbs2_zinb_results_zero_msi)
write.csv(sbs2_results_msi, file = "sbs2_zinb_msi.csv")

sbs10a_zinb_msi <- zeroinfl(SBS10a ~ pgs + age + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Accuracy,
                            data = merge_mutsig_covars_msi, 
                            dist = "negbin")
summary_sbs10a_zinb_msi <- summary(sbs10a_zinb_msi)
sbs10a_zinb_results_msi <- coef(summary_sbs10a_zinb_msi)
sbs10a_zinb_results_count_msi <- sbs10a_zinb_results_msi$count
sbs10a_zinb_results_zero_msi <- sbs10a_zinb_results_msi$zero
sbs10a_results_msi <- rbind(sbs10a_zinb_results_count_msi, sbs10a_zinb_results_zero_msi)
write.csv(sbs10a_results_msi, file = "sbs10a_zinb_msi.csv")

sbs10b_zinb_msi <- zeroinfl(SBS10b ~ pgs + age + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Accuracy,
                            data = merge_mutsig_covars_msi, 
                            dist = "negbin")
summary_sbs10b_zinb_msi <- summary(sbs10b_zinb_msi)
sbs10b_zinb_results_msi <- coef(summary_sbs10b_zinb_msi)
sbs10b_zinb_results_count_msi <- sbs10b_zinb_results_msi$count
sbs10b_zinb_results_zero_msi <- sbs10b_zinb_results_msi$zero
sbs10b_results_msi <- rbind(sbs10b_zinb_results_count_msi, sbs10b_zinb_results_zero_msi)
write.csv(sbs10b_results_msi, file = "sbs10b_zinb_msi.csv")

sbs13_zinb_msi <- zeroinfl(SBS13 ~ pgs + age + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Accuracy,
                           data = merge_mutsig_covars_msi, 
                           dist = "negbin")
summary_sbs13_zinb_msi <- summary(sbs13_zinb_msi)
sbs13_zinb_results_msi <- coef(summary_sbs13_zinb_msi)
sbs13_zinb_results_count_msi <- sbs13_zinb_results_msi$count
sbs13_zinb_results_zero_msi <- sbs13_zinb_results_msi$zero
sbs13_results_msi <- rbind(sbs13_zinb_results_count_msi, sbs13_zinb_results_zero_msi)
write.csv(sbs13_results_msi, file = "sbs13_zinb_msi.csv")

sbs15_zinb_msi <- zeroinfl(SBS15 ~ pgs + age + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Accuracy,
                           data = merge_mutsig_covars_msi, 
                           dist = "negbin")
summary_sbs15_zinb_msi <- summary(sbs15_zinb_msi)
sbs15_zinb_results_msi <- coef(summary_sbs15_zinb_msi)
sbs15_zinb_results_count_msi <- sbs15_zinb_results_msi$count
sbs15_zinb_results_zero_msi <- sbs15_zinb_results_msi$zero
sbs15_results_msi <- rbind(sbs15_zinb_results_count_msi, sbs15_zinb_results_zero_msi)
write.csv(sbs15_results_msi, file = "sbs15_zinb_msi.csv")

sbs40_zinb_msi <- zeroinfl(SBS40 ~ pgs + age + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Accuracy,
                           data = merge_mutsig_covars_msi, 
                           dist = "negbin")
summary_sbs40_zinb_msi <- summary(sbs40_zinb_msi)
sbs40_zinb_results_msi <- coef(summary_sbs40_zinb_msi)
sbs40_zinb_results_count_msi <- sbs40_zinb_results_msi$count
sbs40_zinb_results_zero_msi <- sbs40_zinb_results_msi$zero
sbs40_results_msi <- rbind(sbs40_zinb_results_count_msi, sbs40_zinb_results_zero_msi)
write.csv(sbs40_results_msi, file = "sbs40_zinb_msi.csv")

sbs44_zinb_msi <- zeroinfl(SBS44 ~ pgs + age + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Accuracy,
                           data = merge_mutsig_covars_msi, 
                           dist = "negbin")
summary_sbs44_zinb_msi <- summary(sbs44_zinb_msi)
sbs44_zinb_results_msi <- coef(summary_sbs44_zinb_msi)
sbs44_zinb_results_count_msi <- sbs44_zinb_results_msi$count
sbs44_zinb_results_zero_msi <- sbs44_zinb_results_msi$zero
sbs44_results_msi <- rbind(sbs44_zinb_results_count_msi, sbs44_zinb_results_zero_msi)
write.csv(sbs44_results_msi, file = "sbs44_zinb_msi.csv")

sbs_all_msi <- rbind(sbs2_results_msi,  
                     sbs10a_results_msi, 
                     sbs10b_results_msi, 
                     sbs13_results_msi, 
                     sbs15_results_msi, 
                     sbs40_results_msi,
                     sbs44_results_msi)
write.csv(sbs_all_msi, file = "sbs_zinb_msi.csv")
