### Set working directory and load required packages:
setwd()
library(vroom)
library(tidyverse)
library(coloc)

### Load data:
epo_overall_harmonised <- vroom("Sun_et_al_2023_cis_epo_overall_harmonised.csv")
epo_dnmt3a_harmonised <- vroom("Sun_et_al_2023_cis_epo_dnmt3a_harmonised.csv")
epo_tet2_harmonised <- vroom("Sun_et_al_2023_cis_epo_tet2_harmonised.csv")

### Calculate MAF:
for (j in 1:nrow(epo_overall_harmonised)){
  if(epo_overall_harmonised$eaf.exposure[j] > 0.5){epo_overall_harmonised$MAF1[j] = 1 - epo_overall_harmonised$eaf.exposure[j]} 
  else {epo_overall_harmonised$MAF1[j] = epo_overall_harmonised$eaf.exposure[j]}  
  if(epo_overall_harmonised$eaf.outcome[j] > 0.5){epo_overall_harmonised$MAF2[j] = 1 - epo_overall_harmonised$eaf.outcome[j]} 
  else {epo_overall_harmonised$MAF2[j] = epo_overall_harmonised$eaf.outcome[j]}}

for (j in 1:nrow(epo_dnmt3a_harmonised)){
  if(epo_dnmt3a_harmonised$eaf.exposure[j]>0.5){epo_dnmt3a_harmonised$MAF1[j]=1-epo_dnmt3a_harmonised$eaf.exposure[j]} 
  else {epo_dnmt3a_harmonised$MAF1[j]=epo_dnmt3a_harmonised$eaf.exposure[j]}  
  if(epo_dnmt3a_harmonised$eaf.outcome[j]>0.5){epo_dnmt3a_harmonised$MAF2[j]=1-epo_dnmt3a_harmonised$eaf.outcome[j]} 
  else {epo_dnmt3a_harmonised$MAF2[j]=epo_dnmt3a_harmonised$eaf.outcome[j]}}

for (j in 1:nrow(epo_tet2_harmonised)){
  if(epo_tet2_harmonised$eaf.exposure[j]>0.5){epo_tet2_harmonised$MAF1[j]=1-epo_tet2_harmonised$eaf.exposure[j]} 
  else {epo_tet2_harmonised$MAF1[j]=epo_tet2_harmonised$eaf.exposure[j]}  
  if(epo_tet2_harmonised$eaf.outcome[j]>0.5){epo_tet2_harmonised$MAF2[j]=1-epo_tet2_harmonised$eaf.outcome[j]} 
  else {epo_tet2_harmonised$MAF2[j]=epo_tet2_harmonised$eaf.outcome[j]}}

### Add study sample sizes:
N1 <- 
N_overall <- 368526
N_dnmt3a <- 359088
N_tet2 <- 346787

### Run coloc for epo_overall_harmonised:
coloc_overall <- coloc.abf(dataset1 = list(snp = epo_overall_harmonised$SNP, 
                                           beta = epo_overall_harmonised$beta.exposure, 
                                           varbeta = epo_overall_harmonised$se.exposure^2, 
                                           MAF = epo_overall_harmonised$MAF1,
                                           N = N1, 
                                           type = "cc"),
                           dataset2 = list(snp = epo_overall_harmonised$SNP, 
                                           beta = epo_overall_harmonised$beta.outcome, 
                                           varbeta = epo_overall_harmonised$se.outcome^2, 
                                           MAF = epo_overall_harmonised$MAF2,
                                           N = N_overall, 
                                           type = "cc"),
                           p1 = 1e-04, p2 = 1e-04, p12 = 1e-05)
coloc_overall
coloc_overall_df <- as.data.frame(coloc_overall$summary)
coloc_overall_df

### Run coloc for epo_dnmt3a_harmonised:
coloc_dnmt3a <- coloc.abf(dataset1 = list(snp = epo_dnmt3a_harmonised$SNP, 
                                          beta = epo_dnmt3a_harmonised$beta.exposure, 
                                          varbeta = epo_dnmt3a_harmonised$se.exposure^2, 
                                          MAF = epo_dnmt3a_harmonised$MAF1, 
                                          N = N1, 
                                          type = "cc"),
                          dataset2 = list(snp = epo_dnmt3a_harmonised$SNP, 
                                          beta = epo_dnmt3a_harmonised$beta.outcome, 
                                          varbeta = epo_dnmt3a_harmonised$se.outcome^2, 
                                          MAF = epo_dnmt3a_harmonised$MAF2, 
                                          N = N_dnmt3a, 
                                          type = "cc"),
                          p1 = 1e-04, p2 = 1e-04, p12 = 1e-05)
coloc_dnmt3a
coloc_dnmt3a_df <- as.data.frame(coloc_dnmt3a$summary)
coloc_dnmt3a_df

### Run coloc for epo_tet2_harmonised:
coloc_tet2 <- coloc.abf(dataset1 = list(snp = epo_tet2_harmonised$SNP, 
                                        beta = epo_tet2_harmonised$beta.exposure, 
                                        varbeta = epo_tet2_harmonised$se.exposure^2, 
                                        MAF = epo_tet2_harmonised$MAF1, 
                                        N = N1, 
                                        type = "cc"),
                        dataset2 = list(snp = epo_tet2_harmonised$SNP, 
                                        beta = epo_tet2_harmonised$beta.outcome, 
                                        varbeta = epo_tet2_harmonised$se.outcome^2, 
                                        MAF = epo_tet2_harmonised$MAF2, 
                                        N = N_tet2, 
                                        type = "cc"),
                        p1 = 1e-04, p2 = 1e-04, p12 = 1e-05)
coloc_tet2
coloc_tet2_df <- as.data.frame(coloc_tet2$summary)
coloc_tet2_df

### Calculate alternative coloc scores:
alt_coloc_epo_overall <- coloc_overall_df[6,] / (coloc_overall_df[5,] + coloc_overall_df[6,])
alt_coloc_epo_overall
alt_coloc_epo_dnmt3a <- coloc_dnmt3a_df[6,] / (coloc_dnmt3a_df[5,] + coloc_dnmt3a_df[6,])
alt_coloc_epo_dnmt3a
alt_coloc_epo_tet2 <- coloc_tet2_df[6,] / (coloc_tet2_df[5,] + coloc_tet2_df[6,])
alt_coloc_epo_tet2

### Combine scores into one table:
alt_coloc <- rbind(alt_coloc_epo_overall, alt_coloc_epo_dnmt3a, alt_coloc_epo_tet2)
colnames(alt_coloc) <-  "alt_coloc_score"
alt_coloc

### Combine coloc and alt_coloc scores into one table:
coloc_overall_df_2 <- t(coloc_overall_df)
coloc_dnmt3a_df_2 <- t(coloc_dnmt3a_df)
coloc_tet2_df_2 <- t(coloc_tet2_df)
coloc <- rbind(coloc_overall_df_2, coloc_dnmt3a_df_2, coloc_tet2_df_2)
coloc_alt_coloc <- cbind(coloc, alt_coloc)
write.csv(coloc_alt_coloc, file = "Sun_et_al_2023_cis_epo_ch_coloc_alt_coloc_scores.csv")
