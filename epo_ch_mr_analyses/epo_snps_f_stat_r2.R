### Set working directory and load required packages:
setwd()
library(vroom)
library(tidyverse)

### Load required data:
epo <- vroom("Sun_et_al_2023_cis_epo_overall_harmonised.csv")

### Select relevant columns:
epo_rs1617640 <- epo %>%
  filter(SNP == "rs1617640") %>%
  select(SNP, effect_allele.exposure, other_allele.exposure, beta.exposure, eaf.exposure, se.exposure, pval.exposure)
epo_rs1617640
epo_rs11976235 <- epo %>%
  filter(SNP == "rs11976235") %>%
  select(SNP, effect_allele.exposure, other_allele.exposure, beta.exposure, eaf.exposure, se.exposure, pval.exposure)
epo_rs11976235

### Calculate MAF:
for (j in 1:nrow(epo_rs1617640)) {
  if (epo_rs1617640$eaf.exposure[j] > 0.5) {
    epo_rs1617640$MAF1[j] <- 1 - epo_rs1617640$eaf.exposure[j]
  } else {
    epo_rs1617640$MAF1[j] <- epo_rs1617640$eaf.exposure[j]
  }
}

for (j in 1:nrow(epo_rs11976235)) {
  if (epo_rs11976235$eaf.exposure[j] > 0.5) {
    epo_rs11976235$MAF1[j] <- 1 - epo_rs11976235$eaf.exposure[j]
  } else {
    epo_rs11976235$MAF1[j] <- epo_rs11976235$eaf.exposure[j]
  }
}

### Calculate proporiton of variance explained (R2): 
r2_rs1617640 <- epo_rs1617640 %>%
  mutate(
    numerator = 2 * epo_rs1617640$beta.exposure^2 * epo_rs1617640$MAF1 * (1 - epo_rs1617640$MAF1),
    denominator = numerator + (epo_rs1617640$se.exposure^2 * 2 * 33657 * epo_rs1617640$MAF1 * (1 - epo_rs1617640$MAF1)),
    R2 = numerator / denominator)
r2_rs1617640
r2_rs11976235 <- epo_rs11976235 %>%
  mutate(
    numerator = 2 * epo_rs11976235$beta.exposure^2 * epo_rs11976235$MAF1 * (1 - epo_rs11976235$MAF1),
    denominator = numerator + (epo_rs11976235$se.exposure^2 * 2 * 33657 * epo_rs11976235$MAF1 * (1 - epo_rs11976235$MAF1)),
    R2 = numerator / denominator)
r2_rs11976235

### Use R2 value to calculate the F-statistic:
f_statistic_rs1617640 <- r2_rs1617640 %>%
  mutate(
    F = (R2 * (33657 - 2)) / ((1 - R2) * 1))
f_statistic_rs1617640

f_statistic_rs11976235 <- r2_rs11976235 %>%
  mutate(
    F = (R2 * (33657 - 2)) / ((1 - R2) * 1))
f_statistic_rs11976235

### Display as table:
f_r2_stats_rs1617640 <- f_statistic_rs1617640 %>%
  select(SNP, R2, F)
f_r2_stats_rs11976235 <- f_statistic_rs11976235 %>%
  select(SNP, R2, F)
f_r2 <- rbind(f_r2_stats_rs1617640, f_r2_stats_rs11976235)
f_r2 <- f_r2 %>%
  mutate(R2 = R2*100) %>%
  rename("R2 (%)" = R2)
