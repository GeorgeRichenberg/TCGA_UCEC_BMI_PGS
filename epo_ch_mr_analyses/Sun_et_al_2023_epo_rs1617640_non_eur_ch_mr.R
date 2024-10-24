### Set working directory and load required packages:
setwd()
library(vroom)
library(tidyverse)
library(TwoSampleMR)

### Load data:
epo <- vroom("Sun_et_al_2023_cis_epo_overall_matched.csv")
ch_non_eur <- vroom("epo_common_snp_ch_non_eur.csv")

### Select rs1617640 and define epo as exposure:
colnames(epo)
epo2 <- epo %>%
  filter(rsid == "rs1617640") %>%
  select(c(rsid, CHROM, GENPOS, ALLELE0.epo, ALLELE1.epo, A1FREQ.epo, B.epo, SE.epo, P.epo))
exp_epo <- format_data(epo2,
                       type = "exposure",
                       snps = NULL,
                       header = TRUE,
                       snp_col = "rsid",
                       beta_col = "B.epo",
                       se_col = "SE.epo",
                       eaf_col = "A1FREQ.epo",
                       effect_allele_col = "ALLELE1.epo",
                       other_allele_col = "ALLELE0.epo",
                       pval_col = "P.epo",
                       chr_col = "CHROM",
                       pos_col = "GENPOS")

### Format ch_non_eur:
head(ch_non_eur)
colnames(ch_non_eur)
ch_non_eur <- ch_non_eur %>%
  mutate(rsid = "rs1617640")
ch_non_eur$beta <- log(ch_non_eur$odds_ratio)
ch_non_eur$se <- (((log(ch_non_eur$ci_upper))-(log(ch_non_eur$odds_ratio)))/1.96)

### Divide ch_non_eur by ancestry:
ch_sas <- ch_non_eur %>%
  filter(ancestry == "SAS")
ch_afr <- ch_non_eur %>%
  filter(ancestry == "AFR")
ch_eas <- ch_non_eur %>%
  filter(ancestry == "EAS")

### Define each ch_non_eur as outcome:
out_ch_sas <- format_data(ch_sas,
                          type = "outcome",
                          snps = NULL,
                          header = TRUE,
                          snp_col = "rsid",
                          beta_col = "beta",
                          se_col = "se",
                          eaf_col = "effect_allele_frequency",
                          effect_allele_col = "effect_allele",
                          other_allele_col = "other_allele",
                          pval_col = "p_value",
                          chr_col = "chromosome",
                          pos_col = "base_pair_location")

out_ch_afr <- format_data(ch_afr,
                          type = "outcome",
                          snps = NULL,
                          header = TRUE,
                          snp_col = "rsid",
                          beta_col = "beta",
                          se_col = "se",
                          eaf_col = "effect_allele_frequency",
                          effect_allele_col = "effect_allele",
                          other_allele_col = "other_allele",
                          pval_col = "p_value",
                          chr_col = "chromosome",
                          pos_col = "base_pair_location")

out_ch_eas <- format_data(ch_eas,
                          type = "outcome",
                          snps = NULL,
                          header = TRUE,
                          snp_col = "rsid",
                          beta_col = "beta",
                          se_col = "se",
                          eaf_col = "effect_allele_frequency",
                          effect_allele_col = "effect_allele",
                          other_allele_col = "other_allele",
                          pval_col = "p_value",
                          chr_col = "chromosome",
                          pos_col = "base_pair_location")

### Harmonise epo and each ch_non_eur:
epo_ch_sas_harmonised <- harmonise_data(exposure_dat = exp_epo, 
                                        outcome_dat = out_ch_sas, 
                                        action = 2)
dim(epo_ch_sas_harmonised)
epo_ch_afr_harmonised <- harmonise_data(exposure_dat = exp_epo, 
                                        outcome_dat = out_ch_afr, 
                                        action = 2)
dim(epo_ch_afr_harmonised)
epo_ch_eas_harmonised <- harmonise_data(exposure_dat = exp_epo, 
                                        outcome_dat = out_ch_eas, 
                                        action = 2)
dim(epo_ch_eas_harmonised)

### Run MR:
epo_ch_sas_mr <- generate_odds_ratios(mr(epo_ch_sas_harmonised))
epo_ch_sas_mr
epo_ch_afr_mr <- generate_odds_ratios(mr(epo_ch_afr_harmonised))
epo_ch_afr_mr
epo_ch_eas_mr <- generate_odds_ratios(mr(epo_ch_eas_harmonised))
epo_ch_eas_mr

### Format results:
epo_ch_sas_mr <- epo_ch_sas_mr %>%
  mutate(exposure = "epo_common_rs1617640") %>%
  mutate(outcome = "sas_overall_ch")
epo_ch_afr_mr <- epo_ch_afr_mr %>%
  mutate(exposure = "epo_common_rs1617640") %>%
  mutate(outcome = "afr_overall_ch")
epo_ch_eas_mr <- epo_ch_eas_mr %>%
  mutate(exposure = "epo_common_rs1617640") %>%
  mutate(outcome = "eas_overall_ch")

### Save results:
write_csv(epo_ch_sas_mr, file = "epo_common_rs1617640_sas_ch_mr_res.csv")
write_csv(epo_ch_afr_mr, file = "epo_common_rs1617640_afr_ch_mr_res.csv")
write_csv(epo_ch_eas_mr, file = "epo_common_rs1617640_eas_ch_mr_res.csv")

### Combine results and save:
epo_ch_sas_afr_eas <- rbind(epo_ch_sas_mr, epo_ch_afr_mr,  epo_ch_eas_mr)
write_csv(epo_ch_sas_afr_eas, file = "epo_common_rs1617640_sas_afr_eas_ch_mr_res.csv")