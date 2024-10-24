### Set working directory and load required packages:
setwd()
library(vroom)
library(tidyverse)
library(TwoSampleMR)

### Load data:
epo_overall <- vroom("Sun_et_al_2023_cis_epo_overall_matched.csv")
epo_dnmt3a <- vroom("Sun_et_al_2023_cis_epo_dnmt3a_matched.csv")
epo_tet2 <- vroom("Sun_et_al_2023_cis_epo_tet2_matched.csv")

### Define epo as exposure:
colnames(epo_overall)
exp_epo <- format_data(epo_overall,
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

### Define overall_inc as outcome:
colnames(epo_overall)
out_overall <- format_data(epo_overall,
                           type = "outcome",
                           snps = NULL,
                           header = TRUE,
                           snp_col = "rsid",
                           beta_col = "B.overall",
                           se_col = "SE.overall",
                           eaf_col = "eaf.overall",
                           effect_allele_col = "effect_allele.overall",
                           other_allele_col = "other_allele.overall",
                           pval_col = "P.overall",
                           chr_col = "CHROM",
                           pos_col = "GENPOS")

### Define dnmt3a as outcome:
colnames(epo_dnmt3a)
out_dnmt3a <- format_data(epo_dnmt3a,
                          type = "outcome",
                          snps = NULL,
                          header = TRUE,
                          snp_col = "rsid",
                          beta_col = "B.dnmt3a",
                          se_col = "SE.dnmt3a",
                          eaf_col = "eaf.dnmt3a",
                          effect_allele_col = "effect_allele.dnmt3a",
                          other_allele_col = "other_allele.dnmt3a",
                          pval_col = "P.dnmt3a",
                          chr_col = "CHROM",
                          pos_col = "GENPOS")

### Define tet2 as outcome:
colnames(epo_tet2)
out_tet2 <- format_data(epo_tet2,
                        type = "outcome",
                        snps = NULL,
                        header = TRUE,
                        snp_col = "rsid",
                        beta_col = "B.tet2",
                        se_col = "SE.tet2",
                        eaf_col = "eaf.tet2",
                        effect_allele_col = "effect_allele.tet2",
                        other_allele_col = "other_allele.tet2",
                        pval_col = "P.tet2",
                        chr_col = "CHROM",
                        pos_col = "GENPOS")

### Harmonise epo and overall CH:
epo_overall_harmonised <- harmonise_data(exposure_dat = exp_epo, 
                                         outcome_dat = out_overall, 
                                         action = 2)
dim(epo_overall_harmonised)

### Harmonise epo and dnmt3a CH:
epo_dnmt3a_harmonised <- harmonise_data(exposure_dat = exp_epo, 
                                        outcome_dat = out_dnmt3a, 
                                        action = 2)
dim(epo_dnmt3a_harmonised)

### Harmonise epo and tet2 CH:
epo_tet2_harmonised <- harmonise_data(exposure_dat = exp_epo, 
                                      outcome_dat = out_tet2, 
                                      action = 2)
dim(epo_tet2_harmonised)

### Save harmonised datasets as .csv files:
write_csv(epo_overall_harmonised, file = "Sun_et_al_2023_cis_epo_overall_harmonised.csv")
write_csv(epo_dnmt3a_harmonised, file = "Sun_et_al_2023_cis_epo_dnmt3a_harmonised.csv")
write_csv(epo_tet2_harmonised, file = "Sun_et_al_2023_cis_epo_tet2_harmonised.csv")