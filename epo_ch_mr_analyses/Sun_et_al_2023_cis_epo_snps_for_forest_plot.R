### Set working directory and load required packages:
setwd()
library(tidyverse)
library(vroom)

### Load data:
cis_epo <- vroom("Sun_et_al_2023_cis_epo_overall_matched.csv", col_select = c(SNP, CHROM, GENPOS, ALLELE0.epo, ALLELE1.epo, A1FREQ.epo, B.epo, SE.epo, P.epo, rsid))

### Format epo_rs34164109:
epo_rs34164109 <- epo_rs34164109 %>%
  mutate(ALLELE1.epo = str_replace(ALLELE1.epo, "RUE", ""))

### Filter rs1617640 and rs11976235 from cis_epo:
epo_snps <- cis_epo %>%
  filter(rsid %in% c("rs1617640", "rs11976235"))

### Rename epo data:
epo_snps <- epo_snps %>%
  mutate(phenotype = "EPO") %>%
  rename("ALLELE0" = "ALLELE0.epo" ,
         "ALLELE1" = "ALLELE1.epo" ,
         "A1FREQ" = "A1FREQ.epo",
         "B" = "B.epo",
         "SE" = "SE.epo",
         "P" = "P.epo")

### Re-order row numbers of epo_hgb_chr7:
epo_snps_2 <- epo_snps %>%
  slice(c(2, 1))

### Re-name SNP columns:
epo_snps_2 <- epo_snps_2 %>%
  mutate(variant = c("common variant: rs1617640-A", "rare variant: rs11976235-T"))

### Calculate confidence intervals for combined epo_hgb:
epo_snps_2$lower_CI <- epo_snps_2$B - (epo_snps_2$SE * 1.96)
epo_snps_2$upper_CI <- epo_snps_2$B + (epo_snps_2$SE * 1.96)

### Save output:
write.csv(epo_snps_2, file = "epo_snps.csv", row.names = FALSE)
