### Set working directory and load required packages:
setwd()
library(vroom)
library(tidyverse)
library(TwoSampleMR)

### Load data:
overall <- vroom("overall_summary_stats_cis_epo_region.csv", col_select = c("SNP", "chromosome", "base_pair_location", "other_allele.overall", "effect_allele.overall", "eaf.overall", "B.overall", "SE.overall", "P.overall"))
dnmt3a <- vroom("dnmt3a_summary_stats_cis_epo_region.csv", col_select = c("SNP", "chromosome", "base_pair_location", "other_allele.dnmt3a", "effect_allele.dnmt3a", "eaf.dnmt3a", "B.dnmt3a", "SE.dnmt3a", "P.dnmt3a"))
tet2 <- vroom("tet2_summary_stats_cis_epo_region.csv", col_select = c("SNP", "chromosome", "base_pair_location", "other_allele.tet2", "effect_allele.tet2", "eaf.tet2", "B.tet2", "SE.tet2", "P.tet2"))
  
### Filter rs1617640 and rs11976235 from CH data:
overall_2 <- overall %>%
  filter(SNP %in% c("7_100719675", "7_100722758")) %>%
  mutate(rsid = c("rs1617640", "rs11976235")) %>%
  select(rsid, everything(), -SNP)

dnmt3a_2 <- dnmt3a %>%
  filter(SNP %in% c("7_100719675", "7_100722758")) %>%
  mutate(rsid = c("rs1617640", "rs11976235")) %>%
  select(rsid, everything(), -SNP)

tet2_2 <- tet2 %>%
  filter(SNP %in% c("7_100719675", "7_100722758")) %>%
  mutate(rsid = c("rs1617640", "rs11976235")) %>%
  select(rsid, everything(), -SNP)

### Define overall as outcome:
colnames(overall_2)
out_overall <- format_data(overall_2,
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
                           chr_col = "chromosome",
                           pos_col = "base_pair_location")

### Define dnmt3a as outcome:
colnames(dnmt3a_2)
out_dnmt3a <- format_data(dnmt3a_2,
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
                          chr_col = "chromosome",
                          pos_col = "base_pair_location")

### Defining tet2 as outcome:
colnames(tet2_2)
out_tet2 <- format_data(tet2_2,
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
                        chr_col = "chromosome",
                        pos_col = "base_pair_location")

### Save:
write_csv(out_overall, file = "out_overall.csv")
write_csv(out_dnmt3a, file = "out_dnmt3a.csv")
writecsv(out_tet2, file = "out_tet2.csv")