### Set working directory and load required packages:
setwd()
library(vroom)
library(tidyverse)
library(TwoSampleMR)

### Load data:
epo <- vroom("epo_invnt_metaanalysis_results.txt")

### Filter epo for rs1617640 and rs11976235:
head(epo)
colnames(epo)
epo2 <- epo %>%
  filter(MarkerName %in% c("7:100317298", "7:100320381")) %>%
  select("MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P-value",  "TotalSampleSize")
epo3 <- epo2 %>%
  separate(MarkerName, into = c("Chr", "Pos_37"), sep = ":") %>%
  mutate(Chr = as.numeric(Chr),
         Pos_37 = as.numeric(Pos_37),
         Allele1 = toupper(Allele1),
         Allele2 = toupper(Allele2))

### Liftover from GrCh37 to GrCh38:
epo4 <- epo3 %>%
  mutate(Pos_38 = c(100719675, 100722758)) %>%
  select(-Pos_37) %>%
  select(Chr, Pos_38, everything())

### Add rsIDs to epo4:
epo4 <- epo4 %>%
  mutate(rsid = c("rs1617640", "rs11976235")) %>%
  select(rsid, everything())
         
### Define epo as exposure:
colnames(epo4)
exp_epo <- format_data(epo4,
                       type = "exposure",
                       snps = NULL,
                       header = TRUE,
                       snp_col = "rsid",
                       beta_col = "Effect",
                       se_col = "StdErr",
                       eaf_col = "Freq1",
                       effect_allele_col = "Allele1",
                       other_allele_col = "Allele2",
                       pval_col = "P-value",
                       chr_col = "Chr",
                       pos_col = "Pos_38")

### Save:
write_csv(exp_epo, file = "exp_epo.csv")
