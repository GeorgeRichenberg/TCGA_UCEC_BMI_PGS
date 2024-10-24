### Set working directory and load required packages:
setwd()
library(vroom)
library(tidyverse)

### Load required data:
epo_chr7 <- vroom("discovery_chr7_EPO:P01588:OID20522:v1:Inflammation.gz")

### Filter for SNPs ±100 kb of EPO gene [chr7:100,720,468-100,723,700 in GRCh38]:
dim(epo_chr7)
colnames(epo_chr7)
epo_chr7_2 <- epo_chr7 %>%
  filter(GENPOS >= 100720468 - 100000,  
         GENPOS <= 100723700 + 100000)
dim(epo_chr7_2)

### Format epo_chr7:
epo_chr7_3 <- epo_chr7_2 %>%
  mutate(P = 10^-LOG10P) %>%
  arrange(P)
epo_chr7_4 <- epo_chr7_3 %>%
  mutate(SNP = paste(CHROM, GENPOS, sep = "_")) %>%
  select(SNP,CHROM, GENPOS, ALLELE0, ALLELE1, A1FREQ, BETA, SE, P) %>%
  rename("ALLELE0.epo" = "ALLELE0",
         "ALLELE1.epo" = "ALLELE1",
         "A1FREQ.epo" = "A1FREQ",
         "B.epo" = "BETA",
         "SE.epo" = "SE",
         "P.epo" = "P")
head(epo_chr7_4)
write.csv(epo_chr7_4, file = "Sun_et_al_2023_sum_stats_cis_epo.csv", row.names = FALSE)