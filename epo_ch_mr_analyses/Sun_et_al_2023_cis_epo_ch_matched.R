### Set working directory and load required packages:
setwd()
library(vroom)
library(tidyverse)

### Load required data:
epo <- vroom("Sun_et_al_2023_sum_stats_cis_epo.csv")
overall <- vroom("overall_summary_stats_cis_epo_region.csv")
dnmt3a <- vroom("dnmt3a_summary_stats_cis_epo_region.csv")
tet2 <- vroom("tet2_summary_stats_cis_epo_region.csv")
g1000 <- vroom("all_g1000_38.txt")

### Explore epo and CH data:
dim(epo)
dim(overall)
dim(dnmt3a)
dim(tet2)

### Merge epo and CH data:
epo_overall <- epo %>%
  inner_join(overall, by = "SNP") %>%
  select(SNP, CHROM, GENPOS, ALLELE0.epo, ALLELE1.epo, A1FREQ.epo, B.epo, SE.epo, P.epo, other_allele.overall, effect_allele.overall, eaf.overall, B.overall, SE.overall, P.overall)
dim(epo_overall)
sum(duplicated(epo_overall$SNP))
epo_overall_2 <- epo_overall %>%
  filter((ALLELE0.epo == other_allele.overall & ALLELE1.epo == effect_allele.overall) |
         (ALLELE1.epo == other_allele.overall & ALLELE0.epo == effect_allele.overall))
dim(epo_overall_2)
sum(duplicated(epo_overall_2$SNP))


epo_dnmt3a <- epo %>%
  inner_join(dnmt3a, by = "SNP") %>%
  select(SNP, CHROM, GENPOS, ALLELE0.epo, ALLELE1.epo, A1FREQ.epo, B.epo, SE.epo, P.epo, other_allele.dnmt3a, effect_allele.dnmt3a, eaf.dnmt3a, B.dnmt3a, SE.dnmt3a, P.dnmt3a)
dim(epo_dnmt3a)
sum(duplicated(epo_dnmt3a$SNP))
epo_dnmt3a_2 <- epo_dnmt3a %>%
  filter((ALLELE0.epo == other_allele.dnmt3a & ALLELE1.epo == effect_allele.dnmt3a) |
         (ALLELE1.epo == other_allele.dnmt3a & ALLELE0.epo == effect_allele.dnmt3a))
dim(epo_dnmt3a_2)
sum(duplicated(epo_dnmt3a_2$SNP))


epo_tet2 <- epo %>%
  inner_join(tet2, by = "SNP") %>%
  select(SNP, CHROM, GENPOS, ALLELE0.epo, ALLELE1.epo, A1FREQ.epo, B.epo, SE.epo, P.epo, other_allele.tet2, effect_allele.tet2, eaf.tet2, B.tet2, SE.tet2, P.tet2)
dim(epo_tet2)
sum(duplicated(epo_tet2l$SNP))
epo_tet2_2 <- epo_tet2 %>%
  filter((ALLELE0.epo == other_allele.tet2 & ALLELE1.epo == effect_allele.tet2) |
         (ALLELE1.epo == other_allele.tet2 & ALLELE0.epo == effect_allele.tet2))
dim(epo_tet2_2)
sum(duplicated(epo_tet2_2$SNP))

### Format g1000:
g1000 <- g1000 %>%
  unite("SNP", c(Chr_38, Pos_38), remove = TRUE)
g1000 <- g1000 %>%
  dplyr::select(-ID)

### Generate rsIDs for merged data:
epo_overall_3 <- epo_overall_2 %>%
  inner_join(g1000, by = "SNP") %>%
  dplyr::select(SNP, CHROM, GENPOS, rsid, everything())
dim(epo_overall_3)
write_csv(epo_overall_3, file = "Sun_et_al_2023_cis_epo_overall_matched.csv")

epo_dnmt3a_3 <- epo_dnmt3a_2 %>%
  inner_join(g1000, by = "SNP") %>%
  dplyr::select(SNP, CHROM, GENPOS, rsid, everything())
dim(epo_dnmt3a_3)
write_csv(epo_dnmt3a_3, file = "Sun_et_al_2023_cis_epo_dnmt3a_matched.csv")

epo_tet2_3 <- epo_tet2_2 %>%
  inner_join(g1000, by = "SNP") %>%
  dplyr::select(SNP, CHROM, GENPOS, rsid, everything())
dim(epo_tet2_3)
write_csv(epo_tet2_3, file = "Sun_et_al_2023_cis_epo_tet2_matched.csv")