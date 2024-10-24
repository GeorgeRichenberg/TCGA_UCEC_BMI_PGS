### Set working directory and load required packages:
setwd()
library(vroom)
library(tidyverse)

### Load data:
overall <- vroom("GCST90165261_buildGRCh38.tsv.gz", col_select = c(chromosome, base_pair_location, other_allele, effect_allele, odds_ratio, ci_lower, ci_upper, p_value, effect_allele_frequency))
dnmt3a <- vroom("GCST90165271_buildGRCh38.tsv.gz", col_select = c(chromosome, base_pair_location, other_allele, effect_allele, odds_ratio, ci_lower, ci_upper, p_value, effect_allele_frequency))
tet2 <- vroom("GCST90165281_buildGRCh38.tsv.gz", col_select = c(chromosome, base_pair_location, other_allele, effect_allele, odds_ratio, ci_lower, ci_upper, p_value, effect_allele_frequency))

### OVERALL:
### Filter SNPs ±100 kb of EPO gene [chr7:100,720,468-100,723,700 in GRCh38]:
overall <- overall %>%
  filter(chromosome == 7)
overall <- overall %>%
  filter(base_pair_location >= 100720468 - 100000,  
         base_pair_location <= 100723700 + 100000)
dim(overall)

### Create beta and SE columns:
overall$beta <- log(overall$odds_ratio)
overall$se <- (((log(overall$ci_upper)) - (log(overall$odds_ratio))) / 1.96)

### Format overall:
overall <- overall %>%
  mutate(SNP = paste(chromosome, base_pair_location, sep = "_")) %>%
  select(SNP, chromosome, base_pair_location, other_allele, effect_allele, effect_allele_frequency, beta, se, p_value) %>%
  rename("other_allele.overall" = "other_allele",
         "effect_allele.overall" = "effect_allele",
         "eaf.overall" = "effect_allele_frequency",
         "B.overall" = "beta",
         "SE.overall" = "se",
         "P.overall" = "p_value")
head(overall)
write_csv(overall, file = "overall_summary_stats_cis_epo_region.csv")

### DNMT3A:
### Filter SNPs ±100 kb of EPO gene [chr7:100,720,468-100,723,700 in GRCh38]:
dnmt3a <- dnmt3a %>%
  filter(chromosome == 7)
dnmt3a <- dnmt3a %>%
  filter(base_pair_location >= 100720468 - 100000,  
         base_pair_location <= 100723700 + 100000)
dim(dnmt3a)

### Create beta and SE columns:
dnmt3a$beta <- log(dnmt3a$odds_ratio)
dnmt3a$se <- (((log(dnmt3a$ci_upper)) - (log(dnmt3a$odds_ratio))) / 1.96)

### Format dnmt3a:
dnmt3a <- dnmt3a %>%
  mutate(SNP = paste(chromosome, base_pair_location, sep = "_")) %>%
  select(SNP, chromosome, base_pair_location, other_allele, effect_allele, effect_allele_frequency, beta, se, p_value) %>%
  rename("other_allele.dnmt3a" = "other_allele",
         "effect_allele.dnmt3a" = "effect_allele",
         "eaf.dnmt3a" = "effect_allele_frequency",
         "B.dnmt3a" = "beta",
         "SE.dnmt3a" = "se",
         "P.dnmt3a" = "p_value")
dim(dnmt3a)
write_csv(dnmt3a, file = "dnmt3a_summary_stats_cis_epo_region.csv")

### TET2:
### Filter SNPs ±100 kb of EPO gene [chr7:100,720,468-100,723,700 in GRCh38]:
tet2 <- tet2 %>%
  filter(chromosome == 7)
dim(tet2)
tet2 <- tet2 %>%
  filter(base_pair_location >= 100720468 - 100000,  
         base_pair_location <= 100723700 + 100000)
dim(tet2)

### Create beta and SE columns:
tet2$beta <- log(tet2$odds_ratio)
tet2$se <- (((log(tet2$ci_upper)) - (log(tet2$odds_ratio))) / 1.96)

### Format tet2:
tet2 <- tet2 %>%
  mutate(SNP = paste(chromosome, base_pair_location, sep = "_")) %>%
  select(SNP, chromosome, base_pair_location, other_allele, effect_allele, effect_allele_frequency, beta, se, p_value) %>%
  rename("other_allele.tet2" = "other_allele",
         "effect_allele.tet2" = "effect_allele",
         "eaf.tet2" = "effect_allele_frequency",
         "B.tet2" = "beta",
         "SE.tet2" = "se",
         "P.tet2" = "p_value")
write_csv(tet2, file = "tet2_summary_stats_cis_epo_region.csv")