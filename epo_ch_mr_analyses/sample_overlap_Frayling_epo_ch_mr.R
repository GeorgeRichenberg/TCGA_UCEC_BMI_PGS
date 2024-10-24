### Set working directory and load required packages:
setwd()
library(vroom)
library(tidyverse)
library(TwoSampleMR)

### Load data:
epo_overall <- vroom("epo_overall_harmonised.csv")
epo_dnmt3a <- vroom("epo_dnmt3a_harmonised.csv")
epo_tet2 <- vroom("epo_tet2_harmonised.csv")

### Filter harmonised data for common SNP (rs1617640):
epo_overall_rs1617640 <- epo_overall %>%
  filter(SNP == "rs1617640")
epo_dnmt3a_rs1617640 <- epo_dnmt3a %>%
  filter(SNP == "rs1617640")
epo_tet2_rs1617640 <- epo_tet2 %>%
  filter(SNP == "rs1617640")

### Filter harmonised data for rare variant (rs11976235):
epo_overall_rs11976235 <- epo_overall %>%
  filter(SNP == "rs11976235")
epo_dnmt3a_rs11976235 <- epo_dnmt3a %>%
  filter(SNP == "rs11976235")
epo_tet2_rs11976235 <- epo_tet2 %>%
  filter(SNP == "rs11976235")

### Run MR for rs1617640:
epo_overall_results_rs1617640 <- generate_odds_ratios(mr(epo_overall_rs1617640))
epo_overall_results_rs1617640
epo_overall_results_rs1617640 <- epo_overall_results_rs1617640 %>%
  mutate(outcome = "overall CH",
         exposure = "common variant: rs1617640-A")
epo_overall_results_rs1617640

epo_dnmt3a_results_rs1617640 <- generate_odds_ratios(mr(epo_dnmt3a_rs1617640))
epo_dnmt3a_results_rs1617640
epo_dnmt3a_results_rs1617640 <- epo_dnmt3a_results_rs1617640 %>%
  mutate(outcome = "DNMT3A-mutant CH",
         exposure = "common variant: rs1617640-A")
epo_dnmt3a_results_rs1617640

epo_tet2_results_rs1617640 <- generate_odds_ratios(mr(epo_tet2_rs1617640))
epo_tet2_results_rs1617640
epo_tet2_results_rs1617640 <- epo_tet2_results_rs1617640 %>%
  mutate(outcome = "TET2-mutant CH",
         exposure = "common variant: rs1617640-A")
epo_tet2_results_rs1617640

### Run MR for rs11976235:
epo_overall_results_rs11976235 <- generate_odds_ratios(mr(epo_overall_rs11976235))
epo_overall_results_rs11976235
epo_overall_results_rs11976235 <- epo_overall_results_rs11976235 %>%
  mutate(outcome = "overall CH",
         exposure = "rare variant: rs11976235-T")
epo_overall_results_rs11976235

epo_dnmt3a_results_rs11976235 <- generate_odds_ratios(mr(epo_dnmt3a_rs11976235))
epo_dnmt3a_results_rs11976235
epo_dnmt3a_results_rs11976235 <- epo_dnmt3a_results_rs11976235 %>%
  mutate(outcome = "DNMT3A-mutant CH",
         exposure = "rare variant: rs11976235-T")
epo_dnmt3a_results_rs11976235

epo_tet2_results_rs11976235 <- generate_odds_ratios(mr(epo_tet2_rs11976235))
epo_tet2_results_rs11976235
epo_tet2_results_rs11976235 <- epo_tet2_results_rs11976235 %>%
  mutate(outcome = "TET2-mutant CH",
         exposure = "rare variant: rs11976235-T")
epo_tet2_results_rs11976235

### Combine data:
epo_ch <- rbind(epo_overall_results_rs1617640, 
                epo_dnmt3a_results_rs1617640,
                epo_tet2_results_rs1617640,
                epo_overall_results_rs11976235,
                epo_dnmt3a_results_rs11976235,
                epo_tet2_results_rs11976235)
### Save:
write_csv(epo_ch, file = "Frayling_epo_ch.csv")