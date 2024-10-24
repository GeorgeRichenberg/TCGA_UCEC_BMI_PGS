### Set working directory and load required packages:
setwd()
library(tidyverse)
library(vroom)

### Load data:
epo_overall_results_rs1617640 <- vroom("Sun_et_al_2023_cis_epo_snp_rs1617640_overall_mr_results.csv", col_select = c("outcome", "exposure", "method", "nsnp", "b", "se", "pval", "or", "or_lci95", "or_uci95"))
epo_dnmt3a_results_rs1617640 <- vroom("Sun_et_al_2023_cis_epo_snp_rs1617640_dnmt3a_mr_results.csv", col_select = c("outcome", "exposure", "method", "nsnp", "b", "se", "pval", "or", "or_lci95", "or_uci95"))
epo_tet2_results_rs1617640 <- vroom("Sun_et_al_2023_cis_epo_snp_rs1617640_tet2_mr_results.csv", col_select = c("outcome", "exposure", "method", "nsnp", "b", "se", "pval", "or", "or_lci95", "or_uci95"))

epo_overall_results_rs11976235 <- vroom("Sun_et_al_2023_cis_epo_snp_rs11976235_overall_mr_results.csv", col_select = c("outcome", "exposure", "method", "nsnp", "b", "se", "pval", "or", "or_lci95", "or_uci95"))
epo_dnmt3a_results_rs11976235 <- vroom("Sun_et_al_2023_cis_epo_snp_rs11976235_dnmt3a_mr_results.csv", col_select = c("outcome", "exposure", "method", "nsnp", "b", "se", "pval", "or", "or_lci95", "or_uci95"))
epo_tet2_results_rs11976235 <- vroom("Sun_et_al_2023_cis_epo_snp_rs11976235_tet2_mr_results.csv", col_select = c("outcome", "exposure", "method", "nsnp", "b", "se", "pval", "or", "or_lci95", "or_uci95"))

### Add variant names:
epo_overall_results_rs1617640 <- epo_overall_results_rs1617640 %>%
  mutate(variant = "common variant: rs1617640-A")
epo_overall_results_rs11976235 <- epo_overall_results_rs11976235 %>%
  mutate(variant = "rare variant: rs11976235-T")

epo_dnmt3a_results_rs1617640 <- epo_dnmt3a_results_rs1617640 %>%
  mutate(variant = "common variant: rs1617640-A")
epo_dnmt3a_results_rs11976235 <- epo_dnmt3a_results_rs11976235 %>%
  mutate(variant = "rare variant: rs11976235-T")

epo_tet2_results_rs1617640 <- epo_tet2_results_rs1617640 %>%
  mutate(variant = "common variant: rs1617640-A")
epo_tet2_results_rs11976235 <- epo_tet2_results_rs11976235 %>%
  mutate(variant = "rare variant: rs11976235-T")

### Rename outcome and exposure columns:
epo_overall <- rbind(epo_overall_results_rs1617640, epo_overall_results_rs11976235)
epo_overall <- epo_overall %>%
  mutate(outcome = ifelse(outcome == "outcome", "Overall CH", outcome),
         exposure = ifelse(exposure == "exposure", "EPO", exposure))

epo_dnmt3a <- rbind(epo_dnmt3a_results_rs1617640, epo_dnmt3a_results_rs11976235)
epo_dnmt3a <- epo_dnmt3a %>%
  mutate(outcome = ifelse(outcome == "outcome", "DNMT3A CH", outcome),
         exposure = ifelse(exposure == "exposure", "EPO", exposure))

epo_tet2 <- rbind(epo_tet2_results_rs1617640, epo_tet2_results_rs11976235)
epo_tet2 <- epo_tet2 %>%
  mutate(outcome = ifelse(outcome == "outcome", "TET2 CH", outcome),
         exposure = ifelse(exposure == "exposure", "EPO", exposure))

### Combine data:
epo_ch <- rbind(epo_overall, epo_dnmt3a, epo_tet2)
epo_ch

### Plot forest plot:
epo_ch$outcome <- factor(epo_ch$outcome, levels = c("Overall CH", "DNMT3A CH", "TET2 CH"))
epo_ch$variant <- factor(epo_ch$SNP, levels = c("rare variant: rs11976235-T", "common variant: rs1617640-A"))

plot <- ggplot(data = epo_ch, aes(x = variant, y = or, ymin = or_lci95, ymax = or_uci95)) +
  geom_pointrange(aes(col = variant), shape = 15, size = 2) +
  geom_hline(yintercept = 1, linetype = 2) +
  xlab(NULL) +
  ylab("OR (95% CI)") +
  geom_errorbar(aes(ymin = or_lci95, ymax = or_uci95, col = variant), width = 0.3, cex = 2.0) + 
  facet_wrap(~outcome, strip.position="left", nrow=9, scales = "free_y") +
  theme(
    axis.title.x = element_text(face = "bold", size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_blank(),
    strip.text.y = element_text(hjust=0.5, vjust = 0.5, angle = 90),
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background = element_rect(color = "black", fill = "black"),
    strip.text = element_text(color = "white", face = "bold", size = 20),
    legend.background = element_rect(color = "black", size = 1.0),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 18),
  ) +
  coord_flip() +
  guides(color = guide_legend(reverse = TRUE))
plot

### Save plot:
ggsave(plot, file = "Sun_et_al_2023_epo_ch_forest_plot.pdf")
