### Set working directory and load required packages:
setwd()
library(tidyverse)
library(vroom)

### Load data:
epo_ch_eur <- vroom("Sun_et_al_2023_cis_epo_snp_rs1617640_overall_mr_results.csv", col_select = c("outcome", "exposure", "method", "nsnp", "b", "se", "pval", "or", "or_lci95", "or_uci95"))
epo_ch_non_eur <- vroom("epo_common_rs1617640_sas_afr_eas_ch_mr_res.csv", col_select = c("outcome", "exposure", "method", "nsnp", "b", "se", "pval", "or", "or_lci95", "or_uci95"))

### Combine data:
epo_ch <- rbind(epo_ch_eur, epo_ch_non_eur)

### Add SNP names:
epo_ch <- epo_ch %>%
  mutate(variant = "common variant: rs1617640-A")

### Rename outcome and exposure columns:
epo_ch <- epo_ch %>%
  mutate(outcome = case_when(
    outcome == "outcome" ~ "Overall CH (EUR)\n(EAF = 0.60)",
    outcome == "sas_overall_ch" ~ "Overall CH (SAS)\n(EAF = 0.64)",
    outcome == "afr_overall_ch" ~ "Overall CH (AFR)\n(EAF = 0.67)",
    outcome == "eas_overall_ch" ~ "Overall CH (EAS)\n(EAF = 0.80)",
    TRUE ~ outcome),
    exposure = case_when(
      exposure == "exposure" ~ "EPO",
      exposure == "epo_common_rs1617640" ~ "EPO",
      TRUE ~ exposure))

### Prepare data for plotting:
epo_ch$outcome <- factor(epo_ch$outcome, levels = c(
  "Overall CH (EUR)\n(EAF = 0.60)", 
  "Overall CH (AFR)\n(EAF = 0.67)", 
  "Overall CH (EAS)\n(EAF = 0.80)", 
  "Overall CH (SAS)\n(EAF = 0.64)"))

### Plot forest plot:
plot <- ggplot(data = epo_ch, aes(x = variant, y = or, ymin = or_lci95, ymax = or_uci95)) + 
  geom_pointrange(aes(col = variant), shape = 15, size = 2) + 
  geom_hline(yintercept = 1, linetype = 2) + 
  xlab(NULL) + 
  ylab("OR (95% CI)") + 
  geom_errorbar(aes(ymin = pmax(or_lci95, 0), ymax = pmin(or_uci95, 6), col = variant), width = 0.25, cex = 2) +  
  facet_wrap(~outcome, strip.position = "left", nrow = 9, scales = "free_y") + 
  theme(
    axis.title.x = element_text(face = "bold", size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_blank(),
    strip.text.y = element_text(hjust = 0.5, vjust = 0.5, angle = 90),
    panel.background = element_rect(fill = "white"), # Set background color
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background = element_rect(color = "black", fill = "black"),
    strip.text = element_text(color = "white", face = "bold", size = 20),
    legend.background = element_rect(color = "black", size = 1.0),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 18)
  ) + 
  ylim(0, 6) + 
  scale_color_manual(values = c("darkturquoise")) + 
  coord_flip() + 
  guides(color = guide_legend(reverse = TRUE))
plot

### Save plot:
ggsave(plot, file = "Sun_et_al_2023_epo_ch_non_eur_forest_plot.pdf")

### Format necessary data for figure and save:
epo_ch$`OR (95% CI)` <- paste0(
  round(epo_ch$or, 2), " (epo", 
  round(epo_ch$or_lci95, 2), " - ", 
  round(epo_ch$or_uci95, 2), ") ")
epo_ch$pval <- round(epo_ch$pval, 2)
write_csv(epo_ch, file = "fig3_epo_ch_non_eur_data.csv")