#### Set working directory and load required packages:
setwd()
library(vroom)
library(tidyverse)
library(TwoSampleMR)

### Load data:
sun <- vroom("Sun_epo_ch.csv")
frayling <- vroom("Frayling_epo_ch.csv")

### Format sun:
sun2 <- sun %>%
  mutate(exposure = paste("SUN -", exposure)) %>%
  select(-c(id.exposure, id.outcome, method, nsnp, lo_ci, up_ci))
sun2

### Format frayling:
frayling2 <- frayling %>%
  mutate(exposure = paste("FRAYLING -", exposure)) %>%
  select(-c(id.exposure, id.outcome, method, nsnp, lo_ci, up_ci))
frayling2

### Combine sun2 and frayling2:
epo_ch <- rbind(sun2, frayling2)
epo_ch$outcome <- gsub("-mutant", "", epo_ch$outcome)
epo_ch$outcome <- gsub("overall", "Overall", epo_ch$outcome)

### Format data for plotting:
epo_ch$outcome <- factor(epo_ch$outcome, levels = c("Overall CH", "DNMT3A CH", "TET2 CH"))
epo_ch$exposure <- factor(epo_ch$exposure, levels = c("FRAYLING - rare variant: rs11976235-T", 
                                                      "SUN - rare variant: rs11976235-T", 
                                                      "FRAYLING - common variant: rs1617640-A",
                                                      "SUN - common variant: rs1617640-A"))
epo_ch <- epo_ch %>%
  slice(c(1:3, 7:9, 4:6, 10:12))
write_csv(epo_ch, file = "Sun_Frayling_epo_ch_mr_res_combined.csv")

### Create forest plot:
plot <- ggplot(data = epo_ch, aes(x = exposure, y = or, ymin = or_lci95, ymax = or_uci95, col = exposure)) +
  geom_pointrange(shape = 15, size = 2) +
  geom_hline(yintercept = 1, linetype = 2) +
  xlab(NULL) +
  ylab("OR (95% CI)") +
  geom_errorbar(width = 0.3, cex = 2.0) + 
  facet_wrap(~outcome, strip.position="left", nrow=9, scales = "free_y") +
  scale_color_manual(values = c("SUN - common variant: rs1617640-A" = "dodgerblue",
                                "FRAYLING - common variant: rs1617640-A"  = "lightskyblue",
                                "SUN - rare variant: rs11976235-T" = "darkred",
                                "FRAYLING - rare variant: rs11976235-T" = "indianred")) +
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
    legend.position = "left",
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 18),
  ) +
  coord_flip() +
  guides(color = guide_legend(reverse = TRUE))
plot

### Save:
ggsave(plot, file = "Sun_epo_Frayling_epo_ch_mr_results_forest_plot_compare.pdf")

### Format epo_ch data for forest plot:
epo_ch2 <- epo_ch %>%
  select(-c(b, se)) %>%
  select(outcome, exposure, or, or_lci95, or_uci95, pval) %>%
  mutate(or = round(or, 2),
         or_lci95 = round(or_lci95, 2),
         or_uci95 = round(or_uci95, 2),
         pval = round(pval, 2)) %>%
  unite(or_interval, or_lci95, or_uci95, sep = "-") %>%
  unite("OR (95% CI)", or, or_interval, sep = " ")
write_csv(epo_ch2, "delete_after_use.csv")
