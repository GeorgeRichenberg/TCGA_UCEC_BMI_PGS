### Set working directory and load required packages:
setwd("")
library(tidyverse)
library(vroom)
library(ggpubr)
library(cowplot)

### Load required data:
clinical1 <- read_delim("TCGA.UCEC.sampleMap_UCEC_clinicalMatrix.txt", delim = "\t")
pgs <- read_delim("ucec_bmi_female_primary_std.txt", delim = "\t")
levine <- vroom("levine.csv")

### Plot by histological subtype:

### Format clinical1:
colnames(clinical1)
which(duplicated(clinical1$`_PATIENT`))
clinical2 <- clinical1 %>% 
  distinct(`_PATIENT`, .keep_all = TRUE)
clinical3 <- clinical2 %>% 
  select(`_PATIENT`, histological_type)
clinical3 <- clinical3 %>% 
  rename(tcgaid = `_PATIENT`, 
         subtype = histological_type)

### Join with pgs:
pgs_clinical <- clinical3 %>%
  inner_join(pgs, by = "tcgaid")
dim(pgs_clinical)

### Divide pgs into quintiles:
pgs_clinical2 <- pgs_clinical %>%
  arrange(pgs)
pgs_clinical3 <- pgs_clinical2 %>%
  mutate(quintile = ntile(pgs, 5))

### Rename subtypes:
subtypes <- pgs_clinical3 %>%
  distinct(subtype)
subtypes
pgs_clinical4 <- pgs_clinical3 %>% 
  mutate(subtype = case_when(subtype == "Endometrioid endometrial adenocarcinoma" ~ "Endometrioid",
                             subtype == "Mixed serous and endometrioid" ~ "Mixed",
                             subtype == "Serous endometrial adenocarcinoma" ~ "Serous"))
pgs_clinical4

### Counts per subtype per quintile:
pgs_clinical5 <- pgs_clinical4 %>%
  group_by(quintile, subtype) %>%
  summarise(count_subtype = n())
pgs_clinical5

### Calculate the total count of samples in each quintile:
total_count_per_quintile <- pgs_clinical5 %>%
  group_by(quintile) %>%
  summarise(total_count = sum(count_subtype))
total_count_per_quintile

### Join the two data frames to get the proportion of each subtype per quintile:
histological_proportions <- pgs_clinical5 %>%
  left_join(total_count_per_quintile, by = "quintile") %>%
  mutate(proportion = count_subtype / total_count * 100)
histological_proportions

### Create plot:
histological_subtypes_plot <- ggbarplot(data = histological_proportions,
                                        x = "quintile",
                                        y = "proportion",
                                        fill = "subtype",
                                        label = histological_proportions$count_subtype,
                                        lab.hjust = 0.5,
                                        lab.vjust = 1.2,
                                        lab.size = 5,
                                        xlab = "BMI PGS quintile",
                                        ylab = "Histological subtype proportion (%)") + 
  
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))
histological_subtypes_plot
ggsave(histological_subtypes_plot, file = "histological_subtype_proportions.pdf")

### Plot graph by TCGA subtype

### Combine pgs and levine:
pgs_levine <- pgs %>%
  inner_join(levine, by = "tcgaid")
dim(pgs_levine)
not_pgs_levine <- pgs %>%
  anti_join(levine, by = "tcgaid") %>%
  mutate(subtype = "Not available")
dim(not_pgs_levine)

### Combine pgs_levine and not_pgs_levine:
pgs_levine2 <- bind_rows(pgs_levine, not_pgs_levine)

### Divide pgs into quintiles:
pgs_levine2 <- pgs_levine2 %>%
  arrange(pgs)
pgs_levine2 <- pgs_levine2 %>%
  mutate(quintile = ntile(pgs, 5))

### Plot Proportions of subtypes by quintile:
pgs_levine3 <- pgs_levine2 %>%
  group_by(quintile, subtype) %>%
  summarise(count_subtype = n())
pgs_levine3

### Calculate the total count of samples in each quintile:
total_count_per_quintile <- pgs_levine3 %>%
  group_by(quintile) %>%
  summarise(total_count = sum(count_subtype))
total_count_per_quintile

### Join the two data frames to get the proportion of each subtype per quintile:
proportions <- pgs_levine3 %>%
  left_join(total_count_per_quintile, by = "quintile") %>%
  mutate(proportion = round(count_subtype / total_count * 100)) %>%
  arrange(quintile, subtype == "Notassigned") %>%
  mutate(subtype = ifelse(subtype == "Notassigned", "Not assigned", subtype))

### Create stacked bar plot:
fill_colors <- c("CN high" = "skyblue", 
                 "CN low" = "indianred1", 
                 "MSI" = "limegreen", 
                 "POLE" = "mediumpurple1", 
                 "Not assigned" = "grey60",  
                 "Not available" = "grey80")
proportions$subtype <- factor(proportions$subtype, levels = c("CN high", "CN low", "MSI", "POLE", "Not assigned", "Not available"))
tcga_subtypes_plot <- ggbarplot(proportions, 
                                "quintile", 
                                "proportion",
                                fill = "subtype",
                                label = proportions$count_subtype,
                                lab.hjust = 0.5,
                                lab.vjust = 1.2,
                                lab.size = 5,
                                xlab = "BMI PGS quintile",
                                ylab = "TCGA subtype proportion (%)") +
  scale_fill_manual(values = fill_colors) + 
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))
tcga_subtypes_plot
ggsave(tcga_subtypes_plot, file = "tcga_subtype_proportions.pdf")

### Arrange plots into one frame:
histological_tcga_plots <- ggdraw() +
  draw_plot(histological_subtypes_plot, x = 0, y = 0.0, width = 0.5, height = 1) +
  draw_plot(tcga_subtypes_plot, x = 0.5, y = 0, width = 0.5, height = 1)
histological_tcga_plots
ggsave(histological_tcga_plots, file = "histological_tcga_subtypes_plot.pdf")

### Plot endoemtrioid and serous endometrial cancer vs bmi pgs:
dim(pgs_clinical)
endo_vs_serous <- pgs_clinical %>%
  filter(subtype == c("Serous endometrial adenocarcinoma", "Endometrioid endometrial adenocarcinoma"))
dim(endo_vs_serous)
endo_vs_serous_plot <- ggboxplot(data = endo_vs_serous,
                                        x = "subtype",
                                        y = "pgs",
                                        add = "jitter",
                                        xlab = "Endometrial cancer subtype",
                                        ylab = "BMI PGS") + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))
endo_vs_serous_plot
ggsave(endo_vs_serous_plot, file = "bmi_pgs_with_endo_vs_serous_plot.pdf")
