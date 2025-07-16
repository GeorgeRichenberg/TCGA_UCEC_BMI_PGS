### Set working directory and load required packages:
setwd("")
library(tidyverse)
library(vroom)
library(ggpubr)
library(cowplot)

### Load required data:
clinical1 <- read_delim("TCGA.UCEC.sampleMap_UCEC_clinicalMatrix.txt", delim = "\t")
pgs <- read_delim("ucec_bmi_female_primary_std.txt", delim = "\t")

### Format clinical1:
colnames(clinical1)
which(duplicated(clinical1$`_PATIENT`))
clinical2 <- clinical1 %>% 
  distinct(`_PATIENT`, .keep_all = TRUE)
clinical3 <- clinical2 %>%
  select(`_PATIENT`, CDE_ID_3226963, age_at_initial_pathologic_diagnosis, clinical_stage, height, weight)
clinical3$clinical_stage <- as_factor(clinical3$clinical_stage)
levels(clinical3$clinical_stage)
clinical3 <- clinical3 %>% 
  mutate(stage = recode(clinical_stage, "Stage I"="1", 
                        "Stage II"="2", 
                        "Stage III"="3", 
                        "Stage IV"="4", 
                        "Stage IA"="1", 
                        "Stage IB"="1", 
                        "Stage IC"="1", 
                        "Stage IIA"="2", 
                        "Stage IIB"="2",
                        "Stage IIC"="2",
                        "Stage IIIA"="3", 
                        "Stage IIIB"="3", 
                        "Stage IIIC"="3", 
                        "Stage IIIC1"="3", 
                        "Stage IIIC2"="3", 
                        "Stage IVA"="4",
                        "Stage IVB"="4"))
clinical3$stage <- as.character(clinical3$stage)
clinical3 <- clinical3 %>% 
  rename(tcgaid = `_PATIENT`, 
         age = age_at_initial_pathologic_diagnosis,
         ms_status = CDE_ID_3226963)

### Calculating BMI:
clinical_bmi <- clinical3 %>%
  mutate(height_m =  height / 100)
clinical_bmi <- clinical_bmi %>%
  mutate(height_squared =  height_m * height_m)
clinical_bmi <- clinical_bmi %>%
  mutate(bmi = weight/height_squared)

### Remove cases with BMI<15 and BMI>50:
clinical_bmi2 <- clinical_bmi %>%
  filter(bmi>15) %>%
  filter(bmi<50)
min(clinical_bmi2$bmi)
max(clinical_bmi2$bmi)

### Merge pgs and bmi values:
combined_bmi <-  pgs %>%
  inner_join(clinical_bmi2, by = "tcgaid")
combined_bmi_selected <- combined_bmi %>% 
  select(tcgaid, pgs, bmi)
combined_bmi_distinct <- combined_bmi_selected %>%
  drop_na()
dim(combined_bmi_distinct)

### Categorise low vs high BMI PGS:
low_threshold <- quantile(combined_bmi_distinct$pgs, 0.1)
high_threshold <- quantile(combined_bmi_distinct$pgs, 0.9)

### Density plot for BMI PGS:
pgs_data <- ggplot_build(ggplot(combined_bmi_distinct, aes(x = pgs)) + 
                               geom_density())$data[[1]]
bmi_pgs_density_plot <- ggplot() +
  geom_line(data = pgs_data, aes(x = x, y = y), color = "black") +
  geom_area(data = filter(pgs_data, x > low_threshold & x < high_threshold), aes(x = x, y = y), fill = "gray", alpha = 0.8) +
  geom_ribbon(data = filter(pgs_data, x <= low_threshold), aes(x = x, ymin = 0, ymax = y), fill = "lightblue", alpha = 0.5) +
  geom_ribbon(data = filter(pgs_data, x >= high_threshold), aes(x = x, ymin = 0, ymax = y), fill = "skyblue3", alpha = 0.5) +
  geom_vline(xintercept = low_threshold, linetype = "dashed", color = "lightblue", size = 1.2) +
  geom_vline(xintercept = high_threshold, linetype = "dashed", color = "skyblue3", size = 1.2) +
  geom_vline(xintercept = mean(combined_bmi_distinct$pgs), linetype = "dashed", color = "black", size = 1.2) +
  geom_text(aes(x = mean(combined_bmi_distinct$pgs), y = max(pgs_data$y) * 0.1, label = "Mean BMI PGS"), vjust = -80, hjust = 1.1, color = "black") +
  geom_label(aes(x = low_threshold, y = max(pgs_data$y) * 0.1, label = "UCEC cases with low\n genetically predicted BMI"), 
             vjust = -1, hjust = 1.1, color = "black", size = 3.5, fill = "lightblue") +
  geom_label(aes(x = high_threshold, y = max(pgs_data$y) * 0.1, label = "UCEC cases with high\n genetically predicted BMI"), 
             vjust = -1, hjust = -0.1, color = "black", size = 3.5, fill = "skyblue3") +
  labs(x = "BMI PGS (N=322 UCEC cases)", y = "Density") +
  theme_classic()
bmi_pgs_density_plot
ggsave(bmi_pgs_density_plot, file = "high_vs_low_bmi_pgs_density_plot.pdf")

### Categorise the actual BMI into groups:
combined_bmi_distinct <- combined_bmi_distinct %>%
  mutate(bmi_cat = case_when(
    bmi < 18.5 ~ "Underweight",
    between(bmi, 18.5, 24.9) ~ "Healthy weight",
    between(bmi, 25.0, 29.9) ~ "Overweight",
    between(bmi, 30.0, 39.9) ~ "Obese",
    bmi >= 40.0 ~ "Severely obese"))

### Density plot for actual BMI:
bmi_data <- ggplot_build(ggplot(combined_bmi_distinct, aes(x = bmi)) + 
                               geom_density())$data[[1]]
bmi_density_plot <- ggplot() + 
  geom_line(data = bmi_data, aes(x = x, y = y), color = "black", ) +
  geom_area(data = bmi_data, aes(x = x, y = y), alpha = 0.6, fill ="gray") +
  geom_ribbon(data = subset(bmi_data, x < 18.5), aes(x = x, ymin = 0, ymax = y), fill = "grey60", alpha = 0.5) +
  geom_ribbon(data = subset(bmi_data, x >= 18.5 & x <= 24.9), aes(x = x, ymin = 0, ymax = y), fill = "grey80", alpha = 0.5) +
  geom_ribbon(data = subset(bmi_data, x >= 25.0 & x <= 29.9), aes(x = x, ymin = 0, ymax = y), fill = "lightblue1", alpha = 0.5) +
  geom_ribbon(data = subset(bmi_data, x >= 30.0 & x <= 39.9), aes(x = x, ymin = 0, ymax = y), fill = "skyblue3", alpha = 0.5) +
  geom_ribbon(data = subset(bmi_data, x >= 40), aes(x = x, ymin = 0, ymax = y), fill = "royalblue1", alpha = 0.5) +
  geom_vline(xintercept = mean(combined_bmi_distinct$bmi), linetype = "dashed", color = "black", size = 1.2) +
  geom_text(aes(x = mean(combined_bmi_distinct$bmi), y = max(bmi_data$y) * 0.1, label = "Mean BMI at diagnosis"), 
            vjust = -80, hjust = 1.1, color = "black") +
  labs(x = "BMI at diagnosis (N=322 UCEC cases)", y = "Density") +
  theme_classic()
bmi_density_plot

### Combine plots and save:
bmi_density_combined <- ggdraw() +
  draw_plot(bmi_pgs_density_plot, x = 0, y = 0, width = 0.5, height = 1) +
  draw_plot(bmi_density_plot, x = 0.5, y = 0, width = 0.5, height = 1)
bmi_density_combined
ggsave(bmi_density_combined, file = "tcga_ucec_bmi_pgs_act_bmi_density_plots.pdf")