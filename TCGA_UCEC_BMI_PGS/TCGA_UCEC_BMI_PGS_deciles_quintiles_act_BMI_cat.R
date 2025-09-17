### Set working directory and load required packages:
setwd("")
library(tidyverse)
library(vroom)
library(ggpubr)

### Load required data:
clinical1 <- read_delim("TCGA.UCEC.sampleMap_UCEC_clinicalMatrix.txt", delim = "\t")
pgs <- read_delim("ucec_bmi_female_primary_std.txt", delim = "\t")

### Format clinical1:
colnames(clinical1)
which(duplicated(clinical1$`_PATIENT`))
clinical2 <- clinical1 %>% 
  distinct(`_PATIENT`, .keep_all = TRUE)
clinical3 <- clinical2 %>% 
  select(`_PATIENT`, height, weight)
clinical3 <- clinical3 %>% 
  rename(tcgaid = `_PATIENT`)

### Calculate bmi:
clinical3 <- clinical3 %>%
  mutate(height_m =  height / 100)
clinical3 <- clinical3 %>%
  mutate(height_squared =  height_m * height_m)
clinical3 <- clinical3 %>%
  mutate(bmi = weight/height_squared)

### Remove cases with bmi<15 and bmi>50:
clinical3 <- clinical3 %>%
  filter(bmi>15) %>%
  filter(bmi<50) %>%
  select(-c(height, weight, height_m, height_squared))
min(clinical3$bmi)
max(clinical3$bmi)

### Categorise bmi into five groups:
### (i) underweight (bmi < 18.5) ~ 1
### (ii) healthy weight (bmi 18.5 - 24.9) ~ 2
### (iii) overweight (bmi 25.0 - 29.9) ~ 3
### (iv) obese (bmi 30.0 - 39.9) ~ 4
### (v) severely obese (bmi >= 40.0) ~ 5
clinical4 <- clinical3 %>%
  mutate(bmi = round(bmi, 1)) %>%
  mutate(bmi_cat = case_when(
    bmi < 18.5 ~ "Underweight",
    between(bmi, 18.5, 24.9) ~ "Healthy weight",
    between(bmi, 25.0, 29.9) ~ "Overweight",
    between(bmi, 30.0, 39.9) ~ "Obese",
    bmi >= 40.0 ~ "Severely obese"))
sum(is.na(clinical4$bmi_cat))
clinical4 <- clinical4 %>%
  select(-bmi)

### Merge pgs and clinical3:
pgs_clinical <- pgs %>%
  inner_join(clinical4, by = "tcgaid")
dim(pgs_clinical)
pgs_clinical %>%
  group_by(bmi_cat) %>%
  count()

### Order pgs_clinical by pgs and sort into deciles:
pgs_clinical <- pgs_clinical %>%
  arrange(pgs)
pgs_clinical2 <- pgs_clinical %>%
  mutate(decile = ntile(pgs, 10))

### Plot Proportions of bmi categories by decile:
pgs_clinical3 <- pgs_clinical2 %>%
  group_by(decile, bmi_cat) %>%
  summarise(count_cat = n())
pgs_clinical3

### Calculate the total count of samples in each decile:
total_count_per_decile <- pgs_clinical3 %>%
  group_by(decile) %>%
  summarise(total_count = sum(count_cat))
total_count_per_decile

### Calculate proportion of each subtype per decile:
proportions <- pgs_clinical3 %>%
  left_join(total_count_per_decile, by = "decile") %>%
  mutate(proportion = round(count_cat / total_count * 100))

### Plot pgs decile vs bmi cat:
proportions$bmi_cat <- factor(proportions$bmi_cat,
                              levels = c("Underweight",
                                         "Healthy weight",
                                         "Overweight",
                                         "Obese",
                                         "Severely obese"))
fill_colors <- c("Underweight" = "grey60", 
                 "Healthy weight" = "grey80",
                 "Overweight" = "lightblue1",
                 "Obese" = "skyblue3",
                 "Severely obese" = "royalblue1")
pgs_bmi_cat_dec <- ggbarplot(proportions,
                             x = "decile",
                             y = "proportion",
                             fill = "bmi_cat", 
                             label = proportions$count_cat,
                             lab.hjust = 0.5,
                             lab.vjust = 1.2,
                             lab.size = 5,
                             xlab = "BMI PGS decile",
                             ylab = "Proportion of TCGA UCEC cases by actual BMI category at diagnosis(%)") +
  scale_fill_manual(values = fill_colors) + 
  guides(fill=guide_legend(title="BMI category")) +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  scale_x_continuous(breaks = 1:10, labels = 1:10)
pgs_bmi_cat_dec
ggsave(pgs_bmi_cat_dec, file = "tcga_ucec_prop_bmi_cat_by_bmi_pgs_decile.pdf")

### RE-PLOT WITH BMI PGS QUINTILES:

### Order pgs_clinical by pgs and sort into quintiles:
pgs_clinical <- pgs_clinical %>%
  arrange(pgs)
pgs_clinical2 <- pgs_clinical %>%
  mutate(quintile = ntile(pgs, 5))

### Plot proportions of bmi categories by quintile:
pgs_clinical3 <- pgs_clinical2 %>%
  group_by(quintile, bmi_cat) %>%
  summarise(count_cat = n())
pgs_clinical3

### Calculate the total count of samples in each quintile:
total_count_per_quintile <- pgs_clinical3 %>%
  group_by(quintile) %>%
  summarise(total_count = sum(count_cat))
total_count_per_quintile

### Calculate proportion of each subtype per quintile:
proportions <- pgs_clinical3 %>%
  left_join(total_count_per_quintile, by = "quintile") %>%
  mutate(proportion = round(count_cat / total_count * 100))

### Plot pgs quintile by bmi cat:
proportions$bmi_cat <- factor(proportions$bmi_cat,
                              levels = c("Underweight",
                                         "Healthy weight",
                                         "Overweight",
                                         "Obese",
                                         "Severely obese"))
fill_colors <- c("Underweight" = "grey60", 
                 "Healthy weight" = "grey80",
                 "Overweight" = "lightblue1",
                 "Obese" = "skyblue3",
                 "Severely obese" = "royalblue1")
pgs_bmi_cat_quin <- ggbarplot(proportions,
                              x = "quintile",
                              y = "proportion",
                              fill = "bmi_cat", 
                              label = proportions$count_cat,
                              lab.hjust = 0.5,
                              lab.vjust = 1.2,
                              lab.size = 5,
                              xlab = "BMI PGS quintile",
                              ylab = "Proportion of TCGA UCEC cases by actual BMI category at diagnosis (%)") +
  scale_fill_manual(values = fill_colors) + 
  guides(fill=guide_legend(title="BMI category")) +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))
pgs_bmi_cat_quin
ggsave(pgs_bmi_cat_quin, file = "tcga_ucec_prop_bmi_cat_by_bmi_pgs_quintile.pdf")
