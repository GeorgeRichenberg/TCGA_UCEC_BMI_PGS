### Load required packages:
library(tidyverse)
library(vroom)
library(ggpubr)
library(cowplot)

### Load required data:
clinical1 <- read_delim("TCGA.UCEC.sampleMap_UCEC_clinicalMatrix.txt", delim = "\t")
pgs <- read_delim("ucec_bmi_female_primary_std.txt", delim = "\t")

### Formatting clinical1:
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

### Removing cases with BMI<15 and BMI>50:
clinical_bmi2 <- clinical_bmi %>%
  filter(bmi>15) %>%
  filter(bmi<50)
min(clinical_bmi2$bmi)
max(clinical_bmi2$bmi)

### Merging pgs and BMI value datasets:
combined_bmi <-  pgs %>%
  inner_join(clinical_bmi2, by = "tcgaid")
combined_bmi_selected <- combined_bmi %>% 
  select(tcgaid, pgs, bmi)
combined_bmi_distinct <- combined_bmi_selected %>%
  drop_na()
dim(combined_bmi_distinct)

### Running linear model and plotting BMI PGS vs BMI at diagnosis for cases with 15<BMI<50:
ggscatter(combined_bmi_distinct, 
          x = "pgs",
          y = "bmi", 
          add = "loess",
          cor.coef = TRUE,
          cor.coef.size = 6,
          conf.int = TRUE,
          conf.int.level = 0.95,
          color = "black",
          xlab = "BMI PGS", 
          ylab = "BMI at diagnosis") + font("x.text", size = 12)
ggsave("pgs_vs_15<bmi<50_scatter_plot.pdf")

### Merging PGS and ages:
combined_age <-  pgs %>%
  inner_join(clinical3, by = "tcgaid")
combined_age_selected <- combined_age %>% 
  select(tcgaid, pgs, age)
combined_age_distinct <- combined_age_selected %>%
  drop_na()

### Running linear regression:
pgs_age <- lm(formula = pgs ~ age, data = combined_age_distinct)
summary(pgs_age)

### Merging PGS and stages:
combined_stage <-  pgs %>%
  inner_join(clinical3, by = "tcgaid")
combined_stage_selected <- combined_stage %>% 
  select(tcgaid, pgs, stage)
combined_stage_distinct <- combined_stage_selected %>%
  drop_na()
combined_stage_distinct <- combined_stage_distinct %>%
  arrange(stage)
combined_stage_distinct <- combined_stage_distinct %>%
  mutate(stage = as.numeric(stage))

### Running linear regression:
pgs_stage <- lm(formula = pgs ~ stage, data = combined_stage_distinct)
summary(pgs_stage)

### Removing rows with missing and indeterminate ms status:
dim(clinical3)
clinical_ms <- clinical3 %>%
  drop_na(ms_status)
dim(clinical_ms)
clinical_ms$ms_status <- as_factor(clinical_ms$ms_status)
levels(clinical_ms$ms_status)
which(clinical_ms$ms_status == "Indeterminate")
clinical_ms2 <- clinical_ms %>%
  filter(ms_status != "Indeterminate")
clinical_ms2$ms_status <- as.character(clinical_ms2$ms_status)
clinical_ms3 <- clinical_ms2 %>%
  mutate(ms_status = replace(ms_status, ms_status == "MSS", "MSS")) %>%
  mutate(ms_status = replace(ms_status, ms_status == "MSI-L", "MSI")) %>%
  mutate(ms_status = replace(ms_status, ms_status == "MSI-H", "MSI"))
clinical_ms3$ms_status <- as.character(clinical_ms3$ms_status)

### Merging pgs and ms status:
combined_ms <- pgs %>%
  inner_join(clinical_ms3, by = "tcgaid")
combined_ms_selected <- combined_ms %>% 
  select(tcgaid, pgs, ms_status)
combined_ms_distinct <- combined_ms_selected %>%
  drop_na()
dim(combined_ms_distinct)
combined_ms_distinct <- combined_ms_distinct %>% 
  mutate(ms_status_number = recode(ms_status, 
                                   "MSS"="1", 
                                   "MSI"="2"))
combined_ms_distinct <- combined_ms_distinct %>%
  mutate(ms_status_number = as.numeric(ms_status_number))

### Running linear regression:
pgs_ms_status <- lm(formula = pgs ~ ms_status_number, data = combined_ms_distinct)
summary(pgs_ms_status)

### Plotting results:
age_plot <- ggscatter(combined_age_distinct, 
                      x = "age",
                      y = "pgs", 
                      add = "loess",
                      cor.coef = TRUE,
                      cor.coef.size = 6,
                      conf.int = TRUE,
                      conf.int.level = 0.95,
                      color = "black",
                      xlab = "Age at diagnosis", 
                      ylab = "BMI PGS") + font("x.text", size = 12)
age_plot
ggsave(age_plot, file = "age_vs_pgs_scatter_plot.pdf")
stage_plot <- ggboxplot(combined_stage_distinct, 
                        x = "stage",
                        y = "pgs", 
                        color = "black",
                        xlab = "Cancer stage at diagnosis", 
                        ylab = "BMI PGS") + font("x.text", size = 12)
stage_plot
ggsave(stage_plot, file = "stage_vs_pgs_boxplot.pdf")
ms_plot <- ggboxplot(combined_ms_distinct, 
                     x = "ms_status",
                     y = "pgs",
                     color = "black",
                     xlab = "Microsatellite status", 
                     ylab = "BMI PGS") + font("x.text", size = 12)
ms_plot
ggsave(ms_plot, file = "ms_vs_pgs_boxplot.pdf")

### Combining plots:
ggdraw() +
  draw_plot(ms_plot, x=0.0, y=0, width = 0.5, height = 0.33) +
  draw_plot(stage_plot, x=0, y=0.33, width = 0.5, height = 0.33) +
  draw_plot(age_plot, x=0.0, y=0.66, width = 0.5, height = 0.33)
ggsave(file = "age_stage_ms_plot.pdf")

### Determining median and IQR of BMI:
median(combined_bmi_distinct$bmi)
IQR(combined_bmi_distinct$bmi)

### Determining median and IQR of age:
median(combined_age_distinct$age)
IQR(combined_age_distinct$age)

### Determining percentages of stage and ms_status:
combined_stage_distinct %>% 
  group_by(stage) %>% 
  summarise(percent = 100 * n() / nrow(combined_stage_distinct))

combined_ms_distinct %>% 
  group_by(ms_status) %>% 
  summarise(percent = 100 * n() / nrow(combined_ms_distinct))



