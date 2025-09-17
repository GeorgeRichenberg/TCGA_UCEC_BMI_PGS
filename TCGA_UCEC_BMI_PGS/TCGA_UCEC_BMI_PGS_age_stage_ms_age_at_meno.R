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
  select(`_PATIENT`, CDE_ID_3226963, age_at_initial_pathologic_diagnosis, clinical_stage, menopause_status, height, weight)
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

### Calculate BMI:
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

### Plot bmi pgs vs bmi at diagnosis for cases with 15<BMI<50:
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

### Merge pgs and ages:
combined_age <-  pgs %>%
  inner_join(clinical3, by = "tcgaid")
combined_age_selected <- combined_age %>% 
  select(tcgaid, pgs, age)
combined_age_distinct <- combined_age_selected %>%
  drop_na()

### Run linear regression:
pgs_age <- lm(formula = pgs ~ age, data = combined_age_distinct)
summary(pgs_age)

### Merge pgs and stages:
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

### Run linear regression:
pgs_stage <- lm(formula = pgs ~ stage, data = combined_stage_distinct)
summary(pgs_stage)

### Merge pgs and menopause status:
combined_meno <-  pgs %>%
  inner_join(clinical3, by = "tcgaid")
combined_meno_selected <- combined_meno %>% 
  select(tcgaid, pgs, menopause_status)
count <- combined_meno_selected %>%
  group_by(menopause_status) %>%
  count()
count
combined_meno_distinct <- combined_meno_selected %>%
  drop_na() %>%
  filter(menopause_status != "Indeterminate (neither Pre or Postmenopausal)")
combined_meno_distinct <- combined_meno_distinct %>%
  mutate(menopause_status = replace(menopause_status, menopause_status == "Peri (6-12 months since last menstrual period)", "Peri")) %>%
  mutate(menopause_status = replace(menopause_status, menopause_status == "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)", "Post")) %>%
  mutate(menopause_status = replace(menopause_status, menopause_status == "Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)", "Pre"))
combined_meno_distinct <- combined_meno_distinct %>% 
  mutate(meno_number = recode(menopause_status, 
                              "Peri"="1", 
                              "Post"="2",
                              "Pre"="3"))
combined_meno_distinct <- combined_meno_distinct %>%
  mutate(meno_number = as.numeric(meno_number))

### Run linear regression:
pgs_meno <- lm(formula = pgs ~ meno_number, data = combined_meno_distinct)
summary(pgs_meno)

### Remove rows with missing and indeterminate ms status:
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

### Merge pgs and ms status:
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

### Run linear regression:
pgs_ms_status <- lm(formula = pgs ~ ms_status_number, data = combined_ms_distinct)
summary(pgs_ms_status)

### Create plots:
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

combined_stage_distinct2 <- combined_stage_distinct %>%
  mutate(stage = if_else(stage == "1", "1\n [67% of cases\n (237 of 354)]",
                         if_else(stage == "2", "2\n [8% of cases\n (30 of 354)]",
                                 if_else(stage == "3", "3\n [21% of cases\n (73 of 354)]",
                                         "4\n [4% of cases\n (14 of 354)]"))))
stage_plot <- ggboxplot(combined_stage_distinct2, 
                        x = "stage",
                        y = "pgs", 
                        color = "black",
                        xlab = "Cancer stage at diagnosis", 
                        ylab = "BMI PGS") + font("x.text", size = 12)
stage_plot
ggsave(stage_plot, file = "stage_vs_pgs_boxplot.pdf")

combined_ms_distinct2 <- combined_ms_distinct %>%
  mutate(ms_status = if_else(ms_status == "MSS", "MSS\n [58% of cases\n (204 of 351)]",
                             "MSI\n [42% of cases\n (147 of 351)]"))
ms_plot <- ggboxplot(combined_ms_distinct2, 
                     x = "ms_status",
                     y = "pgs",
                     color = "black",
                     xlab = "Microsatellite status", 
                     ylab = "BMI PGS") + font("x.text", size = 12)
ms_plot
ggsave(ms_plot, file = "ms_vs_pgs_boxplot.pdf")

combined_meno_distinct$menopause_status <- factor(combined_meno_distinct$menopause_status, levels = c("Pre", "Post", "Peri"))
combined_meno_distinct2 <- combined_meno_distinct %>%
  mutate(menopause_status = if_else(menopause_status == "Pre", "Pre\n [5% of cases\n (18 of 323)]",
                                    if_else(menopause_status == "Post", "Post\n [91% of cases\n (293 of 323)]",
                                            "Peri\n [4% of cases\n (12 of 323)]")))
menopause_plot <- ggboxplot(combined_meno_distinct2, 
                            x = "menopause_status",
                            y = "pgs", 
                            color = "black",
                            xlab = "Menopause status at diagnosis", 
                            ylab = "BMI PGS") + font("x.text", size = 12)
menopause_plot
ggsave(menopause_plot, file = "menopause_vs_pgs_boxplot.pdf")

### Combine plots:
ggdraw() +
  draw_plot(ms_plot, x=0.0, y=0, width = 0.5, height = 0.5) +
  draw_plot(stage_plot, x=0.5, y=0.5, width = 0.5, height = 0.5) +
  draw_plot(age_plot, x=0, y=0.5, width = 0.5, height = 0.5) +
  draw_plot(menopause_plot, x=0.5, y=0, width = 0.5, height = 0.5)
ggsave(file = "age_stage_ms_meno_plot.pdf")

### Calculate median and iqr of bmi:
median(combined_bmi_distinct$bmi)
IQR(combined_bmi_distinct$bmi)

### Calculate median and iqr of age:
median(combined_age_distinct$age)
IQR(combined_age_distinct$age)

### Calculate percentages of stage, ms_status and meno:
combined_stage_distinct %>%
  count(stage) %>%
  mutate(percent = n / sum(n) * 100)

combined_ms_distinct %>%
  count(ms_status) %>%
  mutate(percent = n / sum(n) * 100)

combined_meno_distinct %>%
  count(menopause_status) %>%
  mutate(percent = n / sum(n) * 100)

### Menopause status by age:
age_meno_selected <- combined_meno %>% 
  select(tcgaid, age, pgs, menopause_status)
count <- age_meno_selected %>%
  group_by(menopause_status) %>%
  count()
count
age_meno_distinct <- age_meno_selected %>%
  drop_na() %>%
  filter(menopause_status != "Indeterminate (neither Pre or Postmenopausal)")
age_meno_distinct <- age_meno_distinct %>%
  mutate(menopause_status = replace(menopause_status, menopause_status == "Peri (6-12 months since last menstrual period)", "Peri")) %>%
  mutate(menopause_status = replace(menopause_status, menopause_status == "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)", "Post")) %>%
  mutate(menopause_status = replace(menopause_status, menopause_status == "Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)", "Pre"))
age_peri <- age_meno_distinct %>%
  filter(menopause_status == "Peri")
dim(age_peri)
age_post <- age_meno_distinct %>%
  filter(menopause_status == "Post")
dim(age_post)
age_pre <- age_meno_distinct %>%
  filter(menopause_status == "Pre")
dim(age_pre)

### Plot age at diagnosis vs bmi pgs per menopause_status:
age_peri_plot <- ggscatter(age_peri, 
                           x = "age",
                           y = "pgs", 
                           add = "loess",
                           cor.coef = TRUE,
                           cor.coef.size = 6,
                           conf.int = FALSE,
                           conf.int.level = 0.95,
                           color = "black",
                           size = 0.85,
                           xlab = "Age at diagnosis of peri-menopausal UCEC cases (N=12 cases)", 
                           ylab = "BMI PGS") + font("x.text", size = 12)
age_peri_plot

age_post_plot <- ggscatter(age_post, 
                           x = "age",
                           y = "pgs", 
                           add = "loess",
                           cor.coef = TRUE,
                           cor.coef.size = 6,
                           conf.int = FALSE,
                           conf.int.level = 0.95,
                           color = "black",
                           size = 0.85,
                           xlab = "Age at diagnosis of post-menopausal UCEC cases (N=293 cases)", 
                           ylab = "BMI PGS") + font("x.text", size = 12)
age_post_plot

age_pre_plot <- ggscatter(age_pre, 
                          x = "age",
                          y = "pgs", 
                          add = "loess",
                          cor.coef = TRUE,
                          cor.coef.size = 6,
                          conf.int = FALSE,
                          conf.int.level = 0.95,
                          color = "black",
                          size = 0.85,
                          xlab = "Age at diagnosis of pre-menopausal UCEC cases (N=18 cases)", 
                          ylab = "BMI PGS") + font("x.text", size = 12)
age_pre_plot

### Combine plots:
ggdraw() +
  draw_plot(age_post_plot, x=0.0, y=0.66, width = 0.5, height = 0.33) +
  draw_plot(age_pre_plot, x=0, y=0.33, width = 0.5, height = 0.33) +
  draw_plot(age_peri_plot, x=0, y=0, width = 0.5, height = 0.33)
ggsave(file = "age_bmi_pgs_post_pre_peri_meno_plot.pdf")
