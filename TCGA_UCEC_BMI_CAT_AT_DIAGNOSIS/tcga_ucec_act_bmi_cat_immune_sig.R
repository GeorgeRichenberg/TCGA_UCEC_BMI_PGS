### Set working directory and load required packages:
setwd()
library(tidyverse)
library(vroom)
library(broom)
library(ggpubr)
library(cowplot)

### Load required data:
clinical1 <- read_delim("TCGA.UCEC.sampleMap_UCEC_clinicalMatrix.txt", delim = "\t")
cibersort1 <- vroom("TCGA.Kallisto.fullIDs.cibersort.relative.tsv", delim = "\t")
ancestry1 <- read_delim("WashU_PCA_ethnicity_assigned.tsv", delim = "\t")
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
clinical3 <- clinical3 %>% 
  rename(tcgaid = `_PATIENT`, 
         age = age_at_initial_pathologic_diagnosis,
         ms_status = CDE_ID_3226963)
clinical3 <- clinical3 %>% 
  select(tcgaid, age, ms_status, stage, height, weight)

### Calculating BMI:
clinical3 <- clinical3 %>%
  mutate(height_m =  height / 100)
clinical3 <- clinical3 %>%
  mutate(height_squared =  height_m * height_m)
clinical3 <- clinical3 %>%
  mutate(bmi = weight/height_squared)

### Removing cases with BMI<15 and BMI>50:
clinical3 <- clinical3 %>%
  filter(bmi>15) %>%
  filter(bmi<50) %>%
  select(-c(height, weight, height_m, height_squared))
min(clinical3$bmi)
max(clinical3$bmi)

### Removing rows with missing and indeterminate microsatellite status:
dim(clinical3)
sum(is.na(clinical3$ms_status))
clinical4 <- clinical3 %>%
  drop_na(ms_status)
dim(clinical4)
clinical4$ms_status <- as_factor(clinical4$ms_status)
levels(clinical4$ms_status)
which(clinical4$ms_status == "Indeterminate")
clinical5 <- clinical4 %>%
  filter(ms_status != "Indeterminate")
dim(clinical5)

### Re-defining microsatellite status:
clinical5$ms_status <- as.character(clinical5$ms_status)
clinical6 <- clinical5 %>%
  mutate(ms_status = replace(ms_status, ms_status == "MSS", 0)) %>%
  mutate(ms_status = replace(ms_status, ms_status == "MSI-L", 1)) %>%
  mutate(ms_status = replace(ms_status, ms_status == "MSI-H", 1))
clinical6$ms_status <- as.numeric(clinical6$ms_status)

### Categorising BMI into five groups:
### (i) Underweight (BMI < 18.5) ~ 1
### (ii) Healthy weight (BMI 18.5 - 24.9) ~ 1
### (iii) Overweight (BMI 25.0 - 29.9) ~ 2
### (iv) Obese (BMI 30.0 - 39.9) ~ 3
### (v) Severely obese (BMI >= 40.0) ~ 4
clinical6 <- clinical6 %>%
  mutate(bmi = round(bmi, 1)) %>%
  mutate(bmi_cat = case_when(
    between(bmi, 17.0, 24.9) ~ 1,
    between(bmi, 25.0, 29.9) ~ 2,
    between(bmi, 30.0, 39.9) ~ 3,
    bmi >= 40.0 ~ 4))
sum(is.na(clinical6$bmi_cat))
clinical6 <- clinical6 %>%
  select(-bmi)

### Formatting cibersort1:
str_sub(cibersort1$SampleID,5,5) = "-"
str_sub(cibersort1$SampleID,8,8) = "-"
str_sub(cibersort1$SampleID,13,13) = "-"
cibersort1 <- cibersort1 %>%
  filter(CancerType == "UCEC") %>%
  mutate(tcgaid = str_sub(SampleID, 1, 15))
cibersort2 <- cibersort1 %>% 
  select(tcgaid, B.cells.naive:Neutrophils)
cibersort2 <- cibersort2 %>% 
  arrange(tcgaid)
cibersort2 <- cibersort2 %>% 
  separate(tcgaid, into = c("tcgaid1","tcgaid2","tcgaid3","primary"), sep = "-")
cibersort2$primary <- as_factor(cibersort2$primary)
levels(cibersort2$primary)
cibersort2 <- cibersort2 %>%
  filter(str_detect(primary, "11") != TRUE)
cibersort2 <- cibersort2 %>% 
  select(-primary)
cibersort2 <- cibersort2 %>% 
  unite(tcgaid, tcgaid1, tcgaid2, tcgaid3, sep = "-")
which(duplicated(cibersort2$tcgaid))
cibersort2 <- cibersort2 %>% 
  distinct(tcgaid, .keep_all = TRUE)

### Determining immune signatures with a zero count < 90%:
zero_cibersort <- cibersort2 %>%
  summarise_all(list(~ sum(.==0)))
which(zero_cibersort>(0.90*(nrow(cibersort2))))
which(zero_cibersort<(0.90*(nrow(cibersort2))))
cibersort3 <- cibersort2 %>%
  select(-c("T.cells.CD4.naive", "T.cells.gamma.delta", "Eosinophils"))

### Formatting ancestry data:
which(duplicated(ancestry1$Case))
ancestry1 <- ancestry1 %>% 
  arrange(Sample)
ancestry2 <- ancestry1 %>% 
  distinct(Case, .keep_all = TRUE)
ancestry2 <- ancestry2 %>% 
  select(Case, PC1:PC10)
ancestry2 <- ancestry2 %>% 
  rename(tcgaid = Case)

### Merging data :
merge_cibersort_covars <- pgs %>%
  inner_join(clinical6, by = "tcgaid") %>% 
  inner_join(ancestry2, by = "tcgaid") %>% 
  inner_join(cibersort3, by = "tcgaid")
merge_cibersort_covars <- merge_cibersort_covars %>%
  select(-pgs)
dim(merge_cibersort_covars)
merge_cibersort_covars %>%
  group_by(bmi_cat) %>%
  count()

### Running multivariable quasipoisson regression:
quasip_merge_cibersort_covars_result1 <- merge_cibersort_covars %>%
  gather(measure, value, -c(tcgaid, bmi_cat, age, ms_status, stage, PC1:PC10)) %>%
  nest(-measure) %>%
  mutate(fit = map(data, ~ glm(value ~ bmi_cat + age + ms_status + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, quasipoisson(link = "log"), data = .x)),
         tidied = map(fit, tidy)) %>%
  unnest(tidied)
quasip_merge_cibersort_covars_result2 <- quasip_merge_cibersort_covars_result1 %>%
  filter(term == "bmi_cat")
quasip_merge_cibersort_covars_result2 <- quasip_merge_cibersort_covars_result2 %>%
  arrange(p.value)
quasip_merge_cibersort_covars_result2$fdr <- p.adjust(quasip_merge_cibersort_covars_result2$p.value, method = "fdr")
quasip_merge_cibersort_covars_result2 <- quasip_merge_cibersort_covars_result2 %>%
  select(-c(data, fit))
quasip_merge_cibersort_covars_result2
write_csv(quasip_merge_cibersort_covars_result2, file = "act_bmi_pgs_four_cat_quasip_merge_cibersort_covars_result2.csv")

### Creating bar plot of results:
cibersort_plot_data <- quasip_merge_cibersort_covars_result2 %>% 
  arrange(fdr)
cibersort_plot_data$log <- -log10(quasip_merge_cibersort_covars_result2$fdr)
cibersort_plot_data <- cibersort_plot_data %>% 
  arrange(log)
x <- -log10(0.05)
cibersort_bar <- ggbarplot(cibersort_plot_data, 
                           x = "measure", 
                           y = "log", 
                           fill = "statistic",
                           ylab = "-log10(FDR) association with BMI category at diagnosis", 
                           xlab = "Intra-tumor immune cell infiltrates",
                           orientation = "horizontal") + geom_hline(aes(yintercept = x)) +
  labs(fill = "Z-score") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 16)) +
  theme(legend.position = c(0.75, 0.35),
        legend.key.size = unit(0.8, 'cm'), 
        legend.key.height = unit(0.8, 'cm'), 
        legend.key.width = unit(0.8, 'cm'),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))
cibersort_bar
ggsave(cibersort_bar, file = "act_bmi_pgs_four_cat_cibersort_bar.pdf")

### Removing zero counts:
dc_act_non_zero <- merge_cibersort_covars %>%
  filter(Dendritic.cells.activated>0)
dim(merge_cibersort_covars)
dim(dc_act_non_zero)
dc_rest_non_zero <- merge_cibersort_covars %>%
  filter(Dendritic.cells.resting>0)
dim(dc_rest_non_zero)

### Creating scatter plot for non-zero data:
dc_act_non_zero_plot <- ggboxplot(dc_act_non_zero, 
                                  x = "bmi_cat",
                                  y = "Dendritic.cells.activated",
                                  color = "black",
                                  xlab = "BMI category (N=150 cases)", 
                                  ylab = "Non-zero intra-tumor\nactivated dendritic cells") +
  font("x.text", size = 15) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 17))
dc_act_non_zero_plot
ggsave(dc_act_non_zero_plot, file = "act_bmi_pgs_cat_dc_act_non_zero_plot.pdf")

dc_rest_non_zero_plot <- ggboxplot(dc_rest_non_zero, 
                                   x = "bmi_cat",
                                   y = "Dendritic.cells.resting",
                                   color = "black",
                                   xlab = "BMI category (N=126 cases)", 
                                   ylab = "Non-zero intra-tumor\nresting dendritic cells") + 
  font("x.text", size = 15) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 17))
dc_rest_non_zero_plot
ggsave(dc_rest_non_zero_plot, file = "act_bmi_pgs_cat_dc_rest_non_zero_plot.pdf")

### Combining plots:
ggdraw() +
  draw_plot(cibersort_bar, x=0, y=0.5, width = 1, height = 0.5) +
  draw_plot(dc_act_non_zero_plot, x=0, y=0, width = 0.5, height = 0.5) +
  draw_plot(dc_rest_non_zero_plot, x=0.5, y=0, width = 0.5, height = 0.5)
ggsave(file = "act_bmi_pgs_four_cat_cibersort_bar_act_rest_dc_plots.pdf")