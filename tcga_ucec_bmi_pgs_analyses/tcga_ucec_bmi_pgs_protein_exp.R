### Load required packages:
library(tidyverse)
library(vroom)
library(broom)
library(ggpubr)
library(cowplot)

### Load data:
clinical1 <- read_delim("TCGA.UCEC.sampleMap_UCEC_clinicalMatrix.txt", delim = "\t")
rppa1 <- vroom("TCGA.UCEC.sampleMap_RPPA_RBN", delim = "\t")
ancestry1 <- read_delim("WashU_PCA_ethnicity_assigned.tsv", delim = "\t")
pgs <- read_delim("ucec_bmi_female_primary_std.txt", delim = "\t")

### Formatting clinical1:
colnames(clinical1)
which(duplicated(clinical1$`_PATIENT`))
clinical2 <- clinical1 %>% 
  distinct(`_PATIENT`, .keep_all = TRUE)
clinical3 <- clinical2 %>% 
  select(`_PATIENT`, CDE_ID_3226963, age_at_initial_pathologic_diagnosis, clinical_stage)
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
  select(tcgaid, age, ms_status, stage)

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

### Format rppa1:
rppa2 <- rppa1 %>% 
  drop_na()
rppa2 <- rppa2 %>%
  gather(tcgaid, val, 2:ncol(rppa1)) %>%
  spread(Sample_description, val)
rppa2 <- rppa2 %>%
  arrange(tcgaid)
rppa2 <- rppa2 %>%
  separate(tcgaid, into = c("tcgaid1","tcgaid2","tcgaid3","primary"), sep = "-")
rppa2$primary <- as_factor(rppa2$primary)
levels(rppa2$primary)
rppa2 <- rppa2 %>%
  select(-primary)
rppa2 <- rppa2 %>%
  unite(tcgaid, tcgaid1, tcgaid2, tcgaid3, sep = "-")
which(duplicated(rppa2$tcgaid))
rppa2 <- rppa2 %>%
  distinct(tcgaid, .keep_all = TRUE)

### Format ancestry1:
which(duplicated(ancestry1$Case))
ancestry1 <- ancestry1 %>%
  arrange(Sample)
ancestry2 <- ancestry1 %>%
  distinct(Case, .keep_all = TRUE)
ancestry2 <- ancestry2 %>% 
  select(Case, PC1:PC10)
ancestry2 <- ancestry2 %>%
  rename(tcgaid = Case)

### Merge data:
merge_rppa_covars <- pgs %>%
  inner_join(clinical6, by = "tcgaid") %>%
  inner_join(ancestry2, by = "tcgaid") %>%
  inner_join(rppa2, by = "tcgaid")

### Run multivariable regression:
merge_rppa_covars_result1 <- merge_rppa_covars %>%
  gather(measure, value, -c(tcgaid, pgs, age, ms_status, stage, PC1:PC10)) %>%
  nest(-measure) %>%
  mutate(fit = map(data, ~ glm(value ~ pgs + age + ms_status + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, gaussian(link = "identity"), data = .x)),
         tidied = map(fit, tidy)) %>%
  unnest(tidied)
merge_rppa_covars_result2 <- merge_rppa_covars_result1 %>%
  filter(term == "pgs")
merge_rppa_covars_result2 <- merge_rppa_covars_result2 %>%
  arrange(p.value)
merge_rppa_covars_result2$fdr <- p.adjust(merge_rppa_covars_result2$p.value)
write_csv(merge_rppa_covars_result2, file = "merge_rppa_covars_result2.csv")

### Creating graphs:
merge_rppa_covars_result2$log <- -log10(merge_rppa_covars_result2$fdr)
write_csv(merge_rppa_covars_result2, file = "protein_plot_data.csv")
merge_rppa_covars_result2_selected <- merge_rppa_covars_result2[1:50,]
merge_rppa_covars_result2_selected <- merge_rppa_covars_result2_selected %>%
  arrange(log)
x <- -log10(0.05)
protein_exp_plot <- ggbarplot(merge_rppa_covars_result2_selected, 
                              x = "measure", 
                              y = "log", 
                              ylab = "-log10(FDR) association with BMI PGS", 
                              xlab = "Tumor protein levels",
                              orientation = "horizontal") +
  font("x.text", size = 12) +
  geom_hline(aes(yintercept = x))
protein_exp_plot
ggsave(protein_exp_plot, file = "top_50_protein_exp_plot.pdf")
pgs_egfr_plot <- ggscatter(merge_rppa_covars,
                           x = "pgs",
                           y = "EGFRPY1173",
                           cor.coef = TRUE,
                           cor.coef.size = 6,
                           conf.int = TRUE,
                           conf.int.level = 0.95,
                           add = "loess",
                           xlab = "BMI PGS",
                           ylab = "Tumor EGFR expression") +
  font("x.text", size = 12)
pgs_egfr_plot
ggsave(pgs_egfr_plot, file = "pgs_egfr_plot.pdf")

protein_exp_plot_egfr_plot <- ggdraw() +
  draw_plot(protein_exp_plot, x=0, y=0, width = 0.5, height = 1) +
  draw_plot(pgs_egfr_plot, x=0.5, y=0, width = 0.5, height = 1)
protein_exp_plot_egfr_plot
ggsave(protein_exp_plot_egfr_plot, file = "protein_exp_plot_egfr_plot.pdf")

### MSS/MSI separate analyses:
### Selecting MSS cases:
clinical7 <- clinical5 %>%
  filter(ms_status == "MSS")
dim(clinical7)
clinical7$ms_status <- as.character(clinical7$ms_status)
clinical7 <- clinical7 %>%
  mutate(ms_status = replace(ms_status, ms_status == "MSS", 1))

### Merging data for MSS cases:
merge_rppa_covars_mss <- pgs %>% 
  inner_join(clinical7, by = "tcgaid") %>% 
  inner_join(ancestry2, by = "tcgaid") %>% 
  inner_join(rppa2, by = "tcgaid")

### Selecting MSI-H + MSI-L cases:
clinical8 <- clinical5 %>%
  filter(ms_status != "MSS")
dim(clinical8)
clinical8$ms_status <- as.character(clinical8$ms_status)
clinical8 <- clinical8 %>%
  mutate(ms_status = replace(ms_status, ms_status == "MSI-H", 1)) %>%
  mutate(ms_status = replace(ms_status, ms_status == "MSI-L", 1))

### Merging data for MSI-H and MSI-L cases:
merge_rppa_covars_msih_msil <- pgs %>%
  inner_join(clinical8, by = "tcgaid") %>% 
  inner_join(ancestry2, by = "tcgaid") %>% 
  inner_join(rppa2, by = "tcgaid")

### Run multivariable regression for MSS cases:
merge_rppa_covars_result_mss <- merge_rppa_covars_mss %>%
  gather(measure, value, -c(tcgaid, pgs, age, ms_status, stage, PC1:PC10)) %>%
  nest(-measure) %>%
  mutate(fit = map(data, ~ glm(value ~ pgs + age + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, gaussian(link = "identity"), data = .x)),
         tidied = map(fit, tidy)) %>%
  unnest(tidied)
merge_rppa_covars_result_mss2 <- merge_rppa_covars_result_mss %>%
  filter(term == "pgs")
merge_rppa_covars_result_mss2 <- merge_rppa_covars_result_mss2 %>%
  arrange(p.value)
merge_rppa_covars_result_mss2$fdr <- p.adjust(merge_rppa_covars_result_mss2$p.value)
write_csv(merge_rppa_covars_result_mss2, file = "merge_rppa_covars_result_mss2.csv")

### Run multivariable regression for MSI-H and MSI-L cases:
merge_rppa_covars_result_msih_msil <- merge_rppa_covars_msih_msil %>%
  gather(measure, value, -c(tcgaid, pgs, age, ms_status, stage, PC1:PC10)) %>%
  nest(-measure) %>%
  mutate(fit = map(data, ~ glm(value ~ pgs + age + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, gaussian(link = "identity"), data = .x)),
         tidied = map(fit, tidy)) %>%
  unnest(tidied)
merge_rppa_covars_result_msih_msil2 <- merge_rppa_covars_result_msih_msil %>%
  filter(term == "pgs")
merge_rppa_covars_result_msih_msil2 <- merge_rppa_covars_result_msih_msil2 %>%
  arrange(p.value)
merge_rppa_covars_result_msih_msil2$fdr <- p.adjust(merge_rppa_covars_result_msih_msil2$p.value)
write_csv(merge_rppa_covars_result_msih_msil2, file = "merge_rppa_covars_result_msih_msil2.csv")