### Set working directory and load required packages:
setwd()
library(tidyverse)
library(vroom)
library(preprocessCore)
library(RNOmni)
library(CEMiTool)
library(broom)
library(msigdbr)
library(fgsea)
library(ggpubr)
library(irr)
library(cowplot)

### Load required data:
clinical1 <- read_delim("TCGA.UCEC.sampleMap_UCEC_clinicalMatrix.txt", delim = "\t")
gene_exp1 <- vroom("UCEC.uncv2.mRNAseq_raw_counts.txt", delim = "\t")
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
### (ii) Healthy weight (BMI 18.5 - 24.9) ~ 2
### (iii) Overweight (BMI 25.0 - 29.9) ~ 3
### (iv) Obese (BMI 30.0 - 39.9) ~ 4
### (v) Severely obese (BMI >= 40.0) ~ 5
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

### Formatting gene_exp1:
gene_exp2 <- gene_exp1 %>%
  drop_na()
gene_exp3 <- gene_exp2 %>%
  gather(tcgaid, val, 2:ncol(gene_exp2)) %>%
  spread(`HYBRIDIZATION R`, val)
gene_exp3 <- gene_exp3 %>% 
  separate(tcgaid, into = c("tcgaid1","tcgaid2","tcgaid3","primary"), sep = "-")
gene_exp3 <- gene_exp3 %>%
  unite(tcgaid, tcgaid1, tcgaid2, tcgaid3, sep = "-")
gene_exp3$primary <- as_factor(gene_exp3$primary)
levels(gene_exp3$primary)
gene_exp3 <- gene_exp3 %>% 
  filter(str_detect(primary, "11") != TRUE)
which(duplicated(gene_exp3$tcgaid))
which(duplicated(gene_exp3$tcgaid, fromLast = TRUE))
gene_exp3[305,]
gene_exp3[306,]
gene_exp3 <- gene_exp3 %>% 
  select(-primary)
gene_exp3 <- gene_exp3 %>%
  distinct(`tcgaid`, .keep_all = TRUE)
gene_exp4 <- t(gene_exp3)
gene_exp4 <- as.data.frame(gene_exp4[-1,])
gene_exp4 <- gene_exp4 %>% 
  mutate_all(as.numeric)
gene_exp4 <- filter_genes(gene_exp4, pct = 0.90, apply_vst = FALSE)
gene_exp4 <- t(gene_exp4)
gene_exp4 <- as.matrix(gene_exp4)

### Normalising gene_exp:
gene_exp5 <- preprocessCore::normalize.quantiles(as.matrix(gene_exp4))
gene_exp6 <- t(apply(gene_exp5, 1, RankNorm))
colnames(gene_exp6) <- colnames(gene_exp4)
gene_exp6 <- as_tibble(gene_exp6)
gene_exp6 <- gene_exp6 %>% 
  add_column(tcgaid = gene_exp3$tcgaid, .before = 1)
write_csv(gene_exp6, file = "formatted_gene_exp_data.csv")

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

### Merging data:
merge_gene_exp_covars <- pgs %>% 
  inner_join(clinical6, by = "tcgaid") %>%
  inner_join(ancestry2, by = "tcgaid") %>% 
  inner_join(gene_exp6, by = "tcgaid")
merge_gene_exp_covars <- merge_gene_exp_covars %>%
  select(-pgs)
merge_gene_exp_covars %>%
  group_by(bmi_cat) %>%
  count()

### Running multivariable regression:
merge_gene_exp_covars_result1 <- merge_gene_exp_covars %>%
  gather(measure, value, -c(tcgaid, bmi_cat, age, ms_status, stage, PC1:PC10)) %>%
  nest(-measure) %>%
  mutate(fit = map(data, ~ glm(value ~ bmi_cat + age + ms_status + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, gaussian(link = "identity"), data = .x)),
         tidied = map(fit, tidy)) %>%
  unnest(tidied)
merge_gene_exp_covars_result2 <- merge_gene_exp_covars_result1 %>%
  filter(term == "bmi_cat")
merge_gene_exp_covars_result2 <- merge_gene_exp_covars_result2 %>%
  arrange(p.value)
merge_gene_exp_covars_result2$fdr <- p.adjust(merge_gene_exp_covars_result2$p.value, method = "fdr")
merge_gene_exp_covars_result2$measure <- str_remove(merge_gene_exp_covars_result2$measure, "\\|[0-9]*")
merge_gene_exp_covars_result2 <- merge_gene_exp_covars_result2 %>%
  filter(measure != "?")
merge_gene_exp_covars_result2 <- merge_gene_exp_covars_result2 %>% 
  distinct(measure, .keep_all = TRUE)
merge_gene_exp_covars_result2 <- merge_gene_exp_covars_result2 %>%
  select(-c(data, fit))
merge_gene_exp_covars_result2
merge_gene_exp_covars_result2 <- vroom("act_bmi_four_cat_pgs_merge_gene_exp_covars_result2.csv")

### View MSigDB pathways:
hallmark <- msigdbr(species = "human", category = "H")
hallmark <- split(x = hallmark$gene_symbol, f = hallmark$gs_name)
rankstats <- as_vector(abs(merge_gene_exp_covars_result2$statistic))
names(rankstats) <- merge_gene_exp_covars_result2$measure
set.seed(149)
fgsea_hallmark <- fgsea(pathways = hallmark,
                        stats    = rankstats,
                        minSize  = 15,
                        maxSize  = 500,
                        scoreType = "pos")
fgsea_hallmark <- fgsea_hallmark %>%
  arrange(padj)
View(fgsea_hallmark)

### Determining enriched gene-pathway associaitons:
fgsea_hallmark_plot_data <- fgsea_hallmark %>% 
  arrange(-padj)
fgsea_hallmark_plot_data$log <- -log10(fgsea_hallmark_plot_data$padj)
write_csv(fgsea_hallmark_plot_data, file = "act_bmi_pgs_four_cat_fgsea_hallmark_plot_data.csv")
fgsea_hallmark_plot_data <- vroom("act_bmi_pgs_four_cat_fgsea_hallmark_plot_data.csv")
fgsea_hallmark_plot_data$pathway <- gsub("HALLMARK_", "", fgsea_hallmark_plot_data$pathway)
x <- -log10(0.05)
fgsea_hallmark_bar <- ggbarplot(fgsea_hallmark_plot_data, 
                                x = "pathway", 
                                y = "log", 
                                fill = "NES",
                                legend = "right",
                                ylab = "-log10(FDR) association with BMI category", 
                                xlab = "Hallmark pathways",
                                width = 0.8,
                                orientation = "horizontal") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position = c(0.85, 0.175),
        legend.key.size = unit(0.8, 'cm'), 
        legend.key.height = unit(0.8, 'cm'), 
        legend.key.width = unit(0.8, 'cm'),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
  scale_x_discrete(expand = c(0.01, 0.1)) +
  geom_hline(aes(yintercept = x), color = "black")
fgsea_hallmark_bar
ggsave(fgsea_hallmark_bar, file = "act_bmi_pgs_four_cat_fgsea_hallmark_bar.pdf")

### Assign statistics to lead genes from E2F targets hallmark: 
mycv2_lead <- as_tibble(fgsea_hallmark$leadingEdge[1][[1]])
mycv2_lead <- mycv2_lead %>% 
  rename(measure = value)
mycv2_lead_result <- mycv2_lead %>%
  inner_join(merge_gene_exp_covars_result2, by = "measure")
mycv2_lead_result <- mycv2_lead_result %>%
  add_column(sign_gsea = if_else(.$statistic < 0, "Downregulation", "Upregulation"))
mycv2_lead_result <- mycv2_lead_result %>%
  rename("regulation" = sign_gsea)
write_csv(mycv2_lead_result, file = "act_bmi_pgs_four_cat_mycv2_lead_result.csv")

### Plot mycv1 result:
mycv2_bar_gsea <- ggbarplot(mycv2_lead_result, 
                            x = "measure",  
                            y = "statistic",
                            ylab = "Tumor expression Z-score with increasing BMI category", 
                            xlab = "MYC targets V2 leading edge genes (N=39 genes)") +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) +
  scale_x_discrete(expand = c(0.01, 0.1)) +
  rotate()
mycv2_bar_gsea
ggsave(mycv2_bar_gsea, file = "act_bmi_pgs_four_cat_mycv2_bar_gsea.pdf")

### E2F targets:
### Assign statistics to lead genes from E2F targets hallmark: 
e2f_lead <- as_tibble(fgsea_hallmark$leadingEdge[2][[1]])
e2f_lead <- e2f_lead %>% 
  rename(measure = value)
e2f_lead_result <- e2f_lead %>%
  inner_join(merge_gene_exp_covars_result2, by = "measure")
e2f_lead_result <- e2f_lead_result %>%
  add_column(sign_gsea = if_else(.$statistic < 0, "Downregulation", "Upregulation"))
e2f_lead_result <- e2f_lead_result %>%
  rename("regulation" = sign_gsea)
write_csv(e2f_lead_result, file = "act_bmi_pgs_four_cat_e2f_lead_result.csv")

### Plot e2f result:
e2f_bar_gsea <- ggbarplot(e2f_lead_result, 
                          x = "measure",  
                          y = "statistic",
                          ylab = "Tumor expression Z-score with increasing BMI category", 
                          xlab = "E2F targets leading edge genes (N=134 genes)") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  scale_x_discrete(expand = c(0.01, 0.1)) +
  rotate()
e2f_bar_gsea
ggsave(e2f_bar_gsea, file = "act_bmi_pgs_four_cat_e2f_bar_gsea.pdf")

### Assign statistics to lead genes from E2F targets hallmark: 
mycv1_lead <- as_tibble(fgsea_hallmark$leadingEdge[3][[1]])
mycv1_lead <- mycv1_lead %>% 
  rename(measure = value)
mycv1_lead_result <- mycv1_lead %>%
  inner_join(merge_gene_exp_covars_result2, by = "measure")
mycv1_lead_result <- mycv1_lead_result %>%
  add_column(sign_gsea = if_else(.$statistic < 0, "Downregulation", "Upregulation"))
mycv1_lead_result <- mycv1_lead_result %>%
  rename("regulation" = sign_gsea)
write_csv(mycv1_lead_result, file = "act_bmi_pgs_four_cat_mycv1_lead_result.csv")

### Plot mycv1 result:
mycv1_bar_gsea <- ggbarplot(mycv1_lead_result, 
                            x = "measure",  
                            y = "statistic",
                            ylab = "Tumor expression Z-score with increasing BMI category", 
                            xlab = "MYC targets V1 leading edge genes (N=121 genes)") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  scale_x_discrete(expand = c(0.01, 0.1)) +
  rotate()
mycv1_bar_gsea
ggsave(mycv1_bar_gsea, file = "act_bmi_pgs_four_cat_mycv1_bar_gsea.pdf")

### Assign statistics to lead genes from E2F targets hallmark: 
allograft_lead <- as_tibble(fgsea_hallmark$leadingEdge[4][[1]])
allograft_lead <- allograft_lead %>% 
  rename(measure = value)
allograft_lead_result <- allograft_lead %>%
  inner_join(merge_gene_exp_covars_result2, by = "measure")
allograft_lead_result <- allograft_lead_result %>%
  add_column(sign_gsea = if_else(.$statistic < 0, "Downregulation", "Upregulation"))
allograft_lead_result <- allograft_lead_result %>%
  rename("regulation" = sign_gsea)
write_csv(allograft_lead_result, file = "act_bmi_pgs_four_cat_allograft_lead_result.csv")

### Plot allograft result:
allograft_bar_gsea <- ggbarplot(allograft_lead_result, 
                                x = "measure",  
                                y = "statistic",
                                ylab = "Tumor expression Z-score with increasing BMI category", 
                                xlab = "Allograft rejection leading edge genes (N=112 genes)") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  scale_x_discrete(expand = c(0.01, 0.1)) +
  rotate()
allograft_bar_gsea
ggsave(allograft_bar_gsea, file = "act_bmi_pgs_four_cat_allograft_bar_gsea.pdf")

### Assign statistics to lead genes from E2F targets hallmark: 
androgen_lead <- as_tibble(fgsea_hallmark$leadingEdge[5][[1]])
androgen_lead <- androgen_lead %>% 
  rename(measure = value)
androgen_lead_result <- androgen_lead %>%
  inner_join(merge_gene_exp_covars_result2, by = "measure")
androgen_lead_result <- androgen_lead_result %>%
  add_column(sign_gsea = if_else(.$statistic < 0, "Downregulation", "Upregulation"))
androgen_lead_result <- androgen_lead_result %>%
  rename("regulation" = sign_gsea)
write_csv(androgen_lead_result, file = "act_bmi_pgs_four_cat_androgen_lead_result.csv")

### Plot androgen result:
androgen_bar_gsea <- ggbarplot(androgen_lead_result, 
                               x = "measure",  
                               y = "statistic",
                               ylab = "Tumor expression Z-score with increasing BMI category", 
                               xlab = "Androgen response leading edge genes (N=53 genes)") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  scale_x_discrete(expand = c(0.01, 0.1)) +
  rotate()
androgen_bar_gsea
ggsave(androgen_bar_gsea, file = "act_bmi_pgs_four_cat_androgen_bar_gsea.pdf")

### Assign statistics to lead genes from E2F targets hallmark: 
estrogen_lead <- as_tibble(fgsea_hallmark$leadingEdge[6][[1]])
estrogen_lead <- estrogen_lead %>% 
  rename(measure = value)
estrogen_lead_result <- estrogen_lead %>%
  inner_join(merge_gene_exp_covars_result2, by = "measure")
estrogen_lead_result <- estrogen_lead_result %>%
  add_column(sign_gsea = if_else(.$statistic < 0, "Downregulation", "Upregulation"))
estrogen_lead_result <- estrogen_lead_result %>%
  rename("regulation" = sign_gsea)
write_csv(estrogen_lead_result, file = "act_bmi_pgs_four_cat_estrogen_lead_result.csv")

### Plot estrogen result:
estrogen_bar_gsea <- ggbarplot(estrogen_lead_result, 
                               x = "measure",  
                               y = "statistic",
                               ylab = "Tumor expression Z-score with increasing BMI category", 
                               xlab = "Estrogen Response Early leading edge genes (N=81 genes)") +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) +
  scale_x_discrete(expand = c(0.01, 0.1)) +
  rotate()
estrogen_bar_gsea
ggsave(estrogen_bar_gsea, file = "act_bmi_pgs_four_cat_estrogen_bar_gsea.pdf")
