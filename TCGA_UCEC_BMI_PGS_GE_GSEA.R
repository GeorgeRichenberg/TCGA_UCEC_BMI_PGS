### Set working directory and load required packages:
setwd("")
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

### Format clinical1:
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

### Remove rows with missing and indeterminate ms status:
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

### Re-define ms status:
clinical5$ms_status <- as.character(clinical5$ms_status)
clinical6 <- clinical5 %>%
  mutate(ms_status = replace(ms_status, ms_status == "MSS", 0)) %>%
  mutate(ms_status = replace(ms_status, ms_status == "MSI-L", 1)) %>%
  mutate(ms_status = replace(ms_status, ms_status == "MSI-H", 1))
clinical6$ms_status <- as.numeric(clinical6$ms_status)

### Format gene_exp1:
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

### Remove 10% of genes with the lowest expression: 
gene_exp4 <- t(gene_exp3)
gene_exp4 <- as.data.frame(gene_exp4[-1,])
gene_exp4 <- gene_exp4 %>% 
  mutate_all(as.numeric)
gene_exp4 <- filter_genes(gene_exp4, pct = 0.90, apply_vst = FALSE)

### Perform quantile normalisation: 
gene_exp4 <- t(gene_exp4)
gene_exp4 <- as.matrix(gene_exp4)
gene_exp5 <- preprocessCore::normalize.quantiles(as.matrix(gene_exp4))

### Perform inverse normal transformation:
gene_exp6 <- t(apply(gene_exp5, 1, RankNorm))
colnames(gene_exp6) <- colnames(gene_exp4)
gene_exp6 <- as_tibble(gene_exp6)
gene_exp6 <- gene_exp6 %>% 
  add_column(tcgaid = gene_exp3$tcgaid, .before = 1)
write_csv(gene_exp6, file = "gene_exp_data_fig2.csv")

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
merge_gene_exp_covars <- pgs %>% 
  inner_join(clinical6, by = "tcgaid") %>% 
  inner_join(ancestry2, by = "tcgaid") %>% 
  inner_join(gene_exp6, by = "tcgaid")

### Run multivariable regression:
merge_gene_exp_covars_result1 <- merge_gene_exp_covars %>%
  gather(measure, value, -c(tcgaid, pgs, age, ms_status, stage, PC1:PC10)) %>%
  nest(-measure) %>%
  mutate(fit = map(data, ~ glm(value ~ pgs + age + ms_status + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, gaussian(link = "identity"), data = .x)),
         tidied = map(fit, tidy)) %>%
  unnest(tidied)
merge_gene_exp_covars_result2 <- merge_gene_exp_covars_result1 %>%
  filter(term == "pgs")
merge_gene_exp_covars_result2 <- merge_gene_exp_covars_result2 %>%
  arrange(p.value)
merge_gene_exp_covars_result2$fdr <- p.adjust(merge_gene_exp_covars_result2$p.value, method = "fdr")
merge_gene_exp_covars_result2$measure <- str_remove(merge_gene_exp_covars_result2$measure, "\\|[0-9]*")
merge_gene_exp_covars_result2 <- merge_gene_exp_covars_result2 %>%
  filter(measure != "?")
merge_gene_exp_covars_result2 <- merge_gene_exp_covars_result2 %>% 
  distinct(measure, .keep_all = TRUE)
write_csv(merge_gene_exp_covars_result2, file = "merge_gene_exp_covars_result2.csv")

### Run GSEA:
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

### Determine enriched pathways:
fgsea_hallmark_plot_data <- fgsea_hallmark %>% 
  arrange(-padj)
fgsea_hallmark_plot_data$log <- -log10(fgsea_hallmark_plot_data$padj)
write_csv(fgsea_hallmark_plot_data, file = "fgsea_hallmark_plot_data.csv")
fgsea_hallmark_plot_data$pathway <- sub("^HALLMARK_", "", fgsea_hallmark_plot_data$pathway)
x <- -log10(0.05)
fgsea_hallmark_bar <- ggbarplot(fgsea_hallmark_plot_data, 
                                x = "pathway", 
                                y = "log", 
                                fill = "NES",
                                legend = "right",
                                ylab = "-log10(FDR) association with BMI PGS", 
                                xlab = "Hallmark pathways",
                                width = 0.8,
                                orientation = "horizontal") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  theme(legend.position = c(0.85, 0.175),
        legend.key.size = unit(0.8, 'cm'), 
        legend.key.height = unit(0.8, 'cm'), 
        legend.key.width = unit(0.8, 'cm'),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
  scale_x_discrete(expand = c(0.01, 0.1)) +
  geom_hline(aes(yintercept = x), color = "black")
fgsea_hallmark_bar
ggsave(fgsea_hallmark_bar, file = "new_fgsea_hallmark_bar.pdf")

### View IL6-JAK-STAT3 pathway: 
il6_jak_stat3_lead <- as_tibble(fgsea_hallmark$leadingEdge[1][[1]])
il6_jak_stat3_lead <- il6_jak_stat3_lead %>% 
  rename(measure = value)
il6_jak_stat3_lead_result <- il6_jak_stat3_lead %>%
  inner_join(merge_gene_exp_covars_result2, by = "measure")
il6_jak_stat3_lead_result <- il6_jak_stat3_lead_result %>%
  add_column(sign_gsea = if_else(.$statistic < 0, "Downregulation", "Upregulation"))
il6_jak_stat3_lead_result <- il6_jak_stat3_lead_result %>%
  rename("regulation" = sign_gsea)
il6_jak_stat3_lead_result_selected <- select(il6_jak_stat3_lead_result, -c(data, fit))
write_csv(il6_jak_stat3_lead_result_selected, file = "il6_jak_stat3_lead_result.csv")

### Plot IL6-JAK-STAT3 result:
il6_jak_stat3_lead_result$measure <- factor(il6_jak_stat3_lead_result$measure,
                                            levels = rev(unique(il6_jak_stat3_lead_result$measure)))
il6_jak_stat3_bar_gsea <- ggbarplot(il6_jak_stat3_lead_result, 
                                    x = "measure", 
                                    y = "statistic",
                                    ylab = "Tumor expression Z-score with higher BMI PGS", 
                                    xlab = "IL6-JAK-STAT3 signaling leading edge genes",
                                    orientation = "vertical") + 
  rotate() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  scale_x_discrete(expand = c(0.01, 0.1))
il6_jak_stat3_bar_gsea
ggsave(il6_jak_stat3_bar_gsea, file = "il6_jak_stat3_bar_gsea.pdf")
il6_jak_stat3_gsea_plot <- plotEnrichment(hallmark[["HALLMARK_IL6_JAK_STAT3_SIGNALING"]],
                                          rankstats) + labs(title="HALLMARK_IL6_JAK_STAT3_SIGNALING")
il6_jak_stat3_gsea_plot
ggsave(il6_jak_stat3_gsea_plot, file = "il6_jak_stat3_gsea_plot.pdf")

### Plot pgs vs IL6-JAK-STAT3 key genes:
merge_gene_exp_covars2 <- merge_gene_exp_covars %>%
  rename("IL6" = "IL6|3569")
merge_gene_exp_covars2 <- merge_gene_exp_covars2 %>%
  rename("STAT3" = "STAT3|6774")
plot_il6 <- ggscatter(merge_gene_exp_covars2, 
                      x = "pgs",
                      y = "IL6", 
                      add = "loess",
                      cor.coef = TRUE,
                      cor.coef.size = 7.5,
                      conf.int = TRUE,
                      conf.int.level = 0.95,
                      color = "black",
                      size = 1,
                      xlab = "BMI PGS", 
                      ylab = "IL6 tumor gene expression") + 
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))
plot_il6
ggsave(plot_il6, file = "pgs_vs_il6_exp.pdf")
plot_stat3 <- ggscatter(merge_gene_exp_covars2, 
                        x = "pgs",
                        y = "STAT3", 
                        add = "loess",
                        cor.coef = TRUE,
                        cor.coef.size = 7.5,
                        conf.int = TRUE,
                        conf.int.level = 0.95,
                        size = 1,
                        color = "black",
                        xlab = "BMI PGS", 
                        ylab = "STAT3 tumor gene expression") +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))
plot_stat3
ggsave(plot_stat3, file = "pgs_vs_stat3_exp.pdf")
ggdraw() +
  draw_plot(il6_jak_stat3_bar_gsea, x=0, y=0, width = 0.6, height = 1) +
  draw_plot(plot_il6, x=0.6, y=0.5, width = 0.4, height = 0.5) +
  draw_plot(plot_stat3, x=0.6, y=0, width = 0.4, height = 0.5)
ggsave(file = "il6_stat3_expression.pdf")