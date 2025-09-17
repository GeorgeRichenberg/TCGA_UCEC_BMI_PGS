### Set working directory and load required packages:
setwd("")
library(tidyverse)
library(vroom)
library(stringr)
library(ggpubr)
library(gridExtra)

### Load required data:
mut_ind1 <- vroom("mc3_gene_level_UCEC_mc3_gene_level.txt", delim = "\t")
pgs <- read_delim("ucec_bmi_female_primary_std.txt", delim = "\t")

### Format mut_ind1:
mut_ind2 <- mut_ind1 %>%
  subset(sample %in% c("PTEN", 
                       "PIK3CA",
                       "TP53",
                       "ARID1A",
                       "CTNNB1",
                       "PIK3R1",
                       "KRAS", 
                       "KMT2D",
                       "CTCF", 
                       "PPP2R1A",
                       "CHD4",
                       "ZFHX3",
                       "FBXW7"))
mut_ind3 <- mut_ind2 %>%
  drop_na()
mut_ind4 <- mut_ind3 %>%
  gather(tcgaid, val, 2:ncol(mut_ind2)) %>%
  spread(sample, val)
mut_ind5 <- mut_ind4 %>% 
  separate(tcgaid, into = c("tcgaid1","tcgaid2","tcgaid3","primary"), sep = "-")
mut_ind6 <- mut_ind5 %>%
  unite(tcgaid, tcgaid1, tcgaid2, tcgaid3, sep = "-")
mut_ind6$primary <- as_factor(mut_ind6$primary)
levels(mut_ind6$primary)
mut_ind7 <- mut_ind6 %>% 
  select(-primary)
which(duplicated(mut_ind7$tcgaid))
mut_ind8 <- mut_ind7 %>%
  distinct(tcgaid, .keep_all = TRUE)

### Merge data:
merge_mut <- pgs %>%
  inner_join(mut_ind8, by = "tcgaid")
dim(merge_mut)

### Divide pgs into quintiles:
merge_mut <- merge_mut %>%
  arrange(pgs)
merge_mut <- merge_mut %>%
  mutate(quintile = ntile(pgs, 5))
merge_mut %>%
  group_by(quintile) %>%
  count()

### Calculate counts of wild-type and non-silent mutations per quintile per gene:
proportion <- merge_mut %>%
  group_by(quintile) %>%
  summarise(across(ARID1A:ZFHX3, 
                   list(prop.0 = ~sum(. == 0),
                        prop.1 = ~sum(. == 1)),
                   .names = "{.col}_{.fn}"))

### Format proportion data for plotting:
proportion_plot <- proportion %>%
  pivot_longer(cols = -quintile, names_to = c(".value", "mutation_status"), names_sep = "_") %>%
  pivot_longer(cols = -c(quintile, mutation_status), names_to = "gene", values_to = "value")

### Calculate percentages of mutation counts:
proportion_plot <- proportion_plot %>%
  group_by(quintile, gene) %>%
  mutate(percentage = (value / sum(value)) * 100) %>%
  ungroup()

### Split proportion_plot per gene:
proportion_plot_gene <- split(proportion_plot, proportion_plot$gene)
for (gene_name in names(proportion_plot_gene)) {
  cat(gene_name, "\n")
  print(proportion_plot_gene[[gene_name]])}

### Plot mutations per gene:
plot_list <- list()
for (gene_name in names(proportion_plot_gene)) {
  mutation_data <- proportion_plot_gene[[gene_name]]
  mutation_data$mutation_status <- factor(mutation_data$mutation_status, levels = c("prop.1", "prop.0"))
  plot <- ggbarplot(mutation_data, 
                    x = "quintile", 
                    y = "percentage",
                    fill = "mutation_status",
                    palette = "Paired",
                    label = mutation_data$value,
                    lab.pos = "in",
                    lab.vjust = 0.5,
                    lab.hjust = 1.2,
                    lab.size = 6,
                    legend = "null",
                    xlab = "BMI PGS quintile",
                    ylab = paste(gene_name, "mutation proportion (%)")) +
    rotate() +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 13))
  plot_list[[gene_name]] <- plot}

### Orgnaise plots into a palette and save: 
plots <- grid.arrange(grobs = plot_list, ncol = 4)
plots
ggsave(plots, file = "bmi_pgs_mutation_plots.pdf")
