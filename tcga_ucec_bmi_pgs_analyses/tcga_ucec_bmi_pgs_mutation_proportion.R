### Load required packages:
library(tidyverse)
library(vroom)
library(stringr)
library(ggpubr)
library(cowplot)

### Load required data:
mut_ind1 <- vroom("mc3_gene_level_UCEC_mc3_gene_level.txt", delim = "\t")
pgs <- read_delim("ucec_bmi_female_primary_std.txt", delim = "\t")

### Formatting mut_ind1:
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
pgs_mut <- pgs %>% 
  inner_join(mut_ind8, by = "tcgaid")
pgs2 <- pgs_mut %>% 
  select(tcgaid, pgs)

### Dividing pgs dataset into quintiles according to pgs value:
pgs2 <- pgs2 %>%
  arrange(pgs)
pgs3 <- pgs2 %>%
  mutate(quintile = ntile(pgs, 5))

### Separate mut_ind by gene:
mut_ind_arid1a <- mut_ind8 %>%
  select(tcgaid, ARID1A)
mut_ind_chd4 <- mut_ind8 %>%
  select(tcgaid, CHD4)
mut_ind_ctcf <- mut_ind8 %>%
  select(tcgaid, CTCF)
mut_ind_ctnnb1 <- mut_ind8 %>%
  select(tcgaid, CTNNB1)
mut_ind_fbxw7 <- mut_ind8 %>%
  select(tcgaid, FBXW7)
mut_ind_kmt2d <- mut_ind8 %>%
  select(tcgaid, KMT2D)
mut_ind_kras <- mut_ind8 %>%
  select(tcgaid, KRAS)
mut_ind_pik3ca <- mut_ind8 %>%
  select(tcgaid, PIK3CA)
mut_ind_pik3r1 <- mut_ind8 %>%
  select(tcgaid, PIK3R1)
mut_ind_ppp2r1a <- mut_ind8 %>%
  select(tcgaid, PPP2R1A)
mut_ind_pten <- mut_ind8 %>%
  select(tcgaid, PTEN)
mut_ind_tp53 <- mut_ind8 %>%
  select(tcgaid, TP53)
mut_ind_zfhx3 <- mut_ind8 %>%
  select(tcgaid, ZFHX3)

### Merge with pgs3 for plotting:
combined_arid1a <- pgs3 %>%
  inner_join(mut_ind_arid1a, by = "tcgaid")
combined_chd4 <- pgs3 %>%
  inner_join(mut_ind_chd4, by = "tcgaid")
combined_ctcf <- pgs3 %>%
  inner_join(mut_ind_ctcf, by = "tcgaid")
combined_ctnnb1 <- pgs3 %>%
  inner_join(mut_ind_ctnnb1, by = "tcgaid")
combined_fbxw7 <- pgs3 %>%
  inner_join(mut_ind_fbxw7, by = "tcgaid")
combined_kmt2d <- pgs3 %>%
  inner_join(mut_ind_kmt2d, by = "tcgaid")
combined_kras <- pgs3 %>%
  inner_join(mut_ind_kras, by = "tcgaid")
combined_pik3ca <- pgs3 %>%
  inner_join(mut_ind_pik3ca, by = "tcgaid")
combined_pik3r1 <- pgs3 %>%
  inner_join(mut_ind_pik3r1, by = "tcgaid")
combined_ppp2r1a <- pgs3 %>%
  inner_join(mut_ind_ppp2r1a, by = "tcgaid")
combined_pten <- pgs3 %>%
  inner_join(mut_ind_pten, by = "tcgaid")
combined_tp53 <- pgs3 %>%
  inner_join(mut_ind_tp53, by = "tcgaid")
combined_zfhx3 <- pgs3 %>%
  inner_join(mut_ind_zfhx3, by = "tcgaid")

### Calculating proportions of non-silent mutations per pgs quintile:
arid1a <- combined_arid1a %>%
  group_by(quintile, as.character(ARID1A)) %>%
  summarise(n = n()) %>%
  mutate(percent = n / sum(n)*100)
arid1a[,"percent"] <- round(arid1a[,"percent"], digits = 1)
arid1a <- arid1a %>%
  rename("mutation_status" = 'as.character(ARID1A)')
arid1a$mutation_status <- str_replace(arid1a$mutation_status, "0", "wild_type")
arid1a$mutation_status <- str_replace(arid1a$mutation_status, "1", "non_silent_mutation")

chd4 <- combined_chd4 %>%
  group_by(quintile, as.character(CHD4)) %>%
  summarise(n = n()) %>%
  mutate(percent = n / sum(n)*100)
chd4[,"percent"] <-round(chd4[,"percent"], digits = 1)
chd4 <- chd4 %>%
  rename("mutation_status" = 'as.character(CHD4)')
chd4$mutation_status <- str_replace(chd4$mutation_status, "0", "wild_type")
chd4$mutation_status <- str_replace(chd4$mutation_status, "1", "non_silent_mutation")

ctcf <- combined_ctcf %>%
  group_by(quintile, as.character(CTCF)) %>%
  summarise(n = n()) %>%
  mutate(percent = n / sum(n)*100)
ctcf[,"percent"] <-round(ctcf[,"percent"], digits = 1)
ctcf <- ctcf %>%
  rename("mutation_status" = 'as.character(CTCF)')
ctcf$mutation_status <- str_replace(ctcf$mutation_status, "0", "wild_type")
ctcf$mutation_status <- str_replace(ctcf$mutation_status, "1", "non_silent_mutation")

ctnnb1 <- combined_ctnnb1 %>%
  group_by(quintile, as.character(CTNNB1)) %>%
  summarise(n = n()) %>%
  mutate(percent = n / sum(n)*100)
ctnnb1[,"percent"] <-round(ctnnb1[,"percent"], digits = 1)
ctnnb1 <- ctnnb1 %>%
  rename("mutation_status" = 'as.character(CTNNB1)')
ctnnb1$mutation_status <- str_replace(ctnnb1$mutation_status, "0", "wild_type")
ctnnb1$mutation_status <- str_replace(ctnnb1$mutation_status, "1", "non_silent_mutation")

fbxw7 <- combined_fbxw7 %>%
  group_by(quintile, as.character(FBXW7)) %>%
  summarise(n = n()) %>%
  mutate(percent = n / sum(n)*100)
fbxw7[,"percent"] <-round(fbxw7[,"percent"], digits = 1)
fbxw7 <- fbxw7 %>%
  rename("mutation_status" = 'as.character(FBXW7)')
fbxw7$mutation_status <- str_replace(fbxw7$mutation_status, "0", "wild_type")
fbxw7$mutation_status <- str_replace(fbxw7$mutation_status, "1", "non_silent_mutation")

kmt2d <- combined_kmt2d %>%
  group_by(quintile, as.character(KMT2D)) %>%
  summarise(n = n()) %>%
  mutate(percent = n / sum(n)*100)
kmt2d[,"percent"] <-round(kmt2d[,"percent"], digits = 1)
kmt2d <- kmt2d %>%
  rename("mutation_status" = 'as.character(KMT2D)')
kmt2d$mutation_status <- str_replace(kmt2d$mutation_status, "0", "wild_type")
kmt2d$mutation_status <- str_replace(kmt2d$mutation_status, "1", "non_silent_mutation")

kras <- combined_kras %>%
  group_by(quintile, as.character(KRAS)) %>%
  summarise(n = n()) %>%
  mutate(percent = n / sum(n)*100)
kras[,"percent"] <-round(kras[,"percent"], digits = 1)
kras <- kras %>%
  rename("mutation_status" = 'as.character(KRAS)')
kras$mutation_status <- str_replace(kras$mutation_status, "0", "wild_type")
kras$mutation_status <- str_replace(kras$mutation_status, "1", "non_silent_mutation")

pik3ca <- combined_pik3ca %>%
  group_by(quintile, as.character(PIK3CA)) %>%
  summarise(n = n()) %>%
  mutate(percent = n / sum(n)*100)
pik3ca[,"percent"] <- round(pik3ca[,"percent"], digits = 1)
pik3ca <- pik3ca %>%
  rename("mutation_status" = 'as.character(PIK3CA)')
pik3ca$mutation_status <- str_replace(pik3ca$mutation_status, "0", "wild_type")
pik3ca$mutation_status <- str_replace(pik3ca$mutation_status, "1", "non_silent_mutation")

pik3r1 <- combined_pik3r1 %>%
  group_by(quintile, as.character(PIK3R1)) %>%
  summarise(n = n()) %>%
  mutate(percent = n / sum(n)*100)
pik3r1[,"percent"] <- round(pik3r1[,"percent"], digits = 1)
pik3r1 <- pik3r1 %>%
  rename("mutation_status" = 'as.character(PIK3R1)')
pik3r1$mutation_status <- str_replace(pik3r1$mutation_status, "0", "wild_type")
pik3r1$mutation_status <- str_replace(pik3r1$mutation_status, "1", "non_silent_mutation")

ppp2r1a <- combined_ppp2r1a %>%
  group_by(quintile, as.character(PPP2R1A)) %>%
  summarise(n = n()) %>%
  mutate(percent = n / sum(n)*100)
ppp2r1a[,"percent"] <- round(ppp2r1a[,"percent"], digits = 1)
ppp2r1a <- ppp2r1a %>%
  rename("mutation_status" = 'as.character(PPP2R1A)')
ppp2r1a$mutation_status <- str_replace(ppp2r1a$mutation_status, "0", "wild_type")
ppp2r1a$mutation_status <- str_replace(ppp2r1a$mutation_status, "1", "non_silent_mutation")

pten <- combined_pten %>%
  group_by(quintile, as.character(PTEN)) %>%
  summarise(n = n()) %>%
  mutate(percent = n / sum(n)*100)
pten[,"percent"] <- round(pten[,"percent"], digits = 1)
pten <- pten %>%
  rename("mutation_status" = 'as.character(PTEN)')
pten$mutation_status <- str_replace(pten$mutation_status, "0", "wild_type")
pten$mutation_status <- str_replace(pten$mutation_status, "1", "non_silent_mutation")

tp53 <- combined_tp53 %>%
  group_by(quintile, as.character(TP53)) %>%
  summarise(n = n()) %>%
  mutate(percent = n / sum(n)*100)
tp53[,"percent"] <- round(tp53[,"percent"], digits = 1)
tp53 <- tp53 %>%
  rename("mutation_status" = 'as.character(TP53)')
tp53$mutation_status <- str_replace(tp53$mutation_status, "0", "wild_type")
tp53$mutation_status <- str_replace(tp53$mutation_status, "1", "non_silent_mutation")

zfhx3 <- combined_zfhx3 %>%
  group_by(quintile, as.character(ZFHX3)) %>%
  summarise(n = n()) %>%
  mutate(percent = n / sum(n)*100)
zfhx3[,"percent"] <- round(zfhx3[,"percent"], digits = 1)
zfhx3 <- zfhx3 %>%
  rename("mutation_status" = 'as.character(ZFHX3)')
zfhx3$mutation_status <- str_replace(zfhx3$mutation_status, "0", "wild_type")
zfhx3$mutation_status <- str_replace(zfhx3$mutation_status, "1", "non_silent_mutation")

### Plotting proportions of non-silent mutations:
arid1a_plot <- ggbarplot(arid1a, 
                         "quintile", 
                         "percent",
                         fill = "mutation_status", 
                         palette = "Paired",
                         label = arid1a$n,
                         lab.hjust = 1.5,
                         lab.vjust = 0,
                         lab.size = 10,
                         legend = "null",
                         xlab = "BMI PGS quintile",
                         ylab = "ARID1A mutation proportion (%)") +
  rotate() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))
arid1a_plot
chd4_plot <- ggbarplot(chd4, 
                       "quintile", 
                       "percent",
                       fill = "mutation_status", 
                       palette = "Paired",
                       label = chd4$n,
                       lab.hjust = 1.5,
                       lab.vjust = 0,
                       lab.size = 10,
                       legend = "null",
                       xlab = "BMI PGS quintile",
                       ylab = "CHD4 mutation proportion (%)") +
  rotate() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))
chd4_plot
ctcf_plot <- ggbarplot(ctcf, 
                       "quintile", 
                       "percent",
                       fill = "mutation_status", 
                       palette = "Paired",
                       label = ctcf$n,
                       lab.hjust = 1.5,
                       lab.vjust = 0,
                       lab.size = 10,
                       legend = "null",
                       xlab = "BMI PGS quintile",
                       ylab = "CTCF mutation proportion (%)") +
  rotate() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))
ctcf_plot
ctnnb1_plot <- ggbarplot(ctnnb1, 
                         "quintile", 
                         "percent",
                         fill = "mutation_status", 
                         palette = "Paired",
                         label = ctnnb1$n,
                         lab.hjust = 1.5,
                         lab.vjust = 0,
                         lab.size = 10,
                         legend = "null",
                         xlab = "BMI PGS quintile",
                         ylab = "CTNNB1 mutation proportion (%)") +
  rotate() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))
ctnnb1_plot
fbxw7_plot <- ggbarplot(fbxw7, 
                        "quintile", 
                        "percent",
                        fill = "mutation_status", 
                        palette = "Paired",
                        label = fbxw7$n,
                        lab.hjust = 1.5,
                        lab.vjust = 0,
                        lab.size = 10,
                        legend = "null",
                        xlab = "BMI PGS quintile",
                        ylab = "FBXW7 mutation proportion (%)") +
  rotate() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))
fbxw7_plot
kmt2d_plot <- ggbarplot(kmt2d, 
                        "quintile", 
                        "percent",
                        fill = "mutation_status", 
                        palette = "Paired",
                        label = kmt2d$n,
                        lab.hjust = 1.5,
                        lab.vjust = 0,
                        lab.size = 10,
                        legend = "null",
                        xlab = "BMI PGS quintile",
                        ylab = "KMT2D mutation proportion (%)") +
  rotate() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))
kmt2d_plot
kras_plot <- ggbarplot(kras, 
                       "quintile", 
                       "percent",
                       fill = "mutation_status", 
                       palette = "Paired",
                       label = kras$n,
                       lab.hjust = 1.5,
                       lab.vjust = 0,
                       lab.size = 10,
                       legend = "null",
                       xlab = "BMI PGS quintile",
                       ylab = "KRAS mutation proportion (%)") +
  rotate() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))
kras_plot
pik3ca_plot <- ggbarplot(pik3ca, 
                         "quintile", 
                         "percent",
                         fill = "mutation_status", 
                         palette = "Paired",
                         label = pik3ca$n,
                         lab.hjust = 1.5,
                         lab.vjust = 0,
                         lab.size = 10,
                         legend = "null",
                         xlab = "BMI PGS quintile",
                         ylab = "PIK3CA mutation proportion (%)") +
  rotate() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))
pik3ca_plot
pik3r1_plot <- ggbarplot(pik3r1, 
                         "quintile", 
                         "percent",
                         fill = "mutation_status", 
                         palette = "Paired",
                         label = pik3r1$n,
                         lab.hjust = 1.5,
                         lab.vjust = 0,
                         lab.size = 10,
                         legend = "null",
                         xlab = "BMI PGS quintile",
                         ylab = "PIK3R1 mutation proportion (%)") +
  rotate() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))
pik3r1_plot
ppp2r1a_plot <- ggbarplot(ppp2r1a, 
                          "quintile", 
                          "percent",
                          fill = "mutation_status", 
                          palette = "Paired",
                          label = ppp2r1a$n,
                          lab.hjust = 1.5,
                          lab.vjust = 0,
                          lab.size = 10,
                          legend = "null",
                          xlab = "BMI PGS quintile",
                          ylab = "PPP2R1A mutation proportion (%)") +
  rotate() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))
ppp2r1a_plot
pten_plot <- ggbarplot(pten, 
                       "quintile", 
                       "percent",
                       fill = "mutation_status", 
                       palette = "Paired",
                       label = pten$n,
                       lab.hjust = 1.5,
                       lab.vjust = 0,
                       lab.size = 10,
                       legend = "null",
                       xlab = "BMI PGS quintile",
                       ylab = "PTEN mutation proportion (%)") +
  rotate() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))
pten_plot
tp53_plot <- ggbarplot(tp53, 
                       "quintile", 
                       "percent",
                       fill = "mutation_status", 
                       palette = "Paired",
                       label = tp53$n,
                       lab.hjust = 1.5,
                       lab.vjust = 0,
                       lab.size = 10,
                       legend = "null",
                       xlab = "BMI PGS quintile",
                       ylab = "TP53 mutation proportion (%)") +
  rotate() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))
tp53_plot
zfhx3_plot <- ggbarplot(zfhx3, 
                        "quintile", 
                        "percent",
                        fill = "mutation_status", 
                        palette = "Paired",
                        label = zfhx3$n,
                        lab.hjust = 1.5,
                        lab.vjust = 0,
                        lab.size = 10,
                        legend = "null",
                        xlab = "BMI PGS quintile",
                        ylab = "ZFHX3 mutation proportion (%)") +
  rotate() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))
zfhx3_plot

### Arranging plots with cowplot for all TCGA UCEC patients:
ggdraw() +
  draw_plot(arid1a_plot, x=0, y=0, width = 1, height = 0.5)
ggsave(file = "a_non_silent_mutation_proportion.pdf")

ggdraw() +
  draw_plot(chd4_plot, x=0, y=0, width = 1, height = 0.5)
ggsave(file = "b_non_silent_mutation_proportion.pdf")

ggdraw() +
  draw_plot(ctcf_plot, x=0, y=0, width = 1, height = 0.5)
ggsave(file = "c_non_silent_mutation_proportion.pdf")

ggdraw() +
  draw_plot(ctnnb1_plot, x=0, y=0, width = 1, height = 0.5)
ggsave(file = "d_non_silent_mutation_proportion.pdf")

ggdraw() +
  draw_plot(fbxw7_plot, x=0, y=0, width = 1, height = 0.5)
ggsave(file = "e_non_silent_mutation_proportion.pdf")
  
ggdraw() +
  draw_plot(kmt2d_plot, x=0, y=0, width = 1, height = 0.5)
ggsave(file = "f_non_silent_mutation_proportion.pdf")

ggdraw() +
  draw_plot(kras_plot, x=0, y=0, width = 1, height = 0.5)
ggsave(file = "g_non_silent_mutation_proportion.pdf")

ggdraw() +
  draw_plot(pik3ca_plot, x=0, y=0, width = 1, height = 0.5)
ggsave(file = "h_non_silent_mutation_proportion.pdf")

ggdraw() +
  draw_plot(pik3r1_plot, x=0, y=0, width = 1, height = 0.5)
ggsave(file = "i_non_silent_mutation_proportion.pdf")

ggdraw() +
  draw_plot(ppp2r1a_plot, x=0, y=0, width = 1, height = 0.5)
ggsave(file = "j_non_silent_mutation_proportion.pdf")

ggdraw() +
  draw_plot(pten_plot, x=0, y=0, width = 1, height = 0.5)
ggsave(file = "k_non_silent_mutation_proportion.pdf")

ggdraw() +
  draw_plot(tp53_plot, x=0, y=0, width = 1, height = 0.5)
ggsave(file = "l_non_silent_mutation_proportion.pdf")

ggdraw() +
  draw_plot(zfhx3_plot, x=0, y=0, width = 1, height = 0.5)
ggsave(file = "m_non_silent_mutation_proportion.pdf")