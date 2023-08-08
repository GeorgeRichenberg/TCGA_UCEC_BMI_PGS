## Load required packages:
library(tidyverse)
library(vroom)
library(CpGassoc)

### Load required data:
clinical1 <- read_delim("TCGA.UCEC.sampleMap_UCEC_clinicalMatrix.txt", delim = "\t")
dnam1 <- vroom("dnam_for_ewas.csv")
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

### Format dnam1:
dim(dnam1)
sum(rowSums(is.na(dnam1))*100 >80)
dnam1 <- dnam1 %>%
  filter(rowSums(is.na(dnam1))*100 <= 80)
dim(dnam1)
dnam1 <- dnam1 %>%
  gather(tcgaid, val, 2:ncol(dnam1)) %>%
  spread(sample, val)
dnam1 <- dnam1 %>% 
  separate(tcgaid, into = c("tcgaid1","tcgaid2","tcgaid3","primary"), sep = "-")
dnam1 <- dnam1 %>%
  unite(tcgaid, tcgaid1, tcgaid2, tcgaid3, sep = "-")
dnam1$primary <- as_factor(dnam1$primary)
levels(dnam1$primary)
dnam1 <- dnam1 %>% 
  filter(str_detect(primary, "11") != TRUE)
which(duplicated(dnam1$tcgaid))
dnam1 <- dnam1 %>% 
  select(-primary)
dnam1 <- dnam1 %>%
  distinct(`tcgaid`, .keep_all = TRUE)
write_csv(dnam1, file = "dnam_for_ewas.csv")

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
dnam_pgs_covars <- dnam1 %>%
  merge(pgs, by="tcgaid") %>%
  merge(clinical6, by="tcgaid") %>%
  merge(ancestry2, by="tcgaid")

### Transposing DNAm data [row = 450K probe, column = sample]:
betas <- dnam_pgs_covars %>%
  select(-c(pgs, age, ms_status, stage, PC1:PC10))
betas <- betas %>%
  gather(sample, val, 2:ncol(betas)) %>%
  spread(tcgaid, val)
betas <- as.matrix(betas %>%
                   column_to_rownames(var = "sample"))

### Running EWAS:
ewas <- cpg.assoc(beta.val = betas,
                  indep = dnam_pgs_covars$pgs,
                  covariates = dnam_pgs_covars[,c("age", 
                                                  "ms_status", 
                                                  "stage", 
                                                  "PC1", 
                                                  "PC2", 
                                                  "PC3", 
                                                  "PC4", 
                                                  "PC5", 
                                                  "PC6", 
                                                  "PC7", 
                                                  "PC8", 
                                                  "PC9", 
                                                  "PC10")])
capture.output(summary(ewas), file = "bmi_pgs_beta_ewas_summary.txt")

### Format annotation:
load('meffilAnno.RData')
annotation2 <- as_tibble(annotation)
annotation3 <- annotation2 %>%
  select(name, chromosome, position)
annotation4 <- annotation3 %>%
  mutate(chromosome = str_replace(chromosome, "chr", "")) %>%
  mutate(chromosome = str_replace(chromosome, "X", "23")) %>%
  mutate(chromosome = str_replace(chromosome, "Y", "24"))
annotation4$chromosome <- as.numeric(annotation4$chromosome)
annotation5 <- annotation4 %>%
  drop_na()
annotation6 <- annotation5[!grepl("ch.", annotation5$name),]
annotation7 <- annotation6[with(annotation6, order(chromosome, position)),]

### Creating a Manhattan plot of results:
manhattan(ewas, 
          chr = annotation7$chromosome, 
          cpgname = annotation7$name, 
          pos = annotation7$position,
          cpg.labels = "FDR",
          col = c("gray10", "gray60"),
          save.plot = "endometrial_bmi_pgs_beta_ewas_manhattan_plot",
          file.type = "pdf")