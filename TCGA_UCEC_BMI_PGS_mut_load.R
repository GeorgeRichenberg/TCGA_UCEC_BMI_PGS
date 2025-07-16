### Set working directory and load required packages:
setwd("")
library(tidyverse)
library(vroom)
library(broom)
library(DescTools)

### Load required data:
clinical1 <- read_delim("TCGA.UCEC.sampleMap_UCEC_clinicalMatrix.txt", delim = "\t")
mutload1 <- vroom("mutation-load_updated.txt", delim = "\t")
ancestry1 <- read_delim("WashU_PCA_ethnicity_assigned.tsv", delim = "\t")
pgs <- read_delim("ucec_bmi_female_primary_std.txt", delim = "\t")

### Format clinical1:
colnames(clinical1)
which(duplicated(clinical1$`_PATIENT`))
clinical2 <- clinical1 %>% 
  distinct(`_PATIENT`, .keep_all = TRUE)
clinical3 <- clinical2 %>%
  select(`_PATIENT`, CDE_ID_3226963, age_at_initial_pathologic_diagnosis, clinical_stage)
which(is.na(clinical3$clinical_stage))
levels(as_factor(clinical3$clinical_stage))
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

### Re-define ms status:
clinical5$ms_status <- as.character(clinical5$ms_status)
clinical6 <- clinical5 %>%
  mutate(ms_status = replace(ms_status, ms_status == "MSS", 0)) %>%
  mutate(ms_status = replace(ms_status, ms_status == "MSI-L", 1)) %>%
  mutate(ms_status = replace(ms_status, ms_status == "MSI-H", 1))
clinical6$ms_status <- as.numeric(clinical6$ms_status)

### Format mutload1:
mutload2 <- mutload1 %>% 
  filter(Cohort == "UCEC")
mutload3 <- mutload2 %>%
  select(-c(Cohort, Patient_ID))
mutload3 <- mutload3 %>%
  rename(tcgaid = Tumor_Sample_ID)
mutload3 <- mutload3 %>%
  arrange(tcgaid)
mutload3 <- mutload3 %>%
  separate(tcgaid, into = c("tcgaid1","tcgaid2","tcgaid3","primary"), sep = "-")
levels(as_factor(mutload3$primary))
mutload3 <- mutload3 %>%
  filter(str_detect(primary, "06A") != TRUE)
mutload3 <- mutload3 %>%
  select(-primary)
mutload3 <- mutload3 %>%
  unite(tcgaid, tcgaid1, tcgaid2, tcgaid3, sep = "-")
which(duplicated(mutload3$tcgaid))
mutload3 <- mutload3 %>%
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
merge_mutload_covars <- pgs %>% 
  inner_join(clinical6, by = "tcgaid") %>% 
  inner_join(ancestry2, by = "tcgaid") %>% 
  inner_join(mutload3, by = "tcgaid")
dim(merge_mutload_covars)

### Run multivariable regression:
merge_mutload_covars_result1 <- merge_mutload_covars %>%
  gather(measure, value, -c(tcgaid, pgs, age, ms_status, stage, PC1:PC10)) %>%
  nest(-measure) %>%
  mutate(fit = map(data, ~ glm(value ~ pgs + age + ms_status + stage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, quasipoisson(link = "log"), data = .x)),
         tidied = map(fit, tidy)) %>%
  unnest(tidied)
merge_mutload_covars_result2 <- merge_mutload_covars_result1 %>%
  filter(term == "pgs")
merge_mutload_covars_result2 <- merge_mutload_covars_result2 %>%
  arrange(p.value)
merge_mutload_covars_result2
write_csv(merge_mutload_covars_result2, file = "merge_mutload_covars_result2.csv")