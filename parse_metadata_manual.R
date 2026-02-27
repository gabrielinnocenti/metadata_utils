library('dplyr')
# Load your GitHub package
library(vecmatch)

# monoclonal ADMCI
metadata <- read.csv('ADMCI_NED_metadata.csv')
metadata <- metadata[c('sample_name', 'group_test', 'sex', 'age')]
metadata$age <- as.numeric(metadata$age)
metadata$group_test <- as.character(metadata$group_test)
metadata <- metadata %>%
  mutate(
    group_test = case_when(
      grepl("^dementia", group_test) ~ "Dementia",
      grepl("^control", group_test)  ~ "Controls",
      TRUE ~ group_test
    ),
    sex = case_when(
      grepl("^Male", sex)   ~ "M",
      grepl("^Female", sex) ~ "F",
      TRUE ~ sex
    )
  )
metadata['origin'] <- 'Netherlands'
metadata['project'] <- 'ADMCI_NED'
metadata_ADMCI_NED <- metadata

# CRC radiotherapy
metadata <- read.csv('CRC_radiotherapy_metadata.csv', sep=';')
metadata <- metadata[c('sample_name', 'group_test', 'sex', 'age')]

# Filter baseline (group_test == "1") and recode it
metadata <- metadata %>%
  filter(group_test == "1") %>%
  mutate(
    group_test = "CRC"   # simply replace all remaining "1" with "CRC"
  )
metadata
metadata$age <- as.numeric(metadata$age)
metadata$group_test <- as.character(metadata$group_test)
metadata['origin'] <- 'Austria'
metadata['project'] <- 'CRC_radiotherapy'
metadata_crc_radiotherapy <- metadata

# PCa Innsbruck
metadata <- read.csv('PCa_Innsbruck_metadata.csv', sep=';')
colnames(metadata) <- c("sample_id", "sample_name", "group_test", "sex", "age",
                        "Histological_finding_pT_Stadium",
                        "Biochemisches.Rezidiv..1.ja")
metadata <- metadata[c('sample_name', 'group_test', 'sex', 'age')]
metadata$age <- as.numeric(metadata$age)
metadata$group_test <- as.character(metadata$group_test)
metadata['origin'] <- 'Austria'
metadata['project'] <- 'PCa_Innsbruck'
metadata <- metadata %>%
  mutate(
    group_test = case_when(
      grepl("^GR1", group_test) ~ "Low risk prostate cancer",
      grepl("^GR2", group_test) ~ "Medium risk prostate cancer",
      grepl("^GR3", group_test) ~ "High risk prostate cancer",
      grepl("^GR4_matched", group_test) ~ "Controls",
      grepl("^GR4$", group_test) ~ "Controls",
      grepl("Male", sex) ~ "M",
      grepl("Female", sex) ~ "F",
      TRUE ~ group_test
    )
  )
metadata

metadata_pca_innsbruck <- metadata


# CORSA
metadata <- read.csv('CORSA_metadata.csv', sep=';')
metadata <- metadata[c('sample_name', 'group_test', 'sex', 'age')]
metadata$age <- as.numeric(metadata$age)
metadata$group_test <- as.character(metadata$group_test)
metadata['origin'] <- 'Austria'
metadata['project'] <- 'CORSA'
metadata_corsa <- metadata 


# Kiel controls
metadata <- read.csv('Kiel_IBD_metadata.csv', sep=';')
metadata <- metadata[c('sample_name', 'group_test', 'sex', 'age')]
metadata$age <- as.numeric(metadata$age)
metadata$group_test <- as.character(metadata$group_test)
metadata['origin'] <- 'Germany'
metadata['project'] <- 'Kiel_IBD'
metadata_kiel_ibd_controls <- metadata

metadata <- bind_rows(metadata_ADMCI_NED, 
          metadata_kiel_ibd_controls, 
          metadata_corsa, 
          metadata_pca_innsbruck, 
          metadata_crc_radiotherapy)


controls_metadata <- metadata %>%
  filter(group_test == "Controls")


write.table(metadata,
            "metadata_slim.csv", 
            sep = ";",         # semicolon separator
            row.names = FALSE, 
            col.names = TRUE,  # include header
            quote = FALSE,      # quotes around text
            fileEncoding = "UTF-8")

metadata_vecmatch <- metadata %>%
  filter(group_test == "Controls") 
metadata_vecmatch["batch"] <- "db"
metadata_vecmatch <- metadata_vecmatch[c("sample_name","group_test", "sex", "age", "origin", "batch")]

# Simulate data
set.seed(123)  # per riproducibilità

n <- 120  # numero di campioni

# Creazione del dataset
simulated_data <- data.frame(
  sample_name = paste0("Sample_", 1:n),
  group_test = "Controls",
  sex = sample(c("M", "F"), n, replace = TRUE),
  age = round(rnorm(n, mean = 55, sd = 10)), # età media 55, deviazione standard 10
  origin = 'Austria',
  batch = 'simulated' 
)

combined <- rbind(simulated_data, metadata_vecmatch)
combined <- combined %>%
  filter(!is.na(age))  # rimuove righe con NA in age

raincloud(
  data = combined,
  y=age,
  group = batch,
  significance = "t_test",
  sig_label_color = TRUE,
  sig_label_size = 3,
  limits = c(7, 48)
)


formula <- formula(batch ~ age * sex)
gps_matrix <- estimate_gps(formula,
                           data = combined,
                           method = "vglm",
                           reference = "Control"
)
