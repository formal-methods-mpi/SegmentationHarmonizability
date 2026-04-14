library(dplyr)
setwd("/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/scanner_update/Data/")
# load the unadjusted data from the scanner update study 
df_update <- read.csv("df_update.csv")
df_update <- df_update[,2:16]
renaming_table <- data.frame(matrix(nrow = 121, ncol = 2))
colnames(renaming_table) <- c("id", "new_names")
renaming_table$id <- unique(df_update$id)
renaming_table$new_names <- paste0("ID_", sample(1:121))
# rename
df_update_osf <- df_update %>%
  left_join(renaming_table, by = "id") %>%
  mutate(id = new_names) %>%
  select(-new_names)
write_csv(df_update_osf, "df_update_osf.csv")
# load the adjusted data from the scanner update study 
df_update_adj <- read.csv("df_update_adj.csv", col.names = colnames(df_update)[4:10])
df_update_adj[,c("id", "scanner", "segmentation", "SITE")] <- df_update[,c("id", "scanner", "segmentation", "SITE")]

# rename
df_update_adj_osf <- df_update_adj %>%
  left_join(renaming_table, by = "id") %>%
  mutate(id = new_names) %>%
  select(-new_names)
write_csv(df_update_adj_osf, "df_update_adj_osf.csv")

# repeat for the dataset excluding participants with manual editing
# load the unadjusted data from the scanner update study 
df_update_no_correction <- read.csv("df_update_no_correction.csv")
df_update_no_correction <- df_update_no_correction[,2:16]
  # break each scanner-segmentation combination down to a scanner and a segmentation variable
df_update_no_correction$scanner <- ifelse(grepl("VER", df_update_no_correction$SITE), "VERIO", "SKYRA")
df_update_no_correction <- df_update_no_correction %>%
  mutate(segmentation = case_when(
    grepl("samseg", SITE) ~ "samseg",
    grepl("cat", SITE) ~ "cat",
    grepl("recon6", SITE) ~ "recon6",
    grepl("recon7", SITE) ~ "recon7"
  ))
# rename
df_update_no_correction_osf <- df_update_no_correction %>%
  left_join(renaming_table, by = "id") %>%
  mutate(id = new_names) %>%
  select(-new_names)
write_csv(df_update_no_correction_osf, "df_update_no_correction_osf.csv")

# load the adjusted data from the scanner update study 
df_update_no_correction_adj <- read.csv("df_update_no_correction_adj.csv", col.names = colnames(df_update_no_correction)[4:10])
# fill in the missing data from the unadjusted df (these cols must not be included during harmonization) 
df_update_no_correction_adj[,c("id", "scanner", "segmentation", "SITE")] <- df_update_no_correction[,c("id", "scanner", "segmentation", "SITE")]

# rename
df_update_no_correction_adj_osf <- df_update_no_correction_adj %>%
  left_join(renaming_table, by = "id") %>%
  mutate(id = new_names) %>%
  select(-new_names)
write_csv(df_update_no_correction_adj_osf, "df_update_no_correction_adj_osf.csv")

df_update_wmhv <- read.csv("/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/scanner_update/Data/wmhv_update.csv")
df_update_wmhv <- df_update_wmhv[,4:57]
colnames(df_update_wmhv)[1] <- "id"
df_update_wmhv_adj <- read.csv("/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/scanner_update/Data/wmhv_update_adj.csv")
df_update_wmhv$WMHV_adj <- df_update_wmhv_adj$X..WMHV_norm
# rename
df_update_wmhv_osf <- df_update_wmhv %>%
  left_join(renaming_table, by = "id") %>%
  mutate(id = new_names) %>%
  select(-new_names)
write_csv(df_update_wmhv_osf, "df_update_wmhv_osf.csv")

