### this script serves to assemble the different segmentations of the LIFE and the scanner-update datasets

# load required packages
library(tidyverse)
library(lavaan)
library(bestNormalize)
## we start with the data of the scanner update dataset
## we will ave to load the data and harmonize column names to get everything into one dataset
## first we load the samseg data
samseg_update <- read.csv("/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/scanner_update/samseg/results_summary.csv")
# create id, scanner and segmentation col
samseg_update$id <- sub("_.*$", "", samseg_update$participant)
samseg_update$scanner <-  ifelse(grepl("SKYRA", samseg_update$participant), "SKYRA", "VERIO")
samseg_update$segmentation <- "samseg"
# calculate and harmonize the names of the relevant brain regions
# TIV = total intracranial volume
# TGMV = total grey matter volume
# TWMV = total white matter volume
# TCV = total cortical volume
# LVV = lateral ventricle volume
# HCV = hippocampal volume
# AV = amygdala volume
# prefix r = right
# prefix l = left

# calculate TGMV from individual volumes for samseg
gm_structures <- c(
  "Right.Cerebral.Cortex", "Left.Cerebral.Cortex",
  "Left.Cerebellum.Cortex", "Right.Cerebellum.Cortex", 
  "Left.Thalamus", "Right.Thalamus",
  "Left.Caudate", "Right.Caudate",
  "Left.Putamen", "Right.Putamen",
  "Left.Pallidum", "Right.Pallidum",
  "Left.Hippocampus", "Right.Hippocampus",
  "Left.Amygdala", "Right.Amygdala",
  "Left.Accumbens.area", "Right.Accumbens.area",
  "Left.VentralDC", "Right.VentralDC"
)
samseg_update$TGMV <- rowSums(samseg_update[,gm_structures])
# calculate TWMV 
samseg_update$TWMV <- rowSums(samseg_update[,c("Right.Cerebral.White.Matter", "Left.Cerebral.White.Matter", 
                                               "Right.Cerebellum.White.Matter", "Left.Cerebellum.White.Matter")])
# calculate TCV
samseg_update$TCV <- rowSums(samseg_update[,c("Right.Cerebral.Cortex", "Left.Cerebral.Cortex")])
# calculate LVV
samseg_update$LVV <- rowSums(samseg_update[,c("Left.Lateral.Ventricle", "Right.Lateral.Ventricle")])
# calculate HCV
samseg_update$HCV <- rowSums(samseg_update[,c("Left.Hippocampus", "Right.Hippocampus")])
# calculate AV
samseg_update$AV <- rowSums(samseg_update[,c("Left.Amygdala", "Right.Amygdala")])
# rename
samseg_update <- samseg_update %>%
  rename(
    "TIV" = "sbTIV"
  )


## now we load the v 7.4.0 recon-all data
recon7_update <- read.delim("/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/scanner_update/recon/aseg_stats.txt")
# create id, scanner and segmentation col
recon7_update$id <- sub("_.*$", "", recon7_update$Measure.volume)
recon7_update$scanner <-  ifelse(grepl("SKYRA", recon7_update$Measure.volume), "SKYRA", "VERIO")
recon7_update$segmentation <- "recon7"
# calculate and harmonize the names of the relevant brain regions
recon7_update$TWMV <- rowSums(recon7_update[,c("lhCerebralWhiteMatterVol", "rhCerebralWhiteMatterVol", 
                                               "Left.Cerebellum.White.Matter", "Right.Cerebellum.White.Matter")])
recon7_update$LVV <- rowSums(recon7_update[,c("Right.Lateral.Ventricle", "Left.Lateral.Ventricle")])
recon7_update$HCV <- rowSums(recon7_update[,c("Right.Hippocampus", "Left.Hippocampus")])
recon7_update$AV <- rowSums(recon7_update[,c("Right.Amygdala", "Left.Amygdala")])
                                                                                      
#rename
recon7_update <- recon7_update %>%
  rename(
    "TIV" = "EstimatedTotalIntraCranialVol",
    "TGMV" = "TotalGrayVol",
    "TCV" = "CortexVol"
    )

## now we can load the cat 12 data
# load the table according to the neuromorphometrics atlas 
cat_update_csf <- read.csv("/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/scanner_update/cat/ROI_neuromorphometrics_Vcsf.csv")
cat_update_gm <- read.csv("/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/scanner_update/cat/ROI_neuromorphometrics_Vgm.csv")
cat_update_wm <- read.csv("/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/scanner_update/cat/ROI_neuromorphometrics_Vwm.csv")
# add them up
cat_update <- cat_update_gm[,2:ncol(cat_update_gm)] + cat_update_wm[,2:ncol(cat_update_wm)] + cat_update_csf[,2:ncol(cat_update_csf)] 
# create id and scanner cols
cat_update$id <- sub(".*/([^/_]+)_.*", "\\1", cat_update_gm$names)
cat_update$scanner <- ifelse(grepl("SKYRA", cat_update_gm$names), "SKYRA", "VERIO")
# load dfs with global measures
cat_update_global_measures <- read.delim("/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/scanner_update/cat/TIV.txt", col.names = c("path", "TIV", "TGMV", "TWMV", "CSF", "WMHV"), header = F)
cat_update_global_measures$id <- sub(".*/([^/_]+)_.*", "\\1", cat_update_global_measures$path)
cat_update_global_measures$scanner <- ifelse(grepl("SKYRA", cat_update_global_measures$path), "SKYRA", "VERIO")
# merge the two dfs
cat_update <- merge(cat_update, cat_update_global_measures, by = c("id", "scanner"))
# create segmentation col
cat_update$segmentation <- "cat"

# calculate and harmonize the names of the relevant brain regions
# caluclate TCV from individual cortical ROIS
# Step 1: Identify all left/right columns
lr_cols <- grep("^(Left|Right)\\.", names(cat_update), value = TRUE)

# Step 2: Exclude non-cortical columns and subregions (segments, supplementary.motor.cortex)
exclude_patterns <- c("Cerebellum", "White.Matter", "CSF", "Ventricle",
                      "Brain.Stem", "Basal.Forebrain", "Accumbens", "Thalamus", 
                      "Putamen", "Pallidum", "Amygdala", "Hippocampus", "vessel",
                      "Caudate", "Inf.Lat.Vent", "Ventral.DC", "segment", "supplementary.motor.cortex")
cortical_cols <- lr_cols[!grepl(paste(exclude_patterns, collapse = "|"), lr_cols)]
# Step 3: Sum all cortical ROIs row-wise to get total cortical volume
cat_update$TCV <- rowSums(cat_update[, cortical_cols])
# calculate further ROI volumes
cat_update$LVV <- rowSums(cat_update[,c("Left.Lateral.Ventricle", "Right.Lateral.Ventricle")])
cat_update$HCV <- rowSums(cat_update[,c("Left.Hippocampus", "Right.Hippocampus")])
cat_update$AV <- rowSums(cat_update[,c("Left.Amygdala", "Right.Amygdala")])
# all volumes have to be multiplied by 1000 (turn mL into mm^3)
cat_update[,!names(cat_update) %in% c("id", "scanner", "segmentation", "path")] <- 
  cat_update[,!names(cat_update) %in% c("id", "scanner", "segmentation", "path")] * 1000
## lastly we get the segmentation of the original v 6.0.0p recon-all segmentation
recon6_update <- read.delim("/data/pt_nro186_lifeupdate/Data/FREESURFER/FREESURFER_long/aseg_stats.txt")
# create id, scanner and segmentation col
# naming follos a slightl different convention in this table (no "sub-" prefix, . missing in other names)
recon6_update$id <- gsub("\\.", "", paste0("sub-" ,sub("_.*$", "", recon6_update$Measure.volume)))
recon6_update$scanner <- ifelse(grepl("SKYRA", recon6_update$Measure.volume), "SKYRA", "VERIO")
recon6_update$segmentation <- "recon6"

# calculate and harmonize the names of the relevant brain regions
recon6_update$TWMV <- rowSums(recon6_update[,c("lhCerebralWhiteMatterVol", "rhCerebralWhiteMatterVol", 
                                               "Left.Cerebellum.White.Matter", "Right.Cerebellum.White.Matter")])
recon6_update$LVV <- rowSums(recon6_update[, c("Right.Lateral.Ventricle", "Left.Lateral.Ventricle")])
recon6_update$HCV <- rowSums(recon6_update[, c("Right.Hippocampus", "Left.Hippocampus")])
recon6_update$AV  <- rowSums(recon6_update[, c("Right.Amygdala", "Left.Amygdala")])

# rename columns
recon6_update <- recon6_update %>%
  rename(
    "TIV" = "EstimatedTotalIntraCranialVol",
    "TGMV" = "TotalGrayVol",
    "TCV" = "CortexVol" 
  )

# select relevant cols
relevant_cols <- c("id", "scanner", "segmentation", "TIV", "TGMV", "TWMV", "TCV", "LVV", "HCV", "AV")

# bind the dfs together
df_update <- rbind(samseg_update[,relevant_cols], recon6_update[,relevant_cols], recon7_update[,relevant_cols], cat_update[,relevant_cols])
# create a variable that combines the values of scanner and segmentation
df_update$SITE <- as.factor(paste(df_update$scanner, df_update$segmentation, sep = "_"))

# now load the QA data
QA <- read.csv("/data/pt_nro186_lifeupdate/Data/FREESURFER/QA/QA_list_short.csv")
# remove empty rows
QA <- QA[1:121,]
# create id col, add 0 that might got lost at the beginning/end, remove ., add "sub-"
# first add 0 lost at the end
QA$ID <- ifelse(nchar(sub(".*\\.", "", QA$ID)) == 1, paste0(QA$ID, "0"), QA$ID)
QA$id <- ifelse(nchar(QA$ID) != 7, gsub("\\.", "", paste0("sub-" ,QA$ID)), gsub("\\.", "", paste0("sub-", "0", QA$ID)))

# trim trailing spaces
QA$id <- trimws(QA$id)
QA$correction.necessary <- ifelse(is.na(QA$correction.necessary), 0, QA$correction.necessary)
# rename German colnames
QA <- QA %>%
  rename(
    "age" = "Alter_aktuell",
    "gender" = "Geschlecht"
  )
df_update <- df_update %>%
  left_join(QA[,c("id", "correction.necessary", "age", "gender")], by = "id")
# scale age
df_update$age <- (df_update$age - min(df_update$age))/sd(df_update$age)
# square age
df_update$age_squared <- df_update$age^2 
# code gender as numeric
df_update$gender <- ifelse(df_update$gender == "w", 0, 1)

# transform to desired units
df_update[,"TIV"] <- df_update[,"TIV"]/1000000 # mm^3 to dm^3
df_update[,c("TGMV", "TWMV", "TCV")] <- df_update[,c("TGMV", "TWMV", "TCV")]/100000 # mm^3 to 100 cm^3
df_update[,c("HCV", "AV")] <- df_update[,c("HCV", "AV")]/1000 # mm^3 to cm^3

# LVV is not approx. normally distributed
df_update$LVV <- asinh(df_update$LVV)

# produce a df without participants whose scans have been manually edited
df_update_no_correction <- df_update[df_update$correction.necessary == 0,]

# save dfs for harmonization
write.csv(df_update, file = "/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/scanner_update/Data/df_update.csv")
write.csv(df_update_no_correction, file = "/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/scanner_update/Data/df_update_no_correction.csv")

# now we prepare the scanner update WML data

# load segmented data
wmhv_update <- read.csv("/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/scanner_update/wmhv/results_summary.csv")
# create id, scanner and segmentation col
wmhv_update$id <- wmhv_update$participant
wmhv_update$SITE <-  ifelse(grepl("skyra", wmhv_update$scanner), "SKYRA", "VERIO")
bn <- bestNormalize(c(wmhv_update[wmhv_update$SITE== "SKYRA", "Lesions"], wmhv_update[wmhv_update$SITE== "VERIO", "Lesions"]), , allow_orderNorm = F, standardize = F)
# uses Non-Standardized Log_b(x + a) with a = 0 and b = 10
wmhv_update[wmhv_update$SITE== "SKYRA", "WMHV_norm"] <- bn$x.t[1:length(wmhv_update[wmhv_update$SITE== "SKYRA", "Lesions"])]
wmhv_update[wmhv_update$SITE == "VERIO", "WMHV_norm"] <- bn$x.t[(1+length(wmhv_update[wmhv_update$SITE== "SKYRA", "Lesions"])):length(bn$x.t)]

# merge with df_update for age and gender
wmhv_update <- merge(wmhv_update, df_update[!duplicated(df_update$id),c("id", "age", "gender")], by = "id", all.x = T)

# save df for harmonization
write.csv(wmhv_update, file = "/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/scanner_update/Data/wmhv_update.csv")

######################## now we can load the segmented LIFE data ############################
# we start with the samseg data (both long and cross pipeline)
samseg_life <- read.csv("/data/pt_life_freesurfer/samseg/results_summary.csv")
samseg_cross_life <- read.csv("/data/pt_life_freesurfer/samseg/results_cross_summary.csv")
# create segmentation col
samseg_life$segmentation <- "samseg"
samseg_cross_life$segmentation <- "samseg_cross"
# calculate volumes
samseg_life$TGMV <- rowSums(samseg_life[, gm_structures])

samseg_life$TWMV <- rowSums(samseg_life[,c("Right.Cerebral.White.Matter", "Left.Cerebral.White.Matter", 
                                           "Right.Cerebellum.White.Matter", "Left.Cerebellum.White.Matter")])

# calculate TCV
samseg_life$TCV <- rowSums(samseg_life[, c(
  "Right.Cerebral.Cortex", "Left.Cerebral.Cortex"
)])

# calculate LVV
samseg_life$LVV <- rowSums(samseg_life[, c(
  "Left.Lateral.Ventricle", "Right.Lateral.Ventricle"
)])

# calculate HCV
samseg_life$HCV <- rowSums(samseg_life[, c(
  "Left.Hippocampus", "Right.Hippocampus"
)])

# calculate AV
samseg_life$AV <- rowSums(samseg_life[, c(
  "Left.Amygdala", "Right.Amygdala"
)])
samseg_cross_life$TGMV <- rowSums(samseg_cross_life[, gm_structures])

samseg_cross_life$TWMV <- rowSums(samseg_cross_life[,c("Right.Cerebral.White.Matter", "Left.Cerebral.White.Matter", 
                                                       "Right.Cerebellum.White.Matter", "Left.Cerebellum.White.Matter")])
samseg_cross_life$TCV <- rowSums(samseg_cross_life[, c(
  "Right.Cerebral.Cortex", "Left.Cerebral.Cortex"
)])
samseg_cross_life$LVV <- rowSums(samseg_cross_life[, c(
  "Left.Lateral.Ventricle", "Right.Lateral.Ventricle"
)])
samseg_cross_life$HCV <- rowSums(samseg_cross_life[, c(
  "Left.Hippocampus", "Right.Hippocampus"
)])
samseg_cross_life$HCV <- rowSums(samseg_cross_life[, c(
  "Left.Hippocampus", "Right.Hippocampus"
)])
samseg_cross_life$AV <- rowSums(samseg_cross_life[, c(
  "Left.Amygdala", "Right.Amygdala"
)])
# rename
samseg_life <- samseg_life %>%
  rename(
    "TIV" = "sbTIV",
    "id" = "participant",
    "WMHV" = "Lesions"
  )
samseg_cross_life <- samseg_cross_life %>%
  rename(
    "TIV" = "sbTIV",
    "id" = "participant",
    "WMHV" = "Lesions"
  )
# now we continue with the LST WMHV data
lst_life <- read.csv("/data/pt_life_whm/Results/Tables/Datatable_for_LIFE/qa_info_LIFE_WMH_LST.csv")
lst_life <- lst_life[,c("pseudonym", "qa_LST_fu", "vol_long_bl", "vol_long_fu")] %>%
  pivot_longer(
    cols = c(vol_long_bl, vol_long_fu),
    names_to = "timepoint",
    values_to = "WMHV",
    values_drop_na = TRUE
  )
lst_life$timepoint <- ifelse(grepl("fu", lst_life$timepoint), "FU", "BL")
lst_life$segmentation <- "lst"
lst_life <- lst_life %>%
  rename(
    "id" = "pseudonym")
lst_life$wmhv_usable <- ifelse(lst_life$qa_LST_fu != 0, 0, 1)

# now we load the FS v 5.3.0 segmentations of the LIFE data (long and cross piepline)
recon5_life <- read.delim("/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/life/Data/aseg530_stats.txt")
recon5_cross_life <- read.delim("/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/life/Data/aseg530_cross_stats.txt")
# create id, timepoint and segmentation col
# drop everything after the first .
recon5_life$id <- sub("\\..*$", "", recon5_life$Measure.volume)
recon5_cross_life$id <- recon5_cross_life$Measure.volume
# create timepoint variable
recon5_life$timepoint <- ifelse(grepl("fu" ,recon5_life$id), "FU", "BL")
recon5_cross_life$timepoint <- ifelse(grepl("fu" ,recon5_cross_life$id), "FU", "BL")
# drop the fu in the name if present
recon5_life$id <- ifelse(grepl("fu" , recon5_life$id), substr(recon5_life$id, 1, nchar(recon5_life$id) - 3), recon5_life$id)
recon5_cross_life$id <- ifelse(grepl("fu" , recon5_cross_life$id), substr(recon5_cross_life$id, 1, nchar(recon5_cross_life$id) - 3), recon5_cross_life$id)
# calculate and harmonize the names of the relevant brain regions
recon5_life$LVV <- rowSums(recon5_life[, c("Right.Lateral.Ventricle", "Left.Lateral.Ventricle")])
recon5_life$HCV <- rowSums(recon5_life[, c("Right.Hippocampus", "Left.Hippocampus")])
recon5_life$AV <- rowSums(recon5_life[, c("Right.Amygdala", "Left.Amygdala")])

recon5_cross_life$LVV <- rowSums(recon5_cross_life[, c("Right.Lateral.Ventricle", "Left.Lateral.Ventricle")])
recon5_cross_life$HCV <- rowSums(recon5_cross_life[, c("Right.Hippocampus", "Left.Hippocampus")])
recon5_cross_life$AV <- rowSums(recon5_cross_life[, c("Right.Amygdala", "Left.Amygdala")])

# rename
recon5_life <- recon5_life %>%
  rename(
    "TIV" = "EstimatedTotalIntraCranialVol",
    "TGMV" = "TotalGrayVol",
    "TCV" = "CortexVol"
  )

recon5_cross_life <- recon5_cross_life %>%
  rename(
    "TIV" = "EstimatedTotalIntraCranialVol",
    "TGMV" = "TotalGrayVol",
    "TCV" = "CortexVol"
  )

recon5_life$segmentation <- "recon5"
recon5_cross_life$segmentation <- "recon5_cross"
# now load the LIFE QA
QA_LIFE <- read_csv(file = "/data/gh_gr_agingandobesity_share/life_shared/Data/Preprocessed/derivatives/radiological_assessment/radio_assessment_long.csv")
QA_LIFE$radiology_usable <- ifelse(QA_LIFE$med_Befund.Bewertung_Bef != "verwendungsf\xe4hig nein" & !is.na(QA_LIFE$med_Befund.Bewertung_Bef), 1, 0)
# remove files without mrt_pseudonym value or duplicates
QA_LIFE <- QA_LIFE[!is.na(QA_LIFE$mrt_pseudonym),]
QA_LIFE <- QA_LIFE[!duplicated(QA_LIFE[,c("fu", "mrt_pseudonym")]),]
QA_LIFE$timepoint <- ifelse(QA_LIFE$fu == "fu", "FU", "BL")  
colnames(QA_LIFE)[22] <- "id"
# now load the QC data for the GM segmentations
QC <- read.csv("/data/pt_life_freesurfer/Tabular_Data_QA/FreeSurfer/QA/final_qc_freesurfer.csv")
colnames(QC)[1] <- "id"
QC <- QC %>%
  pivot_longer(
    cols = -id,
    names_to = c("timepoint", ".value"),
    names_sep = "_"
  )
# read age and gender for LIFE-Adult dataset
life_age_gender <- read.csv("/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/scanner_update/Data/life_age_gender.csv")
# merge with samseg
samseg_life <- merge(samseg_life, life_age_gender, by = c("id", "timepoint"))
# scale age
samseg_life$AGE <- (samseg_life$AGE - min(samseg_life$AGE, na.rm = T))/sd(samseg_life$AGE, na.rm = T)  
# square age
samseg_life$age_squared <- samseg_life$AGE^2
# now we merge the Life dataframes
relevant_cols <- c("id", "segmentation", "timepoint", "TGMV", "TCV", "LVV", "HCV", "AV")
df_life <- bind_rows(samseg_life[,c(relevant_cols, "AGE", "GENDER", "age_squared")], recon5_life[,relevant_cols], samseg_cross_life[,relevant_cols], recon5_cross_life[,relevant_cols])
df_life_wmhv <- bind_rows(samseg_life[,c("id", "segmentation", "timepoint", "WMHV", "AGE", "GENDER", "age_squared")], lst_life[,c("id", "segmentation", "timepoint", "WMHV", "wmhv_usable")])
df_life <- merge(df_life, QC[,c("id", "corrected", "usable", "timepoint")], by = c("id", "timepoint"), all.x = T)
df_life <- merge(df_life, QA_LIFE[,c("radiology_usable", "id", "timepoint")], by = c("id", "timepoint"),  all.x = T)
# copy lst quality judgement to samseg, copy samseg age and gender to lst
# only keep those that have both samseg and lst
df_life_wmhv <- df_life_wmhv %>%
  group_by(id, timepoint) %>%
  mutate(
    wmhv_usable = wmhv_usable[segmentation == "lst"][1],
    AGE = AGE[segmentation == "samseg"][1],
    age_squared = age_squared[segmentation == "samseg"][1],
    GENDER = GENDER[segmentation == "samseg"][1]
  ) %>%
  filter(n() == 2) %>%
  ungroup()
# copy samseg age and gender to recon
# only keep those that have both samseg and recon
df_life <- df_life %>%
  group_by(id, timepoint) %>%
  mutate(
    AGE = AGE[segmentation == "samseg"][1],
    GENDER = GENDER[segmentation == "samseg"][1],
    age_squared = age_squared[segmentation == "samseg"][1]
  ) %>%
  filter(all(c("samseg", "recon5") %in% segmentation)) %>%
  ungroup()
# apply exclusion criteria
df_life <- df_life[df_life$usable == 1 & df_life$radiology_usable == 1,]
df_life_wmhv <- df_life_wmhv[df_life_wmhv$wmhv_usable == 1,]

# rename segmentation as SITE
df_life <- df_life %>%
  rename(
    "SITE" = "segmentation"
  )
df_life_wmhv <- df_life_wmhv %>%
  rename(
    "SITE" = "segmentation"
  )

# drop rows with NA in age, gender, site
df_life <- df_life[!is.na(df_life$AGE) & !is.na(df_life$GENDER) & !is.na(df_life$SITE),]
df_life_wmhv <- df_life_wmhv[!is.na(df_life_wmhv$AGE) & !is.na(df_life_wmhv$GENDER) & !is.na(df_life_wmhv$SITE),]

# transform to desired units
df_life[,c("TGMV", "TCV")] <- df_life[,c("TGMV", "TCV")]/100000 # mm^3 to 100 cm^3
df_life[,"LVV"] <- df_life[,"LVV"]/10000 # mm^3 to 10 cm^3
df_life[,c("HCV", "AV")] <- df_life[,c("HCV", "AV")]/1000 # mm^3 to cm^3

# LVV is not approx normally distributed
df_life$LVV <- asinh(df_life$LVV)

# WMHV is not approx normally distributed
bn <- bestNormalize(df_life_wmhv$WMHV, allow_orderNorm = F, standardize = F)
# uses Non-Standardized Box Cox Transformation with lambda = -0.4430144  
df_life_wmhv$WMHV_norm <- bn$x.t

# remove empty row
df_life <- df_life[complete.cases(df_life),]

# remove manually edited values
df_life_no_correction <- df_life[df_life$corrected == 0,]

# save files (save cross and long runs separately)
write.csv(df_life[!grepl("cross", df_life$SITE),], file = "/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/life/Data/df_life.csv", row.names = F)
write.csv(df_life_no_correction[!grepl("cross", df_life_no_correction$SITE),], file = "/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/life/Data/df_life_no_correction.csv", row.names = F)
write.csv(df_life[grepl("cross", df_life$SITE),], file = "/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/life/Data/df_life_cross.csv", row.names = F)
write.csv(df_life_wmhv, file = "/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/life/Data/df_life_wmhv.csv", row.names = F)
# save BL and FU values separately for neuroharmonize
write.csv(df_life[df_life$timepoint == "BL" & !grepl("cross", df_life$SITE),], file = "/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/life/Data/df_life_bl.csv", row.names = F)
write.csv(df_life_no_correction[df_life_no_correction$timepoint == "BL" & !grepl("cross", df_life_no_correction$SITE),], file = "/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/life/Data/df_life_no_correction_bl.csv", row.names = F)
write.csv(df_life[df_life$timepoint == "BL" & grepl("cross", df_life$SITE),], file = "/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/life/Data/df_life_cross_bl.csv", row.names = F)
write.csv(df_life_wmhv[df_life_wmhv$timepoint == "BL",], file = "/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/life/Data/df_life_wmhv_bl.csv", row.names = F)

write.csv(df_life[df_life$timepoint == "FU" & !grepl("cross", df_life$SITE),], file = "/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/life/Data/df_life_fu.csv", row.names = F)
write.csv(df_life_no_correction[df_life_no_correction$timepoint == "FU" & !grepl("cross", df_life_no_correction$SITE),], file = "/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/life/Data/df_life_no_correction_fu.csv", row.names = F)
write.csv(df_life[df_life$timepoint == "FU" & grepl("cross", df_life$SITE),], file = "/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/life/Data/df_life_cross_fu.csv", row.names = F)
write.csv(df_life_wmhv[df_life_wmhv$timepoint == "FU",], file = "/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/life/Data/df_life_wmhv_fu.csv", row.names = F)

### calculate change scores for an exploratory analysis
# first turn the dataframes into a wide format
# exclude cross runs while doing so
df_life_wide <- pivot_wider(df_life[!grepl("cross", df_life$SITE),],
            names_from = c(SITE, timepoint),
            values_from = c(TGMV, TCV, LVV, HCV, AV),
            id_cols = id)
df_life_wmhv_wide <- pivot_wider(df_life_wmhv,
                            names_from = c(SITE, timepoint),
                            values_from = c(WMHV_norm),
                            id_cols = id) 

# within each segmentation, calculate the difference from BL to FU
## create list of paired segmentations
pairs <- list(
  TGMV_samseg = c("TGMV_samseg_BL", "TGMV_samseg_FU"),
  TGMV_recon5 = c("TGMV_recon5_BL", "TGMV_recon5_FU"),
  TCV_samseg  = c("TCV_samseg_BL",  "TCV_samseg_FU"),
  TCV_recon5  = c("TCV_recon5_BL",  "TCV_recon5_FU"),
  LVV_samseg = c("LVV_samseg_BL", "LVV_samseg_FU"),
  LVV_recon5 = c("LVV_recon5_BL", "LVV_recon5_FU"),
  HCV_samseg  = c("HCV_samseg_BL",  "HCV_samseg_FU"),
  HCV_recon5  = c("HCV_recon5_BL",  "HCV_recon5_FU"),
  AV_samseg  = c("AV_samseg_BL",  "AV_samseg_FU"),
  AV_recon5  = c("AV_recon5_BL",  "AV_recon5_FU")
)
for (nm in names(pairs)) {
  bl_col <- pairs[[nm]][1]
  fu_col <- pairs[[nm]][2]
  df_life_wide[[paste0("d", nm)]] <- df_life_wide[[fu_col]] - df_life_wide[[bl_col]]
}
pairs_wmhv <- list(
  WMHV_samseg = c("samseg_BL", "samseg_FU"),
  WMHV_lst = c("lst_BL", "lst_FU"))
for (nm in names(pairs_wmhv)) {
  bl_col <- pairs_wmhv[[nm]][1]
  fu_col <- pairs_wmhv[[nm]][2]
  df_life_wmhv_wide[[paste0("d", nm)]] <- df_life_wmhv_wide[[fu_col]] - df_life_wmhv_wide[[bl_col]]
}
# turn them to a somewhat longer format
df_life_wide <- df_life_wide %>%
  pivot_longer(cols = grep("^d", colnames(df_life_wide), value = T),
             names_to = c("variable", "SITE"), 
             names_sep = "_",                    
             values_to = "value") %>%
  dplyr::select(id, variable, SITE, value) %>%
  pivot_wider(
    names_from = "variable",  # each variable becomes a column
    values_from = "value"
  ) %>% drop_na()
df_life_wmhv_wide <- df_life_wmhv_wide %>%
  pivot_longer(cols = grep("^d", colnames(df_life_wmhv_wide), value = T),
               names_to = c("variable", "SITE"), 
               names_sep = "_",                    
               values_to = "value") %>%
  dplyr::select(id, variable, SITE, value) %>%
  pivot_wider(
    names_from = "variable",  # each variable becomes a column
    values_from = "value"
  ) %>% drop_na()

# merge age and gender back in 
df_life_wide <- df_life_wide %>%
  left_join(df_life[df_life$timepoint == "BL",c("id", "AGE", "GENDER")], by = "id") %>% distinct()

df_life_wmhv_wide <- df_life_wmhv_wide %>%
  left_join(df_life_wmhv[df_life_wmhv$timepoint == "BL",c("id", "AGE", "GENDER")], by = "id") %>% distinct()

# save the dfs for harmonization
write.csv(df_life_wide, file = "/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/life/Data/df_life_changes.csv", row.names = F)
write.csv(df_life_wmhv_wide, file = "/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/life/Data/df_life_wmhv_changes.csv", row.names = F)