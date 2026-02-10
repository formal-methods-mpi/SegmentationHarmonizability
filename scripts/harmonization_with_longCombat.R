# load relevant packages
library(longCombat)
# set working directory
setwd("/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/life/Data/")
# load the unadjusted data  
df <- read.csv(file = "df_life.csv")
# load the WMHV data
df_wmhv <- read.csv(file = "df_life_wmhv.csv")
# create a ghost col because longCombat would not work otherwise
df_wmhv$ghost_col <- 1
# load the uncorrected df
df_no_correction <- read.csv(file = "df_life_no_correction.csv")

# apply long combat
df_adj <- longCombat(idvar='id', timevar='timepoint', formula='AGE + GENDER',
                     features = c("TGMV", "TCV", "LVV", "HCV",  "AV"),
                     ranef='(1|id)', batchvar = "SITE", data=df)
df_adj <- df_adj$data_combat
df_wmhv_adj <- longCombat(idvar='id', timevar='timepoint', formula='AGE + GENDER',
                     features = c("WMHV_asinh", "ghost_col"), ranef='(1|id)', batchvar = "SITE", data=df_wmhv)
df_wmhv_adj <- df_wmhv_adj$data_combat
df_no_correction_adj <- longCombat(idvar='id', timevar='timepoint', formula='AGE + GENDER',
                     features = c("TGMV", "TCV", "LVV", "HCV",  "AV"),
                     ranef='(1|id)', batchvar = "SITE", data=df_no_correction)
df_no_correction_adj <- df_no_correction_adj$data_combat
# save dfs
write.csv(df_adj, "df_life_adj.csv")
write.csv(df_wmhv_adj, "df_life_wmhv_adj.csv")
write.csv(df_no_correction_adj, "df_life_no_correction_adj.csv")

################# harmonize cross runs ######################
df <- read.csv(file = "df_life_cross.csv")
# apply long combat
df_adj <- longCombat(idvar='id', timevar='timepoint', formula='AGE + GENDER',
                     features = c("TGMV", "TCV", "LVV", "HCV",  "AV"),
                     ranef='(1|id)', batchvar = "SITE", data=df)
df_adj <- df_adj$data_combat
# save df
write.csv(df_adj, "df_life_cross_adj.csv")


