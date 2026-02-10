# this script runs the analyses for the segmentation harmonization study
# there are 2 datasets (Scanner update & LIFE-Adult) that have to be checked before and after harmonization 
# for each ROI, we want to calculate an ICC, a PD between all batches and an ANOVA / t-test
# for the scanner update dataset, we will additionally calculate SEM-based ICEDs
# we start with the scanner update dataset (section 1)
# 1.1 calculation of ICCs
# 1.2 calculation of PDs
# 1.3 calculation of ANOVAs followed-up by pairwise t-tests of individual batches
# 1.4 calculation of ICED to get ICED-ICCs for specific scanner-segmentation combinations before/after harmonization
# then we continue with the LIFE-Adult dataset (section 2)
# 2.1 calculation of ICCs
# 2.2 calculation of PDs
# 2.3 calculation of ANOVAs

# load required packages
library(tidyverse)
library(lme4)
library(rstatix)
library(glue)
library(lavaan)

setwd("/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/scanner_update/Data/")

############ Section 1: The scanner update dataset #######################
# load the unadjusted data from the scanner update study 
df_update <- read.csv("df_update.csv")
# break each scanner-segmentation combination down to a scanner and a segmentation variable
df_update$scanner <- ifelse(grepl("VER", df_update$SITE), "VERIO", "SKYRA")
df_update <- df_update %>%
  mutate(segmentation = case_when(
    grepl("samseg", SITE) ~ "samseg",
    grepl("cat", SITE) ~ "cat",
    grepl("recon6", SITE) ~ "recon6",
    grepl("recon7", SITE) ~ "recon7"
  ))

# load the adjusted data from the scanner update study 
df_update_adj <- read.csv("summary_df_update_adj.csv", col.names = colnames(df_update)[5:11])
# fill in the missing data from the unadjusted df (these cols must not be included during harmonization) 
df_update_adj[,c("id", "scanner", "segmentation", "SITE")] <- df_update[,c("id", "scanner", "segmentation", "SITE")]

# repeat for the dataset excluding participants with manual editing
# load the unadjusted data from the scanner update study 
df_update_no_correction <- read.csv("df_update_no_correction.csv")
# break each scanner-segmentation combination down to a scanner and a segmentation variable
df_update_no_correction$scanner <- ifelse(grepl("VER", df_update_no_correction$SITE), "VERIO", "SKYRA")
df_update_no_correction <- df_update_no_correction %>%
  mutate(segmentation = case_when(
    grepl("samseg", SITE) ~ "samseg",
    grepl("cat", SITE) ~ "cat",
    grepl("recon6", SITE) ~ "recon6",
    grepl("recon7", SITE) ~ "recon7"
  ))

# load the adjusted data from the scanner update study 
df_update_no_correction_adj <- read.csv("summary_df_update_no_correction_adj.csv", col.names = colnames(df_update_no_correction)[5:11])
# fill in the missing data from the unadjusted df (these cols must not be included during harmonization) 
df_update_no_correction_adj[,c("id", "scanner", "segmentation", "SITE")] <- df_update[,c("id", "scanner", "segmentation", "SITE")]

# store the dfs in a list
dfs <- list(df_update = df_update, df_update_adj = df_update_adj, df_update_no_correction = df_update_no_correction,
            df_update_adj_no_correction = df_update_adj_no_correction)

# some analyses will need the dfs in the wide format
df_update_wide <- df_update %>%
  pivot_wider(
    names_from = c(SITE),
    values_from = c(TIV, TGMV, TWMV, TCV, LVV, HCV, AV),
    id_cols = id  
  )
df_update_adj_wide <- df_update_adj %>%
  pivot_wider(
    names_from = c(scanner, segmentation),
    values_from = c(TIV, TGMV, TWMV, TCV, LVV, HCV, AV),
    id_cols = id  
  )

df_update_no_correction_wide <- df_update_no_correction %>%
  pivot_wider(
    names_from = c(SITE),
    values_from = c(TIV, TGMV, TWMV, TCV, LVV, HCV, AV),
    id_cols = id  
  )
df_update_no_correction_adj_wide <- df_update_no_correction_adj %>%
  pivot_wider(
    names_from = c(scanner, segmentation),
    values_from = c(TIV, TGMV, TWMV, TCV, LVV, HCV, AV),
    id_cols = id  
  )

# store wide dfs in a list
wide_dfs <- list(df_update_wide = df_update_wide, df_update_adj_wide = df_update_adj_wide, 
                 df_update_no_correction_wide = df_update_no_correction_wide, df_update_no_correction_adj_wide = df_update_no_correction_adj_wide)

# produce dfs without NAs for pairise t-tests
df_update_pairwise <- df_update %>%
  group_by(id) %>%
  filter(any(scanner == "SKYRA") & any(scanner == "VERIO")) %>%
  ungroup() %>%
  filter(scanner %in% c("SKYRA","VERIO"))
df_update_adj_pairwise <- df_update_adj %>%
  group_by(id) %>%
  filter(any(scanner == "SKYRA") & any(scanner == "VERIO")) %>%
  ungroup() %>%
  filter(scanner %in% c("SKYRA","VERIO"))
pairwise_dfs <- list(df_update_pairwise = df_update_pairwise, df_update_adj_pairwise = df_update_adj_pairwise)

df_update_no_correction_pairwise <- df_update_no_correction %>%
  group_by(id) %>%
  filter(any(scanner == "SKYRA") & any(scanner == "VERIO")) %>%
  ungroup() %>%
  filter(scanner %in% c("SKYRA","VERIO"))
df_update_no_correction_adj_pairwise <- df_update_no_correction_adj %>%
  group_by(id) %>%
  filter(any(scanner == "SKYRA") & any(scanner == "VERIO")) %>%
  ungroup() %>%
  filter(scanner %in% c("SKYRA","VERIO"))
pairwise_dfs <- list(df_update_pairwise = df_update_pairwise, df_update_adj_pairwise = df_update_adj_pairwise, df_update_no_correction_pairwise =
                       df_update_no_correction_pairwise, df_update_no_correction_adj_pairwise = df_update_no_correction_adj_pairwise)

############ Section 1.1: Calculation of ICCs ####################### 
# prepare a results df to be filled
lme_icc_res <- data.frame(matrix(ncol = 62, nrow = 7))
colnames(lme_icc_res) <- c("ROI", paste0(rep(paste0(rep(paste0(rep(c("id_var", "scanner_var", "segmentation_var", "resdiual_var", "ICC"), each = 3), 
                                       c("", "_lower_95_CI", "_upper_95_CI")), each = 2), c("", "_adj")), times = 2), 
                                       rep(c("", "_no_correction"), each = 30)))
lme_icc_res$ROI <- c("TIV", "TGMV", "TWMV", "TCV", "LVV",  "HCV", "AV")

# formula for bootstrapped CIs for ICC
icc_ci <- function(fit) {
  variances <- as.data.frame(VarCorr(fit))
  id_var = variances[variances$grp == "id", "vcov"]
  scan_var = variances[variances$grp == "scanner", "vcov"]
  seg_var = variances[variances$grp == "segmentation", "vcov"] 
  residual_var = variances[variances$grp == "Residual", "vcov"]
  icc = id_var / (id_var + scan_var + seg_var + residual_var)
  c(id_var, scan_var, seg_var, residual_var, icc)
}

# calculate LMEs
# iterate over outcomes
for (outcome in c("TIV", "TGMV", "TWMV", "TCV", "LVV",  "HCV", "AV")) {
  # iterate over dfs
  for (df_name in names(dfs)){
    df <- dfs[[df_name]] 
    # remove SKYRA-recon6 rows if outcome is TIV (uninformative compared to VERIO-recon6)
    if (outcome == "TIV"){
      df <- df[!(df$scanner == "SKYRA" & df$segmentation == "recon6"),]
    }
    # create outcome specific LME formula
    formula <- as.formula(paste0(outcome, " ~ (1|id) + (1|scanner) + (1|segmentation)"))
    # estimate model and bootstrap confidence intervals
    print(paste0("now doing ", outcome))
    res <- lmer(formula = formula, data = df)
    variances <- as.data.frame(VarCorr(res))
    lme_icc <- (variances[variances$grp == "id", "vcov"])/
      (variances[variances$grp == "scanner", "vcov"] + variances[variances$grp == "id", "vcov"] + 
         variances[variances$grp == "segmentation", "vcov"] + variances[variances$grp == "Residual", "vcov"])
    boot_icc <- bootMer(res, FUN  = icc_ci, nsim = 1000, type = "parametric")
    boot_res <- apply(boot_icc$t, 2, quantile, probs = c(0.025, 0.975))
    # store the results in the appropriate columns (unadjusted/adjusted data)
    if (df_name == "df_update") {
      lme_icc_res[lme_icc_res$ROI == outcome, c("id_var", "scanner_var", "segmentation_var", "resdiual_var")] <- 
        c(variances[variances$grp == "id", "vcov"], variances[variances$grp == "scanner", "vcov"], 
          variances[variances$grp == "segmentation", "vcov"], variances[variances$grp == "Residual", "vcov"])
      lme_icc_res[lme_icc_res$ROI == outcome, "ICC"] <- lme_icc
      lme_icc_res[lme_icc_res$ROI == outcome, paste0(rep(c("id_var", "scanner_var", "segmentation_var", "resdiual_var", "ICC"),
                                     each = 2), c("_lower_95_CI", "_upper_95_CI"))] <- as.vector(boot_res)
    } else if (df_name == "df_update_adj") {
      lme_icc_res[lme_icc_res$ROI == outcome, paste0(c("id_var", "scanner_var", "segmentation_var", "resdiual_var"), "_adj")] <- 
        c(variances[variances$grp == "id", "vcov"], variances[variances$grp == "scanner", "vcov"], 
          variances[variances$grp == "segmentation", "vcov"], variances[variances$grp == "Residual", "vcov"])
      lme_icc_res[lme_icc_res$ROI == outcome, "ICC_adj"] <- lme_icc
      lme_icc_res[lme_icc_res$ROI == outcome, paste0(rep(c("id_var", "scanner_var", "segmentation_var", "resdiual_var", "ICC"),
                                     each = 2), c("_lower_95_CI_adj", "_upper_95_CI_adj"))] <- as.vector(boot_res)
    } else if (df_name == "df_update") {
      lme_icc_res[lme_icc_res$ROI == outcome, paste0(c("id_var", "scanner_var", "segmentation_var", "resdiual_var"), "_no_correction")] <- 
        c(variances[variances$grp == "id", "vcov"], variances[variances$grp == "scanner", "vcov"], 
          variances[variances$grp == "segmentation", "vcov"], variances[variances$grp == "Residual", "vcov"])
      lme_icc_res[lme_icc_res$ROI == outcome, "ICC_no_correction"] <- lme_icc
      lme_icc_res[lme_icc_res$ROI == outcome, paste0(rep(c("id_var", "scanner_var", "segmentation_var", "resdiual_var", "ICC"),
                                     each = 2), c("_lower_95_CI_no_correction", "_upper_95_CI_no_correction"))] <- as.vector(boot_res)
      
    } else {
      lme_icc_res[lme_icc_res$ROI == outcome, paste0(c("id_var", "scanner_var", "segmentation_var", "resdiual_var"), "_adj_no_correction")] <- 
        c(variances[variances$grp == "id", "vcov"], variances[variances$grp == "scanner", "vcov"], 
          variances[variances$grp == "segmentation", "vcov"], variances[variances$grp == "Residual", "vcov"])
      lme_icc_res[lme_icc_res$ROI == outcome, "ICC_adj_no_correction"] <- lme_icc
      lme_icc_res[lme_icc_res$ROI == outcome, paste0(rep(c("id_var", "scanner_var", "segmentation_var", "resdiual_var", "ICC"),
                                     each = 2), c("_lower_95_CI_adj_no_correction", "_upper_95_CI_adj_no_correction"))] <- as.vector(boot_res)
    }
  }
}

############ Section 1.2: Calculation of PDs #######################
# iterate over ROIs
for (outcome in c("TIV", "TGMV", "TWMV", "TCV", "LVV",  "HCV", "AV")) {
  # iterate over datasets
  for (df_name in names(wide_dfs)){
    df <- wide_dfs[[df_name]]
    # drop the SKYRA-recon6 TIV col
    df$TIV_SKYRA_recon6 <- NULL
    # get names of relevant cols
    roi_cols <- grep(outcome, colnames(df), value = TRUE)
    # initialize empty matrix
    mean_pd_mat <- matrix(
      NA_real_,
      nrow = length(roi_cols),
      ncol = length(roi_cols),
      dimnames = list(roi_cols, roi_cols)
    )
    
    # fill upper triangle
    for (i in seq_along(roi_cols)) {
      for (j in seq_along(roi_cols)) {
        if (i < j) {
          # extract all values in the 2 relevant cols
          col1 <- df[[roi_cols[i]]]
          col2 <- df[[roi_cols[j]]]
          # calculate a vector of PDs
          # a negative value indicates that the values of the batch given by the rowname is smaller 
          pd <- (2 * (col1 - col2)) / (col1 + col2)
          # calculate the mean PD between the two batches (NAs exist because 5 participants were not scanned on both scanners)
          mean_pd_mat[i, j] <- mean(pd, na.rm = TRUE)
        }
      }
    }
    # mirror the inverse to lower triangle
    mean_pd_mat[lower.tri(mean_pd_mat)] <- t(mean_pd_mat)[lower.tri(mean_pd_mat)] * (-1)
    # turn into percentage values
    mean_pd_mat <- mean_pd_mat*100
    #calculate mean PD for roi from individual comparisons
    mean_pd <- mean(mean_pd_mat, na.rm = T)
    # save variables to workspace
    if (df_name == "df_update_wide"){
      assign(paste0(outcome, "_mean_pd_mat"), mean_pd_mat)
      assign(paste0(outcome, "_mean_pd"), mean_pd)
    } else if (df_name == "df_update_adj_wide"){
      assign(paste0(outcome, "_mean_pd_mat_adj"), mean_pd_mat)
      assign(paste0(outcome, "_mean_pd_adj"), mean_pd)
    } else if (df_name == "df_update_no_correction_wide"){
      assign(paste0(outcome, "_mean_pd_mat_no_correction"), mean_pd_mat)
      assign(paste0(outcome, "_mean_pd_no_correction"), mean_pd)
    } else {
      assign(paste0(outcome, "_mean_pd_mat_no_correction_adj"), mean_pd_mat)
      assign(paste0(outcome, "_mean_pd_no_correction_adj"), mean_pd)
    }
  }
}

############ Section 1.3: Calculation of ANOVAS and pairwise t-tests #######################
# prepare a results df to be filled
anova_res <- data.frame(matrix(ncol = 25, nrow = 7))
colnames(anova_res) <- c("ROI", paste0(rep(c("DFn",  "DFd",  "F",  "p",  "p<.05",  "ges"), times = 4), 
                                       rep(c("", "_adj", "_no_correction", "_no_correction_adj"), each = 6)))
anova_res$ROI <- c("TIV", "TGMV", "TWMV", "TCV", "LVV",  "HCV", "AV")
# iterate over ROIs
for (outcome in c("TIV", "TGMV", "TWMV", "TCV", "LVV",  "HCV", "AV")) {
  # iterate over datasets
  for (df_name in names(dfs)){
    df <- dfs[[df_name]]
    pairwise_df <- pairwise_dfs[[paste0(df_name, "_pairwise")]]
    # remove SKYRA-recon6 rows if outcome is TIV (uninformative compared to VERIO-recon6)
    if (outcome == "TIV"){
      df <- df[!(df$scanner == "SKYRA" & df$segmentation == "recon6"),]
      pairwise_df <- pairwise_df[!(pairwise_df$scanner == "SKYRA" & pairwise_df$segmentation == "recon6"),]
    }
    res.aov <- anova_test(data = df, dv = outcome, wid = id, within = SITE)
    table <- get_anova_table(res.aov)
    if (df_name  == "df_update"){
      anova_res[anova_res$ROI == outcome, 2:7] <- as.data.frame(table)[1,2:7]
    } else if (df_name  == "df_update_adj"){
      anova_res[anova_res$ROI == outcome, 8:13] <- as.data.frame(table)[1,2:7]
    } else if (df_name  == "df_update_no_correction"){
      anova_res[anova_res$ROI == outcome, 14:19] <- as.data.frame(table)[1,2:7]
    } else {
      anova_res[anova_res$ROI == outcome, 20:25] <- as.data.frame(table)[1,2:7]
    }
    # perform pairwise comparisons of the individual scanner-segmentation comparisons
    formula <- as.formula(paste0(outcome, " ~ SITE"))
    pwc <- pairwise_df %>%
      pairwise_t_test(
        formula = formula, paired = TRUE, id = "id",
        p.adjust.method = "hommel"
      )
    # save pairwise comparison table
    if (df_name == "df_update"){
      assign(paste0(outcome, "_anova_result"), res.aov)
      assign(paste0(outcome, "_anova_result_table"), table)
      assign(paste0(outcome, "_pwc"), pwc)
    } else if (df_name  == "df_update_adj"){
      assign(paste0(outcome, "_anova_result_adj"), res.aov)
      assign(paste0(outcome, "_anova_result_table_adj"), table)
      assign(paste0(outcome, "_pwc_adj"), pwc)
    } else if (df_name  == "df_update_no_correction"){
      assign(paste0(outcome, "_anova_result_no_correction"), res.aov)
      assign(paste0(outcome, "_anova_result_table_no_correction"), table)
      assign(paste0(outcome, "_pwc_no_correction"), pwc)
    } else {
      assign(paste0(outcome, "_anova_result_no_correction_adj"), res.aov)
      assign(paste0(outcome, "_anova_result_table_no_correction_adj"), table)
      assign(paste0(outcome, "_pwc_no_correction_adj"), pwc)
    }
  }
}

############ Section 1.4: ICED #######################
# prepare a results df to be filled
iced_res <- data.frame(matrix(ncol = 43, nrow = 7))
parameters <- c("obsvar", "truevar", "SKYRAvar",  "VERIOvar", "samsegvar", "recon6var", "recon7var", "catvar",
                "VERIOmean",  "recon6mean", "recon7mean", "catmean", "truemean")
col_names <- c(parameters, paste0(rep(c("ICC_SKYRA_", "ICC_VERIO_"), each = 4), rep(c("samseg", "recon6", "recon7", "cat"), times = 2))) 
                
colnames(iced_res) <- c("ROI", paste0(rep(col_names, times = 2), rep(c("", "_adj"), each = 21)))
iced_res$ROI <- c("TIV", "TGMV", "TWMV", "TCV", "LVV",  "HCV", "AV")
# iterate over ROIs
for (outcome in c("TIV", "TGMV", "TWMV", "TCV", "LVV",  "HCV", "AV")) {
  # create colname vectors to adapt the model to the respective ROI
  obsvars <- grep(outcome, colnames(df_update_wide), value = T)
  skyravars <- grep("SKYRA", obsvars, value = T)
  veriovars <- grep("VERIO", obsvars, value = T)
  samsegvars <- grep("samseg", obsvars, value = T)
  recon6vars <- grep("recon6", obsvars, value = T)
  recon7vars <- grep("recon7", obsvars, value = T)
  catvars <- grep("cat", obsvars, value = T)
  # write the lavvan model
  model <- glue("
   ! regressions 
   truth=~1.0*{obsvars[1]}
   truth=~1.0*{obsvars[2]}
   truth=~1.0*{obsvars[3]}
   truth=~1.0*{obsvars[4]}
   truth=~1.0*{obsvars[5]}
   truth=~1.0*{obsvars[6]}
   truth=~1.0*{obsvars[7]}
   truth=~1.0*{obsvars[8]}
   SKYRA=~1.0*{skyravars[1]}
   SKYRA=~1.0*{skyravars[2]}
   SKYRA=~1.0*{skyravars[3]}
   SKYRA=~1.0*{skyravars[4]}
   VERIO=~1.0*{veriovars[1]}
   VERIO=~1.0*{veriovars[2]}
   VERIO=~1.0*{veriovars[3]}
   VERIO=~1.0*{veriovars[4]}
   samseg=~1.0*{samsegvars[1]}
   samseg=~1.0*{samsegvars[2]}
   recon6=~1.0*{recon6vars[1]}
   recon6=~1.0*{recon6vars[2]}
   recon7=~1.0*{recon7vars[1]}
   recon7=~1.0*{recon7vars[2]}
   cat=~1.0*{catvars[1]}
   cat=~1.0*{catvars[2]}
! residuals, variances and covariances
   {obsvars[1]} ~~ obsvar*{obsvars[1]}
   {obsvars[2]} ~~ obsvar*{obsvars[2]}
   {obsvars[3]} ~~ obsvar*{obsvars[3]}
   {obsvars[4]} ~~ obsvar*{obsvars[4]}
   {obsvars[5]} ~~ obsvar*{obsvars[5]}
   {obsvars[6]} ~~ obsvar*{obsvars[6]}
   {obsvars[7]} ~~ obsvar*{obsvars[7]}
   {obsvars[8]} ~~ obsvar*{obsvars[8]}
   truth ~~ truevar*truth
   SKYRA ~~ SKYRAvar*SKYRA
   VERIO ~~ VERIOvar*VERIO
   samseg ~~ samsegvar*samseg
   recon6 ~~ recon6var*recon6
   recon7 ~~ recon7var*recon7
   cat ~~ catvar*cat
   truth ~~ 0.0*SKYRA
   truth ~~ 0.0*VERIO
   truth ~~ 0.0*samseg
   truth ~~ 0.0*recon6
   truth ~~ 0.0*recon7
   truth ~~ 0.0*cat
   SKYRA ~~ 0.0*VERIO
   SKYRA ~~ 0.0*samseg
   SKYRA ~~ 0.0*recon6
   SKYRA ~~ 0.0*recon7
   SKYRA ~~ 0.0*cat
   VERIO ~~ 0.0*samseg
   VERIO ~~ 0.0*recon6
   VERIO ~~ 0.0*recon7
   VERIO ~~ 0.0*cat
   samseg ~~ 0.0*recon6
   samseg ~~ 0.0*recon7
   samseg ~~ 0.0*cat
   recon6 ~~ 0.0*recon7
   recon6 ~~ 0.0*cat
   recon7 ~~ 0.0*cat
! means
   # the means are relative to samseg-SKYRA
   SKYRA~0*1;
   samseg~0*1;
   VERIO~VERIOmean*1
   recon6~recon6mean*1
   recon7~recon7mean*1
   cat~catmean*1
   truth~truemean*1
   {obsvars[1]}~0*1;
   {obsvars[2]}~0*1;
   {obsvars[3]}~0*1;
   {obsvars[4]}~0*1;
   {obsvars[5]}~0*1;
   {obsvars[6]}~0*1;
   {obsvars[7]}~0*1;
   {obsvars[8]}~0*1;
                ")
  # adapt model if the ROI is TIV (remove SKYRA-recon6)
  if (outcome == "TIV"){
    drop_lines <- c(
      "TIV_SKYRA_recon6~0*1;",
      "TIV_SKYRA_recon6 ~~ obsvar*TIV_SKYRA_recon6",
      "SKYRA=~1.0*TIV_SKYRA_recon6",
      "truth=~1.0*TIV_SKYRA_recon6",
      "recon6=~1.0*TIV_SKYRA_recon6")
    pattern <- paste0(
      "(?m)^\\s*(",
      paste0("\\Q", drop_lines, "\\E", collapse = "|"),
      ")\\s*$")
    model <- str_remove_all(model, pattern)
  }
  # iterate over datasets
  for (df_name in names(wide_dfs)){
    df <- wide_dfs[[df_name]]
    # estimate the lavaan model
    fit <- lavaan(model = model, data=df, fixed.x=FALSE, missing="FIML")
    pe <- parameterEstimates(fit)
    pe <- pe[!duplicated(pe$label, na.rm = T), ]
    # set negative variances to 0 (might occur if true variance is very close to 0)
    pe$est[pe$est < 0 & grepl("var", pe$label)] <- 0
    # prepare iced icc list
    iced_iccs <- list(
      SKYRA_samseg = NULL,
      SKYRA_recon6 = NULL,
      SKYRA_recon7 = NULL,
      SKYRA_cat = NULL,
      VERIO_samseg = NULL,
      VERIO_recon6 = NULL,
      VERIO_recon7 = NULL,
      VERIO_cat = NULL
    )
    # calculate ICED-ICC for each batch
    for (scanner in c("SKYRA", "VERIO")) {
      for (segmentation in c("samseg", "recon6", "recon7", "cat")){
       icc <- pe[pe$label == "truevar", "est"] / (pe[pe$label == "truevar", "est"] + pe[pe$label == "obsvar", "est"] +
                                                    pe[pe$label == paste0(scanner, "var"), "est"] + pe[pe$label == paste0(segmentation, "var"), "est"])
       iced_iccs[[paste0(scanner, "_", segmentation)]] <- icc
       }
    }
    # store the results
    if(df_name == "df_update_wide"){
      iced_res[iced_res$ROI == outcome, parameters] <- pe[pe$label %in% parameters, "est"]
      iced_res[iced_res$ROI == outcome, grepl("ICC", colnames(iced_res)) & !grepl("adj|correction", colnames(iced_res))] <- iced_iccs
    } else if(df_name == "df_update_adj_wide"){
      iced_res[iced_res$ROI == outcome, paste0(parameters, "_adj")] <- pe[pe$label %in% parameters, "est"]
      iced_res[iced_res$ROI == outcome, grepl("ICC", colnames(iced_res)) & grepl("adj", colnames(iced_res)) & !grepl("correction", colnames(iced_res))] <- iced_iccs
    } else if(df_name == "df_update_no_correction_wide"){
      iced_res[iced_res$ROI == outcome, paste0(parameters, "_no_correction")] <- pe[pe$label %in% parameters, "est"]
      iced_res[iced_res$ROI == outcome, grepl("ICC", colnames(iced_res)) & !grepl("adj", colnames(iced_res)) & grepl("correction", colnames(iced_res))] <- iced_iccs
    } else {
      iced_res[iced_res$ROI == outcome, paste0(parameters, "_no_correction_adj")] <- pe[pe$label %in% parameters, "est"]
      iced_res[iced_res$ROI == outcome, grepl("ICC", colnames(iced_res)) & grepl("adj", colnames(iced_res))] & grepl("correction", colnames(iced_res)) <- iced_iccs
    }
  }
}



