# this script runs the analyses for the segmentation harmonization study
# there are 2 datasets (Scanner update & LIFE-Adult) that have to be checked before and after harmonization 
# for each ROI, we want to calculate an ICC, a PD between all batches and an ANOVA / t-test
# for the scanner update dataset, we will additionally calculate SEM-based ICEDs
# we start with the scanner update dataset (section 1)
# 1.1 calculation of ICCs
# 1.2 calculation of PDs
# 1.3 calculation of ANOVAs followed-up by pairwise t-tests of individual batches
# 1.4 calculation of ICED to get ICED-ICCs for specific scanner-segmentation combinations before/after harmonization
# 1.5 Bland-Altman plots
# then we continue with the LIFE-Adult dataset (section 2)
# 2.1 calculation of ICCs
# 2.2 calculation of PDs
# 2.3 calculation of t-tests
# 2.4 Bland-Altman plots
# lastly we conduct exploratory analyses (section 3)
# 3.1 normally harmonized BL/FU data
# 3.1.1: Calculation of ICCs
# 3.1.2: Calculation of PDs
# 3.1.3: Calculation of t-tests
# 3.1.4 Bland-Altman plots
# 3.2 immediately harmonized change values
# 3.2.1: Calculation of ICCs
# 3.3 Scanner update WMHV

# load required packages
library(tidyverse) # version 2.0.0
library(lme4) # version 1.1-38
library(rstatix) # version 0.7.3
library(glue) # version 1.8.0
library(lavaan) # version 0.6-21
library(BlandAltmanLeh) # version 0.3.1
library(patchwork) # version 1.3.1
library(showtext) # 0.9-7

set.seed(1848)

############ Section 1: The scanner update dataset #######################
setwd("/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/scanner_update/Data/")
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
df_update_adj <- read.csv("df_update_adj.csv", col.names = colnames(df_update)[5:11])
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
df_update_no_correction_adj <- read.csv("df_update_no_correction_adj.csv", col.names = colnames(df_update_no_correction)[5:11])
# fill in the missing data from the unadjusted df (these cols must not be included during harmonization) 
df_update_no_correction_adj[,c("id", "scanner", "segmentation", "SITE")] <- df_update_no_correction[,c("id", "scanner", "segmentation", "SITE")]

# store the dfs in a list
update_dfs <- list(df_update = df_update, df_update_adj = df_update_adj, df_update_no_correction = df_update_no_correction,
            df_update_no_correction_adj = df_update_no_correction_adj)

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
update_wide_dfs <- list(df_update_wide = df_update_wide, df_update_adj_wide = df_update_adj_wide, 
                 df_update_no_correction_wide = df_update_no_correction_wide, df_update_no_correction_adj_wide = df_update_no_correction_adj_wide)

# produce dfs without NAs for pairise t-tests
update_pairwise_dfs <- map(
  update_dfs,
  ~ .x %>%
    group_by(id) %>%
    filter(any(scanner == "SKYRA") & any(scanner == "VERIO")) %>%
    ungroup())

############ Section 1.1: Calculation of ICCs ####################### 
# prepare a results df to be filled
update_lme_icc_res <- data.frame(matrix(ncol = 61, nrow = 7))
colnames(update_lme_icc_res) <- c("ROI", paste0(rep(paste0(rep(paste0(rep(c("id_var", "scanner_var", "segmentation_var", "resdiual_var", "ICC"), each = 3), 
                                       c("", "_lower_95_CI", "_upper_95_CI")), each = 2), c("", "_no_correction")), times = 2), 
                                       rep(c("", "_adj"), each = 30)))
update_lme_icc_res$ROI <- c("TIV", "TGMV", "TWMV", "TCV", "LVV",  "HCV", "AV")

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
  for (df_name in names(update_dfs)){
    df <- update_dfs[[df_name]] 
    suffix <- substring(df_name, 10)
    # remove SKYRA-recon6 rows if outcome is TIV (uninformative compared to VERIO-recon6)
    if (outcome == "TIV"){
      df <- df[!(df$scanner == "SKYRA" & df$segmentation == "recon6"),]
    }
    # create outcome specific LME formula
    formula <- as.formula(paste0(outcome, " ~ (1|id) + (1|scanner) + (1|segmentation)"))
    # estimate model and bootstrap confidence intervals
    print(paste0("now doing ", outcome))
    res <- lmer(formula = formula, data = df, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
    variances <- as.data.frame(VarCorr(res))
    lme_icc <- (variances[variances$grp == "id", "vcov"])/
      (variances[variances$grp == "scanner", "vcov"] + variances[variances$grp == "id", "vcov"] + 
         variances[variances$grp == "segmentation", "vcov"] + variances[variances$grp == "Residual", "vcov"])
    boot_icc <- bootMer(res, FUN  = icc_ci, nsim = 1000, type = "parametric")
    boot_res <- apply(boot_icc$t, 2, quantile, probs = c(0.025, 0.975))
    # store the results in the appropriate columns (unadjusted/adjusted data)
    update_lme_icc_res[update_lme_icc_res$ROI == outcome, paste0(c("id_var", "scanner_var", "segmentation_var", "resdiual_var"), suffix)] <- 
      c(variances[variances$grp == "id", "vcov"], variances[variances$grp == "scanner", "vcov"], 
        variances[variances$grp == "segmentation", "vcov"], variances[variances$grp == "Residual", "vcov"])
    update_lme_icc_res[update_lme_icc_res$ROI == outcome, paste0("ICC", suffix)] <- lme_icc
    update_lme_icc_res[update_lme_icc_res$ROI == outcome, paste0(paste0(rep(c("id_var", "scanner_var", "segmentation_var", "resdiual_var", "ICC"),
                                                              each = 2), c("_lower_95_CI", "_upper_95_CI")), suffix)] <- as.vector(boot_res)
  }
}

############ Section 1.2: Calculation of PDs #######################
# prepare a list to store the pd_matrices
update_pd_matrices <- list()
# iterate over ROIs
for (outcome in c("TIV", "TGMV", "TWMV", "TCV", "LVV",  "HCV", "AV")) {
  # iterate over datasets
  for (df_name in names(update_wide_dfs)){
    df <- update_wide_dfs[[df_name]]
    suffix <- substr(df_name,10,nchar(df_name)-5)
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
          pd <- (2 * (col1 - col2)) / (abs(col1) + abs(col2))
          # calculate the mean PD between the two batches (NAs exist because 5 participants were not scanned on both scanners)
          mean_pd_mat[i, j] <- mean(pd, na.rm = TRUE)
        }
      }
    }
    # mirror the inverse to lower triangle
    mean_pd_mat[lower.tri(mean_pd_mat)] <- t(mean_pd_mat)[lower.tri(mean_pd_mat)] * (-1)
    # turn into percentage values
    mean_pd_mat <- mean_pd_mat*100
    # save variables to list
    update_pd_matrices[[paste0(outcome, "_mean_pd_mat",suffix)]] <- mean_pd_mat
  }
}

############ Section 1.3: Calculation of ANOVAS and pairwise t-tests #######################
# prepare a results df to be filled
update_anova_res <- data.frame(matrix(ncol = 25, nrow = 7))
colnames(update_anova_res) <- c("ROI", paste0(rep(c("DFn",  "DFd",  "F",  "p",  "p<.05",  "ges"), times = 4), 
                                       rep(c("", "_adj", "_no_correction", "_no_correction_adj"), each = 6)))
update_anova_res$ROI <- c("TIV", "TGMV", "TWMV", "TCV", "LVV",  "HCV", "AV")
#prepare list to store pairwise comparisons
update_pwcs <- list()
# iterate over ROIs
for (outcome in c("TIV", "TGMV", "TWMV", "TCV", "LVV",  "HCV", "AV")) {
  # iterate over datasets
  for (df_name in names(update_dfs)){
    df <- update_dfs[[df_name]]
    suffix <- substring(df_name, 10)
    pairwise_df <- update_pairwise_dfs[[df_name]]
    # remove SKYRA-recon6 rows if outcome is TIV (uninformative compared to VERIO-recon6)
    if (outcome == "TIV"){
      df <- df[!(df$scanner == "SKYRA" & df$segmentation == "recon6"),]
      pairwise_df <- pairwise_df[!(pairwise_df$scanner == "SKYRA" & pairwise_df$segmentation == "recon6"),]
    }
    res.aov <- anova_test(data = df, dv = outcome, wid = id, within = SITE)
    table <- get_anova_table(res.aov)
    update_anova_res[update_anova_res$ROI == outcome, paste0(c("DFn",  "DFd",  "F",  "p",  "p<.05",  "ges"), suffix)] <- as.data.frame(table)[1,2:7]
    # perform pairwise comparisons of the individual scanner-segmentation comparisons
    formula <- as.formula(paste0(outcome, " ~ SITE"))
    pwc <- pairwise_df %>%
      pairwise_t_test(
        formula = formula, paired = TRUE, id = "id",
        p.adjust.method = "hommel"
      )
    update_pwcs[[paste0(outcome, "_pwc", suffix)]] <- pwc
  }
}

############ Section 1.4: ICED #######################
# prepare a results df to be filled
update_iced_res <- data.frame(matrix(ncol = 85, nrow = 7))
parameters <- c("obsvar", "truevar", "SKYRAvar",  "VERIOvar", "samsegvar", "recon6var", "recon7var", "catvar",
                "VERIOmean",  "recon6mean", "recon7mean", "catmean", "truemean")
col_names <- c(parameters, paste0(rep(c("ICC_SKYRA_", "ICC_VERIO_"), each = 4), rep(c("samseg", "recon6", "recon7", "cat"), times = 2))) 
                
colnames(update_iced_res) <- c("ROI", paste0(rep(col_names, times = 4), rep(c("", "_adj", "_no_correction", "_no_correction_adj"), each = 21)))
update_iced_res$ROI <- c("TIV", "TGMV", "TWMV", "TCV", "LVV",  "HCV", "AV")
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
   {obsvars[1]} ~~ start(0.01)*obsvar*{obsvars[1]}
   {obsvars[2]} ~~ start(0.01)*obsvar*{obsvars[2]}
   {obsvars[3]} ~~ start(0.01)*obsvar*{obsvars[3]}
   {obsvars[4]} ~~ start(0.01)*obsvar*{obsvars[4]}
   {obsvars[5]} ~~ start(0.01)*obsvar*{obsvars[5]}
   {obsvars[6]} ~~ start(0.01)*obsvar*{obsvars[6]}
   {obsvars[7]} ~~ start(0.01)*obsvar*{obsvars[7]}
   {obsvars[8]} ~~ start(0.01)*obsvar*{obsvars[8]}
   truth ~~ start(0.1)*truevar*truth
   SKYRA ~~ start(0.0001)*SKYRAvar*SKYRA
   VERIO ~~ start(0.0001)*VERIOvar*VERIO
   samseg ~~ start(0.001)*samsegvar*samseg
   recon6 ~~ start(0.001)*recon6var*recon6
   recon7 ~~ start(0.001)*recon7var*recon7
   cat ~~ start(0.001)*catvar*cat
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
  for (df_name in names(update_wide_dfs)){
    df <- update_wide_dfs[[df_name]]
    suffix <- substr(df_name,10,nchar(df_name)-5)
    # estimate the lavaan model
    fit <- lavaan(model = model, data=df, fixed.x=FALSE, missing="FIML", control = list(iter.max = 5000))
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
    update_iced_res[update_iced_res$ROI == outcome, paste0(parameters, suffix)] <- pe[pe$label %in% parameters, "est"]
    update_iced_res[update_iced_res$ROI == outcome, c(paste0("ICC_SKYRA_", c("samseg", "recon6", "recon7", "cat"), suffix), 
                                        paste0("ICC_VERIO_", c("samseg", "recon6", "recon7", "cat"), suffix))] <- iced_iccs
    # check if non-convergences occured (lavaan just gives starting values in this case)
    if (!lavaan::lavInspect(fit, "converged")) {
      update_iced_res[update_iced_res$ROI == outcome, paste0(parameters, suffix)] <- "uconverged"
      update_iced_res[update_iced_res$ROI == outcome, c(paste0("ICC_SKYRA_", c("samseg", "recon6", "recon7", "cat"), suffix), 
                      paste0("ICC_VERIO_", c("samseg", "recon6", "recon7", "cat"), suffix))] <- "unconverged"
    }
  }
}

############ Section 1.5: Bland-Altman plots #######################
showtext_auto()
update_ba_plots <- list()
update_ba_plot_comparisons <- list()
for (outcome in c("TIV", "TGMV", "TWMV", "TCV", "LVV",  "HCV", "AV")) {
  # iterate over datasets
  for (df_name in names(update_wide_dfs)){
    plot_list <- list()
    df <- update_wide_dfs[[df_name]]
    suffix <- substring(df_name, 10)
    cols <- grep(outcome, colnames(df), value = T)
    # get all combinations of cols
    combs <- combn(cols, 2, simplify = FALSE)
    for (combination in combs){
      plt <- bland.altman.plot(group1 = df[[combination[1]]], group2 = df[[combination[2]]], 
                               graph.sys = "ggplot2", conf.int=.95, pch=19)
      plt <- plt + labs(title = paste0(gsub("_", " ", combination[1]), " vs. ", gsub("_", " ", combination[2]))) + theme_bw(base_family = "garamond") + 
        theme(plot.title = element_text(hjust = 0.5, size = 12))
      plot_list[[paste0(combination[1], " vs. ", combination[2])]] <- plt
    }
    plots <- wrap_plots(plot_list, ncol(6))
    update_ba_plots[[paste0(outcome, suffix)]] <- plots
    update_ba_plots[[paste0(outcome, suffix, "_indiv_plots")]] <- plot_list
  }
  for (i in seq(1, 28, by = 4)) {
    plots_left  <- update_ba_plots[[paste0(outcome, "_wide_indiv_plots")]][i:(i+3)]
    plots_right <- update_ba_plots[[paste0(outcome, "_adj_wide_indiv_plots")]][i:(i+3)]
    left_col  <- wrap_plots(plots_left, ncol = 1)
    right_col <- wrap_plots(plots_right, ncol = 1)
    combined <- wrap_plots(left_col, right_col, ncol = 2) + plot_annotation(title = paste0(outcome, " before (left) and after (right) harmonization"), 
                                                                            theme = theme(plot.title = element_text(hjust = 0.5)))
    update_ba_plot_comparisons[[paste0(outcome, (((i-1)/4)+1))]] <- combined
  }
}
############ Section 2: The LIFE-Adult dataset #######################
# switch working directory
setwd("/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/life/Data/")
# load the unadjusted data from the LIFE_Adult study
df_life <- read.csv(file = "df_life.csv") 
df_life_wmhv <- read.csv(file = "df_life_wmhv.csv") 
df_life_no_correction <- read.csv(file = "df_life_no_correction.csv") 

# load the data adjusted using loongitudinal ComBat
df_life_combat <- read.csv(file = "df_life_adj.csv") 
df_life_wmhv_combat <- read.csv(file = "df_life_wmhv_adj.csv") 
df_life_combat_no_correction <- read.csv(file = "df_life_no_correction_adj.csv") 
# align names with other dfs (delete suffix .combat)
colnames(df_life_combat) <- gsub("\\.combat$", "", colnames(df_life_combat))
colnames(df_life_wmhv_combat) <- gsub("\\.combat$", "", colnames(df_life_wmhv_combat))
colnames(df_life_combat_no_correction) <- gsub("\\.combat$", "", colnames(df_life_combat_no_correction))

# load the data adjusted using neuroHarmonize and their unadjusted counterparts
df_life_nh_bl <- read.csv(file = "df_life_bl_adj.csv")
df_life_nh_fu <- read.csv(file = "df_life_fu_adj.csv")
df_life_bl <- read.csv(file = "df_life_bl.csv")
df_life_fu <- read.csv(file = "df_life_fu.csv")
df_life_nh_bl[,c("id", "timepoint", "SITE", "AGE", "GENDER")] <- df_life_bl[,c("id", "timepoint", "SITE", "AGE", "GENDER")]
df_life_nh_fu[,c("id", "timepoint", "SITE", "AGE", "GENDER")] <- df_life_fu[,c("id", "timepoint", "SITE", "AGE", "GENDER")]
df_life_nh <- rbind(df_life_nh_bl, df_life_nh_fu)

df_life_no_correction_nh_bl <- read.csv(file = "df_life_no_correction_bl_adj.csv")
df_life_no_correction_nh_fu <- read.csv(file = "df_life_no_correction_fu_adj.csv")
df_life_no_correction_bl <- read.csv(file = "df_life_no_correction_bl.csv")
df_life_no_correction_fu <- read.csv(file = "df_life_no_correction_fu.csv")
df_life_no_correction_nh_bl[,c("id", "timepoint", "SITE", "AGE", "GENDER")] <- df_life_no_correction_bl[,c("id", "timepoint", "SITE", "AGE", "GENDER")]
df_life_no_correction_nh_fu[,c("id", "timepoint", "SITE", "AGE", "GENDER")] <- df_life_no_correction_fu[,c("id", "timepoint", "SITE", "AGE", "GENDER")]
df_life_nh_no_correction <- rbind(df_life_no_correction_nh_bl, df_life_no_correction_nh_fu)

df_life_wmhv_nh_bl <- read.csv(file = "df_life_wmhv_bl_adj.csv")
df_life_wmhv_nh_fu <- read.csv(file = "df_life_wmhv_fu_adj.csv")
df_life_wmhv_bl <- read.csv(file = "df_life_wmhv_bl.csv")
df_life_wmhv_fu <- read.csv(file = "df_life_wmhv_fu.csv")
df_life_wmhv_nh_bl[,c("id", "timepoint", "SITE", "AGE", "GENDER")] <- df_life_wmhv_bl[,c("id", "timepoint", "SITE", "AGE", "GENDER")]
df_life_wmhv_nh_fu[,c("id", "timepoint", "SITE", "AGE", "GENDER")] <- df_life_wmhv_fu[,c("id", "timepoint", "SITE", "AGE", "GENDER")]
df_life_wmhv_nh <- rbind(df_life_wmhv_nh_bl, df_life_wmhv_nh_fu)
# align names with other dfs (delete prefix X..)
colnames(df_life_nh) <- gsub("^X\\.\\.", "", colnames(df_life_nh))
colnames(df_life_wmhv_nh) <- gsub("^X\\.\\.", "", colnames(df_life_wmhv_nh))
colnames(df_life_nh_no_correction) <- gsub("^X\\.\\.", "", colnames(df_life_nh_no_correction))

# for a sensitivity analysis, only consider BL/FU values
df_life_nh_bl <- df_life_nh[df_life_nh$timepoint == "BL",]
df_life_combat_bl <- df_life_combat[df_life_combat$timepoint == "BL",]
df_life_wmhv_nh_bl <- df_life_wmhv_nh[df_life_wmhv_nh$timepoint == "BL",]
df_life_wmhv_combat_bl <- df_life_wmhv_combat[df_life_wmhv_combat$timepoint == "BL",]

df_life_nh_fu <- df_life_nh[df_life_nh$timepoint == "FU",]
df_life_combat_fu <- df_life_combat[df_life_combat$timepoint == "FU",]
df_life_wmhv_nh_fu <- df_life_wmhv_nh[df_life_wmhv_nh$timepoint == "FU",]
df_life_wmhv_combat_fu <- df_life_wmhv_combat[df_life_wmhv_combat$timepoint == "FU",]

# for a further sensitivity analysis, load the data from the cross pipeline
df_life_cross_combat <- read.csv(file = "df_life_cross_adj.csv") 
colnames(df_life_cross_combat) <- gsub("\\.combat$", "", colnames(df_life_cross_combat))
df_life_cross_combat$SITE <- gsub("_cross", "", df_life_cross_combat$SITE)
# load the data adjusted using neuroHarmonize and their unadjusted counterparts
df_life_nh_cross_bl <- read.csv(file = "df_life_cross_bl_adj.csv")
df_life_nh_cross_fu <- read.csv(file = "df_life_cross_fu_adj.csv")
df_life_cross_bl <- read.csv(file = "df_life_cross_bl.csv")
df_life_cross_fu <- read.csv(file = "df_life_cross_fu.csv")
df_life_nh_cross_bl[,c("id", "timepoint", "SITE", "AGE", "GENDER")] <- df_life_cross_bl[,c("id", "timepoint", "SITE", "AGE", "GENDER")]
df_life_nh_cross_fu[,c("id", "timepoint", "SITE", "AGE", "GENDER")] <- df_life_cross_fu[,c("id", "timepoint", "SITE", "AGE", "GENDER")]
df_life_cross_nh <- rbind(df_life_nh_cross_bl, df_life_nh_cross_fu)
colnames(df_life_cross_nh) <- gsub("^X\\.\\.", "", colnames(df_life_cross_nh))
df_life_cross_nh$SITE <- gsub("_cross", "", df_life_cross_nh$SITE)

# store the dfs in a list
life_dfs <- list(df_life = df_life, df_life_combat = df_life_combat, df_life_nh = df_life_nh,
                 df_life_no_correction = df_life_no_correction, df_life_combat_no_correction = 
                   df_life_combat_no_correction, df_life_nh_no_correction = df_life_nh_no_correction, 
                 df_life_cross_combat = df_life_cross_combat, df_life_cross_nh = df_life_cross_nh)
life_wmhv_dfs <- list(df_life_wmhv = df_life_wmhv, df_life_wmhv_combat = df_life_wmhv_combat, df_life_wmhv_nh = df_life_wmhv_nh)

# turn dfs into wide format
life_wide_dfs <- lapply(life_dfs, function(x)
  tidyr::pivot_wider(x,
                     names_from = c(SITE, timepoint),
                     values_from = c(TGMV, TCV, LVV, HCV, AV),
                     id_cols = id  
  )
)
life_wide_wmhv_dfs <- lapply(life_wmhv_dfs, function(x)
  tidyr::pivot_wider(x,
                     names_from = c(SITE, timepoint),
                     values_from = c(WMHV_norm),
                     id_cols = id  
  )
)
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
life_wide_dfs <- lapply(life_wide_dfs, function(df) {
  for (nm in names(pairs)) {
    bl_col <- pairs[[nm]][1]
    fu_col <- pairs[[nm]][2]
    df[[paste0("d", nm)]] <- df[[fu_col]] - df[[bl_col]]
  }
  df
})
pairs_wmhv <- list(
  WMHV_samseg = c("samseg_BL", "samseg_FU"),
  WMHV_lst = c("lst_BL", "lst_FU"))
life_wide_wmhv_dfs <- lapply(life_wide_wmhv_dfs, function(df) {
  for (nm in names(pairs_wmhv)) {
    bl_col <- pairs_wmhv[[nm]][1]
    fu_col <- pairs_wmhv[[nm]][2]
    df[[paste0("d", nm)]] <- df[[fu_col]] - df[[bl_col]]
  }
  df
})
# now turn them back into a (somwhat) long format
life_long_dfs <- lapply(life_wide_dfs, function(x)
  tidyr::pivot_longer(x,
                      cols = grep("^d", colnames(life_wide_dfs$df_life), value = T),
                      names_to = c("variable", "segmentation"), 
                      names_sep = "_",                    
                      values_to = "value" 
  ) %>%
    select(id, variable, segmentation, value) %>%
    pivot_wider(
      names_from = "variable",  # each variable becomes a column
      values_from = "value"
    )
)
life_long_wmhv_dfs <- lapply(life_wide_wmhv_dfs, function(x)
  tidyr::pivot_longer(x,
                      cols = c(dWMHV_samseg, dWMHV_lst),
                      names_to = c("variable", "segmentation"), 
                      names_sep = "_",                    
                      values_to = "value" 
  ) %>%
    select(id, variable, segmentation, value) %>%
    pivot_wider(
      names_from = "variable",  # each variable becomes a column
      values_from = "value"
    )
) 
# produce dfs without NAs for pairise t-tests
life_pairwise_dfs <- map(
  life_wide_dfs,
  ~ .x %>%
    select(id, grep("^d", colnames(.x), value = T)) %>%
    pivot_longer(
      cols = grep("^d", colnames(.x), value = T),
      names_to = c("measure", "segmentation"),
      names_sep = "_",
      values_to = "value"
    ) %>% 
    drop_na() %>%
    group_by(id, measure) %>%
    filter(any(segmentation == "samseg") & any(segmentation == "recon5")) %>%
    ungroup()
)
life_pairwise_wmhv_dfs <- map(
  life_wide_wmhv_dfs,
  ~ .x %>%
    select(id, grep("^d", colnames(.x), value = T)) %>%
    pivot_longer(
      cols = grep("^d", colnames(.x), value = T),
      names_to = c("measure", "segmentation"),
      names_sep = "_",
      values_to = "value"
    ) %>% 
    drop_na() %>%
    group_by(id, measure) %>%
    filter(any(segmentation == "samseg") & any(segmentation == "lst")) %>%
    ungroup()
)
############ Section 2.1: Calculation of ICCs ####################### 
# prepare a results df to be filled
life_lme_icc_res <- data.frame(matrix(ncol = 97, nrow = 6))
colnames(life_lme_icc_res) <- c("ROI", paste0(rep(paste0(rep(c("id_var", "segmentation_var", "resdiual_var", "ICC"), each = 3), 
                                                               c("", "_lower_95_CI", "_upper_95_CI")), each = 8), c(paste0(rep(c("", "_combat", "_nh"), each = 2), c("", "_no_correction")), "_cross_combat", "_cross_nh")))
life_lme_icc_res$ROI <- c("dTGMV", "dTCV", "dLVV",  "dHCV", "dAV", "dWMHV")

# formula for bootstrapped CIs for ICC 
icc_ci <- function(fit) {
  variances <- as.data.frame(VarCorr(fit))
  id_var = variances[variances$grp == "id", "vcov"]
  seg_var = variances[variances$grp == "segmentation", "vcov"] 
  residual_var = variances[variances$grp == "Residual", "vcov"]
  icc = id_var / (id_var + seg_var + residual_var)
  c(id_var, seg_var, residual_var, icc)
}

# calculate LMEs
# iterate over no wmhv_outcomes
for (outcome in c("dTGMV", "dTCV", "dLVV",  "dHCV", "dAV")) {
  # iterate over dfs
  for (df_name in names(life_long_dfs)){
    df <- life_long_dfs[[df_name]] 
    # get the suffix for storing results
    suffix <- substring(df_name, 8)
    # create outcome specific LME formula
    formula <- as.formula(paste0(outcome, " ~ (1|id) + (1|segmentation)"))
    # estimate model and bootstrap confidence intervals
    print(paste0("now doing ", outcome))
    res <- lmer(formula = formula, data = df, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
    variances <- as.data.frame(VarCorr(res))
    lme_icc <- (variances[variances$grp == "id", "vcov"])/
      (variances[variances$grp == "id", "vcov"] + 
         variances[variances$grp == "segmentation", "vcov"] + variances[variances$grp == "Residual", "vcov"])
    boot_icc <- bootMer(res, FUN  = icc_ci, nsim = 1000, type = "parametric")
    boot_res <- apply(boot_icc$t, 2, quantile, probs = c(0.025, 0.975))
    # store the results in the appropriate columns (unadjusted/adjusted data)
    life_lme_icc_res[life_lme_icc_res$ROI == outcome, paste0(c("id_var", "segmentation_var", "resdiual_var"), suffix)] <- 
      c(variances[variances$grp == "id", "vcov"],  
        variances[variances$grp == "segmentation", "vcov"], variances[variances$grp == "Residual", "vcov"])
    life_lme_icc_res[life_lme_icc_res$ROI == outcome, paste0("ICC", suffix)] <- lme_icc
    life_lme_icc_res[life_lme_icc_res$ROI == outcome, paste0(paste0(rep(c("id_var", "segmentation_var", "resdiual_var", "ICC"),
                                                       each = 2), c("_lower_95_CI", "_upper_95_CI")), suffix)] <- as.vector(boot_res)
  }
}
# iterate over no wmhv dfs

for (df_name in names(life_long_wmhv_dfs)){
  df <- life_long_wmhv_dfs[[df_name]] 
  # get the suffix for storing results
  suffix <- substring(df_name, 13)
  # estimate model and bootstrap confidence intervals
  res <- lmer(formula = dWMHV ~ (1|id) + (1|segmentation), data = df, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
  variances <- as.data.frame(VarCorr(res))
  lme_icc <- (variances[variances$grp == "id", "vcov"])/
    (variances[variances$grp == "id", "vcov"] + 
       variances[variances$grp == "segmentation", "vcov"] + variances[variances$grp == "Residual", "vcov"])
  boot_icc <- bootMer(res, FUN  = icc_ci, nsim = 1000, type = "parametric")
  boot_res <- apply(boot_icc$t, 2, quantile, probs = c(0.025, 0.975))
  # store the results in the appropriate columns (unadjusted/adjusted data)
  life_lme_icc_res[life_lme_icc_res$ROI == "dWMHV", paste0(c("id_var", "segmentation_var", "resdiual_var"), suffix)] <- 
    c(variances[variances$grp == "id", "vcov"],  
      variances[variances$grp == "segmentation", "vcov"], variances[variances$grp == "Residual", "vcov"])
  life_lme_icc_res[life_lme_icc_res$ROI == "dWMHV", paste0("ICC", suffix)] <- lme_icc
  life_lme_icc_res[life_lme_icc_res$ROI == "dWMHV", paste0(paste0(rep(c("id_var", "segmentation_var", "resdiual_var", "ICC"),
                                                                      each = 2), c("_lower_95_CI", "_upper_95_CI")), suffix)] <- as.vector(boot_res)
}

############ Section 2.2: Calculation of PDs #######################
# prepare a list to store the pd_matrices
life_pd_matrices <- list()
# iterate over ROIs
for (outcome in c("TGMV", "TCV", "LVV",  "HCV", "AV")) {
  # iterate over datasets
  for (df_name in names(life_wide_dfs)){
    df <- life_wide_dfs[[df_name]]
    suffix <- substring(df_name,8)
    # get names of relevant cols
    roi_cols <- grep(paste0("d", outcome), colnames(df), value = TRUE)
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
          pd <- (2 * (col1 - col2)) / (abs(col1) + abs(col2))
          # calculate the mean PD between the two batches 
          mean_pd_mat[i, j] <- mean(pd, na.rm = TRUE)
        }
      }
    }
    # mirror the inverse to lower triangle
    mean_pd_mat[lower.tri(mean_pd_mat)] <- t(mean_pd_mat)[lower.tri(mean_pd_mat)] * (-1)
    # turn into percentage values
    mean_pd_mat <- mean_pd_mat*100
    # save variables to list
    life_pd_matrices[[paste0(outcome, "_mean_pd_mat",suffix)]] <- mean_pd_mat
  }
}
# iterate over WMHV dfs
for (df_name in names(life_wide_wmhv_dfs)){
  df <- life_wide_wmhv_dfs[[df_name]]
  suffix <- substring(df_name,13)
  # get names of relevant cols
  roi_cols <- grep("^d", colnames(df), value = TRUE)
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
        pd <- (2 * (col1 - col2)) / (abs(col1) + abs(col2))
        # calculate the mean PD between the two batches (NAs exist because 5 participants were not scanned on both scanners)
        mean_pd_mat[i, j] <- mean(pd, na.rm = TRUE)
      }
    }
  }
  # mirror the inverse to lower triangle
  mean_pd_mat[lower.tri(mean_pd_mat)] <- t(mean_pd_mat)[lower.tri(mean_pd_mat)] * (-1)
  # turn into percentage values
  mean_pd_mat <- mean_pd_mat*100
  # save variables to list
  life_pd_matrices[[paste0("WMHV_mean_pd_mat",suffix)]] <- mean_pd_mat
}

############ Section 2.3: Calculation of t-tests #######################
#prepare list to store pairwise comparisons
life_pwcs <- list()
# iterate over dfs
for (df_name in names(life_pairwise_dfs)) {
  df <- life_pairwise_dfs[[df_name]]
  suffix <- substring(df_name,8)
  results <- df %>%
    group_by(measure) %>%
    pairwise_t_test(
      value ~ segmentation,
      paired = TRUE,
      p.adjust.method = "none")
  life_pwcs[[paste0("pwc", suffix)]] <- results
}
for (df_name in names(life_pairwise_wmhv_dfs)) {
  df <- life_pairwise_wmhv_dfs[[df_name]]
  suffix <- substring(df_name,8)
  results <- df %>%
    t_test(value ~ segmentation, paired = TRUE)
  life_pwcs[[paste0("pwc", suffix)]] <- results
}

############ Section 2.4: Bland-Altman plots #######################
life_ba_plots <- list()
wrapped_life_ba_plots <- list()
for (outcome in c("dTGMV", "dTCV", "dLVV",  "dHCV", "dAV")) {
  # iterate over datasets
  for (df_name in names(life_wide_dfs)){
    df <- life_wide_dfs[[df_name]]
    suffix <- substring(df_name, 8)
    cols <- grep(outcome, colnames(df), value = T)
    # add a vertical line to show the size of 1 SD
    # determine the x-coordinate (should be left of the smallest value)
    x_coord = min(rowMeans(df[,cols]), na.rm = T)
    # determine the y-coordinates 
    sd_val <- sd(unlist(df[cols]), na.rm = TRUE)
    y_coord_a <- sd_val/2 
    y_coord_b <- -sd_val/2
    plt <- bland.altman.plot(group1 = df[[cols[1]]], group2 = df[[cols[2]]], 
                             graph.sys = "ggplot2", conf.int=.95, pch=19)
    plt <- plt + labs(title = paste0(cols[1], " vs. ", cols[2])) + theme_bw(base_family = "garamond") + theme(plot.title = element_text(family = "garamond", hjust = 0.5, size = 13))
    # add line to idicate size of SD
    plt <- plt + 
      geom_segment(
        x = x_coord, xend = x_coord,
        y = y_coord_a, yend = y_coord_b,
        colour = "blue",
        linewidth = 1
      )
    life_ba_plots[[paste0(outcome, suffix)]] <- plt
  }
  unharm_plot <- life_ba_plots[[outcome]] + labs(title = paste0("\u0394", substring(outcome, 2), " samseg vs. ", "\u0394", substring(outcome, 2), " recon5 unharmonized")) 
  nh_plot <- life_ba_plots[[paste0(outcome, "_nh")]] + labs(title = paste0("\u0394", substring(outcome, 2), " samseg vs. ", "\u0394", substring(outcome, 2), " recon5 neuroHarmonize")) 
  combat_plot <- life_ba_plots[[paste0(outcome, "_combat")]] + labs(title = paste0("\u0394", substring(outcome, 2), " samseg vs. ", "\u0394", substring(outcome, 2), " recon5 longitudinal ComBat")) 
  plots <- wrap_plots(list(unharm_plot, nh_plot, combat_plot), nrow = 3, ) & theme(plot.margin = unit(c(2, 2, 8, 2), units = "mm"))
  wrapped_life_ba_plots[[outcome]] <- plots 
}
############ Section 3: Exploratory analyses #######################  
############ Section 3.1 normally harmonized BL/FU data ################
# collect dfs in lists
tp_life_dfs <- list(df_life_bl = df_life_bl, df_life_combat_bl = df_life_combat_bl, df_life_nh_bl = df_life_nh_bl, 
                    df_life_fu = df_life_fu, df_life_combat_fu = df_life_combat_fu, df_life_nh_fu = df_life_nh_fu)
tp_life_wmhv_dfs <- list(df_life_wmhv_bl = df_life_wmhv_bl, df_life_wmhv_combat_bl = df_life_wmhv_combat_bl, df_life_wmhv_nh_bl = df_life_wmhv_nh_bl,
                         df_life_wmhv_fu = df_life_wmhv_fu, df_life_wmhv_combat_fu = df_life_wmhv_combat_fu, df_life_wmhv_nh_fu = df_life_wmhv_nh_fu)
# rename SITE as segmentation for clarity
tp_life_dfs <- lapply(tp_life_dfs, function(df) {
  names(df)[names(df) == "SITE"] <- "segmentation"
  df
})
tp_life_wmhv_dfs <- lapply(tp_life_wmhv_dfs, function(df) {
  names(df)[names(df) == "SITE"] <- "segmentation"
  df
})
# turn into wide format
tp_life_wide_dfs <- lapply(tp_life_dfs, function(df) {
  tidyr::pivot_wider(df,
                      names_from = c(segmentation),
                      values_from = c(TGMV, TCV, LVV, HCV, AV),
                      id_cols = id)
})
tp_life_wmhv_wide_dfs <- lapply(tp_life_wmhv_dfs, function(df) {
  tidyr::pivot_wider(df,
                     names_from = c(segmentation),
                     values_from = c(WMHV_norm),
                     id_cols = id)
})
# turn into dfs with pairwise complete cases for samseg/recon5
tp_life_long_dfs <- lapply(tp_life_wide_dfs, function(df) {
  tidyr::pivot_longer(df,
                      cols = -c(id),
                      names_to = c("measure", "segmentation"),
                      names_sep = "_",
                      values_to = "value")
})
tp_life_wmhv_long_dfs <- lapply(tp_life_wmhv_wide_dfs, function(df) {
  tidyr::pivot_longer(df,
                      cols = -c(id),
                      names_to = c("segmentation"),
                      values_to = "value")
})
############ Section 3.1.1: Calculation of ICCs ####################### 
# prepare a results df to be filled
tp_life_lme_icc_res <- data.frame(matrix(ncol = 37, nrow = 6))
colnames(tp_life_lme_icc_res) <- c("ROI", paste0(rep(paste0(rep(c("id_var", "segmentation_var", "resdiual_var", "ICC"), each = 3), 
                                                                    c("", "_lower_95_CI", "_upper_95_CI")), each = 3), c("_fu", "_combat_fu", "_nh_fu")))
tp_life_lme_icc_res$ROI <- c("TGMV", "TCV", "LVV",  "HCV", "AV", "WMHV_norm")

# formula for bootstrapped CIs for ICC
icc_ci <- function(fit) {
  variances <- as.data.frame(VarCorr(fit))
  id_var = variances[variances$grp == "id", "vcov"]
  seg_var = variances[variances$grp == "segmentation", "vcov"] 
  residual_var = variances[variances$grp == "Residual", "vcov"]
  icc = id_var / (id_var + seg_var + residual_var)
  c(id_var, seg_var, residual_var, icc)
}

# calculate LMEs
# iterate over no wmhv_outcomes
for (outcome in c("TGMV", "TCV", "LVV",  "HCV", "AV")) {
  # iterate over dfs
  for (df_name in names(tp_life_dfs)){
    df <- tp_life_dfs[[df_name]] 
    # get the suffix for storing results
    suffix <- substring(df_name, 8)
    # create outcome specific LME formula
    formula <- as.formula(paste0(outcome, " ~ (1|id) + (1|segmentation)"))
    # estimate model and bootstrap confidence intervals
    print(paste0("now doing ", outcome))
    res <- lmer(formula = formula, data = df, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
    variances <- as.data.frame(VarCorr(res))
    lme_icc <- (variances[variances$grp == "id", "vcov"])/
      (variances[variances$grp == "id", "vcov"] + 
         variances[variances$grp == "segmentation", "vcov"] + variances[variances$grp == "Residual", "vcov"])
    boot_icc <- bootMer(res, FUN  = icc_ci, nsim = 1000, type = "parametric")
    boot_res <- apply(boot_icc$t, 2, quantile, probs = c(0.025, 0.975))
    # store the results in the appropriate columns (unadjusted/adjusted data)
    tp_life_lme_icc_res[tp_life_lme_icc_res$ROI == outcome, paste0(c("id_var", "segmentation_var", "resdiual_var"), suffix)] <- 
      c(variances[variances$grp == "id", "vcov"],  
        variances[variances$grp == "segmentation", "vcov"], variances[variances$grp == "Residual", "vcov"])
    tp_life_lme_icc_res[tp_life_lme_icc_res$ROI == outcome, paste0("ICC", suffix)] <- lme_icc
    tp_life_lme_icc_res[tp_life_lme_icc_res$ROI == outcome, paste0(paste0(rep(c("id_var", "segmentation_var", "resdiual_var", "ICC"),
                                                                        each = 2), c("_lower_95_CI", "_upper_95_CI")), suffix)] <- as.vector(boot_res)
  }
}
# iterate over wmhv dfs

for (df_name in names(tp_life_wmhv_dfs)){
  df <- tp_life_wmhv_dfs[[df_name]] 
  # get the suffix for storing results
  suffix <- substring(df_name, 13)
  # estimate model and bootstrap confidence intervals
  res <- lmer(formula = WMHV_norm ~ (1|id) + (1|segmentation), data = df, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
  variances <- as.data.frame(VarCorr(res))
  lme_icc <- (variances[variances$grp == "id", "vcov"])/
    (variances[variances$grp == "id", "vcov"] + 
       variances[variances$grp == "segmentation", "vcov"] + variances[variances$grp == "Residual", "vcov"])
  boot_icc <- bootMer(res, FUN  = icc_ci, nsim = 1000, type = "parametric")
  boot_res <- apply(boot_icc$t, 2, quantile, probs = c(0.025, 0.975))
  # store the results in the appropriate columns (unadjusted/adjusted data)
  tp_life_lme_icc_res[tp_life_lme_icc_res$ROI == "WMHV_norm", paste0(c("id_var", "segmentation_var", "resdiual_var"), suffix)] <- 
    c(variances[variances$grp == "id", "vcov"],  
      variances[variances$grp == "segmentation", "vcov"], variances[variances$grp == "Residual", "vcov"])
  tp_life_lme_icc_res[tp_life_lme_icc_res$ROI == "WMHV_norm", paste0("ICC", suffix)] <- lme_icc
  tp_life_lme_icc_res[tp_life_lme_icc_res$ROI == "WMHV_norm", paste0(paste0(rep(c("id_var", "segmentation_var", "resdiual_var", "ICC"),
                                                                      each = 2), c("_lower_95_CI", "_upper_95_CI")), suffix)] <- as.vector(boot_res)
}


############ Section 3.1.2: Calculation of PDs ####################### 
# prepare a list to store the pd_matrices
tp_life_pd_matrices <- list()
# iterate over ROIs
for (outcome in c("TGMV", "TCV", "LVV",  "HCV", "AV")) {
  # iterate over datasets
  for (df_name in names(tp_life_wide_dfs)){
    df <- tp_life_wide_dfs[[df_name]]
    suffix <- substring(df_name,8)
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
          pd <- (2 * (col1 - col2)) / (abs(col1) + abs(col2))
          # calculate the mean PD between the two batches (NAs exist because 5 participants were not scanned on both scanners)
          mean_pd_mat[i, j] <- mean(pd, na.rm = TRUE)
        }
      }
    }
    # mirror the inverse to lower triangle
    mean_pd_mat[lower.tri(mean_pd_mat)] <- t(mean_pd_mat)[lower.tri(mean_pd_mat)] * (-1)
    # turn into percentage values
    mean_pd_mat <- mean_pd_mat*100
    # save variables to list
    tp_life_pd_matrices[[paste0(outcome, "_mean_pd_mat",suffix)]] <- mean_pd_mat
  }
}

# iterate over WMHV dfs
for (df_name in names(tp_life_wmhv_wide_dfs)){
  df <- tp_life_wmhv_wide_dfs[[df_name]]
  suffix <- substring(df_name,13)
  # get names of relevant cols
  roi_cols <- c("samseg", "lst")
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
        pd <- (2 * (col1 - col2)) / (abs(col1) + abs(col2))
        # calculate the mean PD between the two batches (NAs exist because 5 participants were not scanned on both scanners)
        mean_pd_mat[i, j] <- mean(pd, na.rm = TRUE)
      }
    }
  }
  # mirror the inverse to lower triangle
  mean_pd_mat[lower.tri(mean_pd_mat)] <- t(mean_pd_mat)[lower.tri(mean_pd_mat)] * (-1)
  # turn into percentage values
  mean_pd_mat <- mean_pd_mat*100
  # save variables to list
  tp_life_pd_matrices[[paste0("WMHV_mean_pd_mat",suffix)]] <- mean_pd_mat
}


############ Section 3.1.3: Calculation of t-tests ####################### 
#prepare list to store pairwise comparisons
tp_life_pwcs <- list()
# iterate over dfs
for (df_name in names(tp_life_long_dfs)) {
  df <- tp_life_long_dfs[[df_name]]
  suffix <- substring(df_name,8)
  results <- df %>%
    group_by(measure) %>%
    pairwise_t_test(
      value ~ segmentation,
      paired = TRUE,
      p.adjust.method = "none")
  tp_life_pwcs[[paste0("pwc", suffix)]] <- results
}
for (df_name in names(tp_life_wmhv_long_dfs)) {
  df <- tp_life_wmhv_long_dfs[[df_name]]
  suffix <- substring(df_name,8)
  results <- df %>%
    t_test(value ~ segmentation, paired = TRUE)
  tp_life_pwcs[[paste0("pwc", suffix)]] <- results
}

############ Section 3.1.4 Bland-Altman plots ################
tp_life_ba_plots <- list()
tp_wrapped_life_ba_plots <- list()
for (outcome in c("TGMV", "TCV", "LVV",  "HCV", "AV")) {
  # iterate over datasets
  for (df_name in names(tp_life_wide_dfs)){
    plot_list <- list()
    df <- tp_life_wide_dfs[[df_name]]
    suffix <- substring(df_name, 8)
    cols <- grep(outcome, colnames(df), value = T)
    x_coord = min(rowMeans(df[,cols]), na.rm = T)
    sd_val <- sd(unlist(df[cols]), na.rm = TRUE)
    y_coord_a <- sd_val/2 
    y_coord_b <- -sd_val/2
    plt <- bland.altman.plot(group1 = df[[cols[1]]], group2 = df[[cols[2]]], 
                             graph.sys = "ggplot2", conf.int=.95, pch=19)
    plt <- plt + labs(title = paste0(cols[1], " vs. ", cols[2])) + theme_bw(base_family = "garamond") + theme(plot.title = element_text(hjust = 0.5, size = 13))
    plt <- plt + 
      geom_segment(
        x = x_coord, xend = x_coord,
        y = y_coord_a, yend = y_coord_b,
        colour = "blue",
        linewidth = 1
      )
    tp_life_ba_plots[[paste0(outcome, suffix)]] <- plt
  }
  for (tp in c("_bl", "_fu")){
    unharm_plot <- tp_life_ba_plots[[paste0(outcome, tp)]] + labs(title = paste0(outcome, " samseg vs. ", outcome, " recon5 unharmonized")) 
    nh_plot <- tp_life_ba_plots[[paste0(outcome, "_nh", tp)]] + labs(title = paste0(outcome, " samseg vs. ", outcome, " recon5 neuroHarmonize")) 
    combat_plot <- tp_life_ba_plots[[paste0(outcome, "_combat", tp)]] + labs(title = paste0(outcome, " samseg vs. ", outcome, " recon5 longitudinal ComBat")) 
    plots <- wrap_plots(list(unharm_plot, nh_plot, combat_plot),nrow = 3) & theme(plot.margin = unit(c(2, 2, 8, 2), units = "mm"))
    tp_wrapped_life_ba_plots[[paste0(outcome, tp)]] <- plots 
  }
}
############ Section 3.2 immediately harmonized change values ################
# load unadjusted data
change_life <- read_csv(file = "df_life_changes.csv")
wmhv_change_life <- read_csv(file = "df_life_wmhv_changes.csv")
# load adjusted data
change_life_adj <- read_csv(file = "df_life_changes_adj.csv")
colnames(change_life_adj) <- gsub("# ", "", colnames(change_life_adj))
wmhv_change_life_adj <- read_csv(file = "df_life_wmhv_changes_adj.csv")
colnames(wmhv_change_life_adj) <- gsub("# ", "", colnames(wmhv_change_life_adj))
# merge adjusted data with Site info
change_life_adj[,c("segmentation", "id")] <- change_life[,c("SITE", "id")]
wmhv_change_life_adj[,c("segmentation", "id")] <- wmhv_change_life[,c("SITE", "id")]

############ Section 3.2.1: Calculation of ICCs ####################### 
# prepare a results df to be filled
change_life_lme_icc_res <- data.frame(matrix(ncol = 13, nrow = 6))
colnames(change_life_lme_icc_res) <- c("ROI", paste0(rep(c("id_var", "segmentation_var", "resdiual_var", "ICC"), each = 3), 
                                                            c("", "_lower_95_CI", "_upper_95_CI")))
change_life_lme_icc_res$ROI <- c("TGMV", "TCV", "LVV",  "HCV", "AV", "WMHV_norm")

# formula for bootstrapped CIs for ICC
icc_ci <- function(fit) {
  variances <- as.data.frame(VarCorr(fit))
  id_var = variances[variances$grp == "id", "vcov"]
  seg_var = variances[variances$grp == "segmentation", "vcov"] 
  residual_var = variances[variances$grp == "Residual", "vcov"]
  icc = id_var / (id_var + seg_var + residual_var)
  c(id_var, seg_var, residual_var, icc)
}

# calculate LMEs
# iterate over non-wmhv_outcomes
for (outcome in c("TGMV", "TCV", "LVV",  "HCV", "AV")) {
  formula <- as.formula(paste0(outcome, " ~ (1|id) + (1|segmentation)"))
  # estimate model and bootstrap confidence intervals
  print(paste0("now doing ", outcome))
  res <- lmer(formula = formula, data = change_life_adj, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
  variances <- as.data.frame(VarCorr(res))
  lme_icc <- (variances[variances$grp == "id", "vcov"])/
    (variances[variances$grp == "id", "vcov"] + 
       variances[variances$grp == "segmentation", "vcov"] + variances[variances$grp == "Residual", "vcov"])
  boot_icc <- bootMer(res, FUN  = icc_ci, nsim = 1000, type = "parametric")
  boot_res <- apply(boot_icc$t, 2, quantile, probs = c(0.025, 0.975))
  # store the results in the appropriate columns (unadjusted/adjusted data)
  change_life_lme_icc_res[change_life_lme_icc_res$ROI == outcome, c("id_var", "segmentation_var", "resdiual_var")] <- 
    c(variances[variances$grp == "id", "vcov"],  
      variances[variances$grp == "segmentation", "vcov"], variances[variances$grp == "Residual", "vcov"])
  change_life_lme_icc_res[change_life_lme_icc_res$ROI == outcome, "ICC"] <- lme_icc
  change_life_lme_icc_res[change_life_lme_icc_res$ROI == outcome, paste0(rep(c("id_var", "segmentation_var", "resdiual_var", "ICC"),
                                                                            each = 2), c("_lower_95_CI", "_upper_95_CI"))] <- as.vector(boot_res)
}
# now wmhv

res <- lmer(formula = dWMHV ~ (1|id) + (1|segmentation), data = wmhv_change_life_adj, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
variances <- as.data.frame(VarCorr(res))
lme_icc <- (variances[variances$grp == "id", "vcov"])/
  (variances[variances$grp == "id", "vcov"] + 
     variances[variances$grp == "segmentation", "vcov"] + variances[variances$grp == "Residual", "vcov"])
boot_icc <- bootMer(res, FUN  = icc_ci, nsim = 1000, type = "parametric")
boot_res <- apply(boot_icc$t, 2, quantile, probs = c(0.025, 0.975))
# store the results in the appropriate columns (unadjusted/adjusted data)
change_life_lme_icc_res[change_life_lme_icc_res$ROI == "WMHV_norm", c("id_var", "segmentation_var", "resdiual_var")] <- 
  c(variances[variances$grp == "id", "vcov"],  
    variances[variances$grp == "segmentation", "vcov"], variances[variances$grp == "Residual", "vcov"])
change_life_lme_icc_res[change_life_lme_icc_res$ROI == "WMHV_norm", "ICC"] <- lme_icc
change_life_lme_icc_res[change_life_lme_icc_res$ROI == "WMHV_norm", paste0(rep(c("id_var", "segmentation_var", "resdiual_var", "ICC"),
                                                                               each = 2), c("_lower_95_CI", "_upper_95_CI"))] <- as.vector(boot_res)




############ Section 3.3 Scanner update WMHV ################
df_update_wmhv <- read.csv("/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/scanner_update/Data/wmhv_update.csv")
df_update_wmhv_adj <- read.csv("/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/scanner_update/Data/wmhv_update_adj.csv")
df_update_wmhv$WMHV_adj <- df_update_wmhv_adj$X..WMHV_norm
df_update_wmhv_wide <- df_update_wmhv[,c("id", "SITE", "WMHV_norm", "WMHV_adj")] %>%
  pivot_wider(
    names_from = SITE,
    values_from = c(WMHV_norm, WMHV_adj)
  )
############ Section 3.3.1: Calculation of ICCs ####################### 
set.seed(1848)
# add row to update_lme_icc_res
update_lme_icc_res[nrow(update_lme_icc_res) + 1, "ROI"] <- "WMHV"

# formula for bootstrapped CIs for ICC
icc_ci <- function(fit) {
  variances <- as.data.frame(VarCorr(fit))
  id_var = variances[variances$grp == "id", "vcov"]
  scan_var = variances[variances$grp == "scanner", "vcov"]
  residual_var = variances[variances$grp == "Residual", "vcov"]
  icc = id_var / (id_var + scan_var + residual_var)
  c(id_var, scan_var, residual_var, icc)
}

# calculate LMEs
# iterate over outcomes
for (outcome in c("WMHV_norm", "WMHV_adj")) {
  # get suffix
  suffix <- ""
  if (outcome == "WMHV_adj"){
    suffix <- "_adj"
  }
  formula <- as.formula(paste0(outcome, " ~ (1|id) + (1|scanner)"))
  # estimate model and bootstrap confidence intervals
  res <- lmer(formula = formula, data = df_update_wmhv, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
  variances <- as.data.frame(VarCorr(res))
  lme_icc <- (variances[variances$grp == "id", "vcov"])/
    (variances[variances$grp == "scanner", "vcov"] + variances[variances$grp == "id", "vcov"] + 
       variances[variances$grp == "Residual", "vcov"])
  boot_icc <- bootMer(res, FUN  = icc_ci, nsim = 1000, type = "parametric")
  boot_res <- apply(boot_icc$t, 2, quantile, probs = c(0.025, 0.975))
  # store the results in the appropriate columns (unadjusted/adjusted data)
  update_lme_icc_res[update_lme_icc_res$ROI == "WMHV", paste0(c("id_var", "scanner_var", "resdiual_var"), suffix)] <- 
    c(variances[variances$grp == "id", "vcov"], variances[variances$grp == "scanner", "vcov"], 
      variances[variances$grp == "Residual", "vcov"])
  update_lme_icc_res[update_lme_icc_res$ROI == "WMHV", paste0("ICC", suffix)] <- lme_icc
  update_lme_icc_res[update_lme_icc_res$ROI == "WMHV", paste0(paste0(rep(c("id_var", "scanner_var", "resdiual_var", "ICC"),
                                                                          each = 2), c("_lower_95_CI", "_upper_95_CI")), suffix)] <- as.vector(boot_res)
}

############ Section 3.3.2: Calculation of PDs #######################

# iterate over ROIs
for (outcome in c("WMHV_norm", "WMHV_adj")) {
  suffix <- ""
  if (outcome == "WMHV_adj"){
    suffix <- "_adj"
  }  
  # get names of relevant cols
  roi_cols <- grep(outcome, colnames(df_update_wmhv_wide), value = TRUE)
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
        col1 <- df_update_wmhv_wide[[roi_cols[i]]]
        col2 <- df_update_wmhv_wide[[roi_cols[j]]]
        # calculate a vector of PDs
        # a negative value indicates that the values of the batch given by the rowname is smaller 
        pd <- (2 * (col1 - col2)) / (abs(col1) + abs(col2))
        # calculate the mean PD between the two batches (NAs exist because 5 participants were not scanned on both scanners)
        mean_pd_mat[i, j] <- mean(pd, na.rm = TRUE)
      }
    }
  }
  # mirror the inverse to lower triangle
  mean_pd_mat[lower.tri(mean_pd_mat)] <- t(mean_pd_mat)[lower.tri(mean_pd_mat)] * (-1)
  # turn into percentage values
  mean_pd_mat <- mean_pd_mat*100
  # save variables to list
  update_pd_matrices[[paste0(outcome, "_mean_pd_mat",suffix)]] <- mean_pd_mat
}


############ Section 3.3.3: Calculation of pairwise t-tests #######################

# iterate over ROIs
for (outcome in c("WMHV_norm", "WMHV_adj")) {
  suffix <- ""
  if (outcome == "WMHV_adj"){
    suffix <- "_adj"
  }  
  # perform pairwise comparisons 
  formula <- as.formula(paste0(outcome, " ~ SITE"))
  pwc <- df_update_wmhv %>%
    group_by(id) %>%
    filter(n() == 2) %>%
    ungroup() %>%
    pairwise_t_test(
      formula = formula, paired = TRUE, id = "id"
    )
  update_pwcs[[paste0(outcome, "_pwc", suffix)]] <- pwc
}
############ Section 3.3.4 Bland-Altman plots ################
for (outcome in c("WMHV_norm", "WMHV_adj")) {
  cols <- grep(outcome, colnames(df_update_wmhv_wide), value = T)
  # add a vertical line to show the size of 1 SD
  # determine the x-coordinate (should be left of the smallest value)
  x_coord = min(rowMeans(df_update_wmhv_wide[,cols]), na.rm = T)
  # determine the y-coordinates 
  sd_val <- sd(unlist(df_update_wmhv_wide[cols]), na.rm = TRUE)
  y_coord_a <- sd_val/2 
  y_coord_b <- -sd_val/2
  plt <- bland.altman.plot(group1 = df_update_wmhv_wide[[cols[1]]], group2 = df_update_wmhv_wide[[cols[2]]], 
                           graph.sys = "ggplot2", conf.int=.95, pch=19)
  plt <- plt + labs(title = paste0(cols[1], " vs. ", cols[2])) + theme(plot.title = element_text(size = 12), base_family = "garamond")
  # add line to idicate size of SD
  plt <- plt + 
    geom_segment(
      x = x_coord, xend = x_coord,
      y = y_coord_a, yend = y_coord_b,
      colour = "blue",
      linewidth = 1
    )
  update_ba_plots[[paste0(outcome)]] <- plt
}
unharm_plot <- update_ba_plots[["WMHV_norm"]] + labs(title = paste0("WMHV SKYRA vs. WMHV VERIO unharmonized")) 
nh_plot <- update_ba_plots[["WMHV_adj"]] + labs(title = paste0("WMHV SKYRA vs. WMHV VERIO harmonized")) 
plots <- wrap_plots(list(unharm_plot, nh_plot), nrow = 2, )
update_ba_plot_comparisons[["WMHV"]] <- plots 
############ save outcomes #################
outcomes <- list(update_lme_icc_res = update_lme_icc_res, update_pd_matrices = update_pd_matrices, update_anova_res = update_anova_res, 
                 update_pwcs = update_pwcs, update_iced_res = update_iced_res, life_lme_icc_res = life_lme_icc_res, 
                 life_pd_matrices = life_pd_matrices, life_pwcs = life_pwcs, tp_life_lme_icc_res = tp_life_lme_icc_res, 
                 tp_life_pd_matrices = tp_life_pd_matrices, tp_life_pwcs = tp_life_pwcs, change_life_lme_icc_res = change_life_lme_icc_res)
ba_plots <- list(update_ba_plots = update_ba_plots, life_ba_plots = life_ba_plots, update_ba_plot_comparisons = update_ba_plot_comparisons, wrapped_life_ba_plots = wrapped_life_ba_plots,
                 tp_life_ba_plots = tp_life_ba_plots, tp_wrapped_life_ba_plots = tp_wrapped_life_ba_plots)

save(outcomes, file = "/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/Results/outcomes.RData")
save(ba_plots, file = "/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/Results/ba_plots.RData")
