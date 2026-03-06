library(tidyverse)
library(MetBrewer)
library(patchwork)


colors = met.brewer(name="Java", n=5, type="discrete")
colors <- colors[c(1,2,4,5)]
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
# for unscaled age data
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
  select(id, TIV, TGMV, TWMV, TCV, LVV, HCV, AV, SITE) %>%
  left_join(QA[,c("id", "correction.necessary", "age", "gender")], by = "id")

# load the adjusted data from the scanner update study 
df_update_adj <- read.csv("df_update_adj.csv", col.names = colnames(df_update)[2:8])
# fill in the missing data from the unadjusted df (these cols must not be included during harmonization) 
df_update_adj[,c("id", "SITE", "age", "gender")] <- df_update[,c("id", "SITE", "age", "gender")]
df_update$SITE <- gsub("_", " ", df_update$SITE)
df_update_adj$SITE <- gsub("_", " ", df_update_adj$SITE)
df_update$SITE <- factor(df_update$SITE, levels = c("SKYRA recon7", "VERIO recon7", "SKYRA recon6", "VERIO recon6", "SKYRA samseg", "VERIO samseg", "SKYRA cat", "VERIO cat"))  
df_update_adj$SITE <- factor(df_update_adj$SITE, levels = c("SKYRA recon7", "VERIO recon7", "SKYRA recon6", "VERIO recon6", "SKYRA samseg", "VERIO samseg", "SKYRA cat", "VERIO cat"))  
plots <- list()
comparisons <- list()
dfs <- list(df_update = df_update, df_update_adj = df_update_adj)
for (outcome in c("TIV", "TGMV", "TWMV", "TCV", "LVV", "HCV", "AV")){
  for (name in names(dfs)){
    df <- dfs[[name]]
    suffix <- ""
    y_axis_title <- outcome
    if (grepl("adj", name)){
      suffix <- "_adj"
      y_axis_title <- paste0(outcome, " harmonized")
    }
    ordered_ids <- df %>%
      group_by(id) %>%
      summarise(mean_val = mean(.data[[outcome]], na.rm = TRUE)) %>%
      arrange(mean_val) %>%
      pull(id)
    df$id <- factor(df$id, levels = ordered_ids)
    plt <- ggplot(data = df, aes(x = factor(id), y = .data[[outcome]], color = SITE, group = SITE, shape = SITE)) +
      geom_point(size = 2.25, alpha = 0.9, stroke = 0.8) +
      scale_shape_manual(values = c(16, 6, 16, 6, 16, 6, 16, 6)) +
      theme_classic() +
      guides(color = guide_legend(ncol = 2, override.aes = list(size = 4)), shape = guide_legend(ncol = 2)) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), 
            legend.text = element_text(size = 15), legend.title = element_text(size = 17)) +
      ylab(y_axis_title) +
      scale_color_manual(values = rep(colors, each = 2))
    if (grepl("adj", name)){
      plt <- plt + theme(axis.text.y  = element_blank(), axis.ticks.y = element_blank())
    }
    plots[[paste0(outcome, suffix)]] <- plt
  }
  y_range <- range(c(plots[[outcome]]$data[[outcome]], plots[[paste0(outcome, "_adj")]]$data[[outcome]]), na.rm = TRUE)
  plots[[outcome]] <- plots[[outcome]] + coord_cartesian(ylim = y_range)
  plots[[paste0(outcome, "_adj")]] <- plots[[paste0(outcome, "_adj")]] + coord_cartesian(ylim = y_range)
  comparison <- wrap_plots(list(plots[[outcome]], plots[[paste0(outcome, "_adj")]]), nrow = 1) + plot_layout(guides = "collect") & theme(legend.position = "right")
  comparisons[[outcome]] <- comparison
}
legend_only_plot <- ggplot() + 
  theme_void() +
  annotation_custom(cowplot::get_legend(comparisons[[1]]))
combined <- wrap_plots(comparisons[[6]], comparisons[[1]], comparisons[[2]], comparisons[[3]], comparisons[[4]], comparisons[[5]], comparisons[[7]], legend_only_plot, 
  nrow = 4, ncol = 2, guides = "collect") + plot_layout(widths = 1) &
  theme(legend.position = "none", plot.margin = margin(10, 5, 10, 5))
combined
ggsave("/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/Results/update_comparison_plot.pdf", plot = combined, width = 15, height = 8, units = "in")
