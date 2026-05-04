### this script loads the R environment as used in our analyses 

# please not that the analyses were conducted in Linux and issues on other OS might occur
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
# Then load it
library(renv)
# restore from lockfile on github
renv::restore(lockfile = "https://raw.githubusercontent.com/LaurenzLammer/SegmentationHarmonizability/main/renv/renv.lock")