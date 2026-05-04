# Segmentation Harmonizability
This is the repository for "Potential and Limitations of Combat for the Harmonization of Multi-Site Neuroimaging Studies".

# Preprint
The study's preprint can be found [here]().

# Overview

- You can **reproduce** the results, tables and figures based on the travelling participant dataset by rendering the qmd document reproduction.qmd
  - the script downloads the data published on [osf](https://osf.io/n3eg8/overview) and runs all analyses
    -  thus, the script requires an internet connection to run    
  - congruence with package versions used in our analyses can be attained with renv
    - all required documents are in the renv directory of this repository
    - execute the script renv_loader.R first to restore the matching R environment
   

- The data were prepared for analysis with the script: data_assembly.R
- the harmonization was performed with: harmonization_with_longCombat.R & harmonization_with_neuroHarmonize.py
- The data were analysed with the script: analysis.R
- the manuscript was written with manuscript.qmd
- the renaming_for_publication.R script was just used to replace participant IDs with random codes

# Code
This code is written in R with occasional bits in Python.
The manuscript was written in Quarto markdown.
