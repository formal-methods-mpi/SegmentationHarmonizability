from neuroHarmonize import harmonizationLearn, harmonizationApply
import pandas as pd
import numpy as np
import os
# set working directory
os.chdir("/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/scanner_update/Data/")
# iterate over the summary files of the scanner-update dataset
for name in ["", "_no_correction"]:
    data = pd.read_csv('df_update' + name + '.csv')
    volumes = data.loc[:, ["TIV", "TGMV", "TWMV", "TCV", "LVV", "HCV", "AV"]]
    volumes = np.array(volumes)
    covars = data.loc[:, ["SITE", "age", "gender"]]
    # run harmonization and store the adjusted data
    model, data_adj = harmonizationLearn(volumes, covars, smooth_terms=['age'])
    np.savetxt("df_update" + name + '_adj.csv', data_adj, delimiter=",", fmt="%f", header="TIV, TGMV, TWMV, TCV, LVV, HCV, AV")

# switch working directory
os.chdir("/data/pt_life/ResearchProjects/LLammer/intergeneration/segmentation_harmonization/life/Data/")
# iterate over the summary files for the LIFE-Adult study
for name in ["", "_wmhv" , "_no_correction", "_cross"]:
    # read bl data
    bl_data = pd.read_csv('df_life' + name + '_bl.csv')
    if name != "_wmhv":
        bl_volumes = bl_data.loc[:,['TGMV', 'TCV', 'LVV', 'HCV', 'AV']]
    else:
        bl_volumes = bl_data.loc[:, ['WMHV_asinh']]
    # transform to numpy array
    bl_volumes = np.array(bl_volumes)
    # get bl covars
    bl_covars = bl_data.loc[:,["SITE", "AGE", "GENDER"]]
    # read fu volumes
    fu_data = pd.read_csv('df_life' + name + '_fu.csv')
    if name != "_wmhv":
        fu_volumes = fu_data.loc[:,['TGMV', 'TCV', 'LVV', 'HCV', 'AV']]
    else:
        fu_volumes = fu_data.loc[:, ['WMHV_asinh']]
    # transform to numpy array
    fu_volumes = np.array(fu_volumes)
    # get fu covars
    fu_covars = fu_data.loc[:,["SITE", "AGE", "GENDER"]]
    # derive minimum and maximum age across all covariate dfs for smooth term bounds
    age_min = bl_covars["AGE"].min()
    age_max = fu_covars["AGE"].max()
    # run harmonization and store the adjusted data
    # don't use empirical Bayes for WMHV
    if name == "_wmhv":
        bl_model, bl_data_adj = harmonizationLearn(bl_volumes, bl_covars, smooth_terms=['AGE'], eb= False, smooth_term_bounds=(age_min, age_max))
        np.savetxt("df_life" + name + '_bl_adj.csv', bl_data_adj, delimiter=",", fmt="%f", header = "WMHV_asinh")
        # adjust values using bl model
        fu_data_adj = harmonizationApply(fu_volumes, fu_covars, bl_model)
        np.savetxt("df_life" + name + '_fu_adj.csv', fu_data_adj, delimiter=",", fmt="%f", header = "WMHV_asinh")
    else:
        bl_model, bl_data_adj = harmonizationLearn(bl_volumes, bl_covars, smooth_terms=['AGE'], smooth_term_bounds=(age_min, age_max))
        np.savetxt("df_life" + name + '_bl_adj.csv', bl_data_adj, delimiter=",", fmt="%f", header="TGMV, TCV, LVV, HCV, AV")
        # adjust values using bl model
        fu_data_adj = harmonizationApply(fu_volumes, fu_covars, bl_model)
        np.savetxt("df_life" + name + '_fu_adj.csv', fu_data_adj, delimiter=",", fmt="%f", header="TGMV, TCV, LVV, HCV, AV")

## harmonize change scores immediately for exploratory analysis
for name in ["_changes", "_wmhv_changes"]:
    data = pd.read_csv('df_life' + name + '.csv')
    if name != "_wmhv_changes":
        volumes = data.loc[:, ['dTGMV', 'dTCV', 'dLVV', 'dHCV', 'dAV']]
    else:
        volumes = data.loc[:, ['dWMHV']]
    # transform to numpy array and drop empty first col
    volumes = np.array(volumes)
    covars = data.loc[:,["SITE", "AGE", "GENDER"]]
    # run harmonization and store the adjusted data
    if name == "_wmhv_changes":
        model, data_adj = harmonizationLearn(volumes, covars, smooth_terms=['AGE'], eb=False)
        np.savetxt("df_life" + name + '_adj.csv', data_adj, delimiter=",", fmt="%f", header="dWMHV")
    else:
        model, data_adj = harmonizationLearn(volumes, covars, smooth_terms=['AGE'])
        np.savetxt("df_life" + name + '_adj.csv', data_adj, delimiter=",", fmt="%f", header="TGMV, TCV, LVV, HCV, AV")