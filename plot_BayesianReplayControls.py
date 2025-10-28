# %% Imports
import numpy as np
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import h5py
import yaml
import pingouin as pg
from tqdm import tqdm

plt.style.use("plot_style.mplstyle")

# %% Load parameters
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

np.random.seed(params["seed"])

# %% Load all sessions
results_dir = params['path_to_output']+"/controlled_bayesian_replay"
resultsList = os.listdir(results_dir)

#%%
data_list = []

for file_name in tqdm(resultsList):
    if (
        file_name.startswith("controlled_bayesian_replay_")
        and file_name.endswith(".h5")
        and "pv1254" not in file_name # Exclude pv1254
    ):
        h5_file = h5py.File(os.path.join(results_dir, file_name))
        # No FDR
        data_list.append(
                {
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                    "Shuffle": "Position",
                    "Decoding": "Actual",
                    "FDR": 'None',
                    "replayEvents": sum(~np.isnan(h5_file['replayLocs_P'][h5_file['pvalue_P'][()]<0.05])), # num of actual events
                    "replayScores": np.nanmean(h5_file['replayScore_P'][h5_file['pvalue_P'][()]<0.05]),
                    "replayJumpiness": np.nanmean(h5_file['replayJumpiness_P'][h5_file['pvalue_P'][()]<0.05]),
                }
            )
        
        # Bonferroni
        bonferroni_alpha = 0.05/len(h5_file['pvalue_P'][()])
        data_list.append(
                {
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                    "Shuffle": "Position",
                    "Decoding": "Actual",
                    "FDR": 'Bonferroni',
                    "replayEvents": sum(~np.isnan(h5_file['replayLocs_P'][h5_file['pvalue_P'][()]<bonferroni_alpha])), # num of actual events
                    "replayScores": np.nanmean(h5_file['replayScore_P'][h5_file['pvalue_P'][()]<bonferroni_alpha]),
                    "replayJumpiness": np.nanmean(h5_file['replayJumpiness_P'][h5_file['pvalue_P'][()]<bonferroni_alpha]),
                }
            )
        
        # No FDR
        data_list.append(
                {
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                    "Shuffle": "Position",
                    "Decoding": "Control",
                    'FDR':'None',
                    "replayEvents": sum(~np.isnan(h5_file['control_replayLocs_P'][h5_file['control_pvalue_P'][()]<0.05])), # num of actual events
                    "replayScores": np.nanmean(h5_file['control_replayScore_P'][h5_file['control_pvalue_P'][()]<0.05]),
                    "replayJumpiness": np.nanmean(h5_file['control_replayJumpiness_P'][h5_file['control_pvalue_P'][()]<0.05]),
                }
            )
        
        # Bonferroni
        bonferroni_alpha = 0.05/len(h5_file['control_pvalue_P'][()])
        data_list.append(
                {
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                    "Shuffle": "Position",
                    "Decoding": "Control",
                    "FDR": 'Bonferroni',
                    "replayEvents": sum(~np.isnan(h5_file['control_replayLocs_P'][h5_file['control_pvalue_P'][()]<bonferroni_alpha])), # num of actual events
                    "replayScores": np.nanmean(h5_file['control_replayScore_P'][h5_file['control_pvalue_P'][()]<bonferroni_alpha]),
                    "replayJumpiness": np.nanmean(h5_file['control_replayJumpiness_P'][h5_file['control_pvalue_P'][()]<bonferroni_alpha]),
                }
            )

        # No FDR
        data_list.append( 
                {
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                    "Shuffle": "Time",
                    "Decoding": "Actual",
                    'FDR':'None',
                    "replayEvents": sum(~np.isnan(h5_file['replayLocs_T'][h5_file['pvalue_T'][()]<0.05])), # num of actual events
                    "replayScores": np.nanmean(h5_file['replayScore_T'][h5_file['pvalue_T'][()]<0.05]),
                    "replayJumpiness": np.nanmean(h5_file['replayJumpiness_T'][h5_file['pvalue_T'][()]<0.05]),
                }
            )
        
        # Bonferroni
        bonferroni_alpha = 0.05/len(h5_file['pvalue_T'][()])
        data_list.append( 
                {
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                    "Shuffle": "Time",
                    "Decoding": "Actual",
                    'FDR':'Bonferroni',
                    "replayEvents": sum(~np.isnan(h5_file['replayLocs_T'][h5_file['pvalue_T'][()]<bonferroni_alpha])), # num of actual events
                    "replayScores": np.nanmean(h5_file['replayScore_T'][h5_file['pvalue_T'][()]<bonferroni_alpha]),
                    "replayJumpiness": np.nanmean(h5_file['replayJumpiness_T'][h5_file['pvalue_T'][()]<bonferroni_alpha]),
                }
            )
        
        # No FDR
        data_list.append(
                {
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                    "Shuffle": "Time",
                    "Decoding": "Control",
                    'FDR':'None',
                    "replayEvents": sum(~np.isnan(h5_file['control_replayLocs_T'][h5_file['control_pvalue_T'][()]<0.05])), # num of actual events
                    "replayScores": np.nanmean(h5_file['control_replayScore_T'][h5_file['control_pvalue_T'][()]<0.05]),
                    "replayJumpiness": np.nanmean(h5_file['control_replayJumpiness_T'][h5_file['control_pvalue_T'][()]<0.05]),
                }
            )
        
        # Bonferroni
        bonferroni_alpha = 0.05/len(h5_file['control_pvalue_T'][()])
        data_list.append(
                {
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                    "Shuffle": "Time",
                    "Decoding": "Control",
                    'FDR':'Bonferroni',
                    "replayEvents": sum(~np.isnan(h5_file['control_replayLocs_T'][h5_file['control_pvalue_T'][()]<bonferroni_alpha])), # num of actual events
                    "replayScores": np.nanmean(h5_file['control_replayScore_T'][h5_file['control_pvalue_T'][()]<bonferroni_alpha]),
                    "replayJumpiness": np.nanmean(h5_file['control_replayJumpiness_T'][h5_file['control_pvalue_T'][()]<bonferroni_alpha]),
                }
            )

        # No FDR
        data_list.append( 
                {
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                    "Shuffle": "Position_Time",
                    "Decoding": "Actual",
                    'FDR':'None',
                    "replayEvents": sum(~np.isnan(h5_file['replayLocs_PT'][h5_file['pvalue_PT'][()]<0.05])), # num of actual events
                    "replayScores": np.nanmean(h5_file['replayScore_PT'][h5_file['pvalue_PT'][()]<0.05]),
                    "replayJumpiness": np.nanmean(h5_file['replayJumpiness_PT'][h5_file['pvalue_PT'][()]<0.05]),
                }
            )
        
        # Bonferroni
        bonferroni_alpha = 0.05/len(h5_file['pvalue_PT'][()])
        data_list.append( 
                {
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                    "Shuffle": "Position_Time",
                    "Decoding": "Actual",
                    'FDR':'Bonferroni',
                    "replayEvents": sum(~np.isnan(h5_file['replayLocs_PT'][h5_file['pvalue_PT'][()]<bonferroni_alpha])), # num of actual events
                    "replayScores": np.nanmean(h5_file['replayScore_PT'][h5_file['pvalue_PT'][()]<bonferroni_alpha]),
                    "replayJumpiness": np.nanmean(h5_file['replayJumpiness_PT'][h5_file['pvalue_PT'][()]<bonferroni_alpha]),
                }
            )
        
        # No FDR
        data_list.append(
                {
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                    "Shuffle": "Position_Time",
                    "Decoding": "Control",
                    'FDR':'None',
                    "replayEvents": sum(~np.isnan(h5_file['control_replayLocs_PT'][h5_file['control_pvalue_PT'][()]<0.05])), # num of actual events
                    "replayScores": np.nanmean(h5_file['control_replayScore_PT'][h5_file['control_pvalue_PT'][()]<0.05]),
                    "replayJumpiness": np.nanmean(h5_file['control_replayJumpiness_PT'][h5_file['control_pvalue_PT'][()]<0.05]),
                }
            )
        
        # Bonferroni
        bonferroni_alpha = 0.05/len(h5_file['control_pvalue_PT'][()])
        data_list.append(
                {
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                    "Shuffle": "Position_Time",
                    "Decoding": "Control",
                    'FDR':'Bonferroni',
                    "replayEvents": sum(~np.isnan(h5_file['control_replayLocs_PT'][h5_file['control_pvalue_PT'][()]<bonferroni_alpha])), # num of actual events
                    "replayScores": np.nanmean(h5_file['control_replayScore_PT'][h5_file['control_pvalue_PT'][()]<bonferroni_alpha]),
                    "replayJumpiness": np.nanmean(h5_file['control_replayJumpiness_PT'][h5_file['control_pvalue_PT'][()]<bonferroni_alpha]),
                }
            )
        
        # Close files
        h5_file.close()

df = pd.DataFrame(data_list)

# %% Plot results
plt.figure(figsize=(2,1))
sns.barplot(
    data=df.query("Type=='replay' and condition=='LTD1' and Shuffle=='Position'"),
    x='replayEvents',
    y='FDR',
    hue="Decoding",
    legend=True,
    errorbar='se'
)

sns.stripplot(
    data=df.query("Type=='replay' and condition=='LTD1' and Shuffle=='Position'"),
    x='replayEvents',
    y="FDR",
    hue="Decoding",
    legend=False,
    dodge=True,
    size=1
)

plt.title('REM post')
plt.xlabel('N events')
#plt.ylabel('N')
plt.savefig("../../output_REM/bayesian_control.pdf")

# %% Plot results
plt.figure(figsize=(2,1))
sns.barplot(
    data=df.query("Type=='replay' and condition=='LTD1' and Shuffle=='Position'"),
    x='replayScores',
    y="FDR",
    hue="Decoding",
    legend=False,
    errorbar='se'
)

sns.stripplot(
    data=df.query("Type=='replay' and condition=='LTD1' and Shuffle=='Position'"),
    x='replayScores',
    y="FDR",
    hue="Decoding",
    legend=False,
    dodge=True,
    size=1
)

plt.title('REM post')
plt.xlabel('R$^{2}$')
#plt.ylabel('N')
plt.savefig("../../output_REM/bayesian_control_scores.pdf")
# %% STATS
pg.anova(
    data=df.query("Type=='replay' and condition=='LTD1'"),
    dv='replayScores',
    between=['Shuffle','Decoding', 'FDR']
).round(4)
# %%
