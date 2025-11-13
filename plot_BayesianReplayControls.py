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

# %% Plot results: Bonferroni correction for novelty
plt.figure(figsize=(2,1))
sns.barplot(
    data=df.query("Type=='replay' and Decoding=='Actual' and Shuffle=='Position_Time' and condition=='LTD1' or condition=='LTD5'"),
    x='replayEvents',
    hue="condition",
    y='FDR',
    legend=True,
    errorbar='se'
)

sns.stripplot(
    data=df.query("Type=='replay' and Decoding=='Actual' and Shuffle=='Position_Time' and condition=='LTD1' or condition=='LTD5'"),
    x='replayEvents',
    hue="condition",
    y='FDR',
    legend=False,
    dodge=True,
    size=1
)

# sns.pointplot(
#     data=df.query("Type=='replay' and FDR=='None' and Shuffle=='Time'"),
#     y='replayScores',
#     x="Decoding",
#     units='mouse',
#     #estimator=None,
#     #y='Shuffle',
#     legend=False,
#     dodge=True,
#     #size=1
# )

plt.title('REM post')
plt.xlabel('N events')
#plt.ylabel('N')
#plt.xlim(.4,.8)
#plt.savefig("../../output_REM/bayesianReplay_events_shuffleVScontrol.pdf")

#%% STATS
pg.anova(
    data=df.query("Type=='replay' and Decoding=='Actual' and Shuffle=='Position_Time' and condition=='LTD1' or condition=='LTD5'"),
    dv='replayEvents',
    between=['condition', 'FDR']
).round(4)



# %% Plot results: within vs across replay
plt.figure(figsize=(2,1))
sns.barplot(
    data=df.query("Type=='replay' and FDR=='None'"),
    x='replayScores',
    hue="Decoding",
    y='Shuffle',
    legend=False,
    errorbar='se'
)

sns.stripplot(
    data=df.query("Type=='replay' and FDR=='None'"),
    x='replayScores',
    hue="Decoding",
    y='Shuffle',
    legend=False,
    dodge=True,
    size=1
)

plt.title('REM post')
plt.xlabel('R$^{2}$')
#plt.ylabel('N')
plt.savefig("../../output_REM/bayesianReplay_R2_shuffleVScontrol.pdf")
# %% STATS
pg.anova(
    data=df.query("Type=='replay' and FDR=='None' and replayScores>0.5"),
    dv='replayEvents',
    between=['Shuffle','Decoding']
).round(4)

# %% Plot results: LTD1 vs LTD5?
plt.figure(figsize=(2,1))
sns.barplot(
    data=df.query("Type=='replay' and FDR=='None'"),
    x='replayScores',
    hue="Decoding",
    y='Shuffle',
    legend=False,
    errorbar='se'
)

sns.stripplot(
    data=df.query("Type=='replay' and FDR=='None'"),
    x='replayScores',
    hue="Decoding",
    y='Shuffle',
    legend=False,
    dodge=True,
    size=1
)

plt.title('REM post')
plt.xlabel('R$^{2}$')
#plt.ylabel('N')
plt.savefig("../../output_REM/bayesianReplay_R2_shuffleVScontrol.pdf")


# %% Plot results: within vs across replay
plt.figure(figsize=(2,1))
sns.barplot(
    data=df.query("Type=='replay' and FDR=='None' and Shuffle=='Position_Time'"),
    x='replayEvents',
    y="Decoding",
    legend=True,
    errorbar='se'
)

sns.stripplot(
    data=df.query("Type=='replay' and FDR=='None' and Shuffle=='Position_Time'"),
    x='replayEvents',
    y="Decoding",
    legend=False,
    dodge=True,
    size=1
)

plt.title('REM post')
plt.xlabel('N events')
#plt.ylabel('N')
plt.savefig("../../output_REM/bayesianReplay_events_acrossVSwithin.pdf")

#%% STATS

# %% Plot results: within vs across replay
plt.figure(figsize=(2,1))
sns.barplot(
    data=df.query("Type=='replay' and FDR=='None' and Shuffle=='Position_Time'"),
    x='replayScores',
    y="Decoding",
    legend=False,
    errorbar='se'
)

sns.stripplot(
    data=df.query("Type=='replay' and FDR=='None' and Shuffle=='Position_Time'"),
    x='replayScores',
    y="Decoding",
    legend=False,
    dodge=True,
    size=1
)

plt.title('REM post')
plt.xlabel('R$^{2}$')
#plt.ylabel('N')
plt.savefig("../../output_REM/bayesianReplay_R2_acrossVSwithin.pdf")
# %% STATS
pg.anova(
    data=df.query("Type=='replay' and condition=='LTD1'"),
    dv='replayScores',
    between=['Shuffle','Decoding', 'FDR']
).round(4)

#%% STATS

# %% Plot results: noFDR vs Bonferroni
plt.figure(figsize=(2,1))
sns.barplot(
    data=df.query("Type=='replay' and Decoding=='Actual' and Shuffle=='Position_Time'"),
    x='replayEvents',
    y="FDR",
    legend=True,
    errorbar='se'
)

sns.stripplot(
    data=df.query("Type=='replay' and Decoding=='Actual' and Shuffle=='Position_Time'"),
    x='replayEvents',
    y="FDR",
    legend=False,
    dodge=True,
    size=1
)

plt.title('REM post')
plt.xlabel('N events')
#plt.ylabel('N')
plt.savefig("../../output_REM/bayesianReplay_events_BonferroniVSnoFDR.pdf")
# %% STATS

# %% Plot results: noFDR vs Bonferroni
plt.figure(figsize=(2,1))
sns.barplot(
    data=df.query("Type=='replay' and Decoding=='Actual' and Shuffle=='Position_Time'"),
    x='replayScores',
    y="FDR",
    legend=True,
    errorbar='se'
)

sns.stripplot(
    data=df.query("Type=='replay' and Decoding=='Actual' and Shuffle=='Position_Time'"),
    x='replayScores',
    y="FDR",
    legend=False,
    dodge=True,
    size=1
)

plt.title('REM post')
plt.xlabel('N events')
#plt.ylabel('N')
plt.savefig("../../output_REM/bayesianReplay_R2_BonferroniVSnoFDR.pdf")

# %% STATS

# %% Plot results: shuffle type
plt.figure(figsize=(2,1))
sns.barplot(
    data=df.query("Type=='replay' and Decoding=='Actual' and FDR=='None'"),
    x='replayEvents',
    y="Shuffle",
    legend=True,
    errorbar='se'
)

sns.stripplot(
    data=df.query("Type=='replay' and Decoding=='Actual' and FDR=='None'"),
    x='replayEvents',
    y="Shuffle",
    legend=False,
    dodge=True,
    size=1
)

plt.title('REM post')
plt.xlabel('N events')
#plt.ylabel('N')
plt.savefig("../../output_REM/bayesianReplay_events_shuffleType.pdf")
# %% STATS

# %% Plot results: shuffle type
plt.figure(figsize=(2,1))
sns.barplot(
    data=df.query("Type=='replay' and Decoding=='Actual' and FDR=='None'"),
    x='replayScores',
    y="Shuffle",
    legend=True,
    errorbar='se'
)

sns.stripplot(
    data=df.query("Type=='replay' and Decoding=='Actual' and FDR=='None'"),
    x='replayScores',
    y="Shuffle",
    legend=False,
    dodge=True,
    size=1
)

plt.title('REM post')
plt.xlabel('N events')
#plt.ylabel('N')
plt.savefig("../../output_REM/bayesianReplay_R2_shuffleType.pdf")

# %% STATS

# %% Plot results: most conservative caste
plt.figure(figsize=(2,1))
sns.barplot(
    data=df.query("Type=='replay' and Decoding=='Actual' and FDR=='Bonferroni' and Shuffle=='Time'"),
    x='replayEvents',
    y="condition",
    legend=True,
    errorbar='se'
)

sns.stripplot(
    data=df.query("Type=='replay' and Decoding=='Actual' and FDR=='Bonferroni' and Shuffle=='Time'"),
    x='replayEvents',
    y="condition",
    legend=False,
    dodge=True,
    size=1
)

plt.title('REM post')
plt.xlabel('N events')
#plt.ylabel('N')
plt.savefig("../../output_REM/bayesianReplay_events_conservative.pdf")

# %%
