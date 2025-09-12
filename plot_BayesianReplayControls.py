# %% Imports
import numpy as np
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import h5py
import yaml
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
        data_list.append(
                {
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                    "Shuffle": "Position",
                    "Decoding": "Actual",
                    "replayEvents": sum(~np.isnan(h5_file['replayLocs_P'])), # num of actual events
                    "replayScores": np.nanmean(h5_file['replayScore_P']) 
                }
            )
        
        data_list.append(
                {
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                    "Shuffle": "Position",
                    "Decoding": "Control",
                    "replayEvents": sum(~np.isnan(h5_file['control_replayLocs_P'])), # num of actual events
                    "replayScores": np.nanmean(h5_file['control_replayScore_P'])
                }
            )
        
        data_list.append( 
                {
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                    "Shuffle": "Time",
                    "Decoding": "Actual",
                    "replayEvents": sum(~np.isnan(h5_file['replayLocs_T'])), # num of actual events
                    "replayScores": np.nanmean(h5_file['replayScore_T'])
                }
            )
        
        data_list.append(
                {
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                    "Shuffle": "Time",
                    "Decoding": "Control",
                    "replayEvents": sum(~np.isnan(h5_file['control_replayLocs_T'])), # num of actual events
                    "replayScores": np.nanmean(h5_file['control_replayScore_T'])
                }
            )

        data_list.append( 
                {
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                    "Shuffle": "Position_Time",
                    "Decoding": "Actual",
                    "replayEvents": sum(~np.isnan(h5_file['replayLocs_PT'])), # num of actual events
                    "replayScores": np.nanmean(h5_file['replayScore_PT'])
                }
            )
        
        data_list.append(
                {
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                    "Shuffle": "Position_Time",
                    "Decoding": "Control",
                    "replayEvents": sum(~np.isnan(h5_file['control_replayLocs_PT'])), # num of actual events
                    "replayScores": np.nanmean(h5_file['control_replayScore_PT'])
                }
            )
        
        # Close files
        h5_file.close()

df = pd.DataFrame(data_list)

# %% Plot results
plt.figure(figsize=(2,1))
sns.barplot(
    data=df.query("Type=='replay' and condition=='LTD1'"),
    x='replayEvents',
    y="Shuffle",
    hue="Decoding",
    legend=True,
    errorbar='se'
)

sns.stripplot(
    data=df.query("Type=='replay' and condition=='LTD1'"),
    x='replayEvents',
    y="Shuffle",
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
    data=df.query("Type=='replay' and condition=='LTD1'"),
    x='replayScores',
    y="Shuffle",
    hue="Decoding",
    legend=False,
    errorbar='se'
)

sns.stripplot(
    data=df.query("Type=='replay' and condition=='LTD1'"),
    x='replayScores',
    y="Shuffle",
    hue="Decoding",
    legend=False,
    dodge=True,
    size=1
)

plt.title('REM post')
plt.xlabel('R$^{2}$')
#plt.ylabel('N')
plt.savefig("../../output_REM/bayesian_control_scores.pdf")
# %%
