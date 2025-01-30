# %% Imports
import numpy as np
import os
import pandas as pd
import pingouin as pg
import seaborn as sns
import matplotlib.pyplot as plt
import h5py
import yaml
from tqdm import tqdm
from utils.helperFunctions import load_data, extract_seq_score, extract_H

plt.style.use("plot_style.mplstyle")

# %% Load parameters
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

np.random.seed(params["seed"])

#%% Assembly analysis
results_dir = params['path_to_output']+"/PVC"
resultsList = os.listdir(results_dir)
data_list = []

for file_name in tqdm(resultsList):
    if (
        file_name.startswith("PVC_")
        and file_name.endswith(".h5")
        # and "pv1254" not in file_name # Exclude pv1254
    ):
        h5_file = h5py.File(os.path.join(results_dir, file_name))

        data_list.append(
            {
                "mouse": h5_file["mouse"][()].decode("utf-8"),
                "condition": h5_file["condition"][()].decode("utf-8"),
                "state_ref": h5_file["state_ref"][()].decode("utf-8"),
                "state_pred": h5_file["state_pred"][()].decode("utf-8"),
                "PVC":h5_file["PVC"][()]
            }
        )

        # Close files
        h5_file.close()

df = pd.DataFrame(data_list)

#%% Melt data
# melted_df = df.melt(
#     id_vars = ['mouse', 'condition'],
#     value_vars = ['numPreplayedAssemblies', 'numReplayedAssemblies'],
#     var_name = 'Type',
#     value_name = 'Num. assemblies'
# )

plt.figure(figsize=(1,.75))
sns.heatmap(data=df.query("condition=='LTD1'").pivot_table(
    index='state_ref',
    columns='state_pred', 
    values='PVC'
    ),
    order=['REMpre','wake','REMpost'],
    cmap='magma',
    vmin=0,
    vmax=.5,
    rasterized=True,
    cbar_kws={'label':'PV correlation'}
)
#plt.xticks([])
#plt.yticks([])
plt.xlabel('')
plt.ylabel('')
plt.savefig("../../output_REM/PVC_heatmap.pdf")


#%% Plot results
plt.figure(figsize=(.5,.5))
sns.barplot(
    data=df.query("condition=='LTD1' and state_ref=='REMpre'"),
    y='PVC',
    x='state_pred',
    order=['REMpre','wake','REMpost'],
    errorbar='se',
    capsize=.2,
    legend=False
)
sns.stripplot(
    data=df.query("condition=='LTD1' and state_ref=='REMpre'"),
    y='PVC',
    x='state_pred',
    order=['REMpre','wake','REMpost'],
    size=2,
    legend=False
)

plt.title('Preplayed events')
plt.xlabel('')
plt.ylabel('PV\ncorrelation')
#plt.xticks([0,1,2],['REMpre','RUN','REMpost'], rotation=90)
plt.xlim(.5,2.5)
plt.xticks([1,2],['RUN','REMpost'], rotation=90)
plt.ylim(0,.5)
#plt.xlabel('Number of\npre-existing assemblies')
#plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)

plt.savefig("../../output_REM/PVC_REMpre_ref.pdf")

#%%
plt.figure(figsize=(.5,.5))
sns.barplot(
    data=df.query("condition=='LTD1' and state_ref=='wake'"),
    y='PVC',
    x='state_pred',
    order=['wake', 'REMpre','REMpost'],
    errorbar='se',
    capsize=.2,
    legend=False
)
sns.stripplot(
    data=df.query("condition=='LTD1' and state_ref=='wake'"),
    y='PVC',
    x='state_pred',
    order=['wake', 'REMpre','REMpost'],
    size=2,
    legend=False
)

plt.title('RUN events')
plt.xlabel('')
plt.ylabel('PV\ncorrelation')
plt.xticks([1,2],['REMpre','REMpost'], rotation=90)
plt.xlim(.5,2.5)
plt.ylim(0,.5)

#plt.xlabel('Number of\npre-existing assemblies')
#plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)

plt.savefig("../../output_REM/PVC_RUN_ref.pdf")

#%% DESCRIPTIVES
mean_wake = df.query("condition=='LTD1' and state_ref=='REMpre' and state_pred=='wake'")['PVC'].mean()
SEM_wake = df.query("condition=='LTD1' and state_ref=='REMpre' and state_pred=='wake'")['PVC'].sem()
print(f'{mean_wake} +/- {SEM_wake}')

mean_REMpost = df.query("condition=='LTD1' and state_ref=='REMpre' and state_pred=='REMpost'")['PVC'].mean()
SEM_REMpost = df.query("condition=='LTD1' and state_ref=='REMpre' and state_pred=='REMpost'")['PVC'].sem()
print(f'{mean_REMpost} +/- {SEM_REMpost}')

#%% STATS
pg.ttest(
    x=df.query("condition=='LTD1' and state_ref=='REMpre' and state_pred=='wake'")['PVC'],
    y=df.query("condition=='LTD1' and state_ref=='REMpre' and state_pred=='REMpost'")['PVC'],
    paired=True)

#%% Same for wake ref
mean_REMpre = df.query("condition=='LTD1' and state_ref=='wake' and state_pred=='REMpre'")['PVC'].mean()
SEM_REMpre = df.query("condition=='LTD1' and state_ref=='wake' and state_pred=='REMpre'")['PVC'].sem()
print(f'{mean_wake} +/- {SEM_wake}')

mean_REMpost = df.query("condition=='LTD1' and state_ref=='wake' and state_pred=='REMpost'")['PVC'].mean()
SEM_REMpost = df.query("condition=='LTD1' and state_ref=='wake' and state_pred=='REMpost'")['PVC'].sem()
print(f'{mean_REMpost} +/- {SEM_REMpost}')

#%% STATS
pg.ttest(
    x=df.query("condition=='LTD1' and state_ref=='wake' and state_pred=='REMpre'")['PVC'],
    y=df.query("condition=='LTD1' and state_ref=='wake' and state_pred=='REMpost'")['PVC'],
    paired=True,
    )