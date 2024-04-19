# %% Imports
import numpy as np
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import h5py
import yaml
from tqdm import tqdm
from utils.helperFunctions import load_data

plt.style.use("plot_style.mplstyle")

# %% Load parameters
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

np.random.seed(params["seed"])

# %% Select example recording to plot wake posteriors
mouse = 'pv1069'
condition = 'LTD1'

# %% Load data
data_wake = load_data(mouse, condition, 'wake', params)

#%% Load assembly data for that session
with h5py.File(
    os.path.join(
        params['path_to_output'],
        'neuron_selection',
        f"selected_neurons_{condition}_{mouse}.h5"
        ), 'r'
        ) as f:
    selected_neurons = f['place_cells'][()]

with h5py.File(
    os.path.join(
        params["path_to_output"],
        'assembly',
        f"assembly_{condition}_{mouse}.h5"
        ), 'r'
        ) as f:
    sig_cells=f['A_sig_cells'][()]
    weights=f['weights'][()]

#%% Assembly example plot
assembly_ID=15
plt.figure(figsize=(3,1))
plt.subplot(142)
plt.stem(
    weights[assembly_ID],
    orientation='horizontal',
    # linefmt=sig_cells[assembly_ID]
    )
plt.xlabel('Weight')
plt.ylabel('Neuron ID')

plt.subplot(122)
plt.imshow(
    data_wake['binaryData'][:,selected_neurons].T*weights[assembly_ID].T[:,None],
    aspect='auto',
    interpolation='none',
    origin='lower',
    cmap='RdBu',
    vmin=-.1,
    vmax=.1
           )
plt.xlim(340,700)
plt.plot([520,670],[0,0],linewidth=2, color='k')
plt.text(570,-40,'5 s')
plt.colorbar()
plt.axis('off')
plt.savefig("../../output_REM/example_awake_assembly.pdf")

# %% Load all sessions
results_dir = params['path_to_output']+"/assembly"
resultsList = os.listdir(results_dir)

#%%
data_list = []
data_list_num_events = []

for file_name in tqdm(resultsList):
    if (
        file_name.startswith("assembly_")
        and file_name.endswith(".h5")
        # and "pv1254" not in file_name # Exclude pv1254
    ):
        h5_file = h5py.File(os.path.join(results_dir, file_name))
        data = load_data(h5_file["mouse"][0].decode("utf-8"),
                         h5_file["condition"][0].decode("utf-8"),
                         'REMpost',
                         params
                         )
        numFrames = data['binaryData'].shape[0]
        recordingLength = numFrames/params['sampling_frequency']

        data_list_num_events.append(
            {
                "mouse": h5_file["mouse"][0].decode("utf-8"),
                "condition": h5_file["condition"][0].decode("utf-8"),
                "numSigEvents": np.sum(h5_file['post_rem_A_sig'][()]==1),
                "numAssemblies": np.max(h5_file['post_rem_A_ID']),
                "meanReactPerAssembly": np.mean(np.unique(h5_file['post_rem_A_ID'][()],return_counts=True)),
                "meanFreqPerAssembly": np.mean(np.unique(h5_file['post_rem_A_ID'][()],return_counts=True))/recordingLength
            }
        )

        for i in range(len(h5_file["post_rem_react_idx"][()])):
            data_list.append(  # This will create one list entry per cell
                {
                    "eventID": i,
                    "mouse": h5_file["mouse"][0].decode("utf-8"),
                    "condition": h5_file["condition"][0].decode("utf-8"),
                    "replayEventTime": h5_file['post_rem_react_idx'][i]/params['sampling_frequency'],
                    "replayEventID": h5_file['post_rem_A_ID'][i],
                    "replayEventStrength": h5_file['post_rem_A_str'][i],
                    "replayEventSignificance": h5_file['post_rem_A_sig'][i],
                }
            )

        # Close files
        h5_file.close()

df = pd.DataFrame(data_list)
df_numEvents = pd.DataFrame(data_list_num_events)

#%% Plot results
plt.figure(figsize=(1.25,.15))
sns.barplot(
    data=df_numEvents.query("condition=='LTD1'"),
    x='numAssemblies',
    errorbar='se',
    capsize=.2
)
sns.stripplot(
    data=df_numEvents.query("condition=='LTD1'"),
    color='gray',
    x='numAssemblies',
    size=2
)
plt.xlim(0,25)
plt.yticks([])
plt.xlabel('Number of assemblies')
plt.savefig("../../output_REM/REMpost_numAssemblies.pdf")

#%%
mean_numAssemblies = df_numEvents.query("condition=='LTD1'")['numAssemblies'].mean()
SEM_numAssemblies = df_numEvents.query("condition=='LTD1'")['numAssemblies'].sem()
print(f'{mean_numAssemblies} +/- {SEM_numAssemblies}')

#%% Plot frequency
plt.figure(figsize=(1.25,.15))
sns.barplot(
    data=df_numEvents.query("condition=='LTD1'"),
    x='meanFreqPerAssembly',
    errorbar='se',
    capsize=.2
)
sns.stripplot(
    data=df_numEvents.query("condition=='LTD1'"),
    color='gray',
    x='meanFreqPerAssembly',
    size=2
)
#plt.ylim(0,25)
plt.yticks([])
plt.xlabel('Mean reactivation\nfrequency (Hz)')
plt.savefig("../../output_REM/REMpost_meanAssemblyFreq.pdf")

#%%
mean_freqAssemblies = df_numEvents.query("condition=='LTD1'")['meanFreqPerAssembly'].mean()
SEM_freqAssemblies = df_numEvents.query("condition=='LTD1'")['meanFreqPerAssembly'].sem()
print(f'{mean_freqAssemblies} +/- {SEM_freqAssemblies}')

#%%
plt.figure(figsize=(1.25,.75))
sns.histplot(
    data=df.query("condition=='LTD1'"),
    x='replayEventStrength',
    y='replayEventID',
    hue='replayEventID',
    palette='Spectral',
    cbar=False,
    legend=False
)
# plt.title('REM post')
plt.xlabel('Strength (z)')
plt.ylabel('Assembly ID')
plt.xscale('log')
plt.savefig("../../output_REM/REMpost_assemblyStrength.pdf")

# %%
mean_strengthAssemblies = df.query("condition=='LTD1'")['replayEventStrength'].mean()
SEM_strengthAssemblies = df.query("condition=='LTD1'")['replayEventStrength'].sem()
print(f'{mean_strengthAssemblies} +/- {SEM_strengthAssemblies}')
# %%
