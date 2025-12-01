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

# %% Select example recording to plot raw data during wake
mouse = 'pv1069' #TODO pick good example
condition = 'LTD1' #TODO pick good example

# %% Load data
data_wake = load_data(mouse, condition, 'wake', params) #TODO pick good example
    
#%% Load selected neurons
with h5py.File(
    os.path.join(
        params['path_to_output'],
        'neuron_selection',
        f"selected_neurons_{condition}_{mouse}.h5"
        ), 'r'
        ) as f:
    selected_neurons = f['place_cells'][()]

#%% Load tuning curves for that session
with h5py.File(
    os.path.join(
        params["path_to_output"],
        'tuning',
        f"tuning_{condition}_{mouse}.h5"
        ), 'r'
        ) as f:
     selected_peak_loc = f['peak_loc'][()][selected_neurons,0]

#%% Sort neurons
sorting_index = np.argsort(selected_peak_loc)
selected_binary = data_wake['binaryData'][:,selected_neurons]
sorted_binary = selected_binary[:,sorting_index]

# %% Plot posterior probabilities during wake
plt.figure(figsize=(2,1))
plt.subplot(211)
plt.imshow(
    sorted_binary.T,
    interpolation='none',
    aspect='auto',
    # vmin=0,
    rasterized=True,
    vmax=.25,
    origin='lower',
    cmap='gray_r'
)

plt.xlim(0,3000)
plt.xticks(np.arange(0,3000,600),[])
plt.title('Wakefulness')
#plt.ylabel('Position\n(cm)')
# plt.colorbar(label='Posterior probability')

plt.subplot(212)
plt.plot(data_wake['position'][:,0])
plt.xlim(0,3000)
plt.xticks(np.arange(0,3000,600),np.arange(0,3000/60,10))
plt.xlabel('Time (s)')
plt.ylabel('Position (cm)')

plt.savefig("../../output_REM/linear_fit_example.pdf")

# %% Load all sessions
# %% Import place cell data
results_dir = params['path_to_output']+"/raw_linear_replay"
resultsList = os.listdir(results_dir)

#%%
alpha = 0.05

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
                    "replayEvents": sum(~np.isnan(h5_file['replayLocs_P'][h5_file['pvalue_P'][()]<alpha])), # num of actual events
                    "replayScores": np.nanmean(h5_file['replayScore_P'][h5_file['pvalue_P'][()]<alpha]),
                    "replayJumpiness": np.nanmean(h5_file['replayJumpiness_P'][h5_file['pvalue_P'][()]<alpha]),
                    "replayPortion": np.nanmean(h5_file['replayPortion_P'][h5_file['pvalue_P'][()]<alpha]),
                    "replaySlope": np.nanmean(h5_file['replaySlope_P'][h5_file['pvalue_P'][()]<alpha]),
                }
            )
        
        # Bonferroni
        bonferroni_alpha = alpha/len(h5_file['pvalue_P'][()])
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
                    "replayPortion": np.nanmean(h5_file['replayPortion_P'][h5_file['pvalue_P'][()]<bonferroni_alpha]),
                    "replaySlope": np.nanmean(h5_file['replaySlope_P'][h5_file['pvalue_P'][()]<bonferroni_alpha]),
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
                    "replayEvents": sum(~np.isnan(h5_file['control_replayLocs_P'][h5_file['control_pvalue_P'][()]<alpha])), # num of actual events
                    "replayScores": np.nanmean(h5_file['control_replayScore_P'][h5_file['control_pvalue_P'][()]<alpha]),
                    "replayJumpiness": np.nanmean(h5_file['control_replayJumpiness_P'][h5_file['control_pvalue_P'][()]<alpha]),
                    "replayPortion": np.nanmean(h5_file['control_replayPortion_P'][h5_file['control_pvalue_P'][()]<alpha]),
                    "replaySlope": np.nanmean(h5_file['control_replaySlope_P'][h5_file['control_pvalue_P'][()]<alpha]),
                }
            )
        
        # Bonferroni
        bonferroni_alpha = alpha/len(h5_file['control_pvalue_P'][()])
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
                    "replayPortion": np.nanmean(h5_file['control_replayPortion_P'][h5_file['control_pvalue_P'][()]<bonferroni_alpha]),
                    "replaySlope": np.nanmean(h5_file['control_replaySlope_P'][h5_file['control_pvalue_P'][()]<bonferroni_alpha]),
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
                    "replayEvents": sum(~np.isnan(h5_file['replayLocs_T'][h5_file['pvalue_T'][()]<alpha])), # num of actual events
                    "replayScores": np.nanmean(h5_file['replayScore_T'][h5_file['pvalue_T'][()]<alpha]),
                    "replayJumpiness": np.nanmean(h5_file['replayJumpiness_T'][h5_file['pvalue_T'][()]<alpha]),
                    "replayPortion": np.nanmean(h5_file['replayPortion_T'][h5_file['pvalue_P'][()]<alpha]),
                    "replaySlope": np.nanmean(h5_file['replaySlope_T'][h5_file['pvalue_P'][()]<alpha]),
                    
                }
            )
        
        # Bonferroni
        bonferroni_alpha = alpha/len(h5_file['pvalue_T'][()])
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
                    "replayPortion": np.nanmean(h5_file['replayPortion_T'][h5_file['pvalue_P'][()]<bonferroni_alpha]),
                    "replaySlope": np.nanmean(h5_file['replaySlope_T'][h5_file['pvalue_P'][()]<bonferroni_alpha]),
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
                    "replayEvents": sum(~np.isnan(h5_file['control_replayLocs_T'][h5_file['control_pvalue_T'][()]<alpha])), # num of actual events
                    "replayScores": np.nanmean(h5_file['control_replayScore_T'][h5_file['control_pvalue_T'][()]<alpha]),
                    "replayJumpiness": np.nanmean(h5_file['control_replayJumpiness_T'][h5_file['control_pvalue_T'][()]<alpha]),
                    "replayPortion": np.nanmean(h5_file['control_replayPortion_T'][h5_file['control_pvalue_P'][()]<alpha]),
                    "replaySlope": np.nanmean(h5_file['control_replaySlope_T'][h5_file['control_pvalue_P'][()]<alpha]),
                }
            )
        
        # Bonferroni
        bonferroni_alpha = alpha/len(h5_file['control_pvalue_T'][()])
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
                    "replayPortion": np.nanmean(h5_file['control_replayPortion_T'][h5_file['control_pvalue_P'][()]<bonferroni_alpha]),
                    "replaySlope": np.nanmean(h5_file['control_replaySlope_T'][h5_file['control_pvalue_P'][()]<bonferroni_alpha]),
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
                    "replayEvents": sum(~np.isnan(h5_file['replayLocs_PT'][h5_file['pvalue_PT'][()]<alpha])), # num of actual events
                    "replayScores": np.nanmean(h5_file['replayScore_PT'][h5_file['pvalue_PT'][()]<alpha]),
                    "replayJumpiness": np.nanmean(h5_file['replayJumpiness_PT'][h5_file['pvalue_PT'][()]<alpha]),
                    "replayPortion": np.nanmean(h5_file['replayPortion_PT'][h5_file['pvalue_P'][()]<alpha]),
                    "replaySlope": np.nanmean(h5_file['replaySlope_PT'][h5_file['pvalue_P'][()]<alpha]),
                }
            )
        
        # Bonferroni
        bonferroni_alpha = alpha/len(h5_file['pvalue_PT'][()])
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
                    "replayPortion": np.nanmean(h5_file['replayPortion_PT'][h5_file['pvalue_P'][()]<bonferroni_alpha]),
                    "replaySlope": np.nanmean(h5_file['replaySlope_PT'][h5_file['pvalue_P'][()]<bonferroni_alpha]),
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
                    "replayEvents": sum(~np.isnan(h5_file['control_replayLocs_PT'][h5_file['control_pvalue_PT'][()]<alpha])), # num of actual events
                    "replayScores": np.nanmean(h5_file['control_replayScore_PT'][h5_file['control_pvalue_PT'][()]<alpha]),
                    "replayJumpiness": np.nanmean(h5_file['control_replayJumpiness_PT'][h5_file['control_pvalue_PT'][()]<alpha]),
                    "replayPortion": np.nanmean(h5_file['control_replayPortion_PT'][h5_file['control_pvalue_P'][()]<alpha]),
                    "replaySlope": np.nanmean(h5_file['control_replaySlope_PT'][h5_file['control_pvalue_P'][()]<alpha]),
                }
            )
        
        # Bonferroni
        bonferroni_alpha = alpha/len(h5_file['control_pvalue_PT'][()])
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
                    "replayPortion": np.nanmean(h5_file['control_replayPortion_PT'][h5_file['control_pvalue_P'][()]<bonferroni_alpha]),
                    "replaySlope": np.nanmean(h5_file['control_replaySlope_PT'][h5_file['control_pvalue_P'][()]<bonferroni_alpha]),
                }
            )
        
        # Close files
        h5_file.close()

df = pd.DataFrame(data_list)

#%% Control
plt.figure(figsize=(2,1))
sns.barplot(
    data=df.query("Type=='replay' and FDR == 'None'"),
    x='replayEvents',
    hue="Decoding",
    y='Shuffle',
    legend=True,
    errorbar='se'
)

sns.stripplot(
    data=df.query("Type=='replay' and FDR == 'None'"),
    x='replayEvents',
    hue="Decoding",
    y='Shuffle',
    legend=False,
    dodge=True,
    size=1
)












# %% Plot results
plt.figure(figsize=(1,.25))
sns.histplot(
    data=df.query("Type=='replay'"),
    x='replayEventJumpiness',
)
plt.title('REM post')
plt.xlabel('Jumpiness (cm)')
plt.ylabel('N')
plt.savefig("../../output_REM/REMpost_raw_linear_replay_Jumpiness.pdf")

#%%
plt.figure(figsize=(1,.25))
sns.histplot(
    data=df.query("Type=='replay' and replayEventJumpiness>0"),
    x='replayEventScore',
)
# plt.title('REM post')
plt.xlabel('Replay score (R$^{2}$)')
plt.ylabel('N')
plt.savefig("../../output_REM/REMpost_raw_linear_replay_replayScores.pdf")

#%% DESCRIPTIVES
mean_S1 = df.query("Type=='replay' and replayEventJumpiness>0")['replayEventScore'].mean()
SEM_S1 = df.query("Type=='replay' and replayEventJumpiness>0")['replayEventScore'].sem()
print(f'{mean_S1} +/- {SEM_S1}')

#%%
plt.figure(figsize=(1.25,.33))
sns.histplot(
    data=df.query("Type=='replay' and replayEventJumpiness==0"),
    x='replayEventSlope',
)
# plt.title('REM post')
plt.xlabel('Replay slope')
plt.ylabel('N')
plt.savefig("../../output_REM/REMpost_raw_linear_replay_replaySlope.pdf")

#%% Example example replay events
idx = 27 # Top k example
example_idx=df.query("Type=='replay' and replayEventJumpiness>0")['replayEventScore'].sort_values(ascending=False).index[idx]
example_info=df.iloc[example_idx]
eventID=example_info['eventID']
mouse=example_info['mouse']
condition=example_info['condition']

# Load data
data = load_data(mouse, condition, 'REMpost', params)

# Load selected neurons
with h5py.File(
    os.path.join(
        params['path_to_output'],
        'neuron_selection',
        f"selected_neurons_{condition}_{mouse}.h5"
        ), 'r'
        ) as f:
    selected_neurons = f['place_cells'][()]

# Load tuning curves for that session
with h5py.File(
    os.path.join(
        params["path_to_output"],
        'tuning',
        f"tuning_{condition}_{mouse}.h5"
        ), 'r'
        ) as f:
     selected_peak_loc = f['peak_loc'][()][selected_neurons,0]

# Sort neurons
sorting_index = np.argsort(selected_peak_loc)
selected_binary = data['binaryData'][:,selected_neurons]
sorted_binary = selected_binary[:,sorting_index]

# Load replay properties
with h5py.File(
            os.path.join(
                params["path_to_output"],
                "raw_linear_replay",
                f"raw_linear_replay_{condition}_{mouse}_REMpost.h5",
            ),
            "r",
        ) as f:
    replayLocs = f['replay_locs'][()]
    replayScore = f['replay_score'][()]
    replayJumpiness = f['replay_jumpiness'][()]
    replayLength = f['replay_length'][()]
    replaySlope = f['replay_slope'][()]

# Extract main vector
active_neuron = np.argmax(sorted_binary,axis=1).astype('float') #TODO if multiple neurons active, use their average identity?
# Remove silent epochs
active_neuron[np.sum(sorted_binary,axis=1)==0]=np.nan

plt.figure(figsize=(.75,.75))
plt.imshow(sorted_binary.T,
           aspect='auto',
           #vmin=.023,
           #vmax=.035,
           interpolation='none',
           origin='lower',
           cmap='gray_r',
           rasterized=True)

plt.xlim(replayLocs[eventID],replayLocs[eventID]+params['windowSize'])
plt.xticks([replayLocs[eventID], replayLocs[eventID]+params['windowSize']],
          [0,
           params['windowSize']/params['sampling_frequency'],
           ])
plt.yticks([0,sorted_binary.shape[1]],
           [0,100])
plt.xlabel('Time (s)')
plt.ylabel('Decoded\nlocation (cm)')

# plt.plot([replayLocs[eventID], replayLocs[eventID]+int(params['windowSize'])],
#          [42,42],
#          linewidth=2,
#          color='C1')
plt.title(f"R$^{2}$ = {replayScore[eventID].round(2)}")

plt.savefig(f"../../output_REM/example_linear_fit_replay_{mouse}_{condition}_{eventID}.pdf")

# %%
