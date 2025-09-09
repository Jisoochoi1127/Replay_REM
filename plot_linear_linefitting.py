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

#%% Load tuning curves for that session
with h5py.File(
    os.path.join(
        params["path_to_output"],
        'tuning',
        f"tuning_{condition}_{mouse}.h5"
        ), 'r'
        ) as f:
    peak_loc = f['peak_loc'][()]
        
# Load selected neurons
with h5py.File(
    os.path.join(
        params['path_to_output'],
        'neuron_selection',
        "selected_neurons_LTD1_pv1069.h5"
        ), 'r'
        ) as f:
    selected_neurons = f['place_cells'][()]

#%% Select neurons and sort


# %% Plot posterior probabilities during wake
plt.figure(figsize=(2,1))
plt.subplot(211)
plt.imshow(
    sorted_binary.T,
    interpolation='none',
    aspect='auto',
    # vmin=0,
    vmax=.25,
    origin='lower'
)

plt.xlim(0,2500)
plt.xticks(np.arange(0,2500,600),[])
plt.yticks([0,40],[0,100])
#plt.ylabel('Position\n(cm)')
# plt.colorbar(label='Posterior probability')

plt.subplot(212)
plt.plot(data_wake['position'][:,0])
plt.xlim(0,2500)
plt.xticks(np.arange(0,2500,600),np.arange(0,2500/60,10))
plt.xlabel('Time (s)')
plt.ylabel('Position (cm)')
plt.savefig("../../output_REM/linear_fit_example.pdf")

# %% Load all sessions
# %% Import place cell data
results_dir = params['path_to_output']+"/raw_linear_replay"
resultsList = os.listdir(results_dir)

#%%
data_list = []

for file_name in tqdm(resultsList):
    if (
        file_name.startswith("raw_linear_replay_")
        and file_name.endswith(".h5")
        and "pv1254" not in file_name # Exclude pv1254
    ):
        h5_file = h5py.File(os.path.join(results_dir, file_name))
        for i in range(len(h5_file["replay_locs"][()])):
            data_list.append(  # This will create one list entry per cell
                {
                    "eventID": i,
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                    "replayEventTime": h5_file['replay_locs'][i]/params['sampling_frequency'],
                    "replayEventScore": h5_file['replay_score'][i],
                    "replayEventJumpiness": h5_file['replay_jumpiness'][i],
                    "replayEventLength": h5_file['replay_length'][i],
                    "replayEventSlope": h5_file['replay_slope'][i]
                }
            )
        
        # Close files
        h5_file.close()

df = pd.DataFrame(data_list)

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
#TODO adapt to linear replay
# idx = 5 # Pick top examples
# example_idx=df.query("Type=='replay' and replayEventJumpiness>0")['replayEventScore'].sort_values(ascending=False).index[idx]
# example_info=df.iloc[example_idx]
# eventID=example_info['eventID']
# mouse=example_info['mouse']
# condition=example_info['condition']

# with h5py.File(
#                 os.path.join(
#                     params["path_to_output"],
#                     "posterior_probs",
#                     f"posterior_probs_{condition}_{mouse}.h5",
#                 ),
#                 "r",
#             ) as f:
#     REMpost_posterior_probs = f["REMpost_posterior_probs"][()]

# with h5py.File(
#             os.path.join(
#                 params["path_to_output"],
#                 "bayesian_replay",
#                 f"bayesian_replay_{condition}_{mouse}_REMpost.h5",
#             ),
#             "r",
#         ) as f:
#     replayLocs = f['replay_locs'][()]
#     replayScore = f['replay_score'][()]
#     replayJumpiness = f['replay_jumpiness'][()]
#     replayLength = f['replay_length'][()]
#     replaySlope = f['replay_slope'][()]

# max_vec = np.max(REMpost_posterior_probs,axis=1)
# #REMpost_posterior_probs = REMpost_posterior_probs/max_vec[:,None]

# plt.figure(figsize=(.75,.75))
# plt.imshow(REMpost_posterior_probs.T,
#            aspect='auto',
#            vmin=.023,
#            vmax=.035,
#            interpolation='none',
#            origin='lower',
#            rasterized=True)

# plt.xlim(replayLocs[eventID],replayLocs[eventID]+params['windowSize'])
# plt.xticks([replayLocs[eventID], replayLocs[eventID]+params['windowSize']],
#           [0,
#            params['windowSize']/params['sampling_frequency'],
#            ])
# plt.yticks([0,REMpost_posterior_probs.shape[1]],
#            [0,100])
# plt.xlabel('Time (s)')
# plt.ylabel('Decoded\nlocation (cm)')

# # plt.plot([replayLocs[eventID], replayLocs[eventID]+int(params['windowSize'])],
# #          [42,42],
# #          linewidth=2,
# #          color='C1')
# plt.title(f"R$^{2}$ = {replayScore[eventID].round(2)}, {abs(replaySlope[eventID]).round(2)}" + " cm.s$^{-1}$")

# plt.savefig(f"../../output_REM/example_replay_{mouse}_{condition}_{eventID}.pdf")

# # %%
