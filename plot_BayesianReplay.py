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
from pycaan.functions.decoding import bayesian_decode

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

#%% Load tuning curves for that session
with h5py.File(
    os.path.join(
        params["path_to_output"],
        'tuning',
        "tuning_LTD1_pv1069.h5"
        ), 'r'
        ) as f:
    tuning_curves = f['tuning_curves'][()]
    marginal_likelihood = f['marginal_likelihood'][()]
    # Use uniform prior
    occupancy = np.ones(tuning_curves.shape[1])/tuning_curves.shape[1]
        
# Load selected neurons
with h5py.File(
    os.path.join(
        params['path_to_output'],
        'neuron_selection',
        "selected_neurons_LTD1_pv1069.h5"
        ), 'r'
        ) as f:
    selected_neurons = f['place_cells'][()]

#%% Compute posterior probs during wake
posterior_probs, map = bayesian_decode(tuning_curves,
                                       occupancy,
                                       marginal_likelihood,
                                       data_wake['binaryData'])
# %% Plot posterior probabilities during wake
plt.figure(figsize=(2,1))
plt.subplot(211)
plt.imshow(
    posterior_probs.T,
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
plt.savefig("../../output_REM/posterior_probs.pdf")

# %% Load all sessions
# %% Import place cell data
results_dir = params['path_to_output']+"/bayesian_replay"
resultsList = os.listdir(results_dir)

#%%
data_list = []

for file_name in tqdm(resultsList):
    if (
        file_name.startswith("bayesian_replay_")
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
                    "replayEventLength": h5_file['replay_length'][i]
                }
            )
        
        # Close files
        h5_file.close()

df = pd.DataFrame(data_list)

# %% Plot results
sns.histplot(
    data=df,
    x='replayEventJumpiness',
)
plt.title('REM post')
plt.xlabel('Replay jumpiness (cm)')
plt.ylabel('N')
plt.savefig("../../output_REM/REMpost_allEventsJumpiness.pdf")

#%%
sns.histplot(
    data=df.query("Type=='replay' and replayEventJumpiness>0"),
    x='replayEventScore',
)
plt.title('REM post')
plt.xlabel('Replay score (R$^{2}$)')
plt.ylabel('N')
plt.savefig("../../output_REM/REMpost_replayScores.pdf")

#%% Example example replay events
idx = 0 # Pick top examples
example_idx=df.query("Type=='replay' and replayEventJumpiness>0")['replayEventScore'].sort_values(ascending=False).index[idx]
example_info=df.iloc[example_idx]
eventID=example_info['eventID']
mouse=example_info['mouse']
condition=example_info['condition']

# %% Load posterior probabilities during REM
with h5py.File(
                os.path.join(
                    params["path_to_output"],
                    "posterior_probs",
                    f"posterior_probs_{condition}_{mouse}.h5",
                ),
                "r",
            ) as f:
    REMpost_posterior_probs = f["REMpost_posterior_probs"][()]

# %% Load replay info
with h5py.File(
            os.path.join(
                params["path_to_output"],
                "bayesian_replay",
                f"bayesian_replay_{condition}_{mouse}_REMpost.h5",
            ),
            "r",
        ) as f:
    replayLocs = f['replay_locs'][()]
    replayScore = f['replay_score'][()]
    replayJumpiness = f['replay_jumpiness'][()]
    replayLength = f['replay_length'][()]

# %% Plot significant replay revents
plt.figure(figsize=(2,.75))

plt.imshow(REMpost_posterior_probs.T,aspect='auto',vmax=.05,interpolation='none',
           origin='lower')

plt.xlim(replayLocs[eventID]-100,replayLocs[eventID]+100)
plt.xticks([replayLocs[eventID]-100,replayLocs[eventID], replayLocs[eventID]+100],
          [0,int(100/30),int(200/30)])
plt.yticks([0,REMpost_posterior_probs.shape[1]],
           [0,100])
plt.xlabel('Time (s)')
plt.ylabel('Decoded\nlocation (cm)')

plt.plot([replayLocs[eventID], replayLocs[eventID]+int(params['windowSize'])],
         [42,42],
         linewidth=2,
         color='C4')
plt.title(f"R$^{2}$ = {replayScore[eventID].round(2)}")

plt.savefig(f"../../output_REM/example_replay_{mouse}_{condition}_{eventID}.pdf")

# %%
