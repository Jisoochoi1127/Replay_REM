# %% Imports
import os
import numpy as np
import pandas as pd
import yaml
import h5py
import seaborn as sns
import matplotlib.pyplot as plt
from utils.helperFunctions import load_data

plt.style.use("plot_style.mplstyle")

# %% Load params
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

# %%
condition_list = ["LTD1", "LTD5", 'HATD1', 'HATD5']
mouse_list = ["pv1043", "pv1060", "pv1069", "pv1254"]

# %% Replay on LTD1
data_list = []
for mouse in mouse_list:
    for condition in condition_list:
        if os.path.exists(os.path.join(params['path_to_dataset'],
                                    mouse,
                                    condition)):
            # Load assemblies
            with h5py.File(
                os.path.join(
                    params["path_to_output"],
                    'assembly',
                    f'assembly_{condition}_{mouse}.h5'
                ),
                "r",
            ) as assembly_file:
                assembly_sig = assembly_file['post_rem_A_sig'][()]
                assembly_ts = assembly_file['post_rem_react_idx'][()]
                assembly_ts = assembly_ts[assembly_sig==1]

            # Load Bayesian replay
            with h5py.File(
                os.path.join(
                    params["path_to_output"],
                    'bayesian_replay',
                    f'bayesian_replay_{condition}_{mouse}_REMpost.h5'
                ),
                "r",
            ) as bayesianReplay_file:
                bayesian_replay_ts = bayesianReplay_file['replay_locs'][()]

            # Load seqNMF replay
            with h5py.File(
                os.path.join(
                    params["path_to_output"],
                    'seqNMF',
                    f'seqReplayResults_{condition}_{mouse}_wake_REMpost.h5'
                ),
                "r",
            ) as seqNMF_file:
                seqNMF_ts = seqNMF_file['seqReplayLocs'][()] # combined both seq types
                
            # Match events
            REMpost_data = load_data(mouse,condition,'REMpost',params)
            numFrames = REMpost_data['binaryData'].shape[0]

            num_strict_rigid = len(seqNMF_ts)
            rigid_freq = num_strict_rigid/(numFrames/params['sampling_frequency'])

            strict_flexible = []
            flexible_and_rigid = []
            for event_ts in bayesian_replay_ts:
                if len(np.where(np.logical_and(seqNMF_ts>=event_ts-int(params['event_overlap_window']/2), seqNMF_ts<=event_ts+int(params['event_overlap_window']/2)))[0])==0:
                    strict_flexible.append(event_ts)
                # else:
                #     flexible_and_rigid.append(event_ts)

            num_strict_flexible = len(strict_flexible)
            flexible_freq = num_strict_flexible/(numFrames/params['sampling_frequency'])
            num_flexible_and_rigid = len(flexible_and_rigid)

            strict_unstructured = []
            for event_ts in assembly_ts:
                if (len(np.where(np.logical_and(seqNMF_ts>=event_ts-int(params['event_overlap_window']/2), seqNMF_ts<=event_ts+int(params['event_overlap_window']/2)))[0])==0
                    and len(np.where(np.logical_and(bayesian_replay_ts>=event_ts-int(params['event_overlap_window']/2), bayesian_replay_ts<=event_ts+int(params['event_overlap_window']/2)))[0])==0
                ):
                    strict_unstructured.append(event_ts)

            num_strict_unstructured = len(strict_unstructured)
            unstructured_freq = num_strict_unstructured/(numFrames/params['sampling_frequency'])

            # Append data
            data_list.append(
                {
                    'mouse':mouse,
                    'condition':condition,
                    # 'Unstructured_nu':num_strict_unstructured,
                    # 'Flexible':num_strict_flexible,
                    # 'Rigid': num_strict_rigid,
                    'Unstructured':unstructured_freq,
                    'Flexible':flexible_freq,
                    'Rigid': rigid_freq
                }
            )

df = pd.DataFrame(data_list)

#%% Melt data
df = df.melt(
    id_vars=['mouse', 'condition'],
    value_vars=['Unstructured', 'Flexible', 'Rigid'],
    var_name='Type',
    value_name='Event frequency (Hz)')

#%% Plot results
sns.barplot(
    data = df.query("condition=='LTD1'"),
    x = 'Type',
    hue = 'Type',
    y = 'Event frequency (Hz)',
    palette=['C0','C2','C4'],
   errorbar='se',
   capsize=.2,
    # showfliers=False
    )
sns.stripplot(
    data = df.query("condition=='LTD1'"),
    x = 'Type',
    y = 'Event frequency (Hz)',
    color='gray',
    size=2
    )
plt.xticks(rotation=90)
plt.xlabel('')
plt.title('Overlapping events filtered')
#plt.ylabel('Num. specific\nevents')
plt.savefig("../../output_REM/LTD1_cross_validation.pdf")

# %% Plot example overlap
mouse = 'pv1069'
condition = 'LTD1'
REMpost_data = load_data(mouse,condition,'REMpost',params)
numFrames = REMpost_data['binaryData'].shape[0]

with h5py.File(
    os.path.join(
        params["path_to_output"],
        "neuron_selection",
        f"selected_neurons_{condition}_{mouse}.h5",
    ),
    "r",
) as f:
    selected_neurons = f["place_cells"][()]

with h5py.File(
                os.path.join(
                    params["path_to_output"],
                    'assembly',
                    f'assembly_{condition}_{mouse}.h5'
                ),
                "r",
            ) as assembly_file:
                assembly_sig = assembly_file['post_rem_A_sig'][()]
                assembly_ts = assembly_file['post_rem_react_idx'][()]
                assembly_ts = assembly_ts[assembly_sig==1]

# Load Bayesian replay
with h5py.File(
    os.path.join(
        params["path_to_output"],
        'bayesian_replay',
        f'bayesian_replay_{condition}_{mouse}_REMpost.h5'
    ),
    "r",
) as bayesianReplay_file:
    bayesian_replay_ts = bayesianReplay_file['replay_locs'][()]

# Load seqNMF replay
with h5py.File(
    os.path.join(
        params["path_to_output"],
        'seqNMF',
        f'seqReplayResults_{condition}_{mouse}_wake_REMpost.h5'
    ),
    "r",
) as seqNMF_file:
    seqNMF_ts = seqNMF_file['seqReplayLocs'][()] # combined both seq types

# Plot events
plt.figure(figsize=(1.5,3))
plt.subplot(211)
plt.imshow(REMpost_data['binaryData'][:,selected_neurons].T,
           interpolation='none',
           aspect='auto',
           cmap='gray_r')
plt.xlim(1000,1800)
plt.xticks([])
plt.ylabel('Neuron ID')

plt.subplot(413)
plt.eventplot(assembly_ts, lineoffsets=1, color = 'C0')
plt.eventplot(bayesian_replay_ts, lineoffsets=2, color = 'C2')
plt.eventplot(seqNMF_ts, lineoffsets=3, color = 'C4')
plt.xlim(1000,1800)
plt.xticks([1000,1400,1800],[0,30,60])
plt.xlabel('Time (s)')
plt.yticks([1,2,3],['Unstructured','Flexible','Rigid'])
plt.savefig("../../output_REM/LTD1_cross_validation_example.pdf")
# %%
