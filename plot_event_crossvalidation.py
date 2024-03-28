# %% Imports
import os
import numpy as np
import pandas as pd
import yaml
import h5py
from tqdm import tqdm

# %% Load params
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

# %%
condition_list = "LTD1"
mouse_list = ["pv1043", "pv1060", "pv1069", "pv1254"]
data_list = []

# %% Replay on LTD1
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
            num_strict_rigid = len(seqNMF_ts)

            strict_flexible = []
            for event_ts in bayesian_replay_ts:
                if len(np.where(np.logical_and(seqNMF_ts>=event_ts-int(params['event_overlap_window']/2), seqNMF_ts<=event_ts+int(params['event_overlap_window']/2)))[0])==0:
                    strict_flexible.append(event_ts)

            num_strict_flexible = len(strict_flexible)

            strict_unstructured = []
            for event_ts in assembly_ts:
                if (len(np.where(np.logical_and(seqNMF_ts>=event_ts-int(params['event_overlap_window']/2), seqNMF_ts<=event_ts+int(params['event_overlap_window']/2)))[0])==0
                    and len(np.where(np.logical_and(bayesian_replay_ts>=event_ts-int(params['event_overlap_window']/2), bayesian_replay_ts<=event_ts+int(params['event_overlap_window']/2)))[0])==0
                ):
                    strict_unstructured.append(event_ts)

            num_strict_unstructured = len(strict_unstructured)

            # Append data
            data_list.append(
                {
                    'mouse':mouse,
                    'num_strict_unstructured':num_strict_unstructured,
                    'num_strict_flexible':num_strict_flexible,
                    'num_strict_rigid': num_strict_rigid,
                }
            )

df = pd.DataFrame(data_list)
#%% Plot results
