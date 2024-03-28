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
condition = "LTD1"
mouse_list = ["pv1043", "pv1060", "pv1069", "pv1254"]
data_list = []

# %% Replay on LTD1
for mouse in mouse_list:
    # Load assemblies
    with h5py.File(
        os.path.join(
            params["path_to_output"],
            'assembly',
            f'assembly_{condition}_{mouse}'
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
            f'bayesian_replay_{condition}_{mouse}_REMpost'
        ),
        "r",
    ) as bayesianReplay_file:
        bayesian_replay_ts = bayesianReplay_file['replay_locs']

    # Load seqNMF replay
    with h5py.File(
        os.path.join(
            params["path_to_output"],
            'seqNMF',
            f'seqReplayResults_{condition}_{mouse}_wake_REMpost'
        ),
        "r",
    ) as seqNMF_file:
        seqNMF_ts = seqNMF_file['seqReplayLocs'][()].flatten() # combined both seq types
        
    
    # Match events

    # Append data
    data_list.append(
        {
            'mouse':mouse,
            'num_strict_unstructured':num_strict_unstructured,
            'num_strict_flexible':num_strict_flexible,
            'num_strict_rigid':num_strict_rigid,
        }
    )

#%% Plot results
