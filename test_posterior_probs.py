# %% Imports
import h5py
import yaml
import os
from tqdm import tqdm
import itertools
from utils.helperFunctions import load_data
from pycaan.functions.decoding import bayesian_decode, temporal_bayesian_filter
import numpy as np

# %% Load parameters
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

#%%
condition = 'LTD1'
mouse = 'pv1069'

#%%
current_path = os.path.join(
            params["path_to_output"], "tuning", f"tuning_{condition}_{mouse}.h5"
        )
        
with h5py.File(current_path, "r") as f:
    tuning_curves = f["tuning_curves"][()]
    marginal_likelihood = f["marginal_likelihood"][()]
    # Use uniform prior
    occupancy = np.ones(tuning_curves.shape[1]) / tuning_curves.shape[1]

data_REMpost = load_data(mouse, condition, "REMpost", params)

#%%
REMpost_posterior_probs, _ = bayesian_decode(
                tuning_curves[selected_neurons],
                occupancy,
                marginal_likelihood[selected_neurons],
                data_REMpost["binaryData"][:, selected_neurons],
            )