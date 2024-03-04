#%%
from utils.helperFunctions import load_data
from pycaan.functions.decoding import bayesian_decode
import h5py
import os
import yaml

#%% Load parameters
with open('params.yaml','r') as file:
    params = yaml.full_load(file)

#%%
with h5py.File(
    os.path.join(
        params["path_to_output"],
        'tuning',
        "tuning_LTD1_pv1043.h5"
        ), 'r'
        ) as f:
    tuning_curves = f['tuning_curves'][()]
    marginal_likelihood = f['marginal_likelihood'][()]
    occupancy = [] # TODO use uniform prior
        
# Load selected neurons
with h5py.File(
    os.path.join(
        params['path_to_output'],
        'neuron_selection',
        "selected_neurons_LTD1_pv1043.h5"
        ), 'r'
        ) as f:
    selected_neurons = f['place_cells'][()]
# %%
