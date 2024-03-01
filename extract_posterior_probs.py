# %% Imports
import h5py
import yaml
import os
from tqdm import tqdm
import itertools
from utils.helperFunctions import load_data
from pycaan.functions.decoding import baysian_decode

#%% Load parameters
with open('params.yaml','r') as file:
    params = yaml.full_load(file)

#%% Define conditions, dataset
condition_list = ['LTD1','LTD5','HATD1','HATD5']
mouse_list = ['pv1043','pv1060', 'pv1069', 'pv1191', 'pv1192', 'pv1252', 'pv1254']

#%% Same but look at replay between conditions
for condition, mouse, state_ref, state_pred in tqdm(list(itertools.product(condition_list,
                                                           mouse_list)),
                                                           total=len(condition_list)*len(mouse_list)):
    if not os.path.exists(os.path.join(params['path_to_output'],"posterior_probs", f'posterior_probs_{condition}_{mouse}.h5')):
        # Load precomputed tuning curves and accessory data
        with h5py.File(
            os.path.join(
                params["path_to_output"],
                'tuning_curves',
                f"tuning_{condition}_{mouse}.h5"
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
                f'selected_neurons_{condition}_{mouse}.h5'
                ), 'r'
                ) as f:
            selected_neurons = f['place_cells'][()]

        # Load REMpre data
        data_REMpre = load_data(mouse, condition, 'REMpre', params)

        # Load REMpost data
        data_REMpost = load_data(mouse, condition, 'REMpost', params)

        # Extract posteriors on REMpre and save results
        REMpre_posterior_probs, _ = baysian_decode(tuning_curves[selected_neurons], occupancy, marginal_likelihood[selected_neurons], data_REMpre['binaryData'][:,selected_neurons])

        # Extract posterior probs on REMpost and save results
        REMpost_posterior_probs, _ = baysian_decode(tuning_curves[selected_neurons], occupancy, marginal_likelihood[selected_neurons], data_REMpost['binaryData'][:,selected_neurons])

        # Save results
        with h5py.File(
            os.path.join(
                params['path_to_output'],
                'posterior_probs',
                f'posterior_probs_{condition}_{mouse}.h5'
                ), 'r'
                ) as f:
            f.create_dataset('mouse', data=mouse)
            f.create_dataset('condition', data=condition)
            f.create_dataset('REMpre_posterior_probs', data=REMpre_posterior_probs)
            f.create_dataset('REMpost_posterior_probs', data=REMpost_posterior_probs)