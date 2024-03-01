# %% Imports
import h5py
import yaml
import numpy as np
import os
from tqdm import tqdm
import itertools
from utils.helperFunctions import load_data, extract_seqReplay_score

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
        # Load precomputed tuning curves
        
        # Load selected neurons

        # Load REMpre data

        # Load REMpost data
        
        data_REMpre = load_data(mouse, condition, 'REMpre', params)
        data_REMpost = load_data(mouse, condition, 'REMpost', params)



# %% Load tuning curves
for file_name in tqdm(resultsList):
    if file_name.startswith('tuning_') and file_name.endswith('.h5'): # Only include seqResults, not replay results
        h5_file = h5py.File(os.path.join(results_dir,file_name))
        mouse = h5_file['mouse'][()].decode("utf-8")
        condition = h5_file['condition'][()].decode("utf-8")
        marginal_likelihood = h5_file['marginal_likelihood'][()] 
        