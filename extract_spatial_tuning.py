#%% Imports
import numpy as np
import os
from tqdm import tqdm
import itertools

import yaml
import h5py
from utils.helperFunctions import load_data, extract_seq_score
from pycaan.functions.tuning import extract_tuning

#%% Load parameters
with open('params.yaml','r') as file:
    params = yaml.full_load(file)

np.random.seed(params['seed'])

#%% Define conditions, dataset
states_list = ['REMpre', 'wake', 'REMpost']
condition_list = ['LTD1','LTD5','HATD1','HATD5']
mouse_list = ['pv1043','pv1060', 'pv1069', 'pv1191', 'pv1192', 'pv1252', 'pv1254']

#%% First, measure 'sequenceness' in each individuate session
for condition, mouse in tqdm(list(itertools.product(condition_list,
                                                           mouse_list)),
                                                           total=len(condition_list)*len(mouse_list)):
    
    if not os.path.exists(os.path.join(params['path_to_output'],f'tuning_{condition}_{mouse}.h5')):
        try:
            # Load data
            data = load_data(mouse, condition, 'wake', params)

            # Extract tuning info/tuning 
            bins=
            extract_tuning()
            
            with h5py.File(os.path.join(params['path_to_output'],f'tuning_{condition}_{mouse}.h5'),'w') as f:
                f.create_dataset('mouse', data=mouse)
                f.create_dataset('condition', data=condition)
                
        except:
            print('Missing session')