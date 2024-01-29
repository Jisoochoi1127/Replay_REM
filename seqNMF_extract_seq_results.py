#%% Imports
import numpy as np
import os
from tqdm import tqdm
import itertools

import yaml
import h5py
from utils.helperFunctions import load_data, extract_seq_score

#%% Load parameters
with open('params.yaml','r') as file:
    params = yaml.full_load(file)

np.random.seed(params['seed'])

#%% Define conditions, dataset
states_list = ['REMpre', 'wake', 'REMpost']
condition_list = ['LTD1','LTD5','HATD1','HATD5']
mouse_list = ['pv1043','pv1060', 'pv1069', 'pv1191', 'pv1192', 'pv1252', 'pv1254']

#%% First, measure 'sequenceness' in each individuate session
for condition, mouse, state in tqdm(list(itertools.product(condition_list,
                                                           mouse_list,
                                                           states_list)),
                                                           total=len(condition_list)*len(mouse_list)*len(states_list)):
    
    if not os.path.exists(os.path.join(params['path_to_output'], "place_cells", f'seqResults_{condition}_{mouse}_{state}.h5')):
        try: #TODO remove try statement, instead parse existing sessions in advance!
            # Load data
            data = load_data(mouse, condition, state, params)

            # Load selected neurons
            with h5py.File(os.path.join(params['path_to_output'],"neuron_selection", f'selected_neurons_{condition}_{mouse}.h5'),'r') as f:
                selected_neurons = f['place_cells'][()]
            # selected_neurons = np.arange(params['numNeurons']) # If file don't exist, just pick top-k neurons

            # Extract seq score
            seq_scores, seq_pvalues, seq_locs = extract_seq_score(data['binaryData'][:,selected_neurons], params)
            
            with h5py.File(os.path.join(params['path_to_output'], "place_cells", f'seqResults_{condition}_{mouse}_{state}.h5'),'w') as f:
                f.create_dataset('mouse', data=mouse)
                f.create_dataset('condition', data=condition)
                f.create_dataset('state', data=state)
                f.create_dataset('S1_score', data=seq_scores[0])
                f.create_dataset('S2_score', data=seq_scores[1])
                f.create_dataset('S1_pvalue', data=seq_pvalues[0])
                f.create_dataset('S2_pvalue', data=seq_pvalues[1])
                f.create_dataset('S1_numSeqs', data=len(seq_locs[0]))
                f.create_dataset('S2_numSeqs', data=len(seq_locs[1]))
        except:
            print('Missing session')