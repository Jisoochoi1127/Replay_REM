#%% Imports
import os
from tqdm import tqdm
import itertools
import yaml
import h5py
from utils.helperFunctions import load_data, extract_seqReplay_score

#%% Load parameters
with open('params.yaml','r') as file:
    params = yaml.full_load(file)

#%% Define conditions, dataset
states_list = ['REMpre', 'wake', 'REMpost']
condition_list = ['LTD1','LTD5','HATD1','HATD5']
mouse_list = ['pv1043','pv1060', 'pv1069', 'pv1191', 'pv1192', 'pv1252', 'pv1254']

#%% Same but look at replay between conditions
for condition, mouse, state_ref, state_pred in tqdm(list(itertools.product(condition_list,
                                                           mouse_list,
                                                           states_list,
                                                           states_list)),
                                                           total=len(condition_list)*len(mouse_list)*len(states_list)*len(states_list)):
    if not os.path.exists(os.path.join(params['path_to_output'],f'seqReplayResults_{condition}_{mouse}_{state_ref}_{state_pred}.h5')):
        try:
            # Load data for both states
            data_ref = load_data(mouse, condition, state_ref, params)
            data_pred = load_data(mouse, condition, state_pred, params)

            # Extract seq score
            seqReplay_scores, seqReplay_pvalues, seqReplay_locs = extract_seqReplay_score(data_ref['binaryData'], data_pred['binaryData'], params)

            with h5py.File(os.path.join(params['path_to_output'],f'seqReplayResults_{condition}_{mouse}_{state_ref}_{state_pred}.h5'),'w') as f:
                f.create_dataset('mouse', data=mouse)
                f.create_dataset('condition', data=condition)
                f.create_dataset('state_ref', data=state_ref)
                f.create_dataset('state_pred', data=state_pred)
                f.create_dataset('S1_score', data=seqReplay_scores[0])
                f.create_dataset('S2_score', data=seqReplay_scores[1])
                f.create_dataset('S1_pvalue', data=seqReplay_pvalues[0])
                f.create_dataset('S2_pvalue', data=seqReplay_pvalues[1])
                f.create_dataset('S1_numSeqs', data=len(seqReplay_locs[0]))
                f.create_dataset('S2_numSeqs', data=len(seqReplay_locs[1]))
        except:
            print('Missing session')
