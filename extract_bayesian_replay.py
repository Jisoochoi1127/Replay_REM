# %% Imports
import h5py
import yaml
import os
from tqdm import tqdm
import itertools
from utils.bayesian_replay import extract_linear_replay

#%% Load parameters
with open('params.yaml','r') as file:
    params = yaml.full_load(file)

#%% Define conditions, dataset
condition_list = ['LTD1','LTD5','HATD1','HATD5']
mouse_list = ['pv1043','pv1060', 'pv1069', 'pv1191', 'pv1192', 'pv1252', 'pv1254']
state_list = ['REMpre', 'REMpost']

#%% For each mouse, condition
for condition, mouse, state in tqdm(list(itertools.product(condition_list,
                                                           mouse_list,
                                                           state_list)),
                                                           total=len(condition_list)*len(mouse_list)*len(state_list)):
    if not os.path.exists(os.path.join(params['path_to_output'],"bayesian_replay", f'bayesian_replay_{condition}_{mouse}_{state}.h5')):
        # Load precomputed tuning curves and accessory data
        with h5py.File(
            os.path.join(
                params['path_to_output'],
                'posterior_probs',
                f'posterior_probs_{condition}_{mouse}.h5'
                ), 'r'
                ) as f:
            posterior_probs = f['{state}_posterior_probs'][()]
            
        replayLocs, replayScore, replayJumpiness, replayPortion = extract_linear_replay(posterior_probs, params)
        
        # Save results
        with h5py.File(os.path.join(params['path_to_output'],"bayesian_replay", f'bayesian_replay_{condition}_{mouse}_{state}.h5'),'w') as f:
            f.create_dataset('mouse', data=mouse)
            f.create_dataset('condition', data=condition)
            f.create_dataset('replay_locs', data=replayLocs)
            f.create_dataset('replay_score', data=replayScore)
            f.create_dataset('replay_jumpiness', data=replayJumpiness)
            f.create_dataset('replay_length', data=replayPortion)