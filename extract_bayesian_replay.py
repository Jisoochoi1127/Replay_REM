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

results_dir = params['path_to_output']+'/posterior_probs'
resultsList=os.listdir(results_dir)

for file_name in tqdm(resultsList):
    with h5py.File(os.path.join(results_dir, file_name), 'r') as f:
        mouse = f['mouse'][()].decode("utf-8")
        condition = f['condition'][()].decode("utf-8")

        for state in ['REMpre', 'REMpost']:
            if not os.path.exists(os.path.join(params['path_to_output'],"bayesian_replay", f'bayesian_replay_{condition}_{mouse}_{state}.h5')):
            # Load precomputed tuning curves and accessory data
                posterior_probs = f[f'{state}_posterior_probs'][()]
                
                replayLocs, replayScore, replayJumpiness, replayPortion = extract_linear_replay(posterior_probs, params)
            
                # Save results
                with h5py.File(os.path.join(params['path_to_output'],"bayesian_replay", f'bayesian_replay_{condition}_{mouse}_{state}.h5'),'w') as f:
                    f.create_dataset('mouse', data=mouse)
                    f.create_dataset('condition', data=condition)
                    f.create_dataset('replay_locs', data=replayLocs)
                    f.create_dataset('replay_score', data=replayScore)
                    f.create_dataset('replay_jumpiness', data=replayJumpiness)
                    f.create_dataset('replay_length', data=replayPortion)