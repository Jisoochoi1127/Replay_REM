# %% Imports
import h5py
import yaml
import os
from tqdm import tqdm
from utils.bayesian_replay import extract_linear_replay

#%% Load parameters
with open('params.yaml','r') as file:
    params = yaml.full_load(file)

results_dir = params['path_to_output']+'/equal_posterior_probs' # GE 20250130: changed from posterior_probs to equal_posterior_probs
resultsList=os.listdir(results_dir)

#%%
for file_name in tqdm(resultsList):
    if file_name.startswith('posterior'):
        with h5py.File(os.path.join(results_dir, file_name), 'r') as f:
            mouse = f['mouse'][()].decode("utf-8")
            condition = f['condition'][()].decode("utf-8")

            for state in ['REMpre', 'REMpost']:
                if not os.path.exists(os.path.join(params['path_to_output'],"equal_bayesian_replay", f'bayesian_replay_{condition}_{mouse}_{state}.h5')):
                # Load precomputed tuning curves and accessory data
                    posterior_probs = f[f'{state}_posterior_probs'][()]
                    
                    replayLocs, replayScore, replayJumpiness, replayPortion, replaySlope = extract_linear_replay(posterior_probs, params)
                
                    # Save results
                    with h5py.File(os.path.join(params['path_to_output'],"equal_bayesian_replay", f'bayesian_replay_{condition}_{mouse}_{state}.h5'),'w') as f2:
                        f2.create_dataset('mouse', data=mouse)
                        f2.create_dataset('condition', data=condition)
                        f2.create_dataset('replay_locs', data=replayLocs)
                        f2.create_dataset('replay_score', data=replayScore)
                        f2.create_dataset('replay_jumpiness', data=replayJumpiness)
                        f2.create_dataset('replay_length', data=replayPortion)
                        f2.create_dataset('replay_slope', data=replaySlope)