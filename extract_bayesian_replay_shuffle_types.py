# %% Imports
import h5py
import yaml
import os
from tqdm import tqdm
from utils.bayesian_replay import extract_linear_replay_shuffle_types

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
                if not os.path.exists(os.path.join(params['path_to_output'],"controlled_bayesian_replay", f'controlled_bayesian_replay_{condition}_{mouse}_{state}.h5')):
                # Load precomputed tuning curves and accessory data
                    posterior_probs = f[f'{state}_posterior_probs'][()]
                    posterior_control = f[f'control_{state}_posterior_probs'][()]
                    
                    (
                        replayLocs_P,
                        replayScore_P,
                        replayJumpiness_P,
                        replayPortion_P,
                        replaySlope_P,
                        pvalue_P,
                        replayLocs_T,
                        replayScore_T,
                        replayJumpiness_T,
                        replayPortion_T,
                        replaySlope_T,
                        pvalue_T,
                        replayLocs_PT,
                        replayScore_PT,
                        replayJumpiness_PT,
                        replayPortion_PT,
                        replaySlope_PT,
                        pvalue_PT,
                    ) = extract_linear_replay_shuffle_types(posterior_probs, params)
                    
                    (
                        control_replayLocs_P,
                        control_replayScore_P,
                        control_replayJumpiness_P,
                        control_replayPortion_P,
                        control_replaySlope_P,
                        control_pvalue_P,
                        control_replayLocs_T,
                        control_replayScore_T,
                        control_replayJumpiness_T,
                        control_replayPortion_T,
                        control_replaySlope_T,
                        control_pvalue_T,
                        control_replayLocs_PT,
                        control_replayScore_PT,
                        control_replayJumpiness_PT,
                        control_replayPortion_PT,
                        control_replaySlope_PT,
                        control_pvalue_PT,
                    ) = extract_linear_replay_shuffle_types(posterior_control, params)
                    
                    # Save results
                    with h5py.File(os.path.join(params['path_to_output'],"controlled_bayesian_replay", f'controlled_bayesian_replay_{condition}_{mouse}_{state}.h5'),'w') as f2:
                        f2.create_dataset('mouse', data=mouse)
                        f2.create_dataset('condition', data=condition)
                        
                        # Position
                        f2.create_dataset('replayLocs_P', data=replayLocs_P)
                        f2.create_dataset('control_replayLocs_P', data=control_replayLocs_P)
                        f2.create_dataset('replayScore_P', data=replayScore_P)
                        f2.create_dataset('control_replayScore_P', data=control_replayScore_P)
                        f2.create_dataset('replayJumpiness_P', data=replayJumpiness_P)
                        f2.create_dataset('control_replayJumpiness_P', data=control_replayJumpiness_P)
                        f2.create_dataset('replayPortion_P', data=replayPortion_P)
                        f2.create_dataset('control_replayPortion_P', data=control_replayPortion_P)
                        f2.create_dataset('replaySlope_P', data=replaySlope_P)
                        f2.create_dataset('control_replaySlope_P', data=control_replaySlope_P)
                        f2.create_dataset('pvalue_P', data=pvalue_P)
                        f2.create_dataset('control_pvalue_P', data=control_pvalue_P)

                        #Time
                        f2.create_dataset('replayLocs_T', data=replayLocs_T)
                        f2.create_dataset('control_replayLocs_T', data=control_replayLocs_T)
                        f2.create_dataset('replayScore_T', data=replayScore_T)
                        f2.create_dataset('control_replayScore_T', data=control_replayScore_T)
                        f2.create_dataset('replayJumpiness_T', data=replayJumpiness_T)
                        f2.create_dataset('control_replayJumpiness_T', data=control_replayJumpiness_T)
                        f2.create_dataset('replayPortion_T', data=replayPortion_T)
                        f2.create_dataset('control_replayPortion_T', data=control_replayPortion_T)
                        f2.create_dataset('replaySlope_T', data=replaySlope_T)
                        f2.create_dataset('control_replaySlope_T', data=control_replaySlope_T)
                        f2.create_dataset('pvalue_T', data=pvalue_T)
                        f2.create_dataset('control_pvalue_T', data=control_pvalue_T)

                        # Position & time
                        f2.create_dataset('replayLocs_PT', data=replayLocs_PT)
                        f2.create_dataset('control_replayLocs_PT', data=control_replayLocs_PT)
                        f2.create_dataset('replayScore_PT', data=replayScore_PT)
                        f2.create_dataset('control_replayScore_PT', data=control_replayScore_PT)
                        f2.create_dataset('replayJumpiness_PT', data=replayJumpiness_PT)
                        f2.create_dataset('control_replayJumpiness_PT', data=control_replayJumpiness_PT)
                        f2.create_dataset('replayPortion_PT', data=replayPortion_PT)
                        f2.create_dataset('control_replayPortion_PT', data=control_replayPortion_PT)
                        f2.create_dataset('replaySlope_PT', data=replaySlope_PT)
                        f2.create_dataset('control_replaySlope_PT', data=control_replaySlope_PT)
                        f2.create_dataset('pvalue_PT', data=pvalue_PT)
                        f2.create_dataset('control_pvalue_PT', data=control_pvalue_PT)