# %% Imports
import h5py
import yaml
import os
from tqdm import tqdm
import numpy as np
from utils.bayesian_replay import extract_raw_linear_replay
from utils.helperFunctions import load_data

#%% Load parameters
with open('params.yaml','r') as file:
    params = yaml.full_load(file)

results_dir = params['path_to_output']+'/neuron_selection' 
resultsList=os.listdir(results_dir)

#%%
for file_name in tqdm(resultsList):
    if file_name.endswith(".h5"):
        condition, mouse = file_name.split("_")[-2:]
        current_name = condition + "_" + mouse #get task_mouse.h5
        mouse = mouse[:-3] # Remove "".h5"
        
        # Load tuning data
        tuning_file = h5py.File(os.path.join(
            params['path_to_output'],
            "tuning",
            "tuning_"+current_name)
            )
        
        # Load selected neurons
        selected_neurons_file = h5py.File(os.path.join(
            params['path_to_output'],
            "neuron_selection",
            "selected_neurons_"+current_name)
            )
        selected_neurons = selected_neurons_file['place_cells'][()]

        # Select neurons and sort peak_loc
        selected_peak_loc = tuning_file['peak_loc'][()][selected_neurons,0]
        sorting_index = np.argsort(selected_peak_loc)

        # Close files
        tuning_file.close()
        selected_neurons_file.close()

        for state in ['REMpre', 'REMpost']:
            if not os.path.exists(os.path.join(params['path_to_output'],"raw_linear_replay", f'raw_linear_replay_{condition}_{mouse}_{state}.h5')):
                # load recording
                sleep_data = load_data(mouse, condition, state, params)

                # Select and sort
                selected_data = sleep_data['binaryData'][:,selected_neurons]
                sorted_binary = selected_data[:,sorting_index]

                # Extract replay events
                replayLocs, replayScore, replayJumpiness, replayPortion, replaySlope = extract_raw_linear_replay(sorted_binary, params)
            
                # Save results
                with h5py.File(os.path.join(params['path_to_output'],"raw_linear_replay", f'raw_linear_replay_{condition}_{mouse}_{state}.h5'),'w') as f2:
                    f2.create_dataset('mouse', data=mouse)
                    f2.create_dataset('condition', data=condition)
                    f2.create_dataset('replay_locs', data=replayLocs)
                    f2.create_dataset('replay_score', data=replayScore)
                    f2.create_dataset('replay_jumpiness', data=replayJumpiness)
                    f2.create_dataset('replay_length', data=replayPortion)
                    f2.create_dataset('replay_slope', data=replaySlope)