#%%
import numpy as np
import os
import h5py
import yaml
from tqdm import tqdm

#%% Load parameters
with open('params.yaml','r') as file:
    params = yaml.full_load(file)

results_dir = params['path_to_output']+'/tuning'
resultsList=os.listdir(results_dir)

for file_name in tqdm(resultsList):
    if file_name.startswith('tuning_') and file_name.endswith('.h5'): # Only include seqResults, not replay results
        h5_file = h5py.File(os.path.join(results_dir,file_name))
        mouse = h5_file['mouse'][()].decode("utf-8")
        condition = h5_file['condition'][()].decode("utf-8")
        marginal_likelihood = h5_file['marginal_likelihood'][()] 
        
        # Select significant neurons
        p_values = h5_file['p_value'][()]
        indices = np.arange(len(p_values))
        
        sig_info = h5_file['info'][p_values<0.05]
        sig_indices = indices[p_values<0.05]
        nonSig_indices = indices[p_values>0.05]

        # Sort neurons to get top-k neurons
        sorted_sig_info = np.argsort(sig_info*-1) #-1 to find descending idx, instead of ascending
        sorted_info = np.argsort(h5_file['info']) # find ascending idx, for worst place cells
        sorted_activity = np.argsort(marginal_likelihood*-1) # Find the most active neurons
        r_sorted_activity = np.argsort(marginal_likelihood) # Find the least active neurons 
        
        # Save selected neurons
        try:
            place_cells = sig_indices[sorted_sig_info][0:params['numNeurons']]
            non_place_cells = indices[sorted_info][0:params['numNeurons']]
            most_active = indices[sorted_activity][0:params['numNeurons']]
            least_active = indices[r_sorted_activity][0:params['numNeurons']]

            with h5py.File(os.path.join(params['path_to_output'], 'neuron_selection', f'selected_neurons_{condition}_{mouse}.h5'),'w') as f2:
                f2.create_dataset('mouse', data=mouse)
                f2.create_dataset('condition', data=condition)
                f2.create_dataset('place_cells', data=place_cells)
                f2.create_dataset('non_place_cells', data=non_place_cells)
                f2.create_dataset('most_active', data=most_active)
                f2.create_dataset('least_active', data=least_active)
        except:
            print('There might not be enough selected neurons. Try reducing the number in params')

        # Close files
        h5_file.close()