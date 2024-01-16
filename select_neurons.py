#%%
import numpy as np
import os
import h5py
import yaml
from tqdm import tqdm

#%%
results_dir = '../../output_REM/'
resultsList=os.listdir(results_dir)

#%% Load parameters
with open('params.yaml','r') as file:
    params = yaml.full_load(file)

for file_name in tqdm(resultsList):
    if file_name.startswith('tuning_') and file_name.endswith('.h5'): # Only include seqResults, not replay results
        h5_file = h5py.File(os.path.join(results_dir,file_name))
        mouse = h5_file['mouse'][()].decode("utf-8")
        condition = h5_file['condition'][()].decode("utf-8")
        
        # Select significant neurons
        p_values = h5_file['p_value'][()]
        indices = np.arange(len(p_values))
        
        sig_info = h5_file['info'][p_values<0.05]
        sig_indices = indices[p_values<0.05]

        # Sort neurons to get top-k neurons
        sorted_info = np.argsort(sig_info*-1) #-1 to find descending idx, instead of ascending

        # Save selected neurons
        try:
            selected_neurons = sig_indices[sorted_info][0:params['numNeurons']]
            with h5py.File(os.path.join(params['path_to_output'],f'selected_neurons_{condition}_{mouse}.h5'),'w') as f2:
                f2.create_dataset('mouse', data=mouse)
                f2.create_dataset('condition', data=condition)
                f2.create_dataset('selected_neurons', data=selected_neurons)
        except:
            print('There might not be enough significant neurons. Try reducing the number in params')

        # Close files
        h5_file.close()