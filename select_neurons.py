#%%
import numpy as np
import os
import h5py
import yaml
from tqdm import tqdm

#%%
results_dir = '../../output_REM/'
resultsList=os.listdir(results_dir)

for file_name in tqdm(resultsList):
    if file_name.startswith('tuning_') and file_name.endswith('.h5'): # Only include seqResults, not replay results
        h5_file = h5py.File(os.path.join(results_dir,file_name))
        # Remvove non-significant neurons
        
        p_values = h5_file['p_value'][()]
        indices = np.arange(len(p_values))
        significant_neurons = h5_file['info'][()]

        # Sort neurons to get top-k neurons

        # Save selected neurons

        
        # Close files
        h5_file.close()