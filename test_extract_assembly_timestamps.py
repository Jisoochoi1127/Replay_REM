# %% Imports
import os
import numpy as np
import pandas as pd
import yaml
import h5py
from tqdm import tqdm

#%% Load params
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

#%%test
h5_file = h5py.File('/Users/guillaumeetter/Documents/output_REM/assembly/assembly_LTD1_pv1043.h5')

#%% 
results_dir = params['path_to_output']+"/assembly"
resultsList = os.listdir(results_dir)

data_list = []

for file_name in tqdm(resultsList):
    if (
        file_name.startswith("assembly_")
        and file_name.endswith(".h5")
        and "pv1254" not in file_name # Exclude pv1254
    ):
        h5_file = h5py.File(os.path.join(results_dir, file_name))
        for i in range(len(h5_file["replay_locs"][()])):
            data_list.append(  # This will create one list entry per cell
                {
                    "eventID": i,
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    
                }
            )
        
        # Close files
        h5_file.close()

df = pd.DataFrame(data_list)