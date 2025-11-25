# %% Imports
import h5py
import yaml
import os
from utils.helperFunctions import load_data
from utils.bayesian_replay import extract_linear_replay_shuffle_types
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# %% Load parameters
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

#%%
condition = 'LTD1'
mouse = 'pv1069'

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
                posterior_probs = f[f'{state}_posterior_probs'][()]

#%%

(replayLocs_P,
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
 pvalue_P
 ) = extract_linear_replay_shuffle_types(posterior_probs, params)

#%%
plt.imshow(posterior_probs, interpolation=None, aspect='auto')

#%%
nan_locs = np.max(np.isnan(posterior_probs),axis=1)
actual_map = (np.argmax(posterior_probs,axis=1)+params['spatialBinSize']/2)*params['spatialBinSize']
actual_map[np.isnan(nan_locs)] = np.nan # place back the nans

#%%
plt.plot(actual_map)
plt.xlim(0,1000)
# %%
