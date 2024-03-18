# %% Imports
import h5py
import yaml
import os
from utils.helperFunctions import load_data
from pycaan.functions.decoding import bayesian_decode, temporal_bayesian_filter
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# %% Load parameters
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

#%%
condition = 'LTD1'
mouse = 'pv1069'

#%%
current_path = os.path.join(
            params["path_to_output"], "tuning", f"tuning_{condition}_{mouse}.h5"
        )
        
with h5py.File(current_path, "r") as f:
    tuning_curves = f["tuning_curves"][()]
    marginal_likelihood = f["marginal_likelihood"][()]
    # Use uniform prior
    occupancy = np.ones(tuning_curves.shape[1]) / tuning_curves.shape[1]

data_REMpost = load_data(mouse, condition, "REMpost", params)

#%%
REMpost_posterior_probs, _ = bayesian_decode(
                tuning_curves,
                occupancy,
                marginal_likelihood,
                data_REMpost["binaryData"],
            )
#%%
plt.imshow(REMpost_posterior_probs.T,aspect='auto',
           interpolation='none')

#%%
windowSize = params['filtWindowSize']
smoothed_posteriors = np.zeros((REMpost_posterior_probs.shape))*np.nan

#%%
posterior_probs = np.concatenate((
        np.zeros((int(windowSize/2),REMpost_posterior_probs.shape[1])),
        REMpost_posterior_probs,
        np.zeros((int(windowSize/2),REMpost_posterior_probs.shape[1]))
        ))

#%%
currentWindowIdx = np.arange(windowSize)

for i in range(len(smoothed_posteriors)):
    bayesian_step_prob = posterior_probs[currentWindowIdx]
    smoothed_posteriors[i,:] = np.expm1(np.nansum(np.log1p(bayesian_step_prob),axis=0)) # This should be used instead of simple product to avoid numerical underflow
    smoothed_posteriors[i,:] = smoothed_posteriors[i,:]/np.nansum(smoothed_posteriors[i,:]) # Normalize into a probability distribution
    currentWindowIdx+=1 # Step forward



# %%
filtered_posteriors = temporal_bayesian_filter(
                    REMpost_posterior_probs,
                    params['filtWindowSize']
                    )
# %%
