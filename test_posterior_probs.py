#%%
from utils.helperFunctions import load_data
from pycaan.functions.decoding import bayesian_decode
import h5py
import os
import yaml
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("plot_style.mplstyle")

#%% Load parameters
with open('params.yaml','r') as file:
    params = yaml.full_load(file)

#%%
with h5py.File(
    os.path.join(
        params["path_to_output"],
        'tuning',
        "tuning_LTD1_pv1069.h5"
        ), 'r'
        ) as f:
    tuning_curves = f['tuning_curves'][()]
    marginal_likelihood = f['marginal_likelihood'][()]
    # Use uniform prior
    occupancy = np.ones(tuning_curves.shape[1])/tuning_curves.shape[1]
        
# Load selected neurons
with h5py.File(
    os.path.join(
        params['path_to_output'],
        'neuron_selection',
        "selected_neurons_LTD1_pv1069.h5"
        ), 'r'
        ) as f:
    selected_neurons = f['place_cells'][()]


# %%
mouse = 'pv1069'
condition = 'LTD1'
data_wake = load_data(mouse, condition, 'wake', params)

# %%
posterior_probs, map = bayesian_decode(tuning_curves, occupancy, marginal_likelihood, data_wake['binaryData'])
# %%
plt.figure(figsize=(4,2))
plt.subplot(211)
plt.imshow(
    posterior_probs.T,
    interpolation='none',
    aspect='auto',
    # vmin=0,
    vmax=.25,
    origin='lower'
)

plt.xlim(0,2500)
plt.xticks(np.arange(0,2500,600),[])
plt.yticks([0,40],[0,100])
plt.ylabel('Position (cm)')
# plt.colorbar(label='Posterior probability')

plt.subplot(212)
plt.plot(data_wake['position'][:,0])
plt.xlim(0,2500)
plt.xticks(np.arange(0,2500,600),np.arange(0,2500/60,10))
plt.xlabel('Time (s)')
plt.ylabel('Position (cm)')
plt.savefig("../../output_REM/posterior_probs.pdf")

# %%
