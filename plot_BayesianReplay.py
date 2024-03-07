# %% Imports
import numpy as np
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import h5py
import yaml
from tqdm import tqdm
from utils.helperFunctions import load_data
from pycaan.functions.decoding import bayesian_decode

plt.style.use("plot_style.mplstyle")

# %% Load parameters
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

np.random.seed(params["seed"])

# %% Select example recording
mouse = 'pv1069'
condition = 'LTD1'

# %% Load data
data_wake = load_data(mouse, condition, 'wake', params)

#%% Load tuning curves for that session
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

#%% Compute posterior probs during wake
posterior_probs, map = bayesian_decode(tuning_curves,
                                       occupancy,
                                       marginal_likelihood,
                                       data_wake['binaryData'])
# %% Plot posterior probabilities during wake
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

# %% Load posterior probabilities during REM
with h5py.File(
                os.path.join(
                    params["path_to_output"],
                    "posterior_probs",
                    f"posterior_probs_{condition}_{mouse}.h5",
                ),
                "r",
            ) as f:
    REMpost_posterior_probs = f["REMpost_posterior_probs"][()]

# %% Load replay info
with h5py.File(
            os.path.join(
                params["path_to_output"],
                "bayesian_replay",
                f"bayesian_replay_{condition}_{mouse}_REMpost.h5",
            ),
            "r",
        ) as f:
    replayLocs = f['replay_locs'][()]
    replayScore = f['replay_score'][()]
    replayJumpiness = f['replay_jumpiness'][()]
    replayLength = f['replay_length'][()]

# %% Identify top replay events
sortedEventIdx = np.argsort(-1*replayScore) # Reverse to get descending order

# %% Plot significant replay revents
plt.figure(figsize=(4,1))
# Plot posterior probs
plt.imshow(REMpost_posterior_probs.T,aspect='auto',vmax=.05,interpolation='none')

for event in replayLocs:
    plt.fill_between([event, event+int(params['windowSize'])],
            [REMpost_posterior_probs.shape[1],REMpost_posterior_probs.shape[1]],
            facecolor='w',
            alpha=.2,
            )
plt.xlim(4500,5000)

# %% TODO overlay linear fit?
    


# %% Load all sessions
    

# %% Plot results
    



# %% Compute stats