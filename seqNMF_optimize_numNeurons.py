# %% Imports
import numpy as np
import os
import matplotlib.pyplot as plt
import yaml
import h5py
from scipy.stats import skew
from seqnmf import seqnmf
import itertools
from utils.helperFunctions import load_data, extract_seqReplay_score
from tqdm import tqdm

plt.style.use("plot_style.mplstyle")

# %% Params
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

np.random.seed(params["seed"])

# %% Load LT and REM data
condition = "LTD1"
mouse_ID = "pv1069" # "pv1254", "pv1043", "pv1060"
numNeurons_list = [8, 16, 32, 64, 128, 256, 512]

data_LT = load_data(mouse=mouse_ID, condition=condition, state="wake", params=params)
data_REMpost = load_data(
    mouse=mouse_ID, condition=condition, state="REMpost", params=params
)

# Open tuning data
results_dir = params["path_to_output"] + "/tuning"
h5_file = h5py.File(os.path.join(results_dir, f'tuning_{condition}_{mouse_ID}.h5'))

# Select significant neurons
p_values = h5_file["p_value"][()]
indices = np.arange(len(p_values))

sig_info = h5_file["info"][p_values < 0.05]
sig_indices = indices[p_values < 0.05]
nonSig_indices = indices[p_values > 0.05]

# Sort neurons to get top-k neurons
sorted_sig_info = np.argsort(
    sig_info * -1
)  # -1 to find descending idx, instead of ascending

#%%
# Establish list of neurons to use
numNeurons_scores = np.zeros((len(numNeurons_list), params["K"])) # array of size numNeuronslist, numSeqs (usually 2)

for i, numNeurons in enumerate(tqdm(numNeurons_list)):
    params['numNeurons'] = numNeurons
    selected_neurons = sig_indices[sorted_sig_info][0 : numNeurons]
    seqReplay_scores, seqReplay_pvalues, seqReplay_locs, W_ref = extract_seqReplay_score(
        data_LT["binaryData"][:, selected_neurons],
        data_REMpost["binaryData"][:, selected_neurons],
        params,
    )

    numNeurons_scores[i, 0] = len(seqReplay_locs[0])
    numNeurons_scores[i, 1] = len(seqReplay_locs[1])

# %% Save data
with h5py.File(os.path.join(params["path_to_output"], f"{mouse_ID}_optimal_numNeurons.h5"), "w") as f:
    f.create_dataset("numNeurons_scores", data=numNeurons_scores)

# %%
