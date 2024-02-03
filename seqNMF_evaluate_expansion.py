# %% Imports
import numpy as np
import os
import matplotlib.pyplot as plt
import yaml
import h5py
from scipy.stats import skew
from seqnmf import seqnmf
from utils.helperFunctions import load_data
from tqdm import tqdm

plt.style.use("plot_style.mplstyle")

# %% Params
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

np.random.seed(params["seed"])

# %% Load LT and REM data
mouse = 'pv1069'
condition = 'LTD1'
data_LT = load_data(mouse=mouse,
                    condition=condition,
                    state='wake',
                    params=params)
data_REMpost = load_data(mouse=mouse,
                    condition=condition,
                    state='REMpost',
                    params=params)
with h5py.File(os.path.join(params['path_to_output'],"neuron_selection", f'selected_neurons_{condition}_{mouse}.h5'),'r') as f:
    selected_neurons = f['place_cells'][()]

PC_activity = data_LT['binaryData'][:,selected_neurons]

# %% First, let's optimize 'L':
L_list = [25,50,75,125,150,175,200,225,250]
L_scores = np.zeros((len(L_list),params['K']))

for i, L in enumerate(tqdm(L_list)):
    W_LT, H_LT, _, _, _ = seqnmf(
        PC_activity.T,
        K=params['K'],
        L=params['L'],
        Lambda=params['Lambda'],
        max_iter=10
        )
    L_scores[i,:] = skew(H_LT, axis=1)

# %%  Plot results
plt.plot(L_list,L_scores)

# %% Create ~20 different compressed/streched templates by interpolation

# %% For each template, compute matrix multiplication to extract H

# %% Compute