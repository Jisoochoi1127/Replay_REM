# %% Imports
import numpy as np
import os
import matplotlib.pyplot as plt
import yaml
import h5py
from scipy.stats import skew
from seqnmf import seqnmf
from utils.helperFunctions import load_data, extract_seqReplay_score
from tqdm import tqdm

plt.style.use("plot_style.mplstyle")

# %% Params
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

np.random.seed(params["seed"])

# %% Load LT and REM data
mouse = 'pv1060'
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

# %% Create ~20 different compressed/streched templates by interpolation
W_LT, H_LT, _, _, _ = seqnmf(
    PC_activity.T,
    K=params['K'],
    L=params['L'],
    Lambda=params['Lambda'],
    max_iter=10
    )

#%% Test temp
from scipy.ndimage.interpolation import zoom

# %% For each template, compute matrix multiplication to extract H
expansionFactorList = [.5,.75,1,1.5,2]


# %% Compute