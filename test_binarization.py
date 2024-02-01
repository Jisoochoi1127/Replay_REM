# %% Imports
import numpy as np
import yaml
from utils.helperFunctions import load_data
import matplotlib.pyplot as plt
from utils.dataloaders import load_data as pycaan_load
from pycaan.functions.signal_processing import binarize_ca_traces, preprocess_data

plt.style.use("plot_style.mplstyle")

# %% Load parameters
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

np.random.seed(params["seed"])

# %% Load data
mouse = 'pv1069'
condition = 'LTD1'
data = load_data(mouse=mouse,
                 condition=condition,
                 state="wake",
                 params=params)

# %%
cell2plot=300
plt.figure()
plt.plot(data['rawData'][:,cell2plot])
plt.plot(data['binaryData'][:,cell2plot]*-1)

# %% [markdown]
# # Conclusions
# For wakefulness, it seems that the binarization parameters are fine (correct sensitivity for transient detection)
# Question remains: is it also the case for REM? So far I only have access to binary data
