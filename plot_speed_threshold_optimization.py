# %% Imports
import numpy as np
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import h5py
from scipy.ndimage import gaussian_filter1d
import yaml
from tqdm import tqdm
from utils.helperFunctions import load_data
from utils.helperFunctions import extract_seqReplay_score

plt.style.use("plot_style.mplstyle")

# %% Load parameters
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

np.random.seed(params["seed"])

# %% Load example statistics
data = load_data(mouse="pv1069", condition="LTD1", state="wake", params=params)

# %% Plot dist of time spent according to speed threshold
plt.figure(figsize=(4, 4))
speed_list = [2, 3, 4, 5, 6, 7, 8, 9, 10]
for i, speed_i in enumerate(speed_list):
    plt.subplot(3, 3, i + 1)
    plt.hist(data["position"][data["velocity"] > speed_i, 0], rasterized=True)
    plt.title(f"Thresh.: {speed_i} cm.s$^{-1}$")
    plt.axis("off")
plt.savefig("../../output_REM/speed_threshold_optimization.pdf")