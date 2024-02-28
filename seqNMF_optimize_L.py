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
mouse_list = ["pv1043", "pv1060", "pv1069", "pv1254"]
L_list = [25, 50, 75, 125, 150, 175, 200, 225, 250, 300]

for mouse in mouse_list:
    data_LT = load_data(mouse=mouse, condition=condition, state="wake", params=params)
    data_REMpost = load_data(
        mouse=mouse, condition=condition, state="REMpost", params=params
    )

    with h5py.File(
        os.path.join(
            params["path_to_output"],
            "neuron_selection",
            f"selected_neurons_{condition}_{mouse}.h5",
        ),
        "r",
    ) as f:
        selected_neurons = f["place_cells"][()]

    L_scores = np.zeros((len(L_list), params["K"]))
    for i, L in enumerate(tqdm(L_list)):
        params["L"] = L  # Override parameters
        seqReplay_scores, seqReplay_pvalues, seqReplay_locs = extract_seqReplay_score(
            data_LT["binaryData"][:, selected_neurons],
            data_REMpost["binaryData"][:, selected_neurons],
            params,
        )

        L_scores[i, 0] = len(seqReplay_locs[0])
        L_scores[i, 1] = len(seqReplay_locs[1])


    # %% Save data
    with h5py.File(os.path.join(params["path_to_output"], f"{mouse}_optimal_L.h5"), "w") as f:
        f.create_dataset("L_scores", data=L_scores)
