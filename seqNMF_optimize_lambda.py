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
lambda_list = [10e-6,10e-5,10e-4,10e-3,10e-2,10e-1,1,10,10e2,10e3,10e4]

params['maxIters'] = 2
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

    PC_activity = data_LT["binaryData"][:, selected_neurons]

    lambda_scores = np.zeros((len(lambda_list), params["K"]))
    for i, Lambda in enumerate(tqdm(lambda_list)):
        params["Lambda"] = Lambda  # Override parameters
        seqReplay_scores, seqReplay_pvalues, seqReplay_locs, _ = extract_seqReplay_score(
            data_LT["binaryData"][:, selected_neurons],
            data_REMpost["binaryData"][:, selected_neurons],
            params,
        )

        lambda_scores[i, 0] = len(seqReplay_locs[0])
        lambda_scores[i, 1] = len(seqReplay_locs[1])

    with h5py.File(os.path.join(params["path_to_output"], f"{mouse}_optimal_lambda.h5"), "w") as f:
        f.create_dataset("L_scores", data=lambda_scores) #TODO rename to lambda
