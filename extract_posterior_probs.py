# %% Imports
import h5py
import yaml
import os
from tqdm import tqdm
import itertools
from utils.helperFunctions import load_data
from pycaan.functions.decoding import bayesian_decode, temporal_bayesian_filter
import numpy as np

# %% Load parameters
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

# %% Define conditions, dataset
condition_list = ["LTD1", "LTD5", "HATD1", "HATD5", "HATDSwitch"]
mouse_list = ["pv1043", "pv1060", "pv1069", "pv1191", "pv1192", "pv1252", "pv1254"]

# %% For each mouse and condition
for condition, mouse in tqdm(
    list(itertools.product(condition_list, mouse_list)),
    total=len(condition_list) * len(mouse_list),
):
    if not os.path.exists(
        os.path.join(
            params["path_to_output"],
            "equal_posterior_probs",
            f"posterior_probs_{condition}_{mouse}.h5",
        )
    ):
        # Load precomputed tuning curves and accessory data
        current_path = os.path.join(
            params["path_to_output"],
            "equal_tuning",
            f"tuning_{condition}_{mouse}.h5",  # GE 20250130: changed from tuning to equal_tuning
        )
        if os.path.exists(current_path):
            with h5py.File(current_path, "r") as f:
                tuning_curves = f["tuning_curves"][()]
                marginal_likelihood = f["marginal_likelihood"][()]
                # Use uniform prior
                occupancy = np.ones(tuning_curves.shape[1]) / tuning_curves.shape[1]

            # Load selected neurons
            with h5py.File(
                os.path.join(
                    params["path_to_output"],
                    "neuron_selection",
                    f"selected_neurons_{condition}_{mouse}.h5",
                ),
                "r",
            ) as f:
                selected_neurons = f["place_cells"][()]

            # Pick different mouse for control
            keep_going = False
            while (
                keep_going == False
            ):  # ensure we are looking at two distinct mice for control
                control_mouse = np.random.choice(mouse_list)
                control_condition = np.random.choice(condition_list) # necessary if initial condition not available
                control_path = os.path.join(
                    params["path_to_output"],
                    "equal_tuning",
                    f"tuning_{control_condition}_{control_mouse}.h5",  # GE 20250130: changed from tuning to equal_tuning
                )
                if os.path.exists(control_path) and control_mouse != mouse:
                    keep_going = True

            # Load precomputed tuning curves and accessory data
            with h5py.File(control_path, "r") as f:
                control_tuning_curves = f["tuning_curves"][()]
                control_marginal_likelihood = f["marginal_likelihood"][()]
                # Use uniform prior
                control_occupancy = np.ones(control_tuning_curves.shape[1]) / control_tuning_curves.shape[1]

            # Load selected neurons
            with h5py.File(
                os.path.join(
                    params["path_to_output"],
                    "neuron_selection",
                    f"selected_neurons_{control_condition}_{control_mouse}.h5",
                ),
                "r",
            ) as f:
                control_selected_neurons = f["place_cells"][()]

            # Load REMpre data
            data_REMpre = load_data(mouse, condition, "REMpre", params)

            # Load REMpost data
            data_REMpost = load_data(mouse, condition, "REMpost", params)

            # Extract posteriors on REMpre and save results
            REMpre_posterior_probs, _ = bayesian_decode(
                tuning_curves[selected_neurons],
                occupancy,
                marginal_likelihood[selected_neurons],
                data_REMpre["binaryData"][:, selected_neurons],
            )

            # Extract posterior probs on REMpost and save results
            REMpost_posterior_probs, _ = bayesian_decode(
                tuning_curves[selected_neurons],
                occupancy,
                marginal_likelihood[selected_neurons],
                data_REMpost["binaryData"][:, selected_neurons],
            )

            # Decode REM using tuning curves from another mouse
            control_REMpre_posterior_probs, _ = bayesian_decode(
                control_tuning_curves[control_selected_neurons],
                control_occupancy,
                control_marginal_likelihood[control_selected_neurons],
                data_REMpre["binaryData"][:, selected_neurons],
            )

            # Extract posterior probs on REMpost and save results
            control_REMpost_posterior_probs, _ = bayesian_decode(
                control_tuning_curves[control_selected_neurons],
                control_occupancy,
                control_marginal_likelihood[control_selected_neurons],
                data_REMpost["binaryData"][:, selected_neurons],
            )

            if params["filtWindowSize"] > 0:
                REMpre_posterior_probs = temporal_bayesian_filter(
                    REMpre_posterior_probs, params["filtWindowSize"]
                )
                REMpost_posterior_probs = temporal_bayesian_filter(
                    REMpost_posterior_probs, params["filtWindowSize"]
                )
                control_REMpre_posterior_probs = temporal_bayesian_filter(
                    control_REMpre_posterior_probs, params["filtWindowSize"]
                )
                control_REMpost_posterior_probs = temporal_bayesian_filter(
                    control_REMpost_posterior_probs, params["filtWindowSize"]
                )

            # Save results
            with h5py.File(
                os.path.join(
                    params["path_to_output"],
                    "equal_posterior_probs",  # GE 20250130: changed from posterior_probs to equal_posterior_probs
                    f"posterior_probs_{condition}_{mouse}.h5",
                ),
                "w",
            ) as f:
                f.create_dataset("mouse", data=mouse)
                f.create_dataset("condition", data=condition)
                f.create_dataset("REMpre_posterior_probs", data=REMpre_posterior_probs)
                f.create_dataset(
                    "REMpost_posterior_probs", data=REMpost_posterior_probs
                )
                f.create_dataset(
                    "control_REMpre_posterior_probs",
                    data=control_REMpre_posterior_probs,
                )
                f.create_dataset(
                    "control_REMpost_posterior_probs",
                    data=control_REMpost_posterior_probs,
                )
# %%
