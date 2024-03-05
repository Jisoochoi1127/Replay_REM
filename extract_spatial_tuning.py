# %% Imports
import numpy as np
import os
from tqdm import tqdm
import itertools

import yaml
import h5py
from utils.helperFunctions import load_data
from pycaan.functions.tuning import extract_tuning
from pycaan.functions.metrics import (
    extract_firing_properties,
    extract_total_distance_travelled,
)

# %% Load parameters
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

np.random.seed(params["seed"])

# %% Define conditions, dataset
states_list = ["REMpre", "wake", "REMpost"]
condition_list = ["LTD1", "LTD5", "HATD1", "HATD5", "HATDSwitch"]
mouse_list = ["pv1043", "pv1060", "pv1069", "pv1191", "pv1192", "pv1252", "pv1254"]

# %% First, measure 'sequenceness' in each individuate session
for condition, mouse in tqdm(
    list(itertools.product(condition_list, mouse_list)),
    total=len(condition_list) * len(mouse_list),
):
    if (
        not os.path.exists(
            os.path.join(
                params["path_to_output"], "tuning", f"tuning_{condition}_{mouse}.h5"
            )
        )
        or params["overwrite_mode"] == "always"
    ):
        try:
            # Load data
            data = load_data(mouse, condition, "wake", params)

            # Extract tuning info/tuning
            bin_vec = np.arange(
                0, 100 + params["spatialBinSize"], params["spatialBinSize"]
            )
            (
                info,
                p_value,
                occupancy_frames,
                active_frames_in_bin,
                tuning_curves,
                peak_loc,
                peak_val,
            ) = extract_tuning(
                data["binaryData"], data["position"][:, 0], data["running_ts"], bin_vec
            )
            total_distance_travelled = extract_total_distance_travelled(
                data["position"]
            )
            marginal_likelihood, trans_prob = extract_firing_properties(
                data["binaryData"]
            )

            with h5py.File(
                os.path.join(
                    params["path_to_output"],
                    "tuning_curves",
                    f"tuning_{condition}_{mouse}.h5",
                ),
                "w",
            ) as f:
                f.create_dataset("mouse", data=mouse)
                f.create_dataset("condition", data=condition)
                f.create_dataset("info", data=info)
                f.create_dataset("p_value", data=p_value)
                f.create_dataset("occupancy_frames", data=occupancy_frames, dtype=int)
                f.create_dataset(
                    "active_frames_in_bin", data=active_frames_in_bin, dtype=int
                )
                f.create_dataset("tuning_curves", data=tuning_curves)
                f.create_dataset("peak_loc", data=peak_loc)
                f.create_dataset("peak_val", data=peak_val)
                f.create_dataset("bins", data=bin_vec)
                f.create_dataset("marginal_likelihood", data=marginal_likelihood)
                f.create_dataset("trans_prob", data=trans_prob)
                f.create_dataset(
                    "total_distance_travelled", data=float(total_distance_travelled)
                )

        except:
            print("Missing session")
# %%
