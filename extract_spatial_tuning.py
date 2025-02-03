# %% Imports
import numpy as np
import os
from tqdm import tqdm
import itertools

import yaml
import h5py
from utils.helperFunctions import load_data, extract_equal_samples
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


for condition, mouse in tqdm(
        list(itertools.product(condition_list, mouse_list)),
        total=len(condition_list) * len(mouse_list),
    ):
    try:
        if (
            not os.path.exists(
                os.path.join(
                    params["path_to_output"], "tuning", f"tuning_{condition}_{mouse}.h5"
                )
            )
            or params["overwrite_mode"] == "always"
        ):
            if os.path.exists(os.path.join(params['path_to_dataset'],mouse,condition)):
                # Load data
                data = load_data(mouse, condition, "wake", params)

                if params['equalize_sampling']:
                    # Extract equal samples of running periods
                    bin_vec = np.arange(
                    0, 100 + 25, 25 # Divide in fourths
                )

                    timestamps = extract_equal_samples(data['position'][:,0], 
                                    bin_vec,
                                    data['running_ts'])

                else:
                    # Only consider running periods
                    timestamps=data['running_ts']

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
                    data["binaryData"], data["position"][:, 0], timestamps, bin_vec
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
                        "equal_tuning",
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
        print("Could not process this recording")
# %%
