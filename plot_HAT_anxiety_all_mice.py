# %% Imports
import numpy as np
from tqdm import tqdm
import itertools
import yaml
from utils.helperFunctions import load_data
import matplotlib.pyplot as plt

plt.style.use("plot_style.mplstyle")

# %% Load parameters
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

np.random.seed(params["seed"])

# %% Define conditions, dataset
condition_list = ["HATD1", "HATD5"]
mouse_list = ["pv1043", "pv1060", "pv1069", "pv1191", "pv1192", "pv1252", "pv1254"]

# %%
ct = 1
plt.figure(figsize=(len(mouse_list), len(condition_list)))
for condition, mouse in tqdm(
    list(itertools.product(condition_list, mouse_list)),
    total=len(condition_list) * len(mouse_list),
):
    try:
        plt.subplot(3, len(mouse_list), ct)

        data = load_data(mouse=mouse, condition=condition, state="wake", params=params)
        plt.plot(data["position"][:, 0], data["caTime"], "C4", linewidth=0.2)
        plt.plot([0, 100], [0, 0], "k", linewidth=0.6)
        plt.plot([50, 50], [0, 100], "k", linewidth=0.6)
        # plt.plot([0, 0], [0, 15], "k", linewidth=0.6)
        # plt.plot([100, 100], [0, 15], "k", linewidth=0.6)

        plt.xlim(-1, 101)
        # plt.ylim(-80, 20)
        plt.axis("off")
        plt.title(f"{mouse}-{condition}")
        ct += 1
    except:
        print("Missing session")
        plt.axis("off")

plt.savefig("../../output_REM/HAT_exploration_all_mice.pdf")
# %%
