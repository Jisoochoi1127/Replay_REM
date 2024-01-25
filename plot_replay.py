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

plt.style.use("plot_style.mplstyle")
plt.ion()

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
# %%
with h5py.File("../../output_REM/tuning/tuning_LTD1_pv1069.h5") as f:
    info = f["info"][()]
    p_values = f["p_value"][()]
    peak_loc = f["peak_loc"][()][:, 0]
    peak_val = f["peak_val"][()]
    tuning_curves = f["tuning_curves"][()]

# %% Select high information neurons
selected_neurons = np.where((p_values < 0.05) & (peak_val > 0.1))[0]
sorted_index = np.argsort(peak_loc[selected_neurons])

# %% Filter traces
smoothed_curves = np.zeros((len(sorted_index), tuning_curves.shape[1]))
for neuron_i in range(len(sorted_index)):
    smoothed_curves[neuron_i, :] = gaussian_filter1d(
        tuning_curves[selected_neurons[sorted_index]][neuron_i], sigma=1.5
    )
smoothed_curves /= np.max(smoothed_curves, axis=0)
# %% Plot example exploration and activity
plt.figure(figsize=(1, 2.5))
plt.subplot(8, 1, 4)
plt.plot(data["position"][:, 0], data["position"][:, 1], linewidth=0.3)
# plt.axis('equal')
plt.plot([70, 90], [-3, -3], "k", linewidth=2)
plt.text(80, -10, "20 cm", horizontalalignment="center")
plt.xlim(0, 100)
plt.axis("off")

plt.subplot(212)
cmap = matplotlib.cm.get_cmap("plasma")
for i in range(0, len(sorted_index)):
    color = cmap(i / (len(sorted_index) * 1.2))  # Factor to scale color range
    plt.plot(smoothed_curves[i] * 20 + i * 1, c=color, linewidth=0.1, rasterized=True)
plt.xlim(0, 40)
plt.axis("off")
plt.savefig("../../output_REM/example_place_fields.pdf")

# %% Import place cell data
results_dir = "../../output_REM/tuning"
resultsList = os.listdir(results_dir)

data_list = []

tuning_curves = []

for file_name in tqdm(resultsList):
    if (
        file_name.startswith("tuning_")
        and file_name.endswith(".h5")
        and "pv1254" not in file_name
    ):  # Only include seqResults, not replay results
        h5_file = h5py.File(os.path.join(results_dir, file_name))
        for i in range(len(h5_file["info"][()])):
            data_list.append(  # This will create one list entry per cell
                {
                    "cell_ID": i,
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    "info": h5_file["info"][i],
                    "p_value": h5_file["p_value"][i],
                    "peak_loc": h5_file["peak_loc"][i][0],
                    "peak_val": h5_file["peak_val"][i],
                    "marginal_likelihood": h5_file["marginal_likelihood"][i],
                    "trans_prob": h5_file["trans_prob"][i],
                    "total_distance_travelled": h5_file["total_distance_travelled"][()],
                }
            )
            tuning_curves.append(h5_file["tuning_curves"][i])

        # Close files
        h5_file.close()

df = pd.DataFrame(data_list)
tuning_curves = np.array(tuning_curves)
# %% Plot place cell properties
sns.histplot(
    data=df.query("condition == 'LTD1' and p_value<0.05"),
    x="peak_loc",
    color="C2",
    stat="density",
    # bins=10,
    kde=True,
)
plt.title("Place cell centroid locations")
plt.xticks([0, 20, 40], [0, 50, 100])
plt.xlabel("Location (cm)")
plt.savefig("../../output_REM/place_cell_centroids_loc.pdf")

# %%
plt.figure(figsize=(1.5, 2))
plt.subplot(412)
plt.hist(data=df, x="info", bins="auto")
plt.xlim(-0.1, 1)
plt.ylabel("Count")
plt.xlabel("")
plt.xticks([])


plt.subplot(212)
plt.scatter(data=df, x="info", y="p_value")
plt.plot([-0.1, 1], [0.05, 0.05], "C4:")
plt.xlim(-0.1, 1)
plt.xlabel("Info.")
plt.ylabel("p value")
plt.savefig("../../output_REM/info_distribution.pdf")
plt.tight_layout()
