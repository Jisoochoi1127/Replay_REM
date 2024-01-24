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

#%% Load parameters
with open('params.yaml','r') as file:
    params = yaml.full_load(file)

np.random.seed(params['seed'])

# %% Load example statistics
data = load_data(mouse='pv1069', condition='LTD1', state='wake',params=params)

# %%
with h5py.File("../../output_REM/tuning/tuning_LTD1_pv1069.h5") as f:
    info = f["info"][()]
    p_values = f["p_value"][()]
    peak_loc = f["peak_loc"][()][:,0]
    peak_val = f["peak_val"][()] 
    tuning_curves = f["tuning_curves"][()]

# %% Select high information neurons
selected_neurons = np.where((p_values<.05) & (peak_val>.1))[0]
sorted_index = np.argsort(peak_loc[selected_neurons])

# %% Filter traces
smoothed_curves = np.zeros((len(sorted_index),tuning_curves.shape[1]))
for neuron_i in range(len(sorted_index)):
    smoothed_curves[neuron_i,:] = gaussian_filter1d(tuning_curves[selected_neurons[sorted_index]][neuron_i],sigma=1.5)
smoothed_curves /= np.max(smoothed_curves, axis=0) 
# %% Plot example exploration and activity
plt.figure(figsize=(1,3))
plt.subplot(16,1,8)
plt.plot(data['position'][:,0],data['position'][:,1],
         linewidth=.3) 
# plt.axis('equal')
plt.xlim(0,100)
plt.axis('off')

plt.subplot(212)
cmap = matplotlib.cm.get_cmap('viridis')
for i in range(0,len(sorted_index)):
    color = cmap(i/(len(sorted_index)*1.2)) # Factor to scale color range
    plt.plot(smoothed_curves[i]*20+i/1,
            c=color,
            linewidth=.3, rasterized=True)
plt.axis('off')
# plt.imshow(smoothed_curves,
#            aspect='auto',
#            interpolation='none',
        #    cmap='magma') 
# %% Plot activity
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
    data=df.query("p_value<0.05"),
    x="peak_loc",
    y="condition",
    cbar=True,
    cbar_kws={"label": "Number of centroids"},
    cmap="magma",
)
plt.title("Place cell centroid locations")
plt.xlabel("Location (cm)")
plt.show()
# %%
plt.figure()
sns.histplot(
    data=df.query("p_value<0.05 and condition=='HATD1' or condition=='HATD5'"),
    x="peak_loc",
    y="peak_val",
    cbar=True,
    cbar_kws={"label": "Number of centroids"},
)
plt.show()
# %%
plt.figure()
sns.histplot(data=df.query("p_value<0.05"), x="info")
plt.show()
# %%
plt.figure()
plt.scatter(data=df, x="info", y="p_value")
plt.plot([0, 0.5], [0.05, 0.05], "C4:")
plt.xlabel("Info.")
plt.ylabel("p value")
plt.show()

# %%