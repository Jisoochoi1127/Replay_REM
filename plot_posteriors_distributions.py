# %% Imports
import numpy as np
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import h5py
import yaml
from tqdm import tqdm

plt.style.use("plot_style.mplstyle")

# %% Load parameters
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

np.random.seed(params["seed"])

# %% List posteriors results
results_dir = params['path_to_output']+'/posterior_probs'
resultsList=os.listdir(results_dir)

# %% Load data
data_list = []
for file_name in tqdm(resultsList):
    if file_name.startswith('posterior'):
        with h5py.File(os.path.join(results_dir, file_name), 'r') as f:
            mouse = f['mouse'][()].decode("utf-8")
            condition = f['condition'][()].decode("utf-8")

            for state in ['REMpre', 'REMpost']:
                posterior_probs = f[f'{state}_posterior_probs'][()]
                avg_posteriors = np.mean(posterior_probs,axis=0)
                median_posteriors = np.median(posterior_probs,axis=0)
                max_posterior = np.argmax(posterior_probs,axis=1)

                for i, val in enumerate(avg_posteriors):
                    data_list.append(
                    {
                        'Mouse': mouse,
                        'Condition': condition,
                        'State': state,
                        'Location': (i+params['spatialBinSize']/2)*params['spatialBinSize'], # Convert to cm
                        'Average posterior': val,
                        'Median posterior': median_posteriors[i],
                        'MAP': max_posterior[i]
                    }
        )

df = pd.DataFrame(data_list)

# %%
sns.lineplot(
    data=df,
    x='Location',
    y='MAP',
    hue='Condition',
    errorbar='se'
    )
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
# %%
sns.lineplot(
    data=df.query("Condition=='HATD1'"),
    x='Location',
    y='Median posterior',
    hue='Mouse',
    errorbar='se'
    )
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
# %%
