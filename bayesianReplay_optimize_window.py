# %% Imports
import h5py
import yaml
import os
from tqdm import tqdm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from utils.bayesian_replay import extract_linear_replay
plt.style.use("plot_style.mplstyle")
#%% Load parameters
with open('params.yaml','r') as file:
    params = yaml.full_load(file)

results_dir = params['path_to_output']+'/posterior_probs'
resultsList=os.listdir(results_dir)

#%%
state='REMpost'
windowSizeList = [8,16,32,64,128]
replayLocs_dict = []
for i, windowSize in enumerate(tqdm(windowSizeList)):
    params['windowSize']=windowSize
    for file_name in tqdm(resultsList):
        if file_name.startswith('posterior'):
            with h5py.File(os.path.join(results_dir, file_name), 'r') as f:
                mouse = f['mouse'][()].decode("utf-8")
                condition = f['condition'][()].decode("utf-8")

                posterior_probs = f[f'{state}_posterior_probs'][()]

                replayLocs, _, _, _ = extract_linear_replay(posterior_probs, params)
                replayLocs_dict.append(
                    {
                        'mouse': mouse,
                        'condition': condition,
                        'windowSize':windowSize,
                        'numReplayEvents':len(replayLocs)
                    }
                )

#%% Save results
df = pd.DataFrame(replayLocs_dict)
df.to_csv(params['path_to_output']+'/optimal_windowSize_BayesianReplay.csv')

#%% Plot results
sns.lineplot(
    data=df,
    x='windowSize',
    y='numReplayEvents',
    hue='condition',
    errorbar='se',
)
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig(params['path_to_output']+'/optimal_windowSize_BayesianReplay_perCondition.pdf')
# %%
