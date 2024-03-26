# %% Imports
import numpy as np
import os
import pandas as pd
import seaborn as sns
import pingouin as pg
import matplotlib.pyplot as plt
import h5py
import yaml
from tqdm import tqdm
from utils.helperFunctions import load_data
from utils.helperFunctions import extract_seqReplay_score
from utils.helperFunctions import extract_H

plt.style.use("plot_style.mplstyle")

# %% Load parameters
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

np.random.seed(params["seed"])

# Load example recording
mouse = 'pv1069'
condition = 'LTD1'

# %% Representative seqNMF replay example

data_LT = load_data(mouse=mouse,
                    condition=condition,
                    state='wake',
                    params=params)
data_REMpost = load_data(mouse=mouse,
                    condition=condition,
                    state='REMpost',
                    params=params)

with h5py.File(os.path.join(params['path_to_output'],"neuron_selection", f'selected_neurons_{condition}_{mouse}.h5'),'r') as f:
    selected_neurons = f['place_cells'][()]

# %% Run SeqNMF
params['numNeurons']=128
seqReplay_scores, seqReplay_pvalues, seqReplay_locs, W_ref = extract_seqReplay_score(data_LT['binaryData'][:,selected_neurons],
                                                                                     data_REMpost['binaryData'][:,selected_neurons],
                                                                                     params)

# %% Sort sequences
seqSortingIndex = np.argsort(np.argmax(np.concatenate((W_ref[:,0,:],W_ref[:,1,:]),axis=1),axis=1))

# %% Focus on place cells, sort cells according to place field locations
LT_PCs_data = data_LT['binaryData'][:,selected_neurons]
REM_PCs_data = data_REMpost['binaryData'][:,selected_neurons]

# %% Sort by place field (LEGACY)
with h5py.File(os.path.join(params['path_to_output'],"tuning", f'tuning_{condition}_{mouse}.h5'),'r') as f:
    placeFieldLocations = f['peak_loc'][()][:,0]
placeFieldLocations = placeFieldLocations[selected_neurons]
PCsortingIndex = np.argsort(placeFieldLocations)

# Select running epochs only
LT_PCs_data = LT_PCs_data[data_LT['running_ts']]
LT_position = data_LT['position'][data_LT['running_ts'],0]

#%%
H = extract_H(W_ref, REM_PCs_data)
#%%
SF1 = REM_PCs_data*H[0][:,None]
SF2 = REM_PCs_data*-H[1][:,None]
sequence_data = SF1+SF2

#%%
plt.figure()
plt.subplot(121)
plt.imshow(sequence_data[:,seqSortingIndex].T,
           interpolation='none',
           aspect='auto',
           cmap='RdBu',
           vmin=-10,
           vmax=10)
plt.xlim(4875,4975)
plt.axis('off')
plt.colorbar()
plt.title('S1')

plt.subplot(122)
plt.imshow(sequence_data[:,seqSortingIndex].T,
           interpolation='none',
           aspect='auto',
           cmap='RdBu',
           vmin=-10,
           vmax=10)
plt.xlim(5450,5550)
plt.axis('off')
plt.colorbar()
plt.title('S2')
plt.tight_layout()
plt.savefig('../../output_REM/example_seqNMF_replay.pdf')

# %% Extract replay statistics
results_dir = '../../output_REM/seqNMF'
resultsList=os.listdir(results_dir)

data_list = []
for file_name in resultsList:
    if file_name.startswith('seqReplayResults_') and file_name.endswith('.h5'): # Only include seqResults, not replay results
        h5_file = h5py.File(os.path.join(results_dir,file_name))
        data_list.append( #This will create one list entry per cell
                {
                    'state_ref':h5_file['state_ref'][()].decode("utf-8"),
                    'state_pred':h5_file['state_pred'][()].decode("utf-8"),
                    'condition':h5_file['condition'][()].decode("utf-8"),
                    'mouse':h5_file['mouse'][()].decode("utf-8"),
                    'S1_numSeqs':h5_file['S1_numSeqs'][()],
                    'S2_numSeqs':h5_file['S2_numSeqs'][()],
                    'S1_score':h5_file['S1_score'][()],
                    'S2_score':h5_file['S2_score'][()],
                    'S1_pvalue':h5_file['S1_pvalue'][()],
                    'S2_pvalue':h5_file['S2_pvalue'][()],
                }
            )

        # Close files
        h5_file.close()

df_replay = pd.DataFrame(data_list)

#%%
df_replay_stats = df_replay
df_replay_stats=df_replay_stats.melt(id_vars=['state_ref','state_pred','mouse', 'condition'],value_name='numSeqs',value_vars=['S1_numSeqs', 'S2_numSeqs'],var_name='seqType')

# %% 
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df_replay_stats.query("condition == 'LTD1' and state_ref == 'wake' and state_pred == 'REMpost'"),
    x='seqType',
    hue='seqType',
    y='numSeqs',
    palette=(['C0','C4']),
    errorbar='se',
    capsize=.2
)

sns.stripplot(
    data=df_replay_stats.query("condition == 'LTD1' and state_ref == 'wake' and state_pred == 'REMpost'"),
    x='seqType',
    y='numSeqs',
    # hue='seqType',
    color='gray',
    size=2,
    #dodge=True,
    legend=False
)

#plt.xticks([0,1],['S1', 'S2'])
plt.ylabel('Num. significant \nsequences')
plt.xticks([0,1],['S1','S2'])
plt.xlabel('')
plt.ylim(0,22)
plt.xticks(rotation=90)
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/PCs_numSeqReplay_LTD1.pdf')

#%% STATS
pg.kruskal(data=df_replay_stats.query("condition == 'LTD1' and state_ref == 'wake' and state_pred == 'REMpost'"),
         dv='numSeqs',
         between='seqType')
# %%
