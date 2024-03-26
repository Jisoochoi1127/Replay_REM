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
from utils.helperFunctions import extract_seqReplay_score

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

#%% TEMP TEST
from seqnmf import seqnmf
W_ref, _, _, _, _ = seqnmf(
        LT_PCs_data.T,
        K=2,
        L=100,
        Lambda=0.0001,
        max_iter=params['maxIters'],
        lambda_OrthW=1, # 1 enforces event-bases sequences.
        lambda_L1W=1,
        shift=False
        )

plt.figure(figsize=(.5,.75))
plt.subplot(121)
plt.imshow(W_ref[PCsortingIndex,0,:],
           interpolation='none',
           aspect='auto',
           cmap='gray_r',
            vmin=0,
            vmax=.5
           )
plt.axis('off')
plt.title('S$_{1}$')
plt.subplot(122)
plt.imshow(W_ref[PCsortingIndex,1,:],
           interpolation='none',
           aspect='auto',
           cmap='gray_r',
            vmin=0,
            vmax=.5
           )
plt.axis('off')
plt.title('S$_{2}$')

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

# %% Plot sequences sorted by place field location
plt.figure(figsize=(.5,.75))
plt.subplot(121)
plt.imshow(W_ref[PCsortingIndex,0,:],
           interpolation='none',
           aspect='auto',
           cmap='gray_r',
        #    vmin=0,
        #    vmax=.5
           )
plt.axis('off')
plt.title('S$_{1}$')
plt.subplot(122)
plt.imshow(W_ref[PCsortingIndex,1,:],
           interpolation='none',
           aspect='auto',
           cmap='gray_r',
        #    vmin=0,
        #    vmax=.5
           )
plt.axis('off')
plt.title('S$_{2}$')

# %% Plot results (i.e. recording, with location of replayed sequences)
plt.figure(figsize=(3,1.33))
plt.subplot(221)
plt.imshow(LT_PCs_data[:,seqSortingIndex].T,
           interpolation='none',
           aspect='auto',
           cmap='gray_r')
plt.xticks([])
plt.xlim(0,len(LT_PCs_data))
plt.yticks([0,128],[128,0])
plt.ylabel('Neuron ID')
plt.title('Wakefulness')

plt.subplot(2,16,9)
plt.imshow(REM_PCs_data[:,seqSortingIndex].T,
           interpolation='none',
           aspect='auto',
           cmap='gray_r')
plt.xlim(50,100)
plt.axis('off')
# plt.plot([75,90],[128,128],
#          'k',
#          linewidth=2)
# plt.text(65,168,'500 ms',
#          horizontalalignment='center',
#          verticalalignment='bottom')

plt.title('REM post')
    
plt.subplot(427)
plt.plot(-LT_position)
plt.xlim(0,len(LT_position))
plt.yticks([-100,0],[0,100])
plt.ylabel('Position (cm)')
plt.title('Linear track (runs only)')
plt.xticks([])
plt.plot([2000,2600],[-120,-120],
         'k',
         linewidth=2)
plt.text(2300,-180,'20 s',
         horizontalalignment='center',
         verticalalignment='bottom')
plt.savefig('../../output_REM/example_replay_seq2.pdf')

#%%
##Optional
# plt.subplot(428)
plt.figure(figsize=(4,1))
plt.imshow(REM_PCs_data[:,PCsortingIndex].T,
           interpolation='none',
           aspect='auto',
           cmap='gray_r')
for event in seqReplay_locs[0]:
    plt.fill_between([event-int(params['L']/2), event+int(params['L']/2)],
            [params['numNeurons'],params['numNeurons']],
            facecolor='C0',
            alpha=.1,
            )
for event in seqReplay_locs[1]:
    plt.fill_between([event-int(params['L']/2), event+int(params['L']/2)],
            [params['numNeurons'],params['numNeurons']],
            facecolor='C4',
            alpha=.1,
            )
plt.xlim(0,len(REM_PCs_data))
plt.xticks([])
plt.yticks([0,params['numNeurons']],[params['numNeurons'],0])

#%% Replay events, concatenated
replayEvents = np.zeros((0,params['numNeurons']))
for event in seqReplay_locs[0]:
    replayEvents = np.append(replayEvents,
                             REM_PCs_data[event-int(params['L']/2):event+int(params['L']/2)],
                             axis=0)

for event in seqReplay_locs[1]:
    replayEvents = np.append(replayEvents,
                             REM_PCs_data[event-int(params['L']/2):event+int(params['L']/2)],
                             axis=0)
plt.imshow(replayEvents[:,sortingIndex].T,
           interpolation='none',
           aspect='auto',
           cmap='gray_r')
plt.xlim(0,len(replayEvents))
plt.xticks(np.arange(0,len(replayEvents),params['L']),[])
plt.yticks([0,128],[])
plt.plot([0,len(seqReplay_locs[0])*params['L']],
         [0,0],
         'C0',
         linewidth=2)
plt.plot([len(seqReplay_locs[0])*params['L'],len(replayEvents)],
         [0,0],
         'C4',
         linewidth=2)

# plt.plot([6000,6600],[params['numNeurons'],params['numNeurons']],
#          'k',
#          linewidth=2)
# plt.text(6300,params['numNeurons']+80,'20 s',
#          horizontalalignment='center',
#          verticalalignment='bottom')

# plt.savefig('../../output_REM/example_replay_seqNMF.pdf')

# plt.xlim(5000,6000)
#%% Only plot replayed sequences
# Plot additional stats (score, num. sig. sequences)

# %% Extract replay statistics
results_dir = '../../output_REM/replay_results'
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
    x='mouse',
    hue='seqType',
    y='numSeqs',
    palette=(['C0','C4']),
    errorbar='se',
    capsize=.2
)

# sns.stripplot(
#     data=df_replay_stats.query("condition == 'LTD1' and state_ref == 'wake' and state_pred == 'REMpost'"),
#     x='mouse',
#     y='numSeqs',
#     # hue='seqType',
#     # palette=(['C3','C4']),
#     size=2,
#     #dodge=True,
#     legend=True
# )

#plt.xticks([0,1],['S1', 'S2'])
plt.ylabel('Num. significant \nsequences')
plt.xlabel('')
plt.ylim(0,22)
plt.xticks(rotation=90)
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/PCs_numSeqReplay_LTD1.pdf')