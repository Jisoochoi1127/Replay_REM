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

#%% Internal structure analysis
# SeqNMF
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
df_replay_stats=df_replay.melt(id_vars=['state_ref','state_pred','mouse', 'condition'],value_name='numSeqs',value_vars=['S1_numSeqs', 'S2_numSeqs'],var_name='seqType')

#%%
plt.figure(figsize=(1,.75))
sns.heatmap(data=df_replay_stats.query("condition=='LTD1'").pivot_table(
    index='state_ref',
    columns='state_pred', 
    values='numSeqs'
    ),
    cmap='magma',
    vmin=0,
    rasterized=True,
    cbar_kws={'label':'Num. replay events'}
)
#plt.xticks([])
#plt.yticks([])
plt.xlabel('Reference')
plt.ylabel('Target')
plt.savefig('../../output_REM/intStruct_numSeqs_heatmap.pdf')

#%% Plot seqNMF preplay
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df_replay_stats.query("condition=='LTD1' and state_ref=='REMpre' and state_pred=='wake'"),
    x=1,
    y='numSeqs',
    # showfliers=False,
    color='C3',
    errorbar='se',
    capsize=.2
)
sns.barplot(
    data=df_replay_stats.query("condition=='LTD1' and state_ref=='wake' and state_pred=='REMpost'"),
    x=2,
    y='numSeqs',
    # showfliers=False,
    color='C0',
    errorbar='se',
    capsize=.2
)
# sns.stripplot(
#     data=df_replay_stats.query("condition=='LTD1' and state_ref=='REMpre' and state_pred=='wake'"),
#     y='numSeqs',
#     size=1,
#     dodge=True,
#     legend=False
# )
plt.xticks([0,1],['Preplay', 'Replay'],rotation=90)
plt.ylabel('Num. rigid sequences')
plt.xlabel('')
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/seqNMF_preplay_numSeqs.pdf')

#%%



# Assembly patterns


#%% Preplay analysis
results_dir = params['path_to_output']+"/bayesian_replay"
resultsList = os.listdir(results_dir)

#%%
data_list = []
num_event_list = []

for file_name in tqdm(resultsList):
    if (
        file_name.startswith("bayesian_replay_")
        and file_name.endswith(".h5")
        and "pv1254" not in file_name # Exclude pv1254
    ):
        h5_file = h5py.File(os.path.join(results_dir, file_name))
        num_event_list.append(
            {
                "mouse": h5_file["mouse"][()].decode("utf-8"),
                "condition": h5_file["condition"][()].decode("utf-8"),
                "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                "numReplayEvents": len(h5_file["replay_locs"][()])
            }
        )
        replay_freqs = 1/(np.diff(h5_file['replay_locs'][()]/params['sampling_frequency']))
        replay_freqs = np.append(replay_freqs,0)
        for i in range(len(h5_file["replay_locs"][()])):
            data_list.append(  # This will create one list entry per cell
                {
                    "eventID": i,
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                    "replayEventTime": h5_file['replay_locs'][i]/params['sampling_frequency'],
                    "replayFrequency": replay_freqs[i],
                    "replayEventScore": h5_file['replay_score'][i],
                    "replayEventJumpiness": h5_file['replay_jumpiness'][i],
                    "replayEventLength": h5_file['replay_length'][i],
                    "replayEventSlope": abs(h5_file['replay_slope'][i])
                }
            )
        
        # Close files
        h5_file.close()

df = pd.DataFrame(data_list)
df_numEvents = pd.DataFrame(num_event_list)

# %% Plot numEvents preplay vs replay
sns.boxenplot(
    data=df.query("replayEventJumpiness>0"),
    y='replayFrequency',
    x='Type',
    palette=['C3','C0'],
    showfliers=False
    #errorbar='se'
)
# plt.title('REM post')
plt.xlabel('')
plt.yscale('log')
plt.ylabel('Replay frequency (Hz)')
plt.savefig("../../output_REM/REMpre_vs_post_replayFrequency.pdf")

# %% Plot slope preplay vs replay
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df_numEvents,
    y='numReplayEvents',
    x='Type',
    palette=['C3','C0'],
    #showfliers=False
    errorbar='se',
    capsize=.2
    
)
# plt.title('REM post')
plt.xlabel('')
plt.xticks([0,1],['Preplay', 'Replay'],rotation=90)
plt.ylabel('Num. flexible sequences')
plt.savefig("../../output_REM/REMpre_vs_post_numEvents.pdf")


#%%
sns.boxenplot(
    data=df.query("replayEventJumpiness>0"),
    y='replayEventSlope',
    x='Type',
    palette=['C3','C0'],
    showfliers=False
)
plt.title('REM post')
plt.ylabel('Replay event\nspeed (cm.s$^{-1}$)')
plt.xlabel('')
plt.savefig("../../output_REM/REMpre_vs_post_slope.pdf")

#%%
plt.figure(figsize=(1,.25))
sns.histplot(
    data=df.query("Type=='replay' and replayEventJumpiness>0"),
    x='replayEventScore',
)
# plt.title('REM post')
plt.xlabel('Replay score (R$^{2}$)')
plt.ylabel('N')
plt.savefig("../../output_REM/REMpost_replayScores.pdf")

#%% DESCRIPTIVES
mean_S1 = df.query("Type=='replay' and replayEventJumpiness>0")['replayEventScore'].mean()
SEM_S1 = df.query("Type=='replay' and replayEventJumpiness>0")['replayEventScore'].sem()
print(f'{mean_S1} +/- {SEM_S1}')

