# %% Imports
import os
import numpy as np
import pandas as pd
import seaborn as sns
import pingouin as pg
import matplotlib.pyplot as plt
import yaml
import h5py
from tqdm import tqdm

plt.style.use("plot_style.mplstyle")

#%% Load params
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

#%% Assemblies TODO

#%% Bayesian replay
results_dir = params['path_to_output']+"/bayesian_replay"
resultsList = os.listdir(results_dir)

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
        for i in range(len(h5_file["replay_locs"][()])):
            data_list.append(  # This will create one list entry per cell
                {
                    "eventID": i,
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                    "replayEventTime": h5_file['replay_locs'][i]/params['sampling_frequency'],
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

#%% Plot number of events per condition
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df_numEvents.query("Type=='replay' and condition=='LTD1' or condition=='LTD5'"),
    x='condition',
    y='numReplayEvents',
    order=['LTD1', 'LTD5'],
    palette=(['C3','gray']),
    #showfliers=False
    errorbar='se',
    capsize=.2
)

plt.ylabel('Num. replay\nevents')
plt.xticks([0,1],['Novel','Familiar'])
plt.xlabel('')
plt.xticks(rotation=90)
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/novelty_bayesian_numEvents.pdf')

#%% DESCRIPTIVES
mean_novelty = df_numEvents.query("condition=='LTD1' and Type=='replay'")['numReplayEvents'].mean()
SEM_novelty = df_numEvents.query("condition=='LTD1' and Type=='replay'")['numReplayEvents'].sem()
print(f'Novel: {mean_novelty} +/- {SEM_novelty}')

mean_familiar = df_numEvents.query("condition=='LTD5' and Type=='replay'")['numReplayEvents'].mean()
SEM_familiar = df_numEvents.query("condition=='LTD5' and Type=='replay'")['numReplayEvents'].sem()
print(f'Familiar: {mean_familiar} +/- {SEM_familiar}')

#%% STATS
pg.anova(
    data=df_numEvents.query("condition=='LTD1' or condition=='LTD5'"),
    dv='numReplayEvents',
    between='condition',

)

#%% Plot replay score as a function of novelty
plt.figure(figsize=(.75,1))
sns.boxenplot(
    data=df.query("Type=='replay' and condition=='LTD1' or condition=='LTD5'"),
    x='condition',
    y='replayEventScore',
    order=['LTD1', 'LTD5'],
    palette=(['C3','gray']),
    showfliers=False
    #errorbar='se',
    #capsize=.2
)

plt.ylabel('Replay score ($R^{2}$)')
plt.xticks([0,1],['Novel','Familiar'])
plt.xlabel('')
plt.xticks(rotation=90)
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/novelty_bayesian_replay_scores.pdf')

#%% DESCRIPTIVES
mean_novelty = df.query("Type=='replay' and condition=='LTD1'")['replayEventScore'].mean()
SEM_novelty = df.query("Type=='replay' and condition=='LTD1'")['replayEventScore'].sem()
print(f'Novel: {mean_novelty} +/- {SEM_novelty}')

mean_familiar = df.query("Type=='replay' and condition=='LTD5'")['replayEventScore'].mean()
SEM_familiar = df.query("Type=='replay' and condition=='LTD5'")['replayEventScore'].sem()
print(f'Familiar: {mean_familiar} +/- {SEM_familiar}')

#%% STATS
pg.anova(
    data=df.query("Type=='replay' and condition=='LTD1' or condition=='LTD5'"),
    dv='replayEventScore',
    between='condition',
)

#%% Plot jumpiness as a function of novelty
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df.query("Type=='replay' and condition=='LTD1' or condition=='LTD5'"),
    x='condition',
    y='replayEventJumpiness',
    order=['LTD1', 'LTD5'],
    palette=(['C3','gray']),
    #showfliers=False
    errorbar='se',
    capsize=.2
)

plt.ylabel('Max. jumpiness (cm)')
plt.xticks([0,1],['Novel','Familiar'])
plt.xlabel('')
plt.xticks(rotation=90)
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/novelty_bayesian_replay_jumpiness.pdf')

#%% DESCRIPTIVES
mean_novelty = df.query("Type=='replay' and condition=='LTD1'")['replayEventJumpiness'].mean()
SEM_novelty = df.query("Type=='replay' and condition=='LTD1'")['replayEventJumpiness'].sem()
print(f'Novel: {mean_novelty} +/- {SEM_novelty}')

mean_familiar = df.query("Type=='replay' and condition=='LTD5'")['replayEventJumpiness'].mean()
SEM_familiar = df.query("Type=='replay' and condition=='LTD5'")['replayEventJumpiness'].sem()
print(f'Familiar: {mean_familiar} +/- {SEM_familiar}')

#%% STATS
pg.anova(
    data=df.query("Type=='replay' and condition=='LTD1' or condition=='LTD5'"),
    dv='replayEventJumpiness',
    between='condition',
)

#%%
plt.figure(figsize=(.75,1))
sns.boxenplot(
    data=df.query("Type=='replay' and condition=='LTD1' or condition=='LTD5'"),
    x='condition',
    y='replayEventSlope',
    order=['LTD1', 'LTD5'],
    palette=(['C3','gray']),
    showfliers=False
    #errorbar='se',
    #capsize=.2
)

plt.ylabel('Slope (cm.s$^{-1}$)')
plt.xticks([0,1],['Novel','Familiar'])
plt.xlabel('')
plt.xticks(rotation=90)
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/novelty_bayesian_replay_slope.pdf')

#%% STATS
#TODO

#%% SeqNMF
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
df_replay_stats=df_replay.melt(id_vars=['state_ref','state_pred','mouse', 'condition'],value_name='numSeqs',value_vars=['S1_numSeqs', 'S2_numSeqs'],var_name='seqType')
df_replay_scores=df_replay.melt(id_vars=['state_ref','state_pred','mouse', 'condition'],value_name='seqScore',value_vars=['S1_score', 'S2_score'],var_name='seqType')

# %% Plot num. sequences vs novelty
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df_replay_stats.query("condition == 'LTD1' or condition == 'LTD5' and state_ref == 'wake' and state_pred == 'REMpost'"),
    x='condition',
    y='numSeqs',
    palette=(['C3','gray']),
    # showfliers=False
    errorbar='se',
    capsize=.2
)

#plt.xticks([0,1],['S1', 'S2'])
plt.ylabel('Num. significant \nsequences')
plt.xticks([0,1],['Novel','Familiar'])
plt.xlabel('')
plt.xticks(rotation=90)
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/novely_seqnmf_numSeqs.pdf')

#%% STATS
pg.rm_anova(data=df_replay_stats.query("condition == 'LTD1' or condition == 'LTD5' and state_ref == 'wake' and state_pred == 'REMpost'"),
         dv='numSeqs',
         within='condition',
         subject='mouse'
         )

# %% Plot num. sequences vs novelty
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df_replay_scores.query("condition == 'LTD1' or condition == 'LTD5' and state_ref == 'wake' and state_pred == 'REMpost'"),
    x='condition',
    y='seqScore',
    palette=(['C3','gray']),
    # showfliers=False
    errorbar='se',
    capsize=.2
)

#plt.xticks([0,1],['S1', 'S2'])
plt.ylabel('Sequence score (z)')
plt.xticks([0,1],['Novel','Familiar'])
plt.xlabel('')
plt.xticks(rotation=90)
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/novely_seqnmf_seqScore.pdf')

#%% STATS
pg.kruskal(data=df_replay_stats.query("condition == 'LTD1' and state_ref == 'wake' and state_pred == 'REMpost'"),
         dv='numSeqs',
         between='seqType')