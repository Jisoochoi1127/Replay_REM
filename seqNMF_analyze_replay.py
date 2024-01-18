#%% Imports
import numpy as np
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('plot_style.mplstyle')

import pingouin as pg
import h5py

#%%
results_dir = '../../output_REM/'
resultsList=os.listdir(results_dir)

# %% Q1: is there replay during REM?
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
df_replay_stats.loc[df_replay_stats['S1_pvalue']>0.05,'S1_numSeqs']=0
df_replay_stats.loc[df_replay_stats['S2_pvalue']>0.05,'S2_numSeqs']=0
df_replay_stats=df_replay_stats.melt(id_vars=['state_ref','state_pred','mouse', 'condition'],value_name='numSeqs',value_vars=['S1_numSeqs', 'S2_numSeqs'],var_name='seqType')

#%% Plot number of wake sequences replayed during REMpost
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df_replay_stats.query("state_ref == 'wake' and state_pred == 'REMpost'"),
    x='condition',
    # hue='seqType',
    y='numSeqs',
    # palette=(['C0','C4']),
    errorbar='se',
    capsize=.2
)

sns.stripplot(
    data=df_replay_stats.query("state_ref == 'wake' and state_pred == 'REMpost'"),
    x='condition',
    y='numSeqs',
    # hue='seqType',
    # palette=(['C0','C4']),
    size=2,
    #dodge=True,
    legend=True
)

#plt.xticks([0,1],['S1', 'S2'])
plt.ylabel('Num. significant \nsequences')
plt.xlabel('')
plt.xticks(rotation=90)
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/numSeqReplay.pdf')
# %%
